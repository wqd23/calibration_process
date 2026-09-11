# -*- coding:utf-8 -*-
"""L1 parquet cache for the raw-parse layer.

Caches the *raw* output of ``readSci``/``readHK`` (a ``Dict[str, np.ndarray]``)
under ``data/{ver}/l1_cache/`` and re-runs the unchanged post-processing in the
``single_readXX`` functions on every call. Replaces the former ``@cachier`` dill
baseline.

The serializer is a generic traversal: flat ``Dict`` of 1-D arrays are stored
column-per-column; 2-D arrays (e.g. 11B ``wf`` waveform samples) are stored as a
polars fixed-size ``Array`` column. Anything irregular (scalar, >2-D, mismatched
row counts, object dtype) causes a bail so the caller falls back to a fresh parse
and the pipeline never fails because of the cache.
"""
import functools
import hashlib
import inspect
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from addict import Dict

SCHEMA_VER = 1
SCHEMA_VER_PROCESSED = 2


def get_project_root(marker=("justfile", ".gitignore")):
    current = Path(os.getcwd())
    while True:
        if any((current / m).exists() for m in marker):
            return current
        if current == current.parent:
            return current
        current = current.parent


def _cache_root(ver: str):
    return get_project_root() / "data" / ver / "l1_cache"


def _make_key(ver: str, reader: str, kind: str, raw_path: str, parse_kwargs: dict):
    payload = f"{ver}|{reader}|{kind}|{raw_path}|{json.dumps(parse_kwargs, sort_keys=True)}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _extract_parse_kwargs(fn, args, kwargs):
    """Split cache-control flags (``overwrite_cache``) from parse-affecting kwargs."""
    sig = inspect.signature(fn)
    params = list(sig.parameters)
    parse = {k: v for k, v in kwargs.items() if k != "overwrite_cache"}
    if len(args) > 1:
        for name, val in zip(params[1:], args[1:]):
            parse[name] = val
    return parse


def _column_key(kind: str):
    return "events.parquet" if kind == "sci" else "tel.parquet"


# --- generic traversal serializer (auto-adapt to any raw parse) ---------------

def serialize(d: dict):
    """Flatten ``Dict[str, np.ndarray]`` into cache-friendly column values.

    Returns ``(cols, dtypes, rows, shapes)`` where ``cols`` maps each key to a
    1-D array or, for 2-D arrays, a list of row sub-arrays; ``rows`` is the
    (common) number of records. Raises for anything the cache cannot represent.
    """
    cols = {}
    dtypes = {}
    shapes = {}
    rows = None
    for k, v in d.items():
        a = np.asarray(v)
        if a.ndim == 0 or a.ndim > 2:
            raise ValueError(f"column {k!r} has ndim {a.ndim}; cache bails (fallback)")
        if rows is None:
            rows = a.shape[0]
        elif a.shape[0] != rows:
            raise ValueError(f"column {k!r} row count differs; cache bails (fallback)")
        if a.ndim == 2:
            cols[k] = [a[i] for i in range(a.shape[0])]
        else:
            cols[k] = a
        dtypes[k] = str(a.dtype)
        shapes[k] = list(a.shape)
    return cols, dtypes, (rows if rows is not None else 0), shapes


def deserialize(cols: dict):
    """Rebuild an ``addict.Dict`` so attribute access keeps working downstream."""
    return Dict(cols)


# --- polars I/O choke point (the single swappable point) ----------------------

def write_frame(cols: dict, path):
    df = pl.DataFrame({k: pl.Series(k, a) for k, a in cols.items()})
    df.write_parquet(path, compression="zstd")


def read_frame(path, meta=None) -> dict:
    """Read a frame back into ``dict[str, np.ndarray]``.

    Uses ``meta`` (``dtypes``/``shapes``) to restore original dtypes/2-D shapes
    whenever polars promotes a width or nests an array.
    """
    df = pl.read_parquet(path)
    dtypes = (meta or {}).get("dtypes", {})
    shapes = (meta or {}).get("shapes", {})
    out = {}
    for c in df.columns:
        shp = shapes.get(c)
        if shp and len(shp) > 1:
            lst = df[c].to_list()
            dt = np.dtype(dtypes[c]) if c in dtypes else None
            out[c] = np.asarray(lst, dtype=dt)
        else:
            arr = df[c].to_numpy()
            if c in dtypes:
                want = np.dtype(dtypes[c])
                if arr.dtype != want:
                    arr = arr.astype(want)
            out[c] = arr
    return out


# --- index.json (raw file -> cache mapping) -----------------------------------

def _append_index(ver: str, record: dict):
    root = _cache_root(ver)
    index_path = root / "index.json"
    existing = []
    if index_path.exists():
        try:
            existing = json.loads(index_path.read_text())
        except (json.JSONDecodeError, OSError):
            existing = []
    if not isinstance(existing, list):
        existing = []
    replaced = False
    for i, rec in enumerate(existing):
        if rec.get("key") == record.get("key") and rec.get("kind") == record.get("kind"):
            existing[i] = record
            replaced = True
            break
    if not replaced:
        existing.append(record)
    root.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(existing, ensure_ascii=False, indent=2))


# --- decorator ---------------------------------------------------------------

def with_l1_cache(ver: str, reader: str, kind: str):
    """Wrap a raw-parse function with an L1 parquet cache.

    ``kind`` is ``"sci"``/``"tel"`` and selects the parquet file name. The
    wrapped function must return a *fresh* object each call (the downstream
    post-processing mutates the Dict), so a cache hit reconstructs a new Dict.
    """
    if kind not in ("sci", "tel"):
        raise ValueError(f"invalid kind {kind!r}")

    def decorate(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            overwrite = kwargs.pop("overwrite_cache", False)
            raw_path = args[0]
            parse_kwargs = _extract_parse_kwargs(fn, args, kwargs)
            key = _make_key(ver, reader, kind, raw_path, parse_kwargs)

            root = _cache_root(ver)
            cache_dir = root / "cache" / key
            parquet_path = cache_dir / _column_key(kind)
            meta_path = cache_dir / "meta.json"

            if not overwrite and parquet_path.exists() and meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text())
                    if meta.get("schema_ver") == SCHEMA_VER and meta.get("kind") == kind:
                        cols = read_frame(parquet_path, meta)
                        return deserialize(cols)
                except Exception:
                    pass  # corrupt/mismatched cache -> fresh parse below

            result = fn(*args, **kwargs)
            try:
                cols, dtype, rows, shapes = serialize(result)
                cache_dir.mkdir(parents=True, exist_ok=True)
                write_frame(cols, parquet_path)
                meta = {
                    "schema_ver": SCHEMA_VER,
                    "reader": reader,
                    "ver": ver,
                    "kind": kind,
                    "dtypes": dtype,
                    "shapes": shapes,
                    "ncols": len(cols),
                    "rows": rows,
                    "parse_kwargs": parse_kwargs,
                }
                meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2))
                rel_parquet = os.path.relpath(parquet_path, get_project_root())
                record = {
                    "raw_ptr": os.fspath(raw_path),
                    "kind": kind,
                    "reader": reader,
                    "ver": ver,
                    "key": key,
                    "schema_ver": SCHEMA_VER,
                    "parse_kwargs": parse_kwargs,
                }
                record["events"] = rel_parquet if kind == "sci" else None
                record["tel"] = rel_parquet if kind == "tel" else None
                _append_index(ver, record)
            except Exception:
                pass  # never let the cache put the pipeline at risk
            return result

        return wrapper

    return decorate


# --- structured serializer for post-processed output ("four-channel stack") ---
#
# The monolithic readers (04/07/09/03B/05B) return already-processed
# ``(sciExtracted, telExtracted)`` dicts whose channel quantities are lists of 4
# arrays with *different* per-channel lengths ("jagged"), plus shared flat
# arrays and a few empty entries.  To put them in parquet (rectangular tables)
# we go one level up: concatenate each key's 4 channel arrays and record the
# channel boundaries; on read, split them back.  This is lossless and needs no
# change to the scientific processing.

def _classify_section(d: dict):
    """Split a section dict into channel buckets (by sublens), flat buckets (by
    shape) and empty/scalar descriptors. Raises for anything unsupported."""
    channel = {}
    flat = {}
    empties = {}
    for k, v in d.items():
        if isinstance(v, np.ndarray):
            if v.ndim != 1:
                raise ValueError(f"flat key {k!r} has ndim {v.ndim}; cache bails")
            flat.setdefault(v.shape, []).append((k, v))
        elif isinstance(v, list):
            if len(v) == 0:
                empties[k] = {"kind": "empty_list"}
            elif len(v) == 4 and all(isinstance(x, (np.ndarray, list)) for x in v):
                arrs = [np.asarray(x) for x in v]
                if any(a.ndim == 0 for a in arrs):
                    raise ValueError(f"channel key {k!r} has 0-d entries; cache bails")
                sublens = tuple(int(a.shape[0]) for a in arrs)
                channel.setdefault(sublens, []).append((k, arrs))
            else:
                raise ValueError(f"list key {k!r} len={len(v)}; cache bails")
        elif v is None:
            empties[k] = {"kind": "none"}
        else:
            raise ValueError(f"key {k!r} type {type(v).__name__}; cache bails")
    return channel, flat, empties


def serialize_processed(result):
    """Serialize a ``(sci, tel)`` return value into a section descriptor.

    Returns ``(section_d, section_meta)`` where ``section_d`` maps file name ->
    column dict and ``section_meta`` is JSON-able structural metadata.  Raises
    (so the caller falls back to a fresh call) for unsupported structures.
    """
    if not (isinstance(result, tuple) and len(result) == 2):
        raise ValueError("processed result is not a 2-tuple; cache bails")

    section_d = {}
    section_meta = {}
    for name, d in zip(("sci", "tel"), result):
        channel, flat, empties = _classify_section(d)
        frames = []
        for i, (sublens, items) in enumerate(channel.items()):
            cols = {}
            dtypes = {}
            for k, arrs in items:
                cat = np.concatenate(arrs)
                cols[k] = cat
                dtypes[k] = str(cat.dtype)
            cols["__channel__"] = np.repeat(np.arange(4, dtype=np.int8), np.asarray(sublens))
            fname = f"{name}__chan__{i}.parquet"
            section_d[fname] = cols
            frames.append({
                "file": fname,
                "kind": "channel",
                "sublens": [int(x) for x in sublens],
                "keys": [k for k, _ in items],
                "dtypes": dtypes,
            })
        for i, (shape, items) in enumerate(flat.items()):
            cols = {k: v for k, v in items}
            fname = f"{name}__flat__{i}.parquet"
            section_d[fname] = cols
            frames.append({
                "file": fname,
                "kind": "flat",
                "shape": [int(x) for x in shape],
                "keys": [k for k, _ in items],
                "dtypes": {k: str(v.dtype) for k, v in items},
            })
        section_meta[name] = {
            "files": frames,
            "empties": empties,
        }
    return section_d, section_meta


def deserialize_processed(section_d, section_meta):
    """Rebuild the ``(sci, tel)`` dicts from in-memory frames + metadata."""
    out = {}
    for name, meta in section_meta.items():
        d = {}
        for fr in meta.get("files", []):
            cols = section_d[fr["file"]]
            if fr["kind"] == "channel":
                bounds = np.cumsum(fr["sublens"])[:-1]
                for k in fr["keys"]:
                    arr = np.asarray(cols[k])
                    want = np.dtype(fr["dtypes"][k])
                    if arr.dtype != want:
                        arr = arr.astype(want)
                    d[k] = [np.asarray(part) for part in np.split(arr, bounds)]
            else:
                for k in fr["keys"]:
                    arr = np.asarray(cols[k])
                    want = np.dtype(fr["dtypes"][k])
                    if arr.dtype != want:
                        arr = arr.astype(want)
                    d[k] = arr
        for k, desc in meta.get("empties", {}).items():
            kind = desc.get("kind")
            if kind == "empty_list":
                d[k] = []
            elif kind == "none":
                d[k] = None
            else:
                raise ValueError(f"unknown empty kind {kind!r}")
        out[name] = d
    return out["sci"], out["tel"]


def _read_processed(cache_dir: Path, meta: dict):
    section_d = {}
    for name, smeta in meta["sections"].items():
        for fr in smeta["files"]:
            section_d[fr["file"]] = read_frame(cache_dir / fr["file"])
    return deserialize_processed(section_d, meta["sections"])


def with_l1_cache_processed(ver: str, reader: str, kind: str):
    """Wrap a ``(sci, tel)``-returning reader with an L1 parquet cache.

    Caches the *final* processed output (the reader's ``(sci, tel)`` return).  The
    jagged list-of-4 channel quantities are stored as channel-stacked rectangular
    tables; a cache hit rebuilds plain ``dict`` objects (matching the reader's
    return type) so downstream access is unchanged.
    """
    def decorate(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            overwrite = kwargs.pop("overwrite_cache", False)
            raw_path = args[0]
            parse_kwargs = _extract_parse_kwargs(fn, args, kwargs)
            key = _make_key(ver, reader, kind, raw_path, parse_kwargs)

            root = _cache_root(ver)
            cache_dir = root / "cache" / key
            meta_path = cache_dir / "meta.json"

            if not overwrite and meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text())
                    if (meta.get("schema_ver") == SCHEMA_VER_PROCESSED
                            and meta.get("kind") == kind):
                        return _read_processed(cache_dir, meta)
                except Exception:
                    pass  # corrupt/mismatched cache -> fresh call below

            result = fn(*args, **kwargs)
            try:
                section_d, section_meta = serialize_processed(result)
                cache_dir.mkdir(parents=True, exist_ok=True)
                for fname, cols in section_d.items():
                    write_frame(cols, cache_dir / fname)
                meta = {
                    "schema_ver": SCHEMA_VER_PROCESSED,
                    "reader": reader,
                    "ver": ver,
                    "kind": kind,
                    "parse_kwargs": parse_kwargs,
                    "sections": section_meta,
                }
                meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2))
                files = [fname for smeta in section_meta.values() for fname in
                         [fr["file"] for fr in smeta["files"]]]
                _append_index(ver, {
                    "raw_ptr": os.fspath(raw_path),
                    "kind": kind,
                    "reader": reader,
                    "ver": ver,
                    "key": key,
                    "schema_ver": SCHEMA_VER_PROCESSED,
                    "parse_kwargs": parse_kwargs,
                    "files": [os.path.relpath(cache_dir / f, get_project_root())
                              for f in files],
                })
            except Exception:
                pass  # never let the cache put the pipeline at risk
            return result

        return wrapper

    return decorate
