# -*- coding:utf-8 -*-
"""Layered parquet cache for the reader: faithful L1 frames + processed L2.

Layout::

    data/{ver}/l1/{key}/<kind>.parquet + <kind>.meta.json   # faithful frames
    data/{ver}/l2/{key}/<section>.parquet + meta.json       # processed (sci, tel)

``key`` is ``sha256(f"{ver}|{reader}|{raw_path}|{parse_kwargs}")[:16]``.  The L1
layer stores exactly what was decoded from the raw stream (one row per particle
for science, one row per sample for HK/TL); L2 stores the processed ``(sci,
tel)`` output and must reproduce the legacy ``single_readXX`` values.

The serializer is a generic traversal: a flat ``Dict[str, np.ndarray]`` is stored
column-per-column and 2-D arrays are stored as a polars fixed-size ``Array``
column.  Anything irregular causes a bail so the caller falls back to a fresh
decode and the pipeline never fails because of the cache.
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

SCHEMA_VER = 3
SCHEMA_VER_PROCESSED = 4


def get_project_root(marker=("justfile", ".gitignore")):
    current = Path(os.getcwd())
    while True:
        if any((current / m).exists() for m in marker):
            return current
        if current == current.parent:
            return current
        current = current.parent


def _layer_root(ver: str, layer: str):
    return get_project_root() / "data" / ver / layer


def _cache_root(ver: str):
    """L1 frame-cache root (``data/{ver}/l1``)."""
    return _layer_root(ver, "l1")


def _l2_root(ver: str):
    """L2 processed-cache root (``data/{ver}/l2``)."""
    return _layer_root(ver, "l2")


def _make_key(ver: str, reader: str, raw_path: str, parse_kwargs: dict):
    payload = f"{ver}|{reader}|{raw_path}|{json.dumps(parse_kwargs, sort_keys=True)}"
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

def _append_index(ver: str, layer: str, record: dict):
    root = _layer_root(ver, layer)
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


# --- L1 frame cache -----------------------------------------------------------

def _kind_meta_path(cache_dir: Path, kind: str):
    return cache_dir / f"{kind}.meta.json"


def _load_kind_meta(cache_dir: Path, kind: str):
    meta_path = _kind_meta_path(cache_dir, kind)
    parquet_path = cache_dir / f"{kind}.parquet"
    if not (meta_path.exists() and parquet_path.exists()):
        return None
    try:
        meta = json.loads(meta_path.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    if meta.get("schema_ver") != SCHEMA_VER or meta.get("kind") != kind:
        return None
    return meta


def get_l1_frames(ver: str, reader: str, raw_path, decoders: dict,
                  parse_kwargs: dict = None, overwrite: bool = False) -> dict:
    """Return ``{kind: frames}`` for one raw file, using the L1 parquet cache.

    ``decoders`` maps each kind (``sci``/``hk``/``tl``) to a zero-argument
    callable that performs the fresh decode.  On a full hit nothing is decoded;
    a corrupt/missing entry falls back to a fresh decode and rewrites the cache.
    """
    parse_kwargs = dict(parse_kwargs or {})
    key = _make_key(ver, reader, raw_path, parse_kwargs)
    cache_dir = _cache_root(ver) / key

    if not overwrite:
        loaded = {}
        for kind in decoders:
            meta = _load_kind_meta(cache_dir, kind)
            if meta is None:
                loaded = None
                break
            try:
                loaded[kind] = deserialize(read_frame(cache_dir / f"{kind}.parquet", meta))
            except Exception:
                loaded = None
                break
        if loaded is not None:
            return loaded

    result = {kind: decode() for kind, decode in decoders.items()}
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        for kind, frames in result.items():
            cols, dtypes, rows, shapes = serialize(frames)
            write_frame(cols, cache_dir / f"{kind}.parquet")
            meta = {
                "schema_ver": SCHEMA_VER,
                "ver": ver,
                "reader": reader,
                "kind": kind,
                "dtypes": dtypes,
                "shapes": shapes,
                "rows": rows,
                "ncols": len(cols),
                "parse_kwargs": parse_kwargs,
            }
            _kind_meta_path(cache_dir, kind).write_text(
                json.dumps(meta, ensure_ascii=False, indent=2))
            rel_parquet = os.path.relpath(cache_dir / f"{kind}.parquet", get_project_root())
            _append_index(ver, "l1", {
                "raw_ptr": os.fspath(raw_path),
                "kind": kind,
                "reader": reader,
                "ver": ver,
                "key": key,
                "schema_ver": SCHEMA_VER,
                "parse_kwargs": parse_kwargs,
                "file": rel_parquet,
                "events": rel_parquet if kind == "sci" else None,
            })
    except Exception:
        pass  # never let the cache put the pipeline at risk
    return result


def with_l1_cache(ver: str, reader: str, kind: str):
    """Wrap a raw-parse function with the faithful L1 frame cache.

    ``kind`` is ``sci``/``hk``/``tl`` and selects the parquet file name. The
    wrapped function must return a *fresh* object each call (the downstream
    post-processing mutates the Dict), so a cache hit reconstructs a new Dict.
    """
    def decorate(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            overwrite = kwargs.pop("overwrite_cache", False)
            raw_path = args[0]
            parse_kwargs = _extract_parse_kwargs(fn, args, kwargs)
            frames = get_l1_frames(
                ver, reader, raw_path, {kind: lambda: fn(*args, **kwargs)},
                parse_kwargs, overwrite=overwrite,
            )
            return frames[kind]

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


# --- L2 processed cache -------------------------------------------------------

def get_l2_processed(ver: str, reader: str, raw_path, parse_kwargs: dict,
                     process, overwrite: bool = False):
    """Return the processed ``(sci, tel)`` for one raw file, cached under L2.

    ``process`` is a zero-argument callable performing the fresh processing.  A
    cache hit returns rebuilt plain dicts (matching the reader's return type) so
    downstream access is unchanged.
    """
    parse_kwargs = dict(parse_kwargs or {})
    key = _make_key(ver, reader, raw_path, parse_kwargs)
    cache_dir = _l2_root(ver) / key
    meta_path = cache_dir / "meta.json"

    if not overwrite and meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text())
            if meta.get("schema_ver") == SCHEMA_VER_PROCESSED:
                return _read_processed(cache_dir, meta)
        except Exception:
            pass  # corrupt/mismatched cache -> fresh call below

    result = process()
    try:
        section_d, section_meta = serialize_processed(result)
        cache_dir.mkdir(parents=True, exist_ok=True)
        for fname, cols in section_d.items():
            write_frame(cols, cache_dir / fname)
        meta = {
            "schema_ver": SCHEMA_VER_PROCESSED,
            "reader": reader,
            "ver": ver,
            "parse_kwargs": parse_kwargs,
            "sections": section_meta,
        }
        meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2))
        files = [fname for smeta in section_meta.values() for fname in
                 [fr["file"] for fr in smeta["files"]]]
        _append_index(ver, "l2", {
            "raw_ptr": os.fspath(raw_path),
            "kind": parse_kwargs.get("_kind"),
            "reader": reader,
            "ver": ver,
            "key": key,
            "schema_ver": SCHEMA_VER_PROCESSED,
            "parse_kwargs": parse_kwargs,
            "files": [os.path.relpath(cache_dir / f, get_project_root()) for f in files],
        })
    except Exception:
        pass  # never let the cache put the pipeline at risk
    return result


def with_l1_cache_processed(ver: str, reader: str, kind: str):
    """Wrap a ``(sci, tel)``-returning reader with the L2 processed cache."""
    def decorate(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            overwrite = kwargs.pop("overwrite_cache", False)
            raw_path = args[0]
            parse_kwargs = _extract_parse_kwargs(fn, args, kwargs)
            parse_kwargs["_kind"] = kind
            return get_l2_processed(
                ver, reader, raw_path, parse_kwargs,
                lambda: fn(*args, **kwargs), overwrite=overwrite,
            )

        return wrapper

    return decorate
