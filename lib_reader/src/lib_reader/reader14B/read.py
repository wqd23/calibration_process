# -*- coding:utf-8 -*-
"""Reader for the 14B payload (``cali_data/14B15B/03``).

Packet layout, verified against the raw streams (header -> next header
distance, tail check, CRC):

- science: 584-byte ``grid1x_ft_packet`` with head ``1c 1c 22 88`` and tail
  ``cc 11 88 22`` -- 38 events per packet with a 14-byte stride
  (``timestamp2`` 4B + ``data_max`` 2B + ``data_base`` 2B + ``data_sum`` 4B +
  ``data_ccm`` 2B).  Defined in ``yunyao_packet.xml`` shipped in this package.
- HK: 178-byte ``hk_packet`` with head ``1a 2b 3c 4d`` and no tail (the next
  head follows the 2-byte CRC directly).  The 224-byte ``hk_packet`` declared
  in ``yunyao_packet.xml`` does NOT match this data; the stream matches the
  N1 definition in ``readerN1/yingtian_packet.xml`` byte for byte (CRC over
  ``[0, 176)`` stored at 176), so that definition is reused here.

TB (temperature-bias) runs are ``{temp}-{bias}-{idx}.event.dat`` paired with
``{temp}-{bias}-ecu_{n}.hk``; the HK file starts with the bias discharge/ramp
and sometimes long idle epochs, so records are cut to the stable plateau at
the target bias (``hk_bias``, defaulting to the ``{bias}`` code in the file
name) and to the science acquisition utc window (the ground clocks here do
carry a real utc, unlike 12B).  ``amp = data_max - data_base / 4``.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from addict import Dict

from ..packet_parser import parse_grid_data_new
from ..l1_cache import get_l1_frames, get_l2_processed, apply_selection

_XML = str(Path(__file__).with_name("yunyao_packet.xml"))
# 14B HK frames are the 178-byte hk_packet defined for N1/yingtian; the
# 224-byte hk_packet in yunyao_packet.xml does not match this data stream.
_HK_XML = str(Path(__file__).resolve().parent.parent / "readerN1" / "yingtian_packet.xml")

# TB file name: {temp}-{bias}-{idx}.event.dat, bias code in units of 0.1 V
_TB_NAME = re.compile(r"^(m?\d+C)[-_](\d{3})[-_]")

_FT_EVENTS = 38
_FT_STRIDE = 14


def _cache_ver(path) -> str:
    """Version component of ``.../data/<ver>/<root>/...`` (falls back to '14B').

    The same reader serves sibling payloads that share the 14B format (15B), so
    each keeps its own L1/L2 cache under ``data/<ver>/``.
    """
    parts = Path(path).parts
    if "data" in parts:
        i = parts.index("data")
        roots = {"raw_data", "ec_src", "ec_xray"}
        segs = parts[i + 1:]
        for j, seg in enumerate(segs):
            if seg in roots:
                return "/".join(segs[:j]) or "14B"
    return "14B"


def _bias_from_name(path) -> float | None:
    m = _TB_NAME.match(Path(path).stem)
    return int(m.group(2)) / 10.0 if m else None


def _readSci_impl(path, mode="ft"):
    if mode != "ft":
        raise ValueError(f"{path}: 14B science files are ft packets (got mode={mode!r})")
    return Dict(parse_grid_data_new(
        path, xml_file=_XML, data_tag="grid1x_ft_packet", endian="MSB",
        multi_evt=_FT_EVENTS, multi_step=_FT_STRIDE)[0])


def _readHK_impl(path):
    # 14B HK files carry frames shorter than the yunyao declaration; framing
    # and CRC are the N1/yingtian 178-byte definition.
    return Dict(parse_grid_data_new(
        path, xml_file=_HK_XML, data_tag="hk_packet", endian="MSB")[0])


def readSci(path, mode="ft", overwrite_cache=False):
    return get_l1_frames(
        _cache_ver(path), "14b", path,
        {"sci": lambda: _readSci_impl(path, mode=mode)}, {"mode": mode},
        overwrite=overwrite_cache,
    )["sci"]


def readHK(path, overwrite_cache=False):
    return get_l1_frames(
        _cache_ver(path), "14b", path,
        {"hk": lambda: _readHK_impl(path)}, {},
        overwrite=overwrite_cache,
    )["hk"]


def getHK(sciFile):
    """Pair a science file with its HK file.

    ``{temp}-{bias}-{idx}.event.dat`` -> ``{temp}-{bias}-ecu_*.hk``,
    ``{name}_{idx}.event.dat`` -> ``{name}_ecu_*.hk`` (source runs) and
    ``{kV}-ch{n}-{idx}.event.dat`` -> ``{kV}-ch{n}-ecu_*.hk`` (X-ray).  When
    several candidates match, pick the one whose utc range overlaps the
    science utc range the most.
    """
    sciFile = Path(sciFile)
    stem = sciFile.stem
    cands = []
    for sep in ("-", "_"):
        if sep not in stem:
            continue
        prefix = stem.rsplit(sep, 1)[0]
        for pat in (f"{prefix}-ecu_*.hk", f"{prefix}_ecu_*.hk"):
            cands += [f for f in sciFile.parent.glob(pat) if "error" not in f.stem]
    cands = list(dict.fromkeys(cands))
    if not cands:
        raise FileNotFoundError(f"HK file for {sciFile} does not exist.")
    if len(cands) == 1:
        return cands[0]
    sci = readSci(str(sciFile))
    lo, hi = float(np.min(sci.utc)), float(np.max(sci.utc))
    best, best_overlap = None, -1.0
    for f in cands:
        tel = readHK(str(f))
        u = np.asarray(tel.utc_time, dtype=float)
        overlap = max(0.0, min(hi, u.max()) - max(lo, u.min()))
        if overlap > best_overlap:
            best, best_overlap = f, overlap
    return best


def _stable_mask(bias, temp, target, hk_utc, sci_utc):
    """Records on the target-bias plateau during the science acquisition."""
    bm = bias.mean(axis=0)
    if target is not None:
        stable = np.abs(bm - target) < 0.3
    else:
        # fallback (12B style): dominant 0.1 V x 1 C epoch
        bb, tb = np.round(bm, 1), np.round(temp.mean(axis=0))
        keys, counts = np.unique(np.stack([bb, tb], axis=1), axis=0, return_counts=True)
        b0, t0 = keys[counts.argmax()]
        stable = (np.abs(bb - b0) < 0.05) & (np.abs(tb - t0) < 0.6)
    stable &= np.logical_and.reduce(np.abs(bias - bm) < 0.3, axis=0)
    if sci_utc is not None and sci_utc.size and sci_utc.max() > 0:
        in_window = (hk_utc >= sci_utc.min()) & (hk_utc <= sci_utc.max())
        if np.any(stable & in_window):
            stable &= in_window
    return stable


def _single_read14B_impl(path, mode, hk_name, hk_bias, sci_half, overwrite,
                         select=None, ver="14B"):
    sci = readSci(path, mode=mode, overwrite_cache=overwrite)
    if select is not None:
        sci = apply_selection(ver, "14b", path, {"mode": mode},
                              select[0], select[1], sci)
    tel = readHK(str(hk_name), overwrite_cache=overwrite)

    if sci_half is not None:
        n_all = sci.data_max.shape[0]
        half = n_all // 2
        sl = slice(0, half) if sci_half == "first" else slice(half, n_all)
        for k, v in sci.items():
            if np.asarray(v).shape[0] == n_all:
                sci[k] = np.asarray(v)[sl]

    if len(sci.data_max) != len(sci.data_base):
        raise ValueError(f"{path}: data_max.len != data_base.len")
    sci.amp = sci.data_max - sci.data_base / 4.0
    sci["timestampEvt"] = sci.timestamp2

    tel.tempSipm = [np.asarray(tel[f"sipm_temp{i}"], dtype=float) / 100 - 273.15
                    for i in range(4)]
    tel.iMon = [np.asarray(tel[f"sipm_current{i}"], dtype=float) for i in range(4)]
    # The HK monitor reads the regulated SiPM-side voltage (it equals the
    # file-name setpoint), so no series-resistor drop is subtracted; pending
    # hardware confirmation, see docs/intermediate_data.md section 6.1.
    tel.vMon = [np.asarray(tel[f"sipm_voltage{i}"], dtype=float) / 1000 for i in range(4)]
    tel.bias = [tel.vMon[i] for i in range(4)]

    b4 = np.stack([np.asarray(b) for b in tel.bias])
    t4 = np.stack([np.asarray(t) for t in tel.tempSipm])
    target = hk_bias if hk_bias is not None else _bias_from_name(path)
    stable = _stable_mask(b4, t4, target, np.asarray(tel.utc_time, dtype=float),
                          np.asarray(sci.utc, dtype=float))
    if not np.any(stable):
        raise ValueError(f"{path}: no stable-bias HK records found")

    for k in list(tel.keys()):
        v = tel[k]
        if isinstance(v, list) and len(v) == 4:
            tel[k] = [np.asarray(v[i])[stable] for i in range(4)]
        elif isinstance(v, np.ndarray) and v.shape[0] == stable.shape[0]:
            tel[k] = v[stable]

    n = sci.data_max.shape[0]
    for k, v in list(sci.items()):
        v = np.asarray(v)
        if v.ndim == 1 and v.shape[0] == n and k != "channel_n":
            sci[k] = [v[sci.channel_n == i] for i in range(4)]
    return sci, tel


def single_read14B(path: str, config=None, mode="ft", hk_path=None, hk_bias=None,
                   sci_half=None, cache_ver=None, **kwargs):
    overwrite = kwargs.get("overwrite_cache", False)
    select = kwargs.get("select")
    ver = cache_ver or "14B"
    params = {
        "mode": mode,
        "hk_path": None if hk_path is None else str(hk_path),
        "hk_bias": hk_bias,
        "sci_half": sci_half,
    }
    if select is not None:
        params["select"] = select[0]

    def process():
        hk_name = Path(hk_path) if hk_path else getHK(path)
        return _single_read14B_impl(path, mode, hk_name, hk_bias, sci_half,
                                    overwrite, select=select, ver=ver)

    return get_l2_processed(ver, "14b", path, params, process, overwrite=overwrite)
