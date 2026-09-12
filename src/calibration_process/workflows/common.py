# -*- coding:utf-8 -*-
"""Reusable scientific stages shared by the explicit workflows.

All of these call the protected kernel directly (file_lib, util_lib,
lib_plot).  They never re-implement the fit / plot / serialization logic.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import List, Optional

import numpy as np

from .. import file_lib, util_lib as util
from ..config_schema import ManifestEntry
from ..products import FileRunSpec, SingleFitResult, TBPoint, ECPoint
from ..runtime import RuntimeConfig


def _jsonable(obj):
    """Recursively convert numpy scalars/arrays so ``json.dumps`` works."""
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if callable(obj):
        return getattr(obj, "__name__", repr(obj))
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)


def _save_spectrum(fp, path) -> None:
    """Write the 4-channel spectrum as a portable long-format parquet table."""
    import polars as pl

    channel, bin_idx, x, spectrum, spectrum_err = [], [], [], [], []
    for i in range(4):
        xi = np.asarray(fp.x[i])
        channel.append(np.full(len(xi), i, dtype=np.int8))
        bin_idx.append(np.arange(len(xi), dtype=np.int32))
        x.append(xi)
        spectrum.append(np.asarray(fp.spectrum[i]))
        spectrum_err.append(np.asarray(fp.spectrum_err[i]))
    pl.DataFrame({
        "channel": np.concatenate(channel),
        "bin": np.concatenate(bin_idx),
        "x": np.concatenate(x),
        "spectrum": np.concatenate(spectrum),
        "spectrum_err": np.concatenate(spectrum_err),
    }).write_parquet(path, compression="zstd")


# --------------------------------------------------------------------------- #
# FileRunSpec resolution
# --------------------------------------------------------------------------- #
def _selection(version: str, branch: str):
    """Version-defined event selection hook for a branch, or None.

    A version module may expose ``selection(version, branch) -> (selkey, fn)``;
    the hook is written entirely in Python (never YAML) and is applied by the
    reader at L1->L2.
    """
    from .registry import get_workflow

    try:
        wf = get_workflow(version)
    except KeyError:
        return None
    hook = getattr(wf, "selection", None)
    if hook is None:
        return None
    return hook(version, branch)


def single_run_spec(rt: RuntimeConfig, branch: str, m: ManifestEntry) -> FileRunSpec:
    """Resolve the read/bkg/spectrum/fit config for one measurement.

    Mirrors the legacy ``file_config()`` tuple so the kernel sees identical
    settings.  This is where version orchestration differences surface
    (background rotation, 4-channel reconstruction, corr).
    """
    data_dir = rt.data_dir
    sel = _selection(getattr(rt, "version", ""), branch)

    def abspath(rel: str) -> str:
        return str(data_dir / rel)

    def rc(path: str, ending: str, select=sel, **kwargs) -> "file_lib.Read_config":
        params = getattr(rt, "reader_params", lambda _e: {})(ending)
        return file_lib.Read_config(
            path, ending=ending, reader_params=params, select=select,
            version=getattr(rt, "version", ""), **kwargs,
        )

    if branch == "tb":
        pb = rt.payload.tb
        kwarg: dict = {}
        if m.hk_files:
            kwarg["hk_path"] = abspath(m.hk_files[-1])
        if m.metadata.get("sci_half"):
            kwarg["sci_half"] = m.metadata["sci_half"]
        if m.metadata.get("hk_bias") is not None:
            kwarg["hk_bias"] = m.metadata["hk_bias"]
        for key in ("mode", "seg_bias", "quantity"):
            if m.metadata.get(key) is not None:
                kwarg[key] = m.metadata[key]
        read = rc(abspath(m.science_files[-1]), pb.reader, kwarg=kwarg)
        bkg = rc("", "normal", select=None)
        spec = file_lib.Spectrum_config(bin_width=pb.bin_width, adc_max=pb.adc_max)
        fit = file_lib.Fit_config(rt.fit_range(branch, _fit_key(m)), rt.bkg_form(branch, _fit_key(m)))
        return FileRunSpec(read, bkg, spec, fit)

    if branch == "ec_source":
        pb = rt.payload.ec
        src_reader = pb.src_reader or pb.reader
        read = rc(abspath(m.science_files[-1]), src_reader)
        bkg_rel = m.aux_files[-1] if m.aux_files else ""
        bkg = (
            rc(abspath(bkg_rel), src_reader, select=None)
            if bkg_rel
            else rc("", src_reader, select=None)
        )
        spec = file_lib.Spectrum_config(
            corr=rt.corr, bin_width=pb.bin_width, adc_max=pb.adc_max,
            rate_span=getattr(pb, "rate_span", "union"),
        )
        fit = file_lib.Fit_config(rt.fit_range(branch, _fit_key(m)), rt.bkg_form(branch, _fit_key(m)))
        return FileRunSpec(read, bkg, spec, fit)

    if branch == "ec_xray":
        pb = rt.payload.ec
        if pb.xray_single_file:
            # 05B style: one 4-channel file, background = the same file with a
            # rotated (cyclically shifted per channel) time cut
            basename = os.path.basename(m.science_files[-1])
            reader = pb.xray_reader or pb.reader
            read = rc(
                abspath(m.science_files[-1]), reader,
                config_file=pb.xray_config_file or "",
                time_cut=_time_cut(pb, basename),
            )
            bkg = rc(
                abspath(m.science_files[-1]), reader,
                config_file=pb.xray_config_file or "",
                time_cut=_bkg_time_cut(pb, basename),
                select=None,
            )
            spec = file_lib.Spectrum_config(
                corr=rt.corr, bin_width=pb.bin_width, adc_max=pb.adc_max
            )
            fit = file_lib.Fit_config(rt.fit_range(branch, _fit_key(m)), rt.bkg_form(branch, _fit_key(m)))
            return FileRunSpec(read, bkg, spec, fit)

        reads = [
            rc(abspath(f), pb.reader, config_file=pb.xray_config_file or "")
            for f in m.science_files
        ]
        n = pb.channel_count
        rotation = pb.xray_bkg_rotation
        bkg_reads = _rotate_bkg(reads, rotation, n)
        spec = file_lib.Spectrum_config(
            corr=rt.corr, bin_width=pb.bin_width, adc_max=pb.adc_max,
            rate_span=getattr(pb, "rate_span", "union"),
        )
        fit = file_lib.Fit_config(rt.fit_range(branch, _fit_key(m)), rt.bkg_form(branch, _fit_key(m)))
        return FileRunSpec(reads, bkg_reads, spec, fit)

    if branch == "neutron":
        # standalone peak-fit profile: no TB/EC global, identity correction
        pb = rt.payload.neutron
        if pb is None:
            raise ValueError("neutron branch configured but payload.neutron is missing")
        kwarg: dict = {}
        if m.hk_files:
            kwarg["hk_path"] = abspath(m.hk_files[-1])
        for key in ("sci_half", "hk_bias", "mode"):
            if m.metadata.get(key) is not None:
                kwarg[key] = m.metadata[key]
        read = rc(abspath(m.science_files[-1]), pb.reader, kwarg=kwarg)
        bkg = rc("", pb.reader, select=None)
        spec = file_lib.Spectrum_config(bin_width=pb.bin_width, adc_max=pb.adc_max)
        fit = file_lib.Fit_config(rt.fit_range(branch, _fit_key(m)), None)
        return FileRunSpec(read, bkg, spec, fit)

    raise ValueError(f"unknown branch {branch!r}")


def _time_cut(pb, basename: str):
    tc = getattr(pb, "time_cut", None) or {}
    return tc.get(basename)


def _bkg_time_cut(pb, basename: str):
    tc = getattr(pb, "time_cut", None) or {}
    v = tc.get(basename)
    if v is None:
        return None
    # legacy bkg_time_cut rotation: [v[1], v[2], v[3], v[0]]
    return [v[1], v[2], v[3], v[0]]


def _fit_key(m) -> str:
    """Return the fit_range/bkg_form lookup key for a measurement.

    Most versions key fit ranges by the measurement id itself; 11B keys its TB
    fit ranges by the file stem (so ``just list`` shows the basename but the
    range is looked up by stem).  A version workflow may record ``fit_key`` in
    the measurement metadata.
    """
    return m.metadata.get("fit_key", m.id) if m.metadata else m.id


def channel_use(m, ch: int) -> bool:
    """measurement-level then channel-level exclusion (plan section 25)."""
    if not m.use:
        return False
    if m.channels and str(ch) in m.channels:
        return bool(m.channels[str(ch)].get("use", True))
    return True


def _rotate_bkg(reads: List, rotation: str, n: int) -> List:
    if rotation == "circle":
        # ch i background uses channel (i+1) % n  (legacy 03B/07/04/09)
        return [reads[(i + 1) % n] for i in range(n)]
    if rotation == "fixed":
        # legacy 10B/11B/12B: [ch1, ch2, ch0, ch0]
        return [reads[1], reads[2], reads[0], reads[0]]
    raise ValueError(f"unknown xray_bkg_rotation {rotation!r}")


# --------------------------------------------------------------------------- #
# File_operation_05b construction (thin wrappers over the protected kernel)
# --------------------------------------------------------------------------- #
def _dict_4ch_reconstruct(dict_4ch):
    """Rebuild a per-channel dict from a list of 4 single-channel dicts.

    Mirrors legacy ``__dict_4ch_reconstruct``: for each key the value is the
    list of per-channel values, with the i-th channel taken from the i-th
    single-channel dict (falling back to the whole value when it is not
    length-4).
    """
    out = {}
    for key in dict_4ch[0].keys():
        out[key] = []
        for i in range(4):
            if len(dict_4ch[i][key]) == 4:
                out[key].append(dict_4ch[i][key][i])
            else:
                out[key].append(dict_4ch[i][key])
    return out


def _build_fp05b(config, nocache=False) -> "file_lib.File_operation_05b":
    return file_lib.File_operation_05b(config[0].path, *config, nocache=nocache)


def _build_fp03b(config, nocache=False) -> "file_lib.File_operation_05b":
    """Build a 4-channel File_operation_05b by reading each channel file.

    The four per-channel operations are read independently and then merged back
    into a single 4-channel object (exactly the legacy ``__get_fp03B``
    behaviour).
    """
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fps = [
        file_lib.File_operation_05b(
            read_config[i].path, read_config[i], bkg_read_config[i],
            spectrum_config, fit_config, nocache=nocache,
        )
        for i in range(4)
    ]
    sci = _dict_4ch_reconstruct([fps[i].sci for i in range(4)])
    tel = _dict_4ch_reconstruct([fps[i].tel for i in range(4)])
    bkg_sci = _dict_4ch_reconstruct([fps[i].bkg_sci for i in range(4)])
    bkg_tel = _dict_4ch_reconstruct([fps[i].bkg_tel for i in range(4)])
    fp = fps[0]
    fp.sci, fp.tel = sci, tel
    fp.bkg_sci, fp.bkg_tel = bkg_sci, bkg_tel
    return fp


# --------------------------------------------------------------------------- #
# Single-fit stage
# --------------------------------------------------------------------------- #
def build_fit_operation(rt, branch, fc: FileRunSpec, nocache=False) -> object:
    """Construct a File_operation_05b from a resolved spec (protected kernel)."""
    if branch == "ec_xray" and not rt.payload.ec.xray_single_file:
        return _build_fp03b(
            [fc.read_config, fc.bkg_read_config, fc.spectrum_config, fc.fit_config],
            nocache=nocache,
        )
    return _build_fp05b(
        [fc.read_config, fc.bkg_read_config, fc.spectrum_config, fc.fit_config],
        nocache=nocache,
    )


def qa_category(output_dir: str, branch: str) -> str:
    if branch == "tb":
        return "tb"
    if branch == "neutron":
        return "neutron"
    return "ec"


def _fit_subdir(branch: str) -> str:
    return {
        "tb": "TB_fit_result",
        "ec_source": "EC_fit_result",
        "ec_xray": "EC_fit_result",
        "neutron": "NEUTRON_fit_result",
    }.get(branch, "EC_fit_result")


def run_single_fit(
    rt: RuntimeConfig,
    branch: str,
    m: ManifestEntry,
    output_root: Optional[Path] = None,
    nocache: bool = False,
):
    """Run the single-fit stage and save pickle + single-fit figure.

    Mirrors legacy ``process()`` exactly, but file/config selection comes from
    the manifest instead of a directory scan.
    """
    output_root = output_root or rt.data_dir
    category = qa_category(str(output_root), branch)
    fc = single_run_spec(rt, branch, m)
    fp = build_fit_operation(rt, branch, fc, nocache=nocache)
    fp.qa_thresholds = util.load_qa_thresholds(rt.version, category)
    fp.get_spectrum()
    fp.peak_fit()

    from lib_plot import plot

    # the persisted/figure name is the measurement id with its extension
    # stripped: TB/EC-src ids are "<stem>.txt", EC-xray ids are "<energy>"
    stem = os.path.splitext(m.id)[0]
    title = (
        f"{stem}: {np.mean(fp.tel['bias'][0]):.2f}V, "
        f"{np.mean(fp.tel['tempSipm'][0]):.2f}C"
    )
    sub = _fit_subdir(branch)
    fig_dir = output_root / "single_process" / "single_fit_fig"
    save_dir = output_root / "single_process" / sub
    fig_dir.mkdir(parents=True, exist_ok=True)
    save_dir.mkdir(parents=True, exist_ok=True)
    plot.fit_plot(
        fp.spectrum,
        fp.x,
        fp.fit_result,
        title=title,
        bkgForm=fc.fit_config.bkg_form,
        fit_range=fc.fit_config.fit_range,
        save_path=str(fig_dir / f"{stem}.png"),
    )
    # L3: portable fit parameters + spectrum (the pickle stays for compatibility
    # with legacy consumers and for the frozen-oracle regression).
    fit_payload = {"file": str(fp.path), "fit_result": _jsonable(fp.fit_result)}
    if branch == "neutron":
        # the standalone profile does no TB/EC correction, so record the
        # measured temperature alongside the fit for later reference
        fit_payload["temperature"] = [
            float(np.mean(t)) for t in fp.tel.get("tempSipm", [])
        ]
    util.json_save(fit_payload, str(save_dir / f"{stem}.fit.json"))
    _save_spectrum(fp, save_dir / f"{stem}.spectrum.parquet")
    fp.save(str(save_dir / f"{stem}.pickle"))
    return fp


class _LoadedFit:
    """Minimal accessor over a persisted single-fit result."""

    def __init__(self, fit_result, tel):
        self.fit_result = fit_result
        self.tel = tel


def load_single_fp_from_store(rt, branch, m, output_root) -> _LoadedFit:
    """Load a single-fit result: L3 ``fit.json`` + L2 telemetry.

    ``fit_result`` comes from the portable L3 ``fit.json``.  TB telemetry is
    rebuilt from the reader's L2 ``processed.parquet`` cache by re-resolving
    the exact read spec (the reader returns the cached processing); this keeps
    L4 independent of the dill pickle's ``tel``.  Both fall back to the pickle
    when the newer artefacts are absent.
    """
    stem = os.path.splitext(m.id)[0]
    sub = _fit_subdir(branch)
    base = output_root / "single_process" / sub
    # build names by string concat: ids such as "x.event.dat" have an interior
    # dot, so Path.with_suffix would drop ".event"
    fit_json = base / f"{stem}.fit.json"
    pickle_path = base / f"{stem}.pickle"
    if fit_json.exists():
        fit_result = json.loads(fit_json.read_text())["fit_result"]
    else:
        fit_result = util.pickle_load(str(pickle_path))["fit_result"]

    tel = None
    if branch == "tb":
        tel = build_fit_operation(rt, "tb", single_run_spec(rt, "tb", m)).tel
        if tel is None:
            tel = util.pickle_load(str(pickle_path))["tel"]
    return _LoadedFit(fit_result, tel)


# --------------------------------------------------------------------------- #
# Typed intermediate conversion
# --------------------------------------------------------------------------- #
def to_single_fit_result(m: ManifestEntry, fp) -> List[Optional[SingleFitResult]]:
    """Wrap the kernel's per-channel fit_result into typed SingleFitResult."""
    out: List[Optional[SingleFitResult]] = []
    fit_result = getattr(fp, "fit_result", None) or []
    for ch, fr in enumerate(fit_result):
        if fr is None:
            out.append(None)
            continue
        out.append(
            SingleFitResult(
                measurement_id=m.id,
                channel=ch,
                peak_amplitude=fr["a"],
                peak_amplitude_err=fr["a_err"],
                peak_center=fr["b"],
                peak_center_err=fr["b_err"],
                peak_sigma=fr["c"],
                peak_sigma_err=fr["c_err"],
                resolution=fr["resolution"],
                resolution_err=fr["resolution_err"],
                rate=fr["rate"],
                rate_err=fr["rate_err"],
                redchi=fr["redchi"],
                ndf=fr["ndf"],
                success=fr["success"],
                qa_flag=fr["qa_flag"],
                boundary_hit=list(fr.get("boundary_hit", [])),
            )
        )
    return out


# --------------------------------------------------------------------------- #
# TB points
# --------------------------------------------------------------------------- #
def build_tb_points(rt: RuntimeConfig, items: List) -> List[List[TBPoint]]:
    """Build per-channel TBPoint lists from (measurement, fp) pairs.

    Mirrors legacy ``load_data``: temperature/bias are the mean/std of the
    telemetry arrays for that channel; the center is the fitted peak center.
    """
    per_channel: List[List[TBPoint]] = [[] for _ in range(4)]
    tb_cfg = getattr(getattr(rt, "payload", None), "tb", None)
    skip_fail = bool(getattr(tb_cfg, "skip_qa_fail", False))
    for m, fp in items:
        tel_4ch = [
            {k: v[i] for k, v in fp.tel.items() if len(v) == 4} for i in range(4)
        ]
        for ch, fit in enumerate(fp.fit_result):
            if fit is None:
                continue
            if skip_fail and fit.get("qa_flag") == "fail":
                continue
            tel = tel_4ch[ch]
            temp = float(np.average(tel["tempSipm"]))
            temp_err = float(np.std(tel["tempSipm"]))
            bias = float(np.average(tel["bias"]))
            bias_err = float(np.std(tel["bias"]))
            per_channel[ch].append(
                TBPoint(
                    measurement_id=m.id,
                    channel=ch,
                    temperature=temp,
                    temperature_err=temp_err,
                    bias=bias,
                    bias_err=bias_err,
                    peak_center=fit["b"],
                    peak_center_err=fit["b_err"],
                    enabled=channel_use(m, ch),
                )
            )
    return per_channel


def global_tb(rt: RuntimeConfig, per_channel: List[List[TBPoint]],
              result_path: Path) -> List[dict]:
    """Temperature-bias 2D global fit (protected kernel)."""
    from lib_plot import plot

    result_path.mkdir(parents=True, exist_ok=True)
    pb = rt.payload.tb
    channels = getattr(pb, "channels", None) or [0, 1, 2, 3]
    # Keep the historical four-slot result so downstream consumers can index by
    # channel; a channel that is not configured, disabled or has no points stays
    # ``null`` instead of aborting the whole fit.
    result: List[Optional[dict]] = [None] * 4
    for ch in channels:
        points = per_channel[ch] if ch < len(per_channel) else []
        pts = [p for p in points if p.enabled]
        if not pts:
            continue
        data_all = np.array(
            [[p.peak_center, p.peak_center_err, p.temperature, p.temperature_err,
              p.bias, p.bias_err] for p in pts]
        )
        data_all = _apply_bias_filter(data_all, pb.bias_min_filter)
        if data_all.size == 0:
            continue
        center, center_err, temp, temp_err, bias, bias_err = (
            data_all[:, 0], data_all[:, 1], data_all[:, 2],
            data_all[:, 3], data_all[:, 4], data_all[:, 5],
        )
        try:
            if pb.tb_fit_method == "lmfit":
                res = util.temp_bias_lmfit(
                    center, center_err, temp, temp_err, bias, bias_err
                )
            else:
                p0_by_ch = getattr(pb, "tb_fit_p0_by_channel", None) or {}
                p0 = p0_by_ch.get(str(ch), pb.tb_fit_p0)
                res = util.temp_bias_fit_curvefit(
                    center, center_err, temp, bias,
                    p0=p0, maxfev=pb.tb_fit_maxfev,
                )
        except util.FitError as e:
            raise util.FitError(f"failed to do temp bias fit: {e.args[-1]}")
        result[ch] = res
        xy = np.stack([temp, bias], axis=1)
        name = f"temp_bias_fit_{ch}.png"
        plot.fit_err_plot_2d(
            xy,
            center,
            lambda x: util.tempbias2DFunctionInternal(x, *(list(res.values())[:5])),
            ("temp$^\\circ$C", "bias/V", "center"),
            title=f"temp bias fit: channel {ch}",
            save_path=str(result_path / util.headtime(name)),
        )
    util.json_save(
        result,
        str(result_path / util.headtime("temp_bias_fit.json")),
    )
    # Option A: publish the temp-bias reference that EC corr reads, so the new
    # workflow is self-contained on a fresh machine.  Only write it when the
    # path is absent: the historical reference (for the migrated versions) is
    # kept untouched because the frozen-oracle EC was produced against it and it
    # is the source of strict identity.  A truly fresh run (no historical file)
    # gets this produced copy.
    ref = rt.payload.ec.tb_ref_path
    if ref:
        ref_path = rt.data_dir / ref
        if not ref_path.exists():
            ref_path.parent.mkdir(parents=True, exist_ok=True)
            util.json_save(result, str(ref_path))
    return result


def _apply_bias_filter(data_all, bias_min: Optional[float]):
    if bias_min is None:
        return data_all
    return data_all[data_all[:, 4] >= bias_min]


# --------------------------------------------------------------------------- #
# EC points + global fit
# --------------------------------------------------------------------------- #
def build_ec_points(rt: RuntimeConfig, items: List, source_kind: str) -> List[List[ECPoint]]:
    """Build per-channel ECPoint lists from (measurement, fp) pairs.

    Only channels < channel_count are kept (10B/11B legacy 3-channel behavior
    is handled by channel_count).  Energy comes from the runtime energy map.
    """
    n = rt.payload.ec.channel_count
    per_channel: List[List[ECPoint]] = [[] for _ in range(n)]
    for m, fp in items:
        energy = rt.energies[m.id]
        for ch, fit in enumerate(fp.fit_result):
            if ch >= n:
                continue
            if fit is None:
                continue
            per_channel[ch].append(
                ECPoint(
                    measurement_id=m.id,
                    channel=ch,
                    source_kind=source_kind,
                    energy=energy,
                    peak_center=fit["b"],
                    peak_center_err=fit["b_err"],
                    resolution=fit["resolution"],
                    resolution_err=fit["resolution_err"],
                    enabled=channel_use(m, ch),
                )
            )
    return per_channel


def global_ec(rt: RuntimeConfig, src_pts: List[List[ECPoint]], x_pts: List[List[ECPoint]],
              result_path: Path) -> List[dict]:
    """EC energy/resolution global fit (protected kernel).

    Mirrors the legacy ``ec_fit``: points are split at the K-edge, each half
    gets an independent quadratic center fit and a resolution fit, and the
    results are written to ``ec_logs/``.  ``plot.ec_plot`` is called directly
    with the same argument shape it expects.
    """
    from lib_plot import plot

    result_path.mkdir(parents=True, exist_ok=True)
    pb = rt.payload.ec
    n = pb.channel_count
    src = [p for ch in src_pts for p in ch if p.enabled]
    xr = [p for ch in x_pts for p in ch if p.enabled]
    all_pts = src + xr
    all_pts.sort(key=lambda p: p.energy)

    # per-channel arrays: channels may cover different point sets (N1 CLYC gets
    # source anchors only, GAGG gets source + X-ray), unlike the legacy versions
    center, center_err, resolution, resolution_err = [], [], [], []
    en_by_ch = []
    energies_all = None
    for ch in range(n):
        pts = [p for p in all_pts if p.channel == ch]
        en = np.array([p.energy for p in pts])
        en_by_ch.append(en)
        center.append(np.array([p.peak_center for p in pts]))
        center_err.append(np.array([p.peak_center_err for p in pts]))
        resolution.append(np.array([p.resolution for p in pts]))
        resolution_err.append(np.array([p.resolution_err for p in pts]))
        if energies_all is None and en.size:
            energies_all = en

    result = [{} for _ in range(n)]
    for ch in range(n):
        en, c, ce, r, re = en_by_ch[ch], center[ch], center_err[ch], resolution[ch], resolution_err[ch]
        q_low = en < pb.energy_split_low
        q_high = en >= pb.energy_split_high
        try:
            res_low, res_low_err = _resolution_fit(pb.resolution_method, en[q_low], r[q_low], re[q_low])
        except Exception:
            res_low, res_low_err = None, None
        try:
            res_high, res_high_err = _resolution_fit(pb.resolution_method, en[q_high], r[q_high], re[q_high])
        except Exception:
            res_high, res_high_err = None, None
        if _ec_form(pb, ch) == "linear":
            # single unsplit line; mirrored into low/high so the plotting layer
            # (which always reads EC_low/EC_high) keeps working
            ec, ec_err = _center_fit(en, c, ce, deg=1)
            result[ch] = {
                "channel": ch,
                "EC_low": ec,
                "EC_low_err": ec_err,
                "EC_high": ec,
                "EC_high_err": ec_err,
                "resolution_low": res_low,
                "resolution_low_err": res_low_err,
                "resolution_high": res_high,
                "resolution_high_err": res_high_err,
                "ec_form": "linear",
            }
        elif _ec_form(pb, ch) == "quadratic":
            # single unsplit quadratic (N1: no Gd K-edge split)
            ec, ec_err = _center_fit(en, c, ce, deg=2)
            result[ch] = {
                "channel": ch,
                "EC_low": ec,
                "EC_low_err": ec_err,
                "EC_high": ec,
                "EC_high_err": ec_err,
                "resolution_low": res_low,
                "resolution_low_err": res_low_err,
                "resolution_high": res_high,
                "resolution_high_err": res_high_err,
                "ec_form": "quadratic",
            }
        else:
            ec_low, ec_low_err = _center_fit(en[q_low], c[q_low], ce[q_low])
            ec_high, ec_high_err = _center_fit(en[q_high], c[q_high], ce[q_high])
            result[ch] = {
                "channel": ch,
                "EC_low": ec_low,
                "EC_low_err": ec_low_err,
                "EC_high": ec_high,
                "EC_high_err": ec_high_err,
                "resolution_low": res_low,
                "resolution_low_err": res_low_err,
                "resolution_high": res_high,
                "resolution_high_err": res_high_err,
            }
        util.json_save(result[ch], str(result_path / util.headtime(f"ec_coef_sci_ch{ch}.json")))
        save_data = np.array([en, c], dtype=np.float64)
        np.save(str(result_path / util.headtime(f"ec_data_ch{ch}.npy")), arr=save_data)

    src_energy = np.sort(np.array(sorted({p.energy for p in src})))
    x_energy = np.sort(np.array(sorted({p.energy for p in xr})))
    # plot.ec_plot always expects 4 channels; 10B/11B are physically 3 channels
    # and legacy pads channel 3 as a copy of channel 0 (see migration note B8)
    center4 = _pad_to4(center)
    result4 = _pad_to4(result)
    src_result = _pad_group_4ch(_group_4ch(src))
    x_result = _pad_group_4ch(_group_4ch(xr))
    try:
        plot.ec_plot(
            energies_all, center4, result4,
            src_energy, x_energy, src_result, x_result,
            str(result_path), pb.energy_split_low, pb.energy_split_high,
        )
    except Exception as e:
        # a channel with too few points has no resolution model; the fit
        # coefficients are still valid, so do not fail the whole run on the plot
        print(f"WARNING: ec_plot skipped ({e})")
    return result


def _pad_to4(seq):
    if len(seq) >= 4:
        return seq[:4]
    out = list(seq)
    while len(out) < 4:
        out.append(out[0])
    return out


def _group_4ch(points: List[ECPoint]) -> List[list]:
    """Rebuild per-measurement list-of-4-channel fit dicts for ec_plot.

    Each element is a list of dicts ``[{b, b_err, resolution,
    resolution_err}, ...]``, exactly the shape legacy ``ec_fit`` passes to
    ``plot.ec_plot``.  Ordered by energy to make the scatter arrays internally
    consistent.
    """
    by_id = {}
    for p in points:
        by_id.setdefault(p.measurement_id, {})[p.channel] = p
    items = sorted(by_id.items(), key=lambda kv: min(q.energy for q in kv[1].values()))
    out = []
    for _mid, chmap in items:
        fit = []
        for ch in sorted(chmap):
            p = chmap[ch]
            fit.append({
                "b": p.peak_center,
                "b_err": p.peak_center_err,
                "resolution": p.resolution,
                "resolution_err": p.resolution_err,
            })
        out.append(fit)
    return out


def _pad_group_4ch(groups: List[list]) -> List[list]:
    """Pad each per-measurement fit list to 4 channels (ch3 = ch0 copy)."""
    out = []
    for g in groups:
        if len(g) >= 4:
            out.append(g[:4])
        else:
            padded = list(g)
            while len(padded) < 4:
                padded.append(padded[0])
            out.append(padded)
    return out


def _center_fit(energy, center, center_err, deg: int = 2):
    popt, pcov = np.polyfit(center, energy, deg=deg, full=False, cov=True, w=1.0 / center_err)
    perr = np.sqrt(np.diag(pcov))
    return list(popt), list(perr)


def _ec_form(pb, ch: int) -> str:
    """E-C center form for a channel (default: K-edge piecewise quadratic)."""
    return (getattr(pb, "ec_form", None) or {}).get(str(ch), "piecewise_quadratic")


def _resolution_fit(method, energy, resolution, resolution_err):
    if method == "polyfit":
        p0, pcov = util.resolution_polyfit(energy, resolution, resolution_err)
        perr = np.sqrt(np.diag(pcov))
        return list(p0), list(perr)
    if method == "exprfit":
        return util.resolution_ExprFit(energy, resolution, resolution_err)
    if method == "lmfit":
        return util.resolution_lmfit(energy, resolution, resolution_err)
    raise ValueError(f"unknown resolution method {method!r}")



