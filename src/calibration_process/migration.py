# -*- coding:utf-8 -*-
"""One-time migration of the historical JSON configuration to the new YAML.

The legacy config (config.json + data/<ver>/single_process/*.json) is read and
translated into the new configs/ tree.  The legacy JSON is left untouched as
the oracle; the new workflow only reads the YAML.

``old vs new resolved-configuration comparison`` is performed by the
regression tests (see tests/compare_config.py).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import yaml

from . import util_lib as util
from .__init__ import CFG_PATH
from .workflows.versions import v09


def _load_legacy(ver: str) -> dict:
    return util.load_config(CFG_PATH)[ver]


def _sp(ver: str) -> Path:
    return Path("data") / ver / "single_process"


def _w(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)
    print(f"  wrote {path}")


def _rel(p: str) -> str:
    # config paths are absolute-ish under data/<ver>/; store relative to the
    # version data dir, e.g. data/09/raw_data/temp_bias -> raw_data/temp_bias
    if p.startswith("data/"):
        parts = p.split("/")
        # parts = ["data", ver, ...]; drop "data" and the version dir
        return "/".join(parts[2:])
    return p


def migrate_version(ver: str, config_root: Path = None) -> None:
    print(f"[migrate] {ver}")
    cfg = _load_legacy(ver)
    sp = _sp(ver)
    config_root = config_root or (Path("src/calibration_process/configs") / ver)

    fit_range = json.load(open(sp / "fit_range.json"))
    bkg_form = json.load(open(sp / "bkg_form.json"))
    energy = json.load(open(sp / "ec_energy.json"))

    params = _params(ver, cfg)
    _w(config_root / "payload.yaml", _build_payload(ver, cfg, params, energy))
    _w(config_root / "analysis.yaml", _build_analysis(bkg_form, params))
    for branch, fname in (("tb", "fit_range_tb.yaml"),
                          ("ec_source", "fit_range_ec_source.yaml"),
                          ("ec_xray", "fit_range_ec_xray.yaml")):
        sel = {s: fit_range[s] for s in params["selected_ids"][branch] if s in fit_range}
        _w(config_root / fname, {"measurements": sel})


def _build_payload(ver, cfg, params, energy) -> dict:
    tb, ec = params["tb"], params["ec"]
    return {
        "version": ver,
        "tb": {
            "reader": tb["reader"],
            "bin_width": tb["bin_width"],
            "adc_max": tb["adc_max"],
            "science_dir": tb["science_dir"],
            "fit_range_file": "fit_range_tb.yaml",
        },
        "ec": {
            "reader": ec["reader"],
            "bin_width": ec["bin_width"],
            "adc_max": ec["adc_max"],
            "x_path": ec["x_path"],
            "src_path": ec["src_path"],
            "energy_split_low": ec["energy_split_low"],
            "energy_split_high": ec["energy_split_high"],
            "ref_temp": ec["ref_temp"],
            "ref_bias": ec["ref_bias"],
            "tb_ref_path": ec["tb_ref_path"],
            "fit_range_file": "fit_range_ec_xray.yaml",
            "energy_map": dict(energy),
            "resolution_method": ec["resolution_method"],
            "channel_count": ec["channel_count"],
            "xray_drop_energies": ec["xray_drop_energies"],
            "xray_bkg_rotation": ec["xray_bkg_rotation"],
        },
    }


def _build_analysis(bkg_form, params) -> dict:
    branches = {}
    for branch in ("tb", "ec_source", "ec_xray"):
        overrides = {
            k: v for k, v in bkg_form.items()
            if branch != "tb" and k in params["selected_ids"][branch]
        }
        branches[branch] = {
            "default_bkg": params["analysis"][branch]["default_bkg"],
            "background_overrides": overrides,
            "peak_overrides": {},
            "qa": {},
        }
    return branches


def _scan_tb_ids(ver, cfg) -> list:
    full = Path(cfg["tb"]["path"])
    files = [f for f in os.listdir(full) if os.path.splitext(f)[1] == ".txt"]
    for excl in v09.TB_EXCLUDE:
        files.remove(excl)
    return files


def _scan_ec_xray_ids(ver, cfg) -> list:
    full = Path(cfg["ec"]["x_path"])
    x_ch = [f for f in os.listdir(full) if "_ch" in f and "65keV_" not in f]
    return sorted({f.split("_")[0] for f in x_ch})


def _params(ver: str, cfg: dict) -> dict:
    tb = cfg["tb"]
    ec = cfg["ec"]
    return {
        "selected_ids": {
            "tb": _scan_tb_ids(ver, cfg),
            "ec_source": list(v09.SRC_LIST),
            "ec_xray": _scan_ec_xray_ids(ver, cfg),
        },
        "tb": {
            "reader": "09", "bin_width": 6, "adc_max": 65535.0,
            "science_dir": _rel(tb["path"]),
        },
        "ec": {
            "reader": "09", "bin_width": 10, "adc_max": 65535.0,
            "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
            "energy_split_low": 49.0, "energy_split_high": 55.0,
            "ref_temp": 25.0, "ref_bias": 28.5,
            "tb_ref_path": _rel(ec["tb_result_path"]),
            "resolution_method": "polyfit", "channel_count": 4,
            "xray_drop_energies": ["65keV"], "xray_bkg_rotation": "circle",
        },
        "analysis": {
            "tb": {"default_bkg": "lin"},
            "ec_source": {"default_bkg": "lin"},
            "ec_xray": {"default_bkg": None},
        },
    }


def main():
    for ver in ["09"]:
        migrate_version(ver)


if __name__ == "__main__":
    main()
