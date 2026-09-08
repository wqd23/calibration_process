# -*- coding:utf-8 -*-
"""New-architecture payload acceptance (scaffold) and deployment check.

These only touch the new YAML config layer + manifest + Pydantic schemas.
They never read ``config.json`` or instantiate legacy Operation classes.
"""

from __future__ import annotations

import os
from pathlib import Path

import yaml

from . import manifest as man
from .config_schema import AnalysisSchema, PayloadSchema

CONFIG_ROOT = Path("src/calibration_process/configs")
DATA = Path("data")

_OUTPUT_SUBDIRS = (
    "single_process/TB_fit_result",
    "single_process/EC_fit_result",
    "single_process/single_fit_fig",
    "tb_logs",
    "ec_logs",
)


def _cfg_root(ver: str, config_root: Path = None) -> Path:
    return (config_root or CONFIG_ROOT) / ver


def _data_dir(ver: str, data_root: Path = None) -> Path:
    return (data_root or DATA) / ver


def check_version(ver: str, fix: bool = False,
                  config_root: Path = None, data_root: Path = None) -> bool:
    """Validate a version's new config/manifest layer and data links.

    - ``payload.yaml`` / ``analysis.yaml`` must be strict-valid.
    - each ``*_manifest.yaml`` must be strict-valid, have unique ids, and every
      referenced file must exist under ``data/<ver>``.
    - output dirs are reported (created with ``--fix``).
    - returns True when everything is ready.
    """
    root = _cfg_root(ver, config_root)
    ok = True

    print(f"[{ver}] deployment check")
    if not root.is_dir():
        print(f"  MISSING  configs/{ver}/ (run `calib scaffold {ver}`)")
        return False

    for fname, schema in (("payload.yaml", PayloadSchema),
                          ("analysis.yaml", AnalysisSchema)):
        p = root / fname
        if not p.exists():
            print(f"  MISSING  configs/{ver}/{fname}")
            ok = False
            continue
        try:
            schema.model_validate(yaml.safe_load(open(p)))
            print(f"  OK       configs/{ver}/{fname}")
        except Exception as e:
            print(f"  INVALID  configs/{ver}/{fname}: {e}")
            ok = False

    data_dir = _data_dir(ver, data_root)
    raw = data_dir / "raw_data"
    if raw.is_symlink():
        state = "valid" if raw.exists() else "BROKEN"
        print(f"  symlink  data/{ver}/raw_data -> {os.readlink(raw)} ({state})")
        if not raw.exists():
            ok = False
    elif raw.is_dir():
        print(f"  real dir data/{ver}/raw_data (not a symlink, acceptable)")
    else:
        print(f"  ABSENT   data/{ver}/raw_data (run `calib scaffold {ver}` or link data)")
        ok = False

    for branch in ("tb", "ec_source", "ec_xray"):
        mp = root / f"{branch}_manifest.yaml"
        if not mp.exists():
            print(f"  MISSING  configs/{ver}/{branch}_manifest.yaml (run `calib discover`)")
            ok = False
            continue
        try:
            m = man.load_manifest(mp)
            ids = [x.id for x in m.measurements]
            if len(ids) != len(set(ids)):
                print(f"  INVALID  {branch}_manifest.yaml: duplicated ids")
                ok = False
            missing = []
            for x in m.measurements:
                for rel in x.science_files + x.hk_files + x.aux_files:
                    if not (data_dir / rel).exists():
                        missing.append(f"{x.id}:{rel}")
            if missing:
                print(f"  MISSING  {branch} files: {missing[:4]}{' ...' if len(missing) > 4 else ''}")
                ok = False
            else:
                print(f"  OK       {branch}_manifest.yaml ({len(ids)} measurements)")
        except Exception as e:
            print(f"  INVALID  {branch}_manifest.yaml: {e}")
            ok = False

    for sub in _OUTPUT_SUBDIRS:
        d = data_dir / sub
        if d.is_dir():
            print(f"  OK       data/{ver}/{sub}")
        elif fix:
            d.mkdir(parents=True, exist_ok=True)
            print(f"  CREATED  data/{ver}/{sub}")
        else:
            print(f"  NO-DIR   data/{ver}/{sub} (run with --fix)")
            ok = False

    print((f"[{ver}] READY" if ok else f"[{ver}] NOT READY"))
    return ok


# minimal template for a brand-new payload (reader defaults; human fills in)
def _template_payload(ver: str, import_from: str = "") -> dict:
    return {
        "version": ver,
        "tb": {
            "reader": "normal",
            "bin_width": 6,
            "adc_max": 16384.0,
            "science_dir": "raw_data",
            "fit_range_file": "fit_range_tb.yaml",
        },
        "ec": {
            "reader": "normal",
            "bin_width": 4,
            "adc_max": 16384.0,
            "x_path": "raw_data/x_data",
            "src_path": "raw_data/src_data",
            "energy_split_low": 49.0,
            "energy_split_high": 55.0,
            "ref_temp": 25.0,
            "ref_bias": 28.5,
            "tb_ref_path": "single_process/temp_bias_fit.json",
            "fit_range_file": "fit_range_ec_xray.yaml",
            "energy_map": {},
            "resolution_method": "polyfit",
            "channel_count": 4,
        },
    }


def scaffold_version(ver: str, data_dir: str = None,
                     config_root: Path = None, data_root: Path = None) -> None:
    """Create a new payload's config skeleton, data dirs and raw_data link."""
    root = _cfg_root(ver, config_root)
    data_dir_path = _data_dir(ver, data_root)
    root.mkdir(parents=True, exist_ok=True)

    if not (root / "payload.yaml").exists():
        with open(root / "payload.yaml", "w") as f:
            yaml.safe_dump(_template_payload(ver), f, sort_keys=False, allow_unicode=True)
        print(f"  created configs/{ver}/payload.yaml")
    else:
        print(f"  exists   configs/{ver}/payload.yaml (skipping)")

    for fname in ("analysis.yaml", "fit_range_tb.yaml",
                  "fit_range_ec_source.yaml", "fit_range_ec_xray.yaml",
                  "tb_manifest.yaml", "ec_source_manifest.yaml",
                  "ec_xray_manifest.yaml"):
        p = root / fname
        if not p.exists():
            if fname == "analysis.yaml":
                data = {"tb": {"default_bkg": "lin", "background_overrides": {},
                               "peak_overrides": {}, "qa": {}},
                        "ec_source": {"default_bkg": "lin", "background_overrides": {},
                                      "peak_overrides": {}, "qa": {}},
                        "ec_xray": {"default_bkg": "lin", "background_overrides": {},
                                    "peak_overrides": {}, "qa": {}}}
            elif fname == "fit_range_tb.yaml":
                data = {"measurements": {}}
            elif fname in ("fit_range_ec_source.yaml", "fit_range_ec_xray.yaml"):
                data = {"measurements": {}}
            else:  # *_manifest.yaml
                data = {"version": ver,
                        "branch": fname.replace("_manifest.yaml", ""),
                        "measurements": []}
            with open(p, "w") as f:
                yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)
            print(f"  created configs/{ver}/{fname}")

    for sub in _OUTPUT_SUBDIRS:
        (data_dir_path / sub).mkdir(parents=True, exist_ok=True)
        print(f"  created data/{ver}/{sub}")

    link = data_dir_path / "raw_data"
    if data_dir and not os.path.exists(link):
        os.symlink(data_dir, link)
        print(f"  symlinked data/{ver}/raw_data -> {data_dir}")
    elif not os.path.exists(link):
        print(f"  NOTE: run `calib scaffold {ver} <path>` or link data/{ver}/raw_data")

    print("\n  Next steps:")
    print(f"  1. Edit configs/{ver}/payload.yaml (reader/bin_width/adc_max/paths)")
    print(f"  2. `calib discover {ver} <branch>`; review the manifest")
    print(f"  3. Fill fit ranges in configs/{ver}/fit_range_*.yaml")
    print(f"  4. `calib fit {ver} <branch>` and `calib global {ver} <branch>`")
    return
