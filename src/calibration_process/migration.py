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
    tb_out = {
        "reader": tb["reader"],
        "bin_width": tb["bin_width"],
        "adc_max": tb["adc_max"],
        "science_dir": tb["science_dir"],
        "fit_range_file": "fit_range_tb.yaml",
    }
    for k in ("tb_fit_p0", "tb_fit_maxfev", "bias_min_filter", "tb_file_map",
              "tb_fit_method"):
        if tb.get(k) is not None and tb.get(k) != "curvefit":
            tb_out[k] = tb[k]
    ec_out = {
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
        "xray_bkg_rotation": ec.get("xray_bkg_rotation", "circle"),
    }
    for k in ("xray_drop_energies", "xray_drop_old", "xray_require_hk",
              "xray_require_fit_range", "src_bkg_map", "xray_drop_substrs",
              "xray_name_index", "xray_config_file", "time_cut",
              "xray_single_file", "xray_reader", "src_reader"):
        if ec.get(k) is not None:
            ec_out[k] = ec[k]
    return {"version": ver, "tb": tb_out, "ec": ec_out}


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


# per-version parameter table (reader, scientific params, file-selection rules)
def _scan_tb_ids(ver, cfg) -> list:
    if ver == "12B":
        from .workflows.versions.v12B import EXCLUDE
        import json
        tm = json.load(open(cfg["tb"]["file_map"].replace("{ver}", ver)))
        return sorted(
            f"{p['temp_setpoint_C']}C_{p['bias_code']}"
            for p in tm if (p["temp_setpoint_C"], p["bias_code"]) not in EXCLUDE
        )
    full = Path(cfg["tb"]["path"])
    if ver == "05B":
        return [f for f in os.listdir(full) if "rundata" in f and "50C" not in f]
    if ver == "03B":
        return [f for f in os.listdir(full)
                if "rundata" in f and "baseline" not in f and "CI" not in f and "50C" not in f]
    if ver == "10B":
        return [f for f in os.listdir(full) if "observe" in f and "50C_265" not in f]
    if ver == "11B":
        from .workflows.versions.v11B import TB_SUBDIRS, TB_REMOVE
        files = []
        for sub in TB_SUBDIRS:
            for item in Path(full / sub).glob("*_observe*.dat"):
                name = item.name
                if any(k in name for k in ("on", "off")):
                    continue
                files.append(f"{sub}/{name}")
        for rel in TB_REMOVE:
            files.remove(rel)
        files = [f for f in files if "_50_Cs_2" not in f]
        files.sort()
        # fit_range.json is keyed by the file stem (no extension)
        return [Path(f).stem for f in files]
    files = [f for f in os.listdir(full) if os.path.splitext(f)[1] == ".txt"]
    if ver == "09":
        for excl in v09.TB_EXCLUDE:
            files.remove(excl)
    return files


def _scan_ec_source_ids(ver, cfg) -> list:
    if ver in ("12B", "04", "03B", "10B", "11B"):
        from .workflows.versions import v03B, v04, v10B, v11B, v12B
        mod = {"03B": v03B, "04": v04, "10B": v10B, "11B": v11B, "12B": v12B}[ver]
        return list(mod.SRC_LIST)
    if ver == "07":
        full = Path(cfg["ec"]["src_path"])
        return [f for f in os.listdir(full)
                if "src" in f and "_bk_" not in f and "bkg" not in f]
    if ver == "05B":
        full = Path(cfg["ec"]["src_path"])
        return [f for f in os.listdir(full) if "rundata" in f and "bkg" not in f]
    return list(v09.SRC_LIST)


def _scan_ec_xray_ids(ver, cfg, fit_range) -> list:
    if ver == "03B":
        full = Path(cfg["ec"]["x_path"])
        x_ch = [f for f in os.listdir(full) if "_rundata" in f and "CI" not in f]
        return sorted({f.split("_")[1] for f in x_ch})
    if ver == "10B":
        full = Path(cfg["ec"]["x_path"])
        x_ch = [f for f in os.listdir(full) if "_ch" in f and "hk" not in f]
        return sorted({f.split("_")[2] for f in x_ch})
    if ver == "11B":
        full = Path(cfg["ec"]["x_path"])
        x_ch = [f for f in os.listdir(full) if "_ch" in f and "hk" not in f]
        return sorted(e for e in ({f.split("_")[2] for f in x_ch}) if e != "20")
    if ver == "05B":
        full = Path(cfg["ec"]["x_path"])
        return [f for f in os.listdir(full) if "_observe.dat" in f and "XM_22" not in f]
    if ver == "12B":
        full = Path(cfg["ec"]["x_path"])
        x_ch = [f for f in os.listdir(full) if f.endswith(".dat") and "_ch" in f and "old" not in f]
        energies = sorted({f.split("_")[1] for f in x_ch}, key=int)
        return [e for e in energies if e in fit_range]
    if ver == "04":
        full = Path(cfg["ec"]["x_path"])
        x_ch = [f for f in os.listdir(full) if "_ch" in f and "_18p0_" not in f]
        return sorted({f.split("_")[3] for f in x_ch})
    if ver == "07":
        full = Path(cfg["ec"]["x_path"])
        x_ch = [f for f in os.listdir(full)
                if "_ch" in f and not any(s in f for s in ["40p0", "15p0", "12p0", "99p9", "90p1"])]
        return sorted({f.split("_")[2] for f in x_ch})
    full = Path(cfg["ec"]["x_path"])
    x_ch = [f for f in os.listdir(full) if "_ch" in f and "65keV_" not in f]
    return sorted({f.split("_")[0] for f in x_ch})


def _params(ver: str, cfg: dict) -> dict:
    tb = cfg["tb"]
    ec = cfg["ec"]
    return _PARAMS[ver](ver, cfg, tb, ec)


_P09 = {
    "selected_ids": None,  # filled per-version below via scan
    "tb": {"reader": "09", "bin_width": 6, "adc_max": 65535.0},
    "ec": {
        "reader": "09", "bin_width": 10, "adc_max": 65535.0,
        "energy_split_low": 49.0, "energy_split_high": 55.0,
        "ref_temp": 25.0, "ref_bias": 28.5,
        "resolution_method": "polyfit", "channel_count": 4,
        "xray_drop_energies": ["65keV"], "xray_bkg_rotation": "circle",
        "xray_drop_old": False, "xray_require_hk": False,
        "xray_require_fit_range": False,
    },
    "analysis": {"tb": {"default_bkg": "lin"},
                 "ec_source": {"default_bkg": "lin"},
                 "ec_xray": {"default_bkg": None}},
}


def _p09(ver, cfg, tb, ec):
    p = dict(_P09)
    p["tb"] = {**_P09["tb"], "science_dir": _rel(tb["path"])}
    p["ec"] = {**_P09["ec"], "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
               "tb_ref_path": _rel(ec["tb_result_path"])}
    p["selected_ids"] = {"tb": _scan_tb_ids(ver, cfg),
                         "ec_source": _scan_ec_source_ids(ver, cfg),
                         "ec_xray": _scan_ec_xray_ids(ver, cfg, json.load(open(_load(ver, "fit_range.json"))))}
    return p


_P12B = {
    "tb": {"reader": "12b", "bin_width": 6, "adc_max": 16384.0,
           "tb_fit_p0": [-0.02, 0.07, 24.4, -35.0, -1000.0],
           "tb_fit_maxfev": 100000, "bias_min_filter": 27.25,
           "science_dir": "raw_data", "tb_file_map": "single_process/tb_file_map.json"},
    "ec": {
        "reader": "12b", "bin_width": 4, "adc_max": 16384.0,
        "energy_split_low": 49.0, "energy_split_high": 55.0,
        "ref_temp": 25.0, "ref_bias": 28.5,
        "resolution_method": "polyfit", "channel_count": 4,
        "xray_drop_old": True, "xray_require_hk": True, "xray_require_fit_range": True,
        "xray_bkg_rotation": "fixed",
    },
    "analysis": {"tb": {"default_bkg": "lin"},
                 "ec_source": {"default_bkg": "lin"},
                 "ec_xray": {"default_bkg": None}},
}


def _p12b(ver, cfg, tb, ec):
    p = dict(_P12B)
    p["ec"] = {**_P12B["ec"], "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
               "tb_ref_path": _rel(ec["tb_result_path"])}
    p["selected_ids"] = {"tb": _scan_tb_ids(ver, cfg),
                         "ec_source": _scan_ec_source_ids(ver, cfg),
                         "ec_xray": _scan_ec_xray_ids(ver, cfg, json.load(open(_load(ver, "fit_range.json"))))}
    return p


def _load(ver, name):
    return f"data/{ver}/single_process/{name}"


def _p04(ver, cfg, tb, ec):
    p = {
        "tb": {"reader": "04", "bin_width": 6, "adc_max": 65535.0,
               "science_dir": _rel(tb["path"])},
        "ec": {
            "reader": "04", "bin_width": 10, "adc_max": 65535.0,
            "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
            "energy_split_low": 49.0, "energy_split_high": 55.0,
            "ref_temp": 25.0, "ref_bias": 28.5,
            "tb_ref_path": _rel(ec["tb_result_path"]),
            "resolution_method": "exprfit", "channel_count": 4,
            "xray_drop_substrs": ["_18p0_"], "xray_name_index": 3,
            "xray_bkg_rotation": "circle",
        },
        "analysis": {"tb": {"default_bkg": "lin"}, "ec_source": {"default_bkg": "lin"},
                     "ec_xray": {"default_bkg": "lin"}},
    }
    p["selected_ids"] = {"tb": _scan_tb_ids(ver, cfg),
                         "ec_source": _scan_ec_source_ids(ver, cfg),
                         "ec_xray": _scan_ec_xray_ids(ver, cfg, json.load(open(_load(ver, "fit_range.json"))))}
    return p


def _p03B(ver, cfg, tb, ec):
    p = {
        "tb": {"reader": "03b", "bin_width": 6, "adc_max": 16384.0,
               "science_dir": _rel(tb["path"])},
        "ec": {
            "reader": "03b", "src_reader": "03b-src", "bin_width": 10,
            "adc_max": 16384.0,
            "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
            "energy_split_low": 49.0, "energy_split_high": 51.0,
            "ref_temp": 25.0, "ref_bias": 28.5,
            "tb_ref_path": _rel(ec["tb_result_path"]),
            "resolution_method": "lmfit", "channel_count": 4,
            "xray_name_index": 1, "xray_drop_substrs": ["CI"],
            "xray_bkg_rotation": "circle",
        },
        "analysis": {"tb": {"default_bkg": "lin"}, "ec_source": {"default_bkg": "lin"},
                     "ec_xray": {"default_bkg": "lin"}},
    }
    p["selected_ids"] = {"tb": _scan_tb_ids(ver, cfg),
                         "ec_source": _scan_ec_source_ids(ver, cfg),
                         "ec_xray": _scan_ec_xray_ids(ver, cfg, json.load(open(_load(ver, "fit_range.json"))))}
    return p


def _p10B(ver, cfg, tb, ec):
    p = {
        "tb": {"reader": "10b", "bin_width": 6, "adc_max": 16384.0,
               "science_dir": _rel(tb["path"])},
        "ec": {
            "reader": "10b", "bin_width": 4, "adc_max": 16384.0,
            "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
            "energy_split_low": 49.0, "energy_split_high": 55.0,
            "ref_temp": 25.0, "ref_bias": 28.5,
            "tb_ref_path": _rel(ec["tb_result_path"]),
            "resolution_method": "polyfit", "channel_count": 3,
            "xray_name_index": 2, "xray_bkg_rotation": "fixed",
        },
        "analysis": {"tb": {"default_bkg": "lin"}, "ec_source": {"default_bkg": "lin"},
                     "ec_xray": {"default_bkg": "lin"}},
    }
    p["selected_ids"] = {"tb": _scan_tb_ids(ver, cfg),
                         "ec_source": _scan_ec_source_ids(ver, cfg),
                         "ec_xray": _scan_ec_xray_ids(ver, cfg, json.load(open(_load(ver, "fit_range.json"))))}
    return p


def _p11B(ver, cfg, tb, ec):
    p = {
        "tb": {"reader": "11b", "bin_width": 6, "adc_max": 16384.0,
               "science_dir": _rel(tb["path"]), "tb_fit_method": "lmfit"},
        "ec": {
            "reader": "11b", "bin_width": 4, "adc_max": 16384.0,
            "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
            "energy_split_low": 49.0, "energy_split_high": 55.0,
            "ref_temp": 25.0, "ref_bias": 28.5,
            "tb_ref_path": _rel(ec["tb_result_path"]),
            "resolution_method": "polyfit", "channel_count": 3,
            "xray_name_index": 2, "xray_bkg_rotation": "fixed",
        },
        "analysis": {"tb": {"default_bkg": "lin"}, "ec_source": {"default_bkg": "lin"},
                     "ec_xray": {"default_bkg": "lin"}},
    }
    p["selected_ids"] = {"tb": _scan_tb_ids(ver, cfg),
                         "ec_source": _scan_ec_source_ids(ver, cfg),
                         "ec_xray": _scan_ec_xray_ids(ver, cfg, json.load(open(_load(ver, "fit_range.json"))))}
    return p


def _p05B(ver, cfg, tb, ec):
    import json as _j
    tc = _j.load(open(_load(ver, "time_cut.json")))
    p = {
        "tb": {"reader": "normal", "bin_width": 6, "adc_max": 16384.0,
               "science_dir": _rel(tb["path"])},
        "ec": {
            "reader": "normal", "xray_reader": "xray", "bin_width": 10,
            "adc_max": 16384.0,
            "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
            "energy_split_low": 49.0, "energy_split_high": 52.0,
            "ref_temp": 25.0, "ref_bias": 28.5,
            "tb_ref_path": _rel(ec["tb_result_path"]),
            "resolution_method": "polyfit", "channel_count": 4,
            "xray_single_file": True,
            # keep the full path so it matches the legacy config_file string
            "xray_config_file": ec["x_config"],
            "time_cut": tc,
        },
        "analysis": {"tb": {"default_bkg": "lin"}, "ec_source": {"default_bkg": "lin"},
                     "ec_xray": {"default_bkg": "lin"}},
    }
    p["selected_ids"] = {"tb": _scan_tb_ids(ver, cfg),
                         "ec_source": _scan_ec_source_ids(ver, cfg),
                         "ec_xray": _scan_ec_xray_ids(ver, cfg, _j.load(open(_load(ver, "fit_range.json"))))}
    return p


def _p07(ver, cfg, tb, ec):
    p = {
        "tb": {"reader": "07", "bin_width": 6, "adc_max": 65535.0,
               "science_dir": _rel(tb["path"])},
        "ec": {
            "reader": "07", "bin_width": 10, "adc_max": 65535.0,
            "x_path": _rel(ec["x_path"]), "src_path": _rel(ec["src_path"]),
            "energy_split_low": 49.0, "energy_split_high": 55.0,
            "ref_temp": 25.0, "ref_bias": 28.5,
            "tb_ref_path": _rel(ec["tb_result_path"]),
            "resolution_method": "polyfit", "channel_count": 4,
            "xray_drop_substrs": ["40p0", "15p0", "12p0", "99p9", "90p1"],
            "xray_name_index": 2, "xray_bkg_rotation": "circle",
        },
        "analysis": {"tb": {"default_bkg": "lin"}, "ec_source": {"default_bkg": "lin"},
                     "ec_xray": {"default_bkg": "lin"}},
    }
    p["selected_ids"] = {"tb": _scan_tb_ids(ver, cfg),
                         "ec_source": _scan_ec_source_ids(ver, cfg),
                         "ec_xray": _scan_ec_xray_ids(ver, cfg, json.load(open(_load(ver, "fit_range.json"))))}
    return p


_PARAMS = {"09": _p09, "12B": _p12b, "04": _p04, "07": _p07, "05B": _p05B,
           "03B": _p03B, "10B": _p10B, "11B": _p11B}


def main():
    for ver in ["03B", "04", "05B", "07", "10B", "11B", "09", "12B"]:
        migrate_version(ver)


if __name__ == "__main__":
    main()
