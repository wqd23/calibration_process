# -*- coding:utf-8 -*-
"""Example: reuse intermediate data + inject a custom pipeline step.

Two separate demos:
  A) Read the on-disk intermediate products (single-fit pickle, EC npy, TB/EC
     json) for extra analysis, from this project or another one.
  B) Build a custom pipeline: run single-fit(s), build TB points, run a CUSTOM
     step on the typed points, then hand them to the existing global fit.

Run:  python scripts/extra_analysis_example.py 09
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from calibration_process import pipeline  # noqa: E402
from calibration_process import manifest as man  # noqa: E402
from calibration_process.workflows import common as stages  # noqa: E402

VER = sys.argv[1] if len(sys.argv) > 1 else "09"


# --------------------------------------------------------------------------- #
# A) Read intermediate products
# --------------------------------------------------------------------------- #
def demo_read_intermediate():
    import dill
    import json

    print("\n== A) reading intermediate products ==")
    out_root = Path(".new") / VER

    # single-fit pickle (dill).  Full spectrum / x / tel / fit_result / config.
    pkl = out_root / "single_process/EC_fit_result/15keV.pickle"
    if pkl.exists():
        d = dill.load(open(pkl, "rb"))
        ch0 = d["fit_result"][0]
        print(f"  pickle {pkl.name}: center ch0 = {ch0['b']:.2f} +/- {ch0['b_err']:.2f}, "
              f"resolution = {ch0['resolution']:.4f}")
        print(f"    keys: {list(d.keys())} | spectrum[0].shape={d['spectrum'][0].shape}")
    else:
        print(f"  pickle {pkl.name}: not present")

    # EC data arrays (npz-free, plain .npy): [energy, center_per_channel]
    for npy in sorted((out_root / "ec_logs").glob("*ec_data_ch0.npy"))[-1:]:
        arr = np.load(npy)
        print(f"  npy {npy.name}: shape={arr.shape}, energy[0..3]={arr[0][:4]}")
    for j in sorted((out_root / "ec_logs").glob("*ec_coef_sci_ch0.json"))[-1:]:
        c = json.load(open(j))
        print(f"  json {j.name}: EC_low={c['EC_low']} EC_high={c['EC_high']}")

    # TB/EC global (plain JSON)
    for j in sorted((out_root / "tb_logs").glob("*temp_bias_fit.json"))[-1:]:
        tb = json.load(open(j))
        print(f"  tb json {j.name}: ch0 G0={tb[0]['G0']:.4e} V0={tb[0]['V0']:.3f}")


# --------------------------------------------------------------------------- #
# B) Custom pipeline: single fit -> points -> CUSTOM step -> global fit
# --------------------------------------------------------------------------- #
def demo_custom_pipeline():
    print("\n== B) custom pipeline (insert a step) ==")
    rt = pipeline.load_rt(VER)
    manifest = man.load_manifest(pipeline._manifest_path(VER, "tb"))
    out_root = Path(".new") / VER

    # step 1: ensure single fits exist, then load them from the store
    for m in man.filtered_measurements(manifest)[:3]:
        stages.run_single_fit(rt, "tb", m, out_root)
    items = [
        (m, stages.load_single_fp_from_store(rt, "tb", m, out_root))
        for m in man.filtered_measurements(manifest)
    ]
    per_channel = stages.build_tb_points(rt, items)
    print(f"  built points: per-channel counts = {[len(p) for p in per_channel]}")

    # ---- CUSTOM step: e.g. only keep points with temperature < 35 degC ----
    custom = []
    for ch, pts in enumerate(per_channel):
        kept = [p for p in pts if p.temperature < 35.0]
        print(f"    ch{ch}: kept {len(kept)}/{len(pts)} points (T<35C)")
        custom.append(kept)

    # step 3: hand the transformed points to the existing global fit
    result_path = out_root / "tb_logs_custom"
    res = stages.global_tb(rt, custom, result_path)
    print(f"  custom global_tb done: {len(res)} channels fit (result_path={result_path})")
    return res


if __name__ == "__main__":
    demo_read_intermediate()
    demo_custom_pipeline()
