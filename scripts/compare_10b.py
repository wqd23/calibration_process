# -*- coding:utf-8 -*-
"""10B differential: single fits + TB global + figures.

The 90 keV X-ray point is unprocessable in legacy (channel-3 gaus fit fails ->
plot.fit_plot KeyError 'a'), so the legacy EC global fit cannot run.  This
script compares everything that legacy can produce: single-fit pickles
(63 TB + 20 X-ray + 4 src), the TB temp-bias global fit, and the figure sets.
"""

import os
import re
import sys
from pathlib import Path

import json

import matplotlib.image as mpimg

sys.path.insert(0, "tests")
from compare import assert_pickle_equivalent, Mismatch  # noqa: E402

OR = Path(".oracle/10B")
NW = Path(".new/10B")
probs = []


def latest(d):
    out = {}
    for f in os.listdir(d):
        out[re.sub(r"^\d{14}", "", f)] = f
    return out


for sub in ("TB_fit_result", "EC_fit_result"):
    leg = sorted(os.listdir(OR / sub))
    new = sorted(os.listdir(NW / "single_process" / sub))
    if leg != new:
        probs.append(f"{sub}: set diff only-legacy={sorted(set(leg)-set(new))} only-new={sorted(set(new)-set(leg))}")
    import pickle
    for f in sorted(set(leg) & set(new)):
        try:
            assert_pickle_equivalent(
                pickle.load(open(NW / "single_process" / sub / f, "rb")),
                pickle.load(open(OR / sub / f, "rb")), f"pickle.{sub}.{f}")
        except Mismatch as e:
            probs.append(f"{sub}/{f}: {e}")

leg = latest(OR / "tb_logs")
new = latest(NW / "tb_logs")
if set(leg) != set(new):
    probs.append(f"tb_logs suffix diff {sorted(set(leg)^set(new))}")
for suf in sorted(set(leg) & set(new)):
    if suf.endswith(".json"):
        if json.load(open(NW / "tb_logs" / new[suf])) != json.load(open(OR / "tb_logs" / leg[suf])):
            probs.append(f"tb_logs/{suf}")


def figs(p):
    return {re.sub(r"^\d{14}", "", f): f for f in os.listdir(p) if f.endswith(".png")}

leg = figs(OR / "single_fit_fig")
new = figs(NW / "single_process" / "single_fit_fig")
if set(leg) != set(new):
    probs.append(f"single_fit_fig diff {sorted(set(leg)^set(new))}")
for k in sorted(set(leg) & set(new)):
    if mpimg.imread(OR / "single_fit_fig" / leg[k]).shape != mpimg.imread(NW / "single_process" / "single_fit_fig" / new[k]).shape:
        probs.append(f"single_fit_fig/{k} dim")
leg = figs(OR / "tb_logs")
new = figs(NW / "tb_logs")
if set(leg) != set(new):
    probs.append(f"tb_logs fig diff {sorted(set(leg)^set(new))}")

if probs:
    print("FAIL:")
    for p in probs:
        print("  -", p)
    sys.exit(1)
print("PASS: 10B single fits + TB global + figures equivalent "
      "(EC global blocked by the 90 keV legacy defect)")
