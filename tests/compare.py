# -*- coding:utf-8 -*-
"""Differential (legacy vs new) comparison helpers.

Default comparison is **strict / exact**.  A tolerance is only allowed for a
specific field that is demonstrably a meaningless float reorder within the
scientific kernel; those exceptions are listed in ``TOLERANCE_EXCEPTIONS``
with the recorded reason.
"""

from __future__ import annotations

import os
from typing import Dict, Tuple

import numpy as np


class Mismatch(Exception):
    pass


# Per-field tolerances.  Keyed by ("kind", field) with allowed rel/abs tol.
# Kept empty unless a field is proven to be a meaningless float change.
TOLERANCE_EXCEPTIONS: Dict[Tuple[str, str], Dict[str, float]] = {}


def tol_for(kind: str, field: str):
    return TOLERANCE_EXCEPTIONS.get((kind, field), {"rtol": 0.0, "atol": 0.0})


def _is_lambda(v) -> bool:
    return callable(v) and not isinstance(v, (np.ndarray,))


def assert_float_equal(a, b, kind, field, exact=True):
    if np.isnan(a) and np.isnan(b):
        return
    t = tol_for(kind, field)
    rtol = t["rtol"]
    atol = t["atol"]
    if not np.isclose(a, b, rtol=rtol, atol=atol):
        raise Mismatch(
            f"[{kind}.{field}] new={a!r} legacy={b!r} "
            f"(rtol={rtol}, atol={atol})"
        )


def assert_array_equal(a, b, kind, field):
    # list-of-arrays (per channel, possibly different lengths): compare channel
    # by channel at the C level (each channel array is homogeneous)
    if isinstance(a, (list,)) and a and isinstance(a[0], np.ndarray):
        if len(a) != len(b):
            raise Mismatch(f"[{kind}.{field}] channel count {len(a)} != {len(b)}")
        for i, (ca, cb) in enumerate(zip(a, b)):
            assert_array_equal(ca, cb, kind, f"{field}.ch{i}")
        return
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        raise Mismatch(f"[{kind}.{field}] shape {a.shape} != {b.shape}")
    if not np.array_equal(a, b):
        t = tol_for(kind, field)
        if t["rtol"] == 0 and t["atol"] == 0:
            raise Mismatch(f"[{kind}.{field}] arrays differ (exact)")
        if not np.allclose(a, b, rtol=t["rtol"], atol=t["atol"]):
            raise Mismatch(f"[{kind}.{field}] arrays differ within tol")


def _compare_scalar(a, b, kind, field):
    if isinstance(a, float) or isinstance(b, float) or isinstance(a, np.floating):
        assert_float_equal(float(a), float(b), kind, field)
    elif isinstance(a, (np.ndarray, list)):
        assert_array_equal(np.asarray(a), np.asarray(b), kind, field)
    else:
        pass  # ints/bools: exact


def compare_fit_result(new_fr, legacy_fr, kind):
    """Compare two single-fit per-channel lists (each None or dict)."""
    if len(new_fr) != len(legacy_fr):
        raise Mismatch(f"[{kind}] channel count {len(new_fr)} != {len(legacy_fr)}")
    for ch, (nf, lf) in enumerate(zip(new_fr, legacy_fr)):
        if nf is None and lf is None:
            continue
        if nf is None or lf is None:
            raise Mismatch(f"[{kind}] channel {ch}: None mismatch")
        for k in nf:
            if k == "bkg":
                continue  # bkg holds a lambda; compared separately
            vn, vl = nf.get(k), lf.get(k)
            if isinstance(vn, (np.ndarray, list)) or isinstance(vl, (np.ndarray, list)):
                assert_array_equal(vn, vl, kind, f"fit.ch{ch}.{k}")
            elif isinstance(vn, float) or isinstance(vl, float):
                assert_float_equal(float(vn), float(vl), kind, f"fit.ch{ch}.{k}")
            elif isinstance(vn, dict) or isinstance(vl, dict):
                for kk in vn:
                    _compare_scalar(vn[kk], vl.get(kk), kind, f"fit.ch{ch}.{k}.{kk}")
            else:
                if vn != vl:
                    raise Mismatch(f"[{kind}] fit.ch{ch}.{k}: {vn!r} != {vl!r}")


def assert_pickle_equivalent(new_dict, legacy_dict, kind):
    if new_dict["file"] != legacy_dict["file"]:
        raise Mismatch(f"[{kind}] file field {new_dict['file']!r} != {legacy_dict['file']!r}")
    compare_fit_result(new_dict["fit_result"], legacy_dict["fit_result"], kind)
    for i, (a, b) in enumerate(zip(new_dict["spectrum"], legacy_dict["spectrum"])):
        assert_array_equal(a, b, kind, f"spectrum.{i}")
    for i, (a, b) in enumerate(zip(new_dict["x"], legacy_dict["x"])):
        assert_array_equal(a, b, kind, f"x.{i}")
    for key in legacy_dict["tel"]:
        a, b = new_dict["tel"].get(key), legacy_dict["tel"][key]
        assert_array_equal(a, b, kind, f"tel.{key}")


def assert_json_equivalent(new, legacy, kind):
    if type(new) is not type(legacy):
        raise Mismatch(f"[{kind}] json type {type(new)} != {type(legacy)}")
    if isinstance(new, list):
        if len(new) != len(legacy):
            raise Mismatch(f"[{kind}] list len {len(new)} != {len(legacy)}")
        for i, (a, b) in enumerate(zip(new, legacy)):
            assert_json_equivalent(a, b, f"{kind}[{i}]")
        return
    if isinstance(new, dict):
        for k in new:
            assert_json_equivalent(new[k], legacy.get(k), f"{kind}.{k}")
        return
    if isinstance(new, float) or isinstance(legacy, float):
        assert_float_equal(float(new), float(legacy), kind.split(".")[0], ".")
        return
    if isinstance(new, (list, np.ndarray)):
        assert_array_equal(new, legacy, kind.split(".")[0], ".")


def assert_npy_equivalent(new, legacy, kind):
    assert_array_equal(new, legacy, "npy", kind)


def compare_dir(legacy_dir, new_dir, expected_exts):
    """Compare the file set of two output dirs (names + count + arrays)."""
    issues = []
    leg = {f for f in os.listdir(legacy_dir) if f.endswith(expected_exts)}
    new = {f for f in os.listdir(new_dir) if f.endswith(expected_exts)}
    if leg != new:
        issues.append(f"file set differs: only-legacy={sorted(leg - new)}, only-new={sorted(new - leg)}")
    for f in sorted(leg & new):
        lp = os.path.join(legacy_dir, f)
        np_ = os.path.join(new_dir, f)
        if f.endswith(".json"):
            import json
            a = json.load(open(lp))
            b = json.load(open(np_))
            try:
                assert_json_equivalent(b, a, f)
            except Mismatch as e:
                issues.append(f"{f}: {e}")
        elif f.endswith(".npy"):
            try:
                assert_npy_equivalent(np.load(np_), np.load(lp), f)
            except Mismatch as e:
                issues.append(f"{f}: {e}")
        elif f.endswith(".pickle"):
            import pickle
            try:
                assert_pickle_equivalent(pickle.load(open(np_, "rb")), pickle.load(open(lp, "rb")), "pickle")
            except Mismatch as e:
                issues.append(f"{f}: {e}")
    return issues
