# -*- coding:utf-8 -*-
"""
deployment completeness check for each GRID payload
----------
Validates, for a single payload version, that every path referenced by
config.json is in place before any processing starts:

- input paths (raw data dirs, fit range / energy / bkg form config files,
  TB fit result consumed by EC) must exist, otherwise processing aborts at
  a much later and more confusing point
- output dirs (normally created by `just init`) are reported and can be
  created on the fly with --fix
- the raw_data symlink is reported for information (broken link surfaces
  as missing input paths above)

Exits non-zero when any required input is missing.
"""

import argparse
import os
import sys

from . import util_lib as util
from .__init__ import CFG_PATH

cfg = util.load_config(CFG_PATH)

# config keys pointing to input files/dirs that must exist before processing
INPUT_KEYS = (
    "path",
    "fit_range",
    "tb_result_path",
    "energy",
    "bkg_form",
    "x_path",
    "src_path",
    "time_cut",
    "x_config",
)
# config keys pointing to output dirs, created by `just init`
OUTPUT_KEYS = ("save_path", "save_fig_path", "result_path")


def _check_entries(ver: str, section: str, fix: bool) -> bool:
    """Check one section ("tb"/"ec") of cfg[ver]. Returns True if all OK."""
    all_ok = True
    for key, path in cfg[ver][section].items():
        if not isinstance(path, str):
            continue
        if key in INPUT_KEYS:
            if os.path.exists(path):
                print(f"  OK       {section}.{key} -> {path}")
            else:
                print(f"  MISSING  {section}.{key} -> {path}")
                all_ok = False
        elif key in OUTPUT_KEYS:
            if os.path.isdir(path):
                print(f"  OK       {section}.{key} -> {path}")
            elif fix:
                os.makedirs(path, exist_ok=True)
                print(f"  CREATED  {section}.{key} -> {path}")
            else:
                print(f"  NO-DIR   {section}.{key} -> {path} (run with --fix to create)")
                all_ok = False
        else:
            # unknown key: still report existence for visibility
            status = "OK      " if os.path.exists(path) else "MISSING "
            print(f"  {status} {section}.{key} -> {path}")
    return all_ok


def _report_raw_data_link(ver: str) -> None:
    link = f"data/{ver}/raw_data"
    if os.path.islink(link):
        target = os.readlink(link)
        state = "valid" if os.path.exists(link) else "BROKEN"
        print(f"  symlink  {link} -> {target} ({state})")
    elif os.path.isdir(link):
        print(f"  real dir {link} (not a symlink, acceptable)")
    else:
        print(f"  ABSENT   {link} (run `just init {ver} <path-to-data>`)")


def check_version(ver: str, fix: bool = False) -> bool:
    if ver not in cfg:
        print(f"unknown version {ver!r}, available: {', '.join(cfg.keys())}")
        return False
    print(f"[{ver}] deployment check")
    print("raw_data link:")
    _report_raw_data_link(ver)
    ok = True
    for section in ("tb", "ec"):
        print(f"{section} paths:")
        ok = _check_entries(ver, section, fix) and ok
    if ok:
        print(f"[{ver}] READY")
    else:
        print(f"[{ver}] NOT READY: fix the MISSING/NO-DIR items above")
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("ver", help="payload version, e.g. 03B")
    parser.add_argument(
        "--fix", action="store_true", help="create missing output dirs"
    )
    args = parser.parse_args()
    sys.exit(0 if check_version(args.ver, fix=args.fix) else 1)


if __name__ == "__main__":
    main()
