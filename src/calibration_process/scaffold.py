# -*- coding:utf-8 -*-
"""
scaffold a new payload version: generate config entry, directory structure,
and lib_reader skeleton.

Usage:
    python3 -m calibration_process.scaffold <ver> [--data-dir <path>]
"""

import argparse
import json
import os
import sys
from pathlib import Path

from .__init__ import CFG_PATH


# minimal config entry for a new version (inherits _defaults on load)
_TEMPLATE = {
    "tb": {
        "path": "data/{ver}/raw_data/tb_data"
    },
    "ec": {
        "tb_result_path": "data/{ver}/single_process/PLACEHOLDER_temp_bias_fit.json",
        "x_path": "data/{ver}/raw_data/x_data",
        "src_path": "data/{ver}/raw_data/src_data"
    }
}

# reader skeleton template
_READER_INIT = '''# -*- coding:utf-8 -*-
"""
Reader for GRID {ver} payload data.
"""
'''


def scaffold_version(ver: str, data_dir: str = None) -> None:
    """Generate config entry, directory structure, and reader skeleton."""

    # 1. config entry
    cfg_path = Path(CFG_PATH)
    with open(cfg_path, "r") as f:
        raw = json.load(f)

    if ver in raw:
        print(f"  WARNING: {ver} already exists in config.json, skipping config entry")
    else:
        # insert before closing brace, after last version entry
        entry = json.loads(json.dumps(_TEMPLATE).replace("{ver}", ver))
        # find insertion point: last version key before }
        keys = [k for k in raw.keys() if k != "_defaults"]
        if keys:
            # insert after last version
            raw[ver] = entry
        else:
            raw[ver] = entry
        with open(cfg_path, "w") as f:
            json.dump(raw, f, indent=4, ensure_ascii=False)
        print(f"  added {ver} to config.json")

    # 2. directory structure
    dirs = [
        f"data/{ver}/single_process/TB_fit_result",
        f"data/{ver}/single_process/EC_fit_result",
        f"data/{ver}/single_process/single_fit_fig",
        f"data/{ver}/tb_logs",
        f"data/{ver}/ec_logs",
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
        print(f"  created {d}")

    # 3. raw_data symlink (if data_dir provided)
    link_path = f"data/{ver}/raw_data"
    if data_dir and not os.path.exists(link_path):
        os.symlink(data_dir, link_path)
        print(f"  symlinked {link_path} -> {data_dir}")
    elif not os.path.exists(link_path):
        print(f"  NOTE: run `just init {ver} <path>` to symlink raw data")

    # 4. placeholder config files
    placeholder_dir = Path(f"data/{ver}/single_process")
    for fname in ("fit_range.json", "bkg_form.json", "ec_energy.json"):
        fpath = placeholder_dir / fname
        if not fpath.exists():
            with open(fpath, "w") as f:
                json.dump({}, f)
            print(f"  created placeholder {fpath}")

    # 5. reader skeleton
    reader_dir = Path(f"lib_reader/src/lib_reader/reader{ver}")
    if not reader_dir.exists():
        reader_dir.mkdir(parents=True)
        with open(reader_dir / "__init__.py", "w") as f:
            f.write(_READER_INIT.format(ver=ver))
        print(f"  created reader skeleton at {reader_dir}")
        print(f"  NOTE: implement reader functions in {reader_dir}/read.py")
    else:
        print(f"  reader {reader_dir} already exists, skipping")

    # 6. QA thresholds
    qa_path = placeholder_dir / "qa_thresholds.json"
    if not qa_path.exists():
        with open(qa_path, "w") as f:
            json.dump({
                "tb": {"redchi_warn": 2.0, "redchi_fail": 5.0},
                "ec": {"redchi_warn": 5.0, "redchi_fail": 50.0}
            }, f, indent=2)
        print(f"  created {qa_path}")

    print(f"\n  DONE. Next steps:")
    print(f"  1. Add raw data: just init {ver} <path-to-raw-data>")
    print(f"  2. Configure fit_range/bkg_form/ec_energy in data/{ver}/single_process/")
    print(f"  3. Implement reader in lib_reader/src/lib_reader/reader{ver}/")
    print(f"  4. Run: just check {ver}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("ver", help="payload version, e.g. 12B")
    parser.add_argument("--data-dir", help="raw data directory to symlink", default=None)
    args = parser.parse_args()
    scaffold_version(args.ver, args.data_dir)


if __name__ == "__main__":
    main()
