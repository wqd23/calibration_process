# -*- coding:utf-8 -*-
"""The new explicit-workflow command-line interface.

``calib discover VERSION BRANCH``  -> generate a manifest draft
``calib list     VERSION BRANCH``  -> show confirmed measurements
``calib fit-one  VERSION BRANCH ID``-> single fit one measurement
``calib fit      VERSION BRANCH``  -> single fit a whole branch
``calib global   VERSION BRANCH``  -> global fit (tb / ec)
``calib all      VERSION``         -> full formal workflow
``calib config   VERSION BRANCH ID``-> show resolved config for one measurement
"""

from __future__ import annotations

import argparse
import sys

from . import pipeline


def cmd_discover(args) -> int:
    n = pipeline.discover(args.version, args.branch)
    print(f"discovered {n} measurements for {args.version}/{args.branch}")
    return 0


def cmd_list(args) -> int:
    for i, bid in enumerate(pipeline.list_measurements(args.version, args.branch)):
        print(i, bid)
    return 0


def cmd_fit_one(args) -> int:
    pipeline.single_fit(args.version, args.branch, args.id, nocache=args.nocache,
                        out=args.out)
    return 0


def cmd_fit(args) -> int:
    n = pipeline.fit_branch(args.version, args.branch, nocache=args.nocache, out=args.out)
    print(f"single fit {n} measurements for {args.version}/{args.branch}")
    return 0


def cmd_global(args) -> int:
    if args.branch == "tb":
        pipeline.global_tb(args.version, out=args.out)
    elif args.branch in ("ec", "ec-src", "ec-xray", "ec_source", "ec_xray"):
        pipeline.global_ec(args.version, out=args.out)
    else:
        raise SystemExit(f"unknown global branch {args.branch!r}")
    return 0


def cmd_all(args) -> int:
    pipeline.all_version(args.version, nocache=args.nocache, out=args.out)
    return 0


def cmd_config(args) -> int:
    from .workflows import common as stages
    rt = pipeline.load_rt(args.version)
    branch = pipeline.branch_for_manifest(args.branch)
    from . import manifest as man
    manifest = man.load_manifest(pipeline._manifest_path(args.version, branch))
    for m in man.filtered_measurements(manifest):
        if m.id == args.id:
            fc = stages.single_run_spec(rt, branch, m)
            print(f"version: {args.version}")
            print(f"branch:  {branch}")
            print(f"id:      {m.id}")
            print(f"read:    {fc.read_config}")
            print(f"bkg:     {fc.bkg_read_config}")
            print(f"spectrum:{fc.spectrum_config}")
            print(f"fit:     {fc.fit_config}")
            return 0
    raise SystemExit(f"measurement {args.id!r} not in {args.version}/{branch} manifest")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="calib")
    sub = p.add_subparsers(dest="cmd", required=True)

    def add(name, fn, help_, extra=None):
        sp = sub.add_parser(name, help=help_)
        sp.add_argument("version")
        sp.add_argument("branch")
        if extra:
            extra(sp)
        sp.set_defaults(func=fn)
        return sp

    add("discover", cmd_discover, "scan dirs -> manifest draft")
    add("list", cmd_list, "list confirmed measurements")
    sp = sub.add_parser("fit-one", help="single fit one measurement")
    sp.add_argument("version")
    sp.add_argument("branch")
    sp.add_argument("id")
    sp.add_argument("--nocache", action="store_true")
    sp.add_argument("-o", "--out")
    sp.set_defaults(func=cmd_fit_one)
    sp = sub.add_parser("fit", help="single fit a whole branch")
    sp.add_argument("version")
    sp.add_argument("branch")
    sp.add_argument("--nocache", action="store_true")
    sp.add_argument("-o", "--out")
    sp.set_defaults(func=cmd_fit)
    sp = sub.add_parser("global", help="global fit")
    sp.add_argument("version")
    sp.add_argument("branch")
    sp.add_argument("-o", "--out")
    sp.set_defaults(func=cmd_global)
    sp = sub.add_parser("all", help="full formal workflow (no discover/preview)")
    sp.add_argument("version")
    sp.add_argument("--nocache", action="store_true")
    sp.add_argument("-o", "--out")
    sp.set_defaults(func=cmd_all)
    sp = sub.add_parser("config", help="show resolved config")
    sp.add_argument("version")
    sp.add_argument("branch")
    sp.add_argument("id")
    sp.set_defaults(func=cmd_config)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
