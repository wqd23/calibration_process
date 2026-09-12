# -*- coding:utf-8 -*-
"""Merge per-source temperature-bias fit results by channel.

Each input JSON is the ``global_tb`` output: a four-slot list (index = channel)
whose entries are ``null`` or a ``{G0, k, V0, b, c, ...}`` dict.  A source that
only covers some channels (e.g. Am-241 on ch1/2, Na-22 on ch0/3) contributes
just those::

    python scripts/merge_tb.py -o merged.json \\
        --from am241_tb.json:1,2 --from na22_tb.json:0,3

The result is a four-slot list suitable for ``payload.ec.tb_ref_path``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

_SLOTS = 4


def _parse_source(spec: str):
    path, sep, channels = spec.rpartition(":")
    if not sep or not path:
        raise ValueError(f"source must be 'FILE:ch1,ch2', got {spec!r}")
    try:
        chans = [int(c) for c in channels.split(",") if c != ""]
    except ValueError as e:
        raise ValueError(f"bad channel list in {spec!r}: {e}") from e
    return path, chans


def merge_tb(sources):
    """Merge ``[(path, channels), ...]`` into a four-slot list."""
    out = [None] * _SLOTS
    for path, channels in sources:
        data = json.loads(Path(path).read_text())
        for ch in channels:
            if not 0 <= ch < _SLOTS:
                raise ValueError(f"channel {ch} out of range 0..{_SLOTS - 1}")
            if out[ch] is not None:
                raise ValueError(f"channel {ch} provided by more than one source")
            out[ch] = data[ch]
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--out", required=True, help="merged JSON path")
    p.add_argument("--from", dest="sources", action="append", required=True,
                   metavar="FILE:ch1,ch2", help="a source file and its channels")
    args = p.parse_args(argv)

    merged = merge_tb([_parse_source(s) for s in args.sources])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(merged, ensure_ascii=False))

    filled = [i for i, v in enumerate(merged) if v is not None]
    print(f"wrote {out} (channels {filled})")
    if len(filled) != _SLOTS:
        print(f"WARNING: channels {sorted(set(range(_SLOTS)) - set(filled))} are null")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
