# -*- coding:utf-8 -*-
"""Generate the git reader golden (offline, one-time).

For each sample it commits a small truncated real raw sample under
``tests/golden/reader/<sample>/raw/`` plus the frozen ``(sci, tel)`` output
(``expected.npz`` + ``structure.json``).  At test time
``tests/test_reader_golden.py`` runs the *new* unified reader on the committed
sample and asserts the output equals the frozen output, with no ``raw_data`` /
``.oracle`` dependency.

The committed expected output is the legacy-validated ``(sci, tel)`` (the same
values the frozen oracle produces, byte-for-byte).  Re-running this script
re-freezes the *current* reader output: only do it intentionally, after the
oracle regression (``tests/test_pipeline_run.py``) has confirmed the reader is
still legacy-identical, otherwise the legacy baseline would be overwritten by
current code.

For 03B/05B the HK / timeline / config sibling files (derived from the rundata
name) are committed alongside the truncated rundata.  The 03B feature sample is
a post-time-cut window so the per-file ``cutFileRef`` cut leaves events.

Usage (data-dependent, run once):

    uv run python scripts/gen_reader_golden.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "tests" / "golden" / "reader"

from lib_reader.reader07.frame_adapter import single_read07, single_read04  # noqa: E402
from lib_reader.reader05.frame_adapter import (  # noqa: E402
    single_read05b_normal,
    single_read05b_xray,
    single_read03b,
    src_read03b,
)

D05B = "data/05B/raw_data/test"
D05X = "data/05B/raw_data/X光机实验-天格"
D03S = "data/03B/raw_data/20210504_source_03B"
D03X = "data/03B/raw_data/20210429_Xray_03B"

# sample -> (legacy reader, rundata name, files[(src, dest, offset, size)], args)
#   size None -> copy to EOF; offset is a byte offset (window for 03B_xray).
SAMPLES = {
    "07_hex": (
        single_read07,
        "bnu_Am241_10cm_15min_220423180307_COM3-Data.txt",
        [("data/07/raw_data/北师大上胶后补标定/bnu_Am241_10cm_15min_220423180307_COM3-Data.txt",
          "bnu_Am241_10cm_15min_220423180307_COM3-Data.txt", 0, 150000)],
        (),
    ),
    "04_hex": (
        single_read04,
        "210501124621_COM7_tb_-20C_27p0V_4m_5cm-Data.txt",
        [("data/04/raw_data/20210501_tempbias_Am241_GRID04/210501124621_COM7_tb_-20C_27p0V_4m_5cm-Data.txt",
          "210501124621_COM7_tb_-20C_27p0V_4m_5cm-Data.txt", 0, 150000)],
        (),
    ),
    "05B_normal": (
        single_read05b_normal,
        "TB_0C_275_rundata2022-04-11-02-29-33.dat",
        [(f"{D05B}/TB_0C_275_rundata2022-04-11-02-29-33.dat",
          "TB_0C_275_rundata2022-04-11-02-29-33.dat", 0, 60000),
         (f"{D05B}/TB_0C_275_HK2022-04-11-02-29-33.dat",
          "TB_0C_275_HK2022-04-11-02-29-33.dat", 0, 30000),
         (f"{D05B}/TB_0C_275_TimeLine2022-04-11-02-29-33.dat",
          "TB_0C_275_TimeLine2022-04-11-02-29-33.dat", 0, 30000),
         (f"{D05B}/TB_0C_275_SciConfig2022-04-11-02-29-33.json",
          "TB_0C_275_SciConfig2022-04-11-02-29-33.json", 0, None)],
        (),
    ),
    "05B_xray": (
        single_read05b_xray,
        "XM_100_rundata2022-04-25-13-46-36.dat",
        [(f"{D05X}/XM_100_rundata2022-04-25-13-46-36.dat",
          "XM_100_rundata2022-04-25-13-46-36.dat", 0, 60000),
         (f"{D05X}/XM_100_HK2022-04-25-13-46-36.dat",
          "XM_100_HK2022-04-25-13-46-36.dat", 0, 30000),
         (f"{D05X}/XM_100_TimeLine2022-04-25-13-46-36.dat",
          "XM_100_TimeLine2022-04-25-13-46-36.dat", 0, 30000),
         (f"{D05X}/XM_100_SciConfig2022-04-25-13-46-36.json",
          "XM_100_SciConfig2022-04-25-13-46-36.json", 0, None)],
        ("",),
    ),
    "03B_src": (
        src_read03b,
        "src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat",
        [(f"{D03S}/src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat",
          "src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat", 0, 250000),
         (f"{D03S}/src_Cs137_12m_10cm_HK2021-05-05-12-12-27.dat",
          "src_Cs137_12m_10cm_HK2021-05-05-12-12-27.dat", 0, 20000),
         (f"{D03S}/src_Cs137_12m_10cm_TimeLine2021-05-05-12-12-27.dat",
          "src_Cs137_12m_10cm_TimeLine2021-05-05-12-12-27.dat", 0, 20000),
         (f"{D03S}/src_Cs137_12m_10cm_scienceConfig2021-05-05-12-12-27.json",
          "src_Cs137_12m_10cm_scienceConfig2021-05-05-12-12-27.json", 0, None)],
        ("",),
    ),
    "03B_xray": (
        single_read03b,
        "jly_18p0_ch0_30s_rundata2021-04-29-15-21-50.dat",
        [(f"{D03X}/jly_18p0_ch0_30s_rundata2021-04-29-15-21-50.dat",
          "jly_18p0_ch0_30s_rundata2021-04-29-15-21-50.dat", 2600000, 120000),
         (f"{D03X}/jly_18p0_ch0_30s_HK2021-04-29-15-21-50.dat",
          "jly_18p0_ch0_30s_HK2021-04-29-15-21-50.dat", 0, None),
         (f"{D03X}/jly_18p0_ch0_30s_TimeLine2021-04-29-15-21-50.dat",
          "jly_18p0_ch0_30s_TimeLine2021-04-29-15-21-50.dat", 0, None),
         (f"{D03X}/jly_18p0_ch0_30s_scienceConfig2021-04-29-15-21-50.json",
          "jly_18p0_ch0_30s_scienceConfig2021-04-29-15-21-50.json", 0, None)],
        ("",),
    ),
}


def flatten(sci, tel):
    flat = {}
    structure = {}
    for section, d in (("sci", sci), ("tel", tel)):
        structure[section] = {}
        for k, v in d.items():
            if isinstance(v, list):
                if len(v) == 0:
                    structure[section][k] = "empty"
                else:
                    structure[section][k] = "list4"
                    for i, arr in enumerate(v):
                        flat[f"{section}.{k}.{i}"] = np.asarray(arr)
            else:
                structure[section][k] = "flat"
                flat[f"{section}.{k}"] = np.asarray(v)
    return flat, structure


def main():
    for sample, (reader, rundata_name, files, args) in SAMPLES.items():
        raw_dir = OUT / sample / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)

        for src, name, offset, size in files:
            data = (ROOT / src).read_bytes()
            chunk = data[offset:offset + size] if size is not None else data[offset:]
            (raw_dir / name).write_bytes(chunk)

        sci, tel = reader(str(raw_dir / rundata_name), *args)
        flat, structure = flatten(sci, tel)
        np.savez(OUT / sample / "expected.npz", **flat)
        (OUT / sample / "structure.json").write_text(json.dumps(structure, indent=2))

        n_sci = sum(len(np.asarray(v)) for v in sci["amp"])
        n_tel = len(np.asarray(tel["timestamp"]))
        print(f"{sample}: sci events={n_sci}, tel records={n_tel}")


if __name__ == "__main__":
    main()
