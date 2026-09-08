# -*- coding:utf-8 -*-
"""Explicit workflow for version 11B.

- TB: files gathered by globbing four temperature-bias subdirectories for
  ``*_observe*.dat``, excluding any path containing ``on``/``off``, removing a
  hardcoded set of points and dropping ``_50_Cs_2``; reader ``11b``, bin_width 6.
  The temp-bias global fit uses the **lmfit** fit (which needs temp/bias errors).
- EC source: 4 hardcoded sources (Ba133 has no background).
- EC X-ray: per-channel ``observe_*_ch{n}.dat`` grouped by tube energy (name
  token 2), dropping energy ``20``; ``fixed`` [ch1,ch2,ch0,ch0] rotation.
- EC is physically 3 channels (channel_count=3) padded to 4 downstream.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

TB_SUBDIRS = [
    "温度-偏压实验-20~-10℃",
    "温度偏压0～10摄氏度",
    "温度偏压20～30摄氏度",
    "温度偏压40～50摄氏度",
]
TB_REMOVE = [
    "温度偏压0～10摄氏度/247_0_Cs137_27.5_observe_1.dat",
    "温度偏压0～10摄氏度/004_10_Cs137_26.5_1_observe.dat",
    "温度偏压40～50摄氏度/039_Cs_40_26.5_observe.dat",
    "温度-偏压实验-20~-10℃/232_Cs_-10_observe.dat",
    "温度-偏压实验-20~-10℃/233_Cs_-10_26.5_observe.dat",
    "温度偏压0～10摄氏度/013_10_Cs137_29_observe .dat",  # trailing space
]
SRC_LIST = [
    "182_Ba133_20min_observe.dat",
    "187_Cs137_ch012_3min_observe.dat",
    "189_Eu152_2min_observe.dat",
    "191_Co60_2min_observe.dat",
]
SRC_BKG = [
    "",
    "188_Cs137_ch3_15min_observe.dat",
    "190_Eu152_15min_ch3_observe.dat",
    "192_Co60_15min_ch3_observe.dat",
]


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    assert version == "11B", "v11B workflow used for non-11B version"
    if branch == "tb":
        return _enumerate_tb(rt, data_dir)
    if branch == "ec_source":
        return _enumerate_ec_source(rt, data_dir)
    if branch == "ec_xray":
        return _enumerate_ec_xray(rt, data_dir)
    raise ValueError(f"unknown branch {branch!r}")


def _enumerate_tb(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    sci_dir = rt.payload.tb.science_dir
    full = data_dir / sci_dir

    def absfile(sub, name):
        return f"{sci_dir}/{sub}/{name}"

    files = []
    for sub in TB_SUBDIRS:
        for item in Path(full / sub).glob("*_observe*.dat"):
            name = item.name
            if any(k in name for k in ("on", "off")):
                continue
            files.append(absfile(sub, name))
    for rel in TB_REMOVE:
        files.remove(f"{sci_dir}/{rel}")
    files = [f for f in files if "_50_Cs_2" not in f]
    files.sort()
    out = []
    for sci in files:
        name = os.path.basename(sci)
        # id = basename (matches `just list`), fit ranges are keyed by stem
        out.append({"id": name, "branch": "tb",
                    "science_files": [sci], "hk_files": [], "aux_files": [],
                    "metadata": {"fit_key": Path(name).stem}, "use": True})
    return out


def _enumerate_ec_source(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    src_dir = rt.payload.ec.src_path
    out = []
    for f, b in zip(SRC_LIST, SRC_BKG):
        aux = [src_dir + "/" + b] if b else []
        out.append({"id": f, "branch": "ec_source",
                    "science_files": [f"{src_dir}/{f}"], "hk_files": [],
                    "aux_files": aux,
                    "metadata": {"energy": rt.energies.get(f)}, "use": True})
    return out


def _enumerate_ec_xray(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    x_dir = rt.payload.ec.x_path
    full = data_dir / x_dir
    x_ch = [f for f in os.listdir(full) if "_ch" in f and "hk" not in f]
    x_list = sorted({f.split("_")[2] for f in x_ch})
    x_list = [e for e in x_list if e != "20"]
    out = []
    for energy in x_list:
        ch_files = []
        for i in range(4):
            matched = [f for f in x_ch if f"{energy}_ch{i}" in f]
            ch_files.append(f"{x_dir}/{matched[0]}")
        out.append({"id": energy, "branch": "ec_xray", "science_files": ch_files,
                    "hk_files": [], "aux_files": [],
                    "metadata": {"energy": rt.energies.get(energy)}, "use": True})
    return out
