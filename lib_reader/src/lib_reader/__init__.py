from pathlib import Path
import os
from .reader05.frame_adapter import (
    src_read03b,
    single_read03b,
    single_read05b_normal,
    single_read05b_xray,
)
from .reader07.read import single_read07, single_read04, single_read09
from .reader10.read import single_read10
from .reader11.read import single_read11
from .reader12.read import single_read12
from .readerN1.read import single_readN1
from .select import EventTable


def get_project_root() -> Path:
    current_path = Path(os.getcwd())
    # 逐级向上查找包含特定标记的目录（例如.git、.project_root或src）
    while True:
        if (current_path / "justfile").exists() or (current_path / ".gitignore").exists():
            return current_path
        # 到达文件系统根目录仍未找到，抛出异常或返回当前路径
        if current_path == current_path.parent:
            return current_path  # 或 raise FileNotFoundError
        current_path = current_path.parent


PROJECT_ROOT = get_project_root()

# Registry of the reader handlers, keyed by the ``ending`` strings used by the
# payload config (see calibration_process/configs/{ver}/reader.yaml).  All
# readers cache internally: faithful L1 frames under data/{ver}/l1/ and the
# processed (sci, tel) output under data/{ver}/l2/.
READERS = {
    "normal": single_read05b_normal,
    "xray": single_read05b_xray,
    "03b": single_read03b,
    "03b-src": src_read03b,
    "04": single_read04,
    "07": single_read07,
    "09": single_read09,
    "10b": single_read10,
    "11b": single_read11,
    "12b": single_read12,
    "n1": single_readN1,
    "n1wf": single_readN1,
}


def read_frames(path, ver, kind="sci", **kwargs):
    """Engine-dispatched faithful L1 decode (one row per particle / sample).

    Returns the raw frame columns for ``kind`` (``sci`` / ``hk`` / ``tl``).
    The per-version payload readers know which frames exist; ``kwargs`` carries
    the reader-specific switches (e.g. 03B ``feature_mode``, or the B ``ending``
    for HK/timeline files).
    """
    if ver in ("10B", "11B", "12B"):
        from . import reader10, reader11, reader12

        mod = {"10B": reader10.read, "11B": reader11.read, "12B": reader12.read}[ver]
        if kind == "sci":
            return mod.readSci(path, **kwargs)
        if kind == "hk":
            return mod.readHK(path)
        raise ValueError(f"{ver}: unknown L1 kind {kind!r}")
    if ver in ("04", "07", "09"):
        from .reader07 import frame_adapter as fa
        from .frame_io import load_hex_text

        buf = load_hex_text(path)
        if kind == "sci":
            return fa._read_sci_frames(buf)
        if kind == "tl":
            return fa._read_tel_frames(buf)
        raise ValueError(f"{ver}: unknown L1 kind {kind!r}")
    if ver == "05B":
        from .reader05 import frame_adapter as fa

        if kind == "sci":
            return fa._decode_sci_l1(path, feature_mode=False, no_udp=True)[0]
        if kind == "hk":
            return fa._decode_hk_l1(path, kwargs.get("ending", "normal"))
        if kind == "tl":
            return fa._decode_tl_l1(path, kwargs.get("ending", "normal"))
        raise ValueError(f"{ver}: unknown L1 kind {kind!r}")
    if ver == "03B":
        from .reader05 import frame_adapter as fa

        if kind == "sci":
            return fa._decode_sci_l1(path, kwargs.get("feature_mode", False), no_udp=False)[0]
        if kind == "hk":
            return fa._decode_hk_l1(path, "03b")
        if kind == "tl":
            return fa._decode_tl_l1(path, "03b")
        raise ValueError(f"{ver}: unknown L1 kind {kind!r}")
    if ver.startswith("N1"):
        from .readerN1 import read as n1

        if kind == "sci":
            return n1._readSci_impl(path, kwargs.get("mode", "ft"))
        if kind == "hk":
            return n1._readHK_impl(path)
        raise ValueError(f"{ver}: unknown L1 kind {kind!r}")
    raise ValueError(f"unknown version {ver!r}")


__all__ = [
    "single_read03b", "src_read03b", "single_read04", "single_read05b_normal",
    "single_read05b_xray", "single_read07", "single_read09", "single_read10",
    "single_read11", "single_read12", "single_readN1", "READERS", "read_frames",
    "EventTable",
]
