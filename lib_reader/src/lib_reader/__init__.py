from pathlib import Path
import os
from .reader05.frame_adapter import (
    src_read03b as _src_read03b,
    single_read03b as _single_read03b,
    single_read05b_normal as _single_read05b_normal,
    single_read05b_xray as _single_read05b_xray,
)
from .reader07.read import single_read07, single_read04, single_read09
from .reader10.read import single_read10
from .reader11.read import single_read11
from .reader12.read import single_read12
from .l1_cache import with_l1_cache_processed


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

# Pipeline entry points for the 03B/05B readers: L2 processed cache of the
# final (sci, tel) output.  The 04/07/09 readers cache L1 frames + L2 output
# internally (see reader07/frame_adapter.py).
src_read03b = with_l1_cache_processed(ver="03B", reader="03b-src", kind="single")(_src_read03b)
single_read03b = with_l1_cache_processed(ver="03B", reader="03b", kind="single")(_single_read03b)
single_read05b_normal = with_l1_cache_processed(ver="05B", reader="normal", kind="single")(_single_read05b_normal)
single_read05b_xray = with_l1_cache_processed(ver="05B", reader="xray", kind="single")(_single_read05b_xray)

__all__ = [
    "single_read03b", "src_read03b", "single_read04", "single_read05b_normal",
    "single_read05b_xray", "single_read07", "single_read09", "single_read10",
    "single_read11", "single_read12",
]
