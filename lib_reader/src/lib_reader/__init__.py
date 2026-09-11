from pathlib import Path
import os
from .reader05.frame_adapter import (
    src_read03b as _src_read03b,
    single_read03b as _single_read03b,
    single_read05b_normal as _single_read05b_normal,
    single_read05b_xray as _single_read05b_xray,
)
from .reader07.read import single_read07 as _single_read07
from .reader04.read import single_read04 as _single_read04
from .reader10.read import single_read10 as single_read10
from .reader11.read import single_read11 as single_read11
from .reader12.read import single_read12 as single_read12
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

# Pipeline entry points: L1 parquet cache of the final processed (sci, tel)
# output, stored as channel-stacked parquet tables.
src_read03b = with_l1_cache_processed(ver="03B", reader="03b-src", kind="single")(_src_read03b)
single_read03b = with_l1_cache_processed(ver="03B", reader="03b", kind="single")(_single_read03b)
single_read05b_normal = with_l1_cache_processed(ver="05B", reader="normal", kind="single")(_single_read05b_normal)
single_read05b_xray = with_l1_cache_processed(ver="05B", reader="xray", kind="single")(_single_read05b_xray)
single_read07 = with_l1_cache_processed(ver="07", reader="07", kind="single")(_single_read07)
single_read04 = with_l1_cache_processed(ver="04", reader="04", kind="single")(_single_read04)
single_read09 = with_l1_cache_processed(ver="09", reader="07", kind="single")(_single_read07)
