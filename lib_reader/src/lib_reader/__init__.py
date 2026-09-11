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

# All readers cache internally: faithful L1 frames under data/{ver}/l1/ and the
# processed (sci, tel) output under data/{ver}/l2/.

__all__ = [
    "single_read03b", "src_read03b", "single_read04", "single_read05b_normal",
    "single_read05b_xray", "single_read07", "single_read09", "single_read10",
    "single_read11", "single_read12",
]
