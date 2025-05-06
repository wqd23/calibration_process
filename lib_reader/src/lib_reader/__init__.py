from cachier import cachier
from pathlib import Path
import os
from .reader05.version_lib import (
    src_read03b as original_src_read03b,
    single_read03b as original_single_read03b,
    single_read05b_normal as original_single_read05b_normal,
    single_read05b_xray as original_single_read05b_xray,
)
from .reader07.read import single_read07 as original_single_read07
from .reader04.read import single_read04 as original_single_read04
from .reader10.read import single_read10 as original_single_read10
from .reader11.read import single_read11 as original_single_read11

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
CACHE_DIR = PROJECT_ROOT / ".cache"

src_read03b = cachier(cache_dir=CACHE_DIR)(original_src_read03b)
single_read03b = cachier(cache_dir=CACHE_DIR)(original_single_read03b)
single_read05b_normal = cachier(cache_dir=CACHE_DIR)(original_single_read05b_normal)
single_read05b_xray = cachier(cache_dir=CACHE_DIR)(original_single_read05b_xray)
single_read07 = cachier(cache_dir=CACHE_DIR)(original_single_read07)
single_read04 = cachier(cache_dir=CACHE_DIR)(original_single_read04)
single_read10 = cachier(cache_dir=CACHE_DIR / "10B", separate_files=True)(
    original_single_read10
)
single_read11 = cachier(cache_dir=CACHE_DIR / "11B", separate_files=True)(
    original_single_read11
)