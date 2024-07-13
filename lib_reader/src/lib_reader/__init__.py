from cachier import cachier
from pathlib import Path

from .reader05.version_lib import (
    src_read03b as original_src_read03b,
    single_read03b as original_single_read03b,
    single_read05b_normal as original_single_read05b_normal,
    single_read05b_xray as original_single_read05b_xray,
)
from .reader07.read import single_read07 as original_single_read07
from .reader04.read import single_read04 as original_single_read04
from .reader10.read import single_read10 as original_single_read10

CACHE_DIR = Path("./.cache")

src_read03b = cachier(cache_dir=CACHE_DIR)(original_src_read03b)
single_read03b = cachier(cache_dir=CACHE_DIR)(original_single_read03b)
single_read05b_normal = cachier(cache_dir=CACHE_DIR)(original_single_read05b_normal)
single_read05b_xray = cachier(cache_dir=CACHE_DIR)(original_single_read05b_xray)
single_read07 = cachier(cache_dir=CACHE_DIR)(original_single_read07)
single_read04 = cachier(cache_dir=CACHE_DIR)(original_single_read04)
single_read10 = cachier(cache_dir=CACHE_DIR / "10B", separate_files=True)(
    original_single_read10
)
