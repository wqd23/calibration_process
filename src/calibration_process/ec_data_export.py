from . import util_lib as util
import numpy as np
from collections import namedtuple, defaultdict
from pathlib import Path

ver = ["03B", "04", "05B", "07"]

path = [Path(f"data/{v}/single_process/EC_fit_result") for v in ver]
ec_data = defaultdict(lambda: [[], [], [], []])
