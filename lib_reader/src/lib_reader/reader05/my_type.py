from typing import List, Tuple, Dict, Generator, Any, Optional, Union, Callable
from numpy.typing import NDArray
import numpy as np
from collections import namedtuple

# 1 dimension float data， [float]
Float1D = NDArray[np.float64]
# 4 channel float data, [float]*4
Float_4channel = List[float]
# 4 channel float array, [Float1D]*4
Float_array_4channel = List[Float1D]