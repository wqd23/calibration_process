from .. import util_lib as util
import os
import numpy as np
import matplotlib.pyplot as plt
from ..operation import TB_operation_03B
# plt.ion()
from .__init__ import tb_op

data_all = tb_op.load_data()
tb_op.temp_bias_fit(data_all)