from ..operation import TB_operation_07B
from .. import util_lib as util
tb_op = TB_operation_07B(**{
    "path": "data/07/raw_data/北师大正样温度偏压标定数据-20211124",
    "fit_range": "data/07/single_process/fit_range.json",
    "save_path": "data/07/single_process/TB_fit_result",
    "save_fig_path": "data/07/single_process/single_fit_fig",
    "result_path": "data/07/tb_logs"
})