from ..operation import TB_operation_04, EC_operation_04

tb_op = TB_operation_04(**{
    'path': "data/04/raw_data/20210501_tempbias_Am241_GRID04",
    "fit_range": "data/04/single_process/fit_range.json",
    "save_path": "data/04/single_process/TB_fit_result",
    "save_fig_path": "data/04/single_process/single_fit_fig",
    "result_path": "data/04/tb_logs"
})