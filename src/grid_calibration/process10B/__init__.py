from ..operation import TB_operation_10B, EC_operation_10B
tb_op = TB_operation_10B(**{
    "path": "data/10B/raw_data/tb_data",
    "fit_range": "data/10B/single_process/fit_range.json",
    "save_path": "data/10B/single_process/TB_fit_result",
    "save_fig_path": "data/10B/single_process/single_fit_fig",
    "result_path": "data/10B/tb_logs"
})
ec_op = EC_operation_10B(**{
    'tb_result_path' : 'data/10B/tb_logs/20240629150801_temp_bias_fit.json',
    'fit_range' : "data/10B/single_process/fit_range.json",
    'energy' : 'data/10B/single_process/ec_energy.json',
    'bkg_form' : 'data/10B/single_process/bkg_form.json',
    'x_path' : "data/10B/raw_data/x_data",
    'src_path' : 'data/10B/raw_data/src_data',
    'save_path' : "data/10B/single_process/EC_fit_result",
    'save_fig_path' : "data/10B/single_process/single_fit_fig",
    'result_path' : "data/10B/ec_logs"
})