from ..operation import TB_operation_03B, EC_operation_03B
tb_op = TB_operation_03B(**{
    'path': "data/03B/raw_data/20210501_tempbias_03B",
    "fit_range": "data/03B/single_process/fit_range.json",
    "save_path": "data/03B/single_process/TB_fit_result",
    "save_fig_path": "data/03B/single_process/single_fit_fig",
    "result_path": "data/03B/tb_logs"
})
ec_op = EC_operation_03B(**{
    'tb_result_path' : 'data/03B/tb_logs/20230812062110_temp_bias_fit.json',
    'fit_range' : "data/03B/single_process/fit_range.json",
    'energy' : 'data/03B/single_process/ec_energy.json',
    'bkg_form' : 'data/03B/single_process/bkg_form.json',
    'x_path' : "data/03B/raw_data/20210429_Xray_03B",
    'src_path' : 'data/03B/raw_data/20210504_source_03B',
    'save_path' : "data/03B/single_process/EC_fit_result",
    'save_fig_path' : "data/03B/single_process/single_fit_fig",
    'result_path' : "data/03B/ec_logs"
})