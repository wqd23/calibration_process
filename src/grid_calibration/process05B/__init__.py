from ..operation import TB_operation_05B, EC_operation_05B
tb_op = TB_operation_05B(**{
    'path': "data/05B/raw_data/温度偏压实验",
    "fit_range": "data/05B/single_process/fit_range.json",
    "save_path": "data/05B/single_process/TB_fit_result",
    "save_fig_path": "data/05B/single_process/single_fit_fig",
    "result_path": "data/05B/tb_logs"
})
ec_op = EC_operation_05B(**{
    'tb_result_path' : 'data/05B/single_process/20230903214245_temp_bias_fit.json',
    'fit_range' : "data/05B/single_process/fit_range.json",
    "time_cut": "data/05B/single_process/time_cut.json",
    'energy' : 'data/05B/single_process/ec_energy.json',
    'bkg_form' : 'data/05B/single_process/bkg_form.json',
    'x_path' : "data/05B/raw_data/X光机实验-天格_reset",
    'src_path' : 'data/05B/raw_data/放射源实验',
    'save_path' : "data/05B/single_process/EC_fit_result",
    'save_fig_path' : "data/05B/single_process/single_fit_fig",
    'result_path' : "data/05B/ec_logs",
    "x_config": "data/05B/raw_data/X光机实验-天格_reset/config.json"
})