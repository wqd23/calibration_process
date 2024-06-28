from ..operation import EC_operation_04, TB_operation_04

tb_op = TB_operation_04(**{
    'path': "data/04/raw_data/20210501_tempbias_Am241_GRID04",
    "fit_range": "data/04/single_process/fit_range.json",
    "save_path": "data/04/single_process/TB_fit_result",
    "save_fig_path": "data/04/single_process/single_fit_fig",
    "result_path": "data/04/tb_logs"
})
ec_op = EC_operation_04(**{
    'tb_result_path' : 'data/04/single_process/20231120140952_temp_bias_fit.json',
    'fit_range' : "data/04/single_process/fit_range.json",
    'energy' : 'data/04/single_process/ec_energy.json',
    'bkg_form' : 'data/04/single_process/bkg_form.json',
    'x_path' : "data/04/raw_data/20210429_Xray_GRID04/data_GRID",
    'src_path' : 'data/04/raw_data/20210504_source_GRID04',
    'save_path' : "data/04/single_process/EC_fit_result",
    'save_fig_path' : "data/04/single_process/single_fit_fig",
    'result_path' : "data/04/ec_logs"
})