from ..operation import EC_operation_07B, TB_operation_07B

tb_op = TB_operation_07B(
    **{
        "path": "data/07/raw_data/北师大正样温度偏压标定数据-20211124",
        "fit_range": "data/07/single_process/fit_range.json",
        "save_path": "data/07/single_process/TB_fit_result",
        "save_fig_path": "data/07/single_process/single_fit_fig",
        "result_path": "data/07/tb_logs",
    }
)
ec_op = EC_operation_07B(
    **{
        "tb_result_path": "data/07/single_process/20231120143925_temp_bias_fit.json",
        "fit_range": "data/07/single_process/fit_range.json",
        "energy": "data/07/single_process/ec_energy.json",
        "bkg_form": "data/07/single_process/bkg_form.json",
        "x_path": "data/07/raw_data/X光机标定实验/正样",
        "src_path": "data/07/raw_data/北师大豁免源标定数据-20211127",
        "save_path": "data/07/single_process/EC_fit_result",
        "save_fig_path": "data/07/single_process/single_fit_fig",
        "result_path": "data/07/ec_logs",
    }
)
