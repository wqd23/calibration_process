default:
  just --list
  @echo "Available ver: 03B 05B 07"
run any:
  python3 -m grid_calibration.{{any}}
tb ver *flags:
  python3 -m grid_calibration.process{{ver}}.temp_bias {{flags}}
tbfit ver:
  python3 -m grid_calibration.process{{ver}}.temp_fit
ec ver *flags:
  python3 -m grid_calibration.process{{ver}}.ec_preprocess {{flags}}
ecfit ver:
  python3 -m grid_calibration.process{{ver}}.ec_fit

init ver path:
  ln -s {{path}} ./data/{{ver}}/raw_data
  mkdir ./data/{{ver}}/ec_logs
  mkdir ./data/{{ver}}/tb_logs
  mkdir ./data/{{ver}}/single_process/TB_fit_result
  mkdir ./data/{{ver}}/single_process/EC_fit_result
  mkdir ./data/{{ver}}/single_process/single_fit_fig
