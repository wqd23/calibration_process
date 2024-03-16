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
