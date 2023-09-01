default:
  just --list
  @echo "Available ver: 03B 05B 07"

tb ver *flags:
  python3 -m grid_calibration.process{{ver}}.temp_bias {{flags}}
tbfit ver:
  python3 -m grid_calibration.process{{ver}}.temp_fit
ec ver *flags:
  python3 -m grid_calibration.process{{ver}}.ec_preprocess {{flags}}
ecfit ver:
  python3 -m grid_calibration.process{{ver}}.ec_fit
