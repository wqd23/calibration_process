default:
  just --list
  @echo "Available ver: $(python3 -c "import json; print(' '.join(json.load(open('src/calibration_process/config.json')).keys()))")"
run any:
  python3 -m calibration_process.{{any}}
tb ver *flags:
  python3 -m calibration_process.process {{ver}} tb {{flags}}
tbfit ver:
  python3 -m calibration_process.process {{ver}} tbfit
ec ver *flags:
  python3 -m calibration_process.process {{ver}} ec {{flags}}
ecfit ver:
  python3 -m calibration_process.process {{ver}} ecfit

check ver *flags:
  python3 -m calibration_process.check {{ver}} {{flags}}

new-payload ver *flags:
  python3 -m calibration_process.scaffold {{ver}} {{flags}}

test ver n:
  @just tb {{ver}} list
  @just tb {{ver}} run {{n}}
  @just ec {{ver}} x list
  @just ec {{ver}} x run {{n}}
  @just ec {{ver}} src list
  @just ec {{ver}} src run {{n}}

all ver:
  @just tb {{ver}} run all
  @just ec {{ver}} src run all
  @just ec {{ver}} x run all

cover:
  coverage run --source src/ -m pytest tests/ && coverage report -m

init ver path:
  -ln -s {{path}} ./data/{{ver}}/raw_data
  mkdir ./data/{{ver}}/ec_logs
  mkdir ./data/{{ver}}/tb_logs
  mkdir -p ./data/{{ver}}/single_process/TB_fit_result
  mkdir -p ./data/{{ver}}/single_process/EC_fit_result
  mkdir -p ./data/{{ver}}/single_process/single_fit_fig

