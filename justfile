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

# --- explicit-workflow (calib) thin wrappers --------------------------------
calib *args:
  python3 -m calibration_process.cli {{args}}

# new explicit workflow: full formal run for one version
new-all ver:
  python3 -m calibration_process.cli all {{ver}}

# freeze a legacy oracle snapshot: copy data/{{ver}} outputs into .oracle/{{ver}}
oracle ver:
  @mkdir -p .oracle/{{ver}}
  @cp -r data/{{ver}}/single_process/TB_fit_result .oracle/{{ver}}/TB_fit_result
  @cp -r data/{{ver}}/single_process/EC_fit_result .oracle/{{ver}}/EC_fit_result
  @cp -r data/{{ver}}/single_process/single_fit_fig .oracle/{{ver}}/single_fit_fig
  @cp -r data/{{ver}}/tb_logs .oracle/{{ver}}/tb_logs
  @cp -r data/{{ver}}/ec_logs .oracle/{{ver}}/ec_logs

# differential comparison: legacy oracle vs new output
compare ver:
  python3 scripts/compare_full.py {{ver}}

init ver path:
  -ln -s {{path}} ./data/{{ver}}/raw_data
  mkdir ./data/{{ver}}/ec_logs
  mkdir ./data/{{ver}}/tb_logs
  mkdir -p ./data/{{ver}}/single_process/TB_fit_result
  mkdir -p ./data/{{ver}}/single_process/EC_fit_result
  mkdir -p ./data/{{ver}}/single_process/single_fit_fig

