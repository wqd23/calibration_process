py := ".venv/bin/python"

default:
  just --list

# --- explicit workflow (calib) thin wrappers --------------------------------
calib *args:
  {{py}} -m calibration_process.cli {{args}}

# full formal workflow for one version (no discover/preview)
all ver:
  {{py}} -m calibration_process.cli all {{ver}}

# single-fit one measurement / whole branch
fit-one ver branch id:
  {{py}} -m calibration_process.cli fit-one {{ver}} {{branch}} {{id}}
fit ver branch:
  {{py}} -m calibration_process.cli fit {{ver}} {{branch}}

# global fit (tb / ec)
global ver branch:
  {{py}} -m calibration_process.cli global {{ver}} {{branch}}

# discovery / manifest / config introspection
discover ver branch:
  {{py}} -m calibration_process.cli discover {{ver}} {{branch}}
list ver branch:
  {{py}} -m calibration_process.cli list {{ver}} {{branch}}
config ver branch id:
  {{py}} -m calibration_process.cli config {{ver}} {{branch}} {{id}}

# deployment check / new-payload scaffold
check ver *flags:
  {{py}} -m calibration_process.cli check {{ver}} {{flags}}
new-payload ver *flags:
  {{py}} -m calibration_process.cli scaffold {{ver}} {{flags}}

# differential regression: new vs frozen legacy oracle
compare ver:
  {{py}} scripts/compare_full.py {{ver}}
compare-10b:
  {{py}} scripts/compare_10b.py

# freeze a legacy oracle snapshot: copy data/{{ver}} outputs into .oracle/{{ver}}
oracle ver:
  @mkdir -p .oracle/{{ver}}
  @cp -r data/{{ver}}/single_process/TB_fit_result .oracle/{{ver}}/TB_fit_result
  @cp -r data/{{ver}}/single_process/EC_fit_result .oracle/{{ver}}/EC_fit_result
  @cp -r data/{{ver}}/single_process/single_fit_fig .oracle/{{ver}}/single_fit_fig
  @cp -r data/{{ver}}/tb_logs .oracle/{{ver}}/tb_logs
  @cp -r data/{{ver}}/ec_logs .oracle/{{ver}}/ec_logs

# initialize a version data dir: link raw data + create output dirs
init ver path:
  -ln -s {{path}} ./data/{{ver}}/raw_data
  mkdir -p ./data/{{ver}}/single_process/TB_fit_result
  mkdir -p ./data/{{ver}}/single_process/EC_fit_result
  mkdir -p ./data/{{ver}}/single_process/single_fit_fig
  mkdir -p ./data/{{ver}}/tb_logs
  mkdir -p ./data/{{ver}}/ec_logs

cover:
  {{py}} -m coverage run --source src/ -m pytest tests/ && {{py}} -m coverage report -m
