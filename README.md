# grid-calibration

对原有标定代码进行重构，代码仓库中去除了原始数据与结果数据，完整结构如下
```
.
├── README.md
├── data
│   ├── 03B
│   │   ├── ec_logs
│   │   ├── single_process
│   │   └── tb_logs
│   └── 07
├── lib_reader
│   ├── README.md
│   ├── pyproject.toml
│   └── src
│       └── lib_reader
├── pyproject.toml
├── requirements-dev.lock
├── requirements.lock
└── src
    └── grid_calibration
        ├── __init__.py
        ├── file_lib.py
        ├── main.py
        ├── operation.py
        ├── process03B
        └── util_lib.py
```