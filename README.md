# 仓库结构

对原有标定代码进行重构，代码仓库中去除了原始数据与结果数据，完整结构如下
```
.
├── README.md
├── data
│   ├── 03B
│   │   ├── ec_logs
│   │   ├── raw_data
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
其中`./data/03B/raw_data`指向实际03B标定数据，clone后需要手动完成。

# lib_reader

`lib_reader`为单独 python module 用于定义数据文件读取接口，对各版本载荷实现读取函数。
```python
# lib_reader/src/lib_reader/__init__.py
from .reader05.version_lib import src_read03b, single_read03b, single_read05b_normal, single_read05b_xray
from .reader07.read import single_read07
```
可以看出lib_reader 目前提供 03B、05B与07的读取函数。

# grid_calibration

此模块为实际处理模块，完成实际通过调用数据读取模块(`lib_reader`)与能谱统计、曲线拟合等工具完成工作。
下为`grid_calibration`代码结构：
```
src/grid_calibration
├── __init__.py
├── file_lib.py
├── main.py
├── operation.py
├── process03B
│   ├── __init__.py
│   ├── ec_fit.py
│   ├── ec_preprocess.py
│   ├── single.py
│   ├── temp_bias.py
│   └── temp_fit.py
├── process07
│   ├── __init__.py
│   └── single.py
└── util_lib.py
```
其中`main.py`为待完成的通用接口，暂无了解必要。

## file_lib

对于标定实验，单次实验的处理模式总是相似的，包括
1. 数据读取(实验数据及可能存在的背景数据)
2. 能谱统计
3. 全能峰拟合

因而在`file_lib.py`中定义`class File_operation_05b`用于对单次实验的数据进行处理。使用方法可以参考`src/grid_calibration/process03B/single.py`

## operation

对于完整标定数据处理，每代卫星载荷之间存在一定差距，但是在通过`file_lib`对单次数据处理进行抽象后，所需的就只是为`file_lib`提供合适的参数，并予以调用。

`operation.py`完成了这一工作。可以通过`TB_operation`与`EC_operation`对一代载荷的标定数据处理过程进行封装.
以03B为例，`TB_operation_03B`提供了以下接口
```
.files = list of all file to process
.file_config(file) = (read_config, bkg_read_config, spectrum_config, fit_config) of file
.load_data() = data_all (all temp bias data processed)
.temp_bias_fit(data_all) fit for all data
```

## 03B代码使用方法

当前03B已经实现了命令行接口调用，调用方式如下：
在仓库根目录下有以下命令可供执行
```bash
python3 -m grid_calibration.process03B.temp_bias -h      # 查看temp_bias帮助
python3 -m grid_calibration.process03B.temp_bias list    # 列出所有温度偏压实验数据文件
python3 -m grid_calibration.process03B.temp_bias run {n} # 处理{n}号文件
python3 -m grid_calibration.process03B.temp_bias run all # 处理所有文件

python3 -m grid_calibration.process03B.temp_fit          # 使用temp_bias处理结果完成温度偏压实验处理
```
```bash
python3 -m grid_calibration.process03B.ec_preprocess -h         # 查看ec_preprocess帮助
python3 -m grid_calibration.process03B.ec_preprocess x -h       # 查看x光机实验数据处理帮助
python3 -m grid_calibration.process03B.ec_preprocess src -h     # 查看放射源实验数据处理帮助
python3 -m grid_calibration.process03B.ec_preprocess x/src ...  # 与temp_bias有同样参数

python3 -m grid_calibration.process03B.ec_fit                   # 使用ec_preprocess 结果完成能量响应
```

# TODO
- [ ] 命令行接口重写，批量生成
- [ ] operation 拆分
