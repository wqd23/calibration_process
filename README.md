# 如何使用

本仓库使用[rye](https://rye-up.com/)管理python环境，并使用[just](https://github.com/casey/just)展示使用方法。但是实际使用的python版本与依赖可以在`.python-version`与`pyproject.toml`中查看并手动配置；下文中介绍的just命令也可以在`justfile`中查看其实际命令。

## 初次使用

初次使用时需要先完成python环境配置与标定数据软链接及目录创建，请注意`just init`中`ver`为载荷版本号，使用`just`命令查看可选版本，`path`为载荷标定数据文件路径。

```bash
rye sync                 # 可能耗时较长，请耐心等待
just init {ver} {path}
```

## 数据处理

对特定载荷的标定数据进行处理
```bash
just tb {ver} run all   #处理该载荷全部温度偏压实验数据
just tbfit {ver}        #拟合温度偏压响应
just ec {ver} x run all   #处理该载荷全部X光机实验数据
just ec {ver} src run all   #处理该载荷全部放射源实验数据
just ecfit {ver}        #拟合E-C关系
```
要重新处理特定一次实验数据，先通过
```bash
just tb {ver} list
just ec {ver} x list
just ec {ver} src list
```
得到文件对应编号`idx`，例如
```
0 src_Na22_20m_10cm_rundata2021-05-05-15-29-56.dat
1 src_Am241_5m_10cm_rundata2021-05-05-15-15-48.dat
2 src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat
```
再使用下列命令处理对应数据
```bash
just tb {ver} run {idx}
just ec {ver} x run {idx}
just ec {ver} src run {idx}
```

# 仓库结构

对原有标定代码进行重构，代码仓库主体结构如下
```
.
├── README.md
├── data
│   ├── {ver}
│   │   ├── ec_logs
│   │   ├── raw_data
│   │   ├── single_process
│   │   └── tb_logs
├── lib_reader
├── lib_plot
├── pyproject.toml
├── requirements-dev.lock
├── requirements.lock
└── src
    └── grid_calibration
        ├── cmd.py
        ├── file_lib.py
        ├── __init__.py
        ├── main.py
        ├── operation.py
        ├── process03B
        ├── process04
        ├── process05B
        ├── process07
        └── util_lib.py
```
其中`./data/{ver}/raw_data`指向实际标定数据。

# TODO
- [ ] 处理数据包中，没有4 channel的数据
- [ ] operation 拆分