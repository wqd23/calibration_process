# 数据准备与目录约定

## 目录结构

```
data/{ver}/
├── raw_data -> /path/to/actual/data       # 软链，指向原始标定数据
├── single_process/
│   ├── fit_range.json                      # 每通道拟合区间配置
│   ├── bkg_form.json                       # 每通道本底模型配置
│   ├── ec_energy.json                      # 文件名 → 能量映射
│   ├── qa_thresholds.json                  # QA 阈值配置
│   ├── TB_fit_result/                      # TB 单谱拟合产物（pickle）
│   ├── EC_fit_result/                      # EC 单谱拟合产物（pickle）
│   └── single_fit_fig/                     # 拟合图（png）
├── tb_logs/                                # TB 二维面拟合产物
└── ec_logs/                                # E-C 关系拟合产物
```

`just init {ver} {path}` 会自动创建目录结构并软链 `raw_data`。`just check {ver}` 会验证所有路径是否就绪。

## 配置文件

### fit_range.json

每通道的拟合区间。格式为 4 元素数组（对应 ch0–ch3），每个元素是 `[low, high]`（ADC 道址范围），`null` 表示跳过该通道。

参考：`data/03B/single_process/fit_range.json`

### bkg_form.json

每通道的本底模型。可选值：`"lin"`（线性）、`"quad"`（二次）、`"exp"`（指数）、`"gaus"`（高斯）、`null`（无本底）。

参考：`data/03B/single_process/bkg_form.json`

### ec_energy.json

文件名到光子能量（keV）的映射。键为 pickle 文件名（不含扩展名），值为能量值（标量或每通道数组）。

参考：`data/07/single_process/ec_energy.json`

### qa_thresholds.json

拟合质量阈值。格式：
```json
{
  "tb": {"redchi_warn": 2.0, "redchi_fail": 5.0},
  "ec": {"redchi_warn": 5.0, "redchi_fail": 50.0}
}
```

文件不存在时使用默认值。新版本接入后建议先跑一遍数据再根据实际分布调整阈值。

## 新版本接入流程

```bash
just new-payload {ver}          # 生成 config 条目 + 目录 + reader 骨架
just init {ver} {data-path}     # 软链数据
# 编辑 data/{ver}/single_process/ 下的 fit_range / bkg_form / ec_energy
# 实现 lib_reader/src/lib_reader/reader{ver}/
just check {ver}                # 验证配置
just all {ver}                  # 跑数据
just tbfit {ver}                # TB 二维面拟合
just ecfit {ver}                # E-C 关系拟合
```

配置文件的具体格式可参考已有版本（如 `data/03B/` 或 `data/07/`）。
