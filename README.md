# GRID 标定数据处理

本仓库用于 GRID 载荷标定数据的处理、拟合与质量检验。

## 快速上手

**前置依赖：** [uv](https://docs.astral.sh/uv/)（Python 环境管理）和 [just](https://github.com/casey/just)（命令运行器）。
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh   # 安装 uv
cargo install just                                  # 安装 just（或 brew install just）
```

**部署：**
```bash
git clone <repo-url> && cd calibration_process
uv sync                              # 安装依赖
just init {ver} {data-path}          # 软链数据 + 创建目录（ver 用 just 查看）
just check {ver}                     # 验证部署
```

**跑数据：**
```bash
just all {ver}                       # 处理全部 TB + EC 数据
just tbfit {ver}                     # TB 二维面拟合
just ecfit {ver}                     # E-C 关系拟合
```

## 文档

详细文档在 [docs/](docs/) 目录（清单与写作约定见 [docs/README.md](docs/README.md)）：

| 文档 | 内容 |
|------|------|
| [docs/deploy.md](docs/deploy.md) | 部署指南：新机器上装环境、挂数据、验证部署 |
| [docs/data.md](docs/data.md) | 数据准备与目录约定：配置文件格式、新版本载荷接入流程与方法论 |
| [docs/results.md](docs/results.md) | 结果产物与 QA 指标说明 |
| [docs/12B_13B/data.md](docs/12B_13B/data.md) | 12B/13B 类载荷的数据说明（点位对照表、各文件的特殊情况） |
| [docs/N1/data.md](docs/N1/data.md) | N1（GRIDN1）的数据说明（GAGG/CLYC 双数据集、扫描分段、坏点） |

每类载荷的数据说明单独放一个目录（如 `docs/12B_13B/`），新增载荷时仿照添加。

## 命令速查

| 命令 | 说明 |
|------|------|
| `just` | 查看所有命令和可用版本 |
| `just init {ver} {path}` | 初始化版本：软链数据 + 建目录 |
| `just check {ver}` | 验证部署完整性（`--fix` 自动创建缺失目录） |
| `just all {ver}` | 处理该版本全部数据 |
| `just tb {ver} run {idx}` | 处理单个 TB 文件 |
| `just ec {ver} x run {idx}` | 处理单个 EC X光机文件 |
| `just ec {ver} src run {idx}` | 处理单个 EC 放射源文件 |
| `just tbfit {ver}` | TB 二维面拟合 |
| `just ecfit {ver}` | E-C 关系拟合 |
| `just new-payload {ver}` | 为新载荷生成配置和目录骨架 |

## 仓库结构

```
.
├── README.md
├── docs/                           # 详细文档（清单见 docs/README.md）
│   ├── README.md                   # 文档目录说明与清单
│   ├── deploy.md                   # 部署指南
│   ├── data.md                     # 数据准备与目录约定、接入方法论
│   ├── results.md                  # 结果产物与 QA 指标
│   ├── 12B_13B/                    # 12B/13B 类载荷的数据说明（每类载荷一个目录）
│   │   └── data.md                 # 点位对照表、各数据文件的特殊情况
│   └── N1/                         # N1（GRIDN1）的数据说明
│       └── data.md
├── data/{ver}/                     # 各版本数据目录
│   ├── raw_data -> /path/to/data   # 原始数据软链
│   ├── single_process/             # 配置文件和拟合产物
│   ├── tb_logs/                    # TB 面拟合产物
│   └── ec_logs/                    # E-C 拟合产物
├── lib_reader/                     # 各版本数据读取库（workspace 子包）
├── lib_plot/                       # 绘图库（workspace 子包）
├── src/calibration_process/
│   ├── process.py                  # CLI 入口（惰性加载）
│   ├── cmd.py                      # 命令行包装
│   ├── operation.py                # TB/EC 操作类
│   ├── file_lib.py                 # 单文件读取与拟合
│   ├── fitting.py                  # 通用峰型拟合
│   ├── util_lib.py                 # 工具函数（拟合、缓存、阈值）
│   ├── check.py                    # 部署校验
│   ├── scaffold.py                 # 新版本脚手架
│   └── config.json                 # 各版本配置（模板化）
├── pyproject.toml
├── uv.lock
└── justfile
```

## TODO

- [x] 12B 简易结果
- [x] 12B 完整结果（TB: 备份目录全 54 点位温偏面拟合，适用偏压 ≥27.5V 内残差 <3.5%；EC: 与其他载荷相同的 K 边拆段二次拟合（EC_low/EC_high 输出 schema 一致），锚点为 Am241/Na22/Cs137/Co60（双高斯拟合 1332 keV）+ X光机 20-100 kV（低管压点的能量按管压赋值、有已知系统偏差，仅作参考；见 docs/12B_13B/data.md））
- [x] N1 温偏（TB: GAGG/CLYC 双数据集分别二维面拟合，残差 std 1.4-2.1%，复用 gridN_cali 的包定义与拟合区间；见 docs/N1/data.md；EC 未接入）
- [ ] 13B 完整结果
- [ ] N1 EC（放射源 260326/260327、计量院标定 260129/260202）
- [ ] 10B, 11B 塑闪结果

## Contributor

- wqd
- 丘智勇
- Lee
- yhl
- xsq(Haruka Kujo)
