# GRID 标定数据处理

本仓库用于 GRID 载荷标定数据的处理、拟合与质量检验。

## 快速上手

**前置依赖：** [uv](https://docs.astral.sh/uv/)（Python 环境管理）和 [just](https://github.com/casey/just)（命令运行器）。Python 版本由 uv 按 `.python-version` 自动管理，无需手动安装。

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh   # 安装 uv
cargo install just                                  # 安装 just（或 brew install just）

git clone <repo-url> && cd calibration_process
uv sync                                             # 安装依赖到 .venv
```

**接入某个版本并跑数据：**

```bash
just init {ver} {data-path}   # 软链数据 + 建目录（仅首次接入）
just check {ver}              # 验证部署完整性
just all {ver}                # 处理全部 TB + EC 数据
just tbfit {ver}              # TB 二维面拟合
just ecfit {ver}              # E-C 关系拟合
```

`{ver}` 为载荷版本号（如 `12B`），完整列表可用 `just`（无参数）查看。

## 部署（Cat 服务器）

修订日期：2026-08-29。部署目标为本机服务器 **Cat**（`hostname: cat`）。以下命令在首次接入时运行即可一次性初始化全部版本，原始数据会软链到 `data/{ver}/raw_data`，不复制数据。

```bash
# 一次性初始化 12B 及之前全部载荷（路径内嵌，直接复制粘贴运行）
just init 03B /home/wqd/cali_data/03B
just init 04  /home/wqd/cali_data/04
just init 05B /home/wqd/cali_data/05B
just init 07  /home/wqd/cali_data/07
just init 09  /home/wqd/cali_data/09
just init 10B /home/wqd/cali_data/10B
just init 11B /home/wqd/cali_data/11B
just init 12B "/home/wqd/cali_data/12B、13B/data"   # 目录名含顿号，需加引号
```

`just init {ver} {path}` 把 `{path}` 软链为 `data/{ver}/raw_data`，并创建输出目录（`tb_logs`、`ec_logs`、`single_process/{TB_fit_result,EC_fit_result,single_fit_fig}`）。各版本数据路径在 `src/calibration_process/config.json` 里以 `data/{ver}/raw_data/...` 相对引用。`just check {ver}` 会逐项验证 raw_data 软链、输入配置与输出目录是否就绪（全 OK 输出 `[{ver}] READY`，缺失可用 `--fix` 自动创建输出目录）。已知的运维坑（如 cachier 缓存锁）见根目录 [AGENTS.md](AGENTS.md)。

## 文档

仓库里的三类文档各自有明确面向的对象，不要混放：

- **[AGENTS.md](AGENTS.md)（面向 AI 助手）**：工程约定（代码风格、常用命令）、**接入新载荷的方法论与坑**、验证方式。
- **[docs/](docs/)（面向人，深入某领域）**：数据准备与目录约定、结果产物与 QA 指标、12B/13B 数据说明。清单与写作约定见 [docs/README.md](docs/README.md)。
- **本 README（面向首次接触者）**：装环境、跑数据、部署入口、命令速查。

| 文档                                        | 内容                                                           |
| ------------------------------------------- | -------------------------------------------------------------- |
| [AGENTS.md](AGENTS.md)                       | 面向 AI 助手的工程约定：命令、代码风格、接入新载荷的方法论与坑 |
| [docs/data.md](docs/data.md)                 | 数据准备与目录约定：配置文件格式、新版本接入步骤骨架           |
| [docs/results.md](docs/results.md)           | 结果产物与 QA 指标说明                                         |
| [docs/12B_13B/data.md](docs/12B_13B/data.md) | 12B/13B 类载荷的数据说明（点位对照表、各文件的特殊情况）       |

每类载荷的数据说明单独放一个目录（如 `docs/12B_13B/`），新增载荷时仿照添加。

## 命令速查

| 命令                            | 说明                                         |
| ------------------------------- | -------------------------------------------- |
| `just`                        | 查看所有命令和可用版本                       |
| `just init {ver} {path}`      | 初始化版本：软链数据 + 建目录                |
| `just check {ver}`            | 验证部署完整性（`--fix` 自动创建缺失目录） |
| `just all {ver}`              | 处理该版本全部数据                           |
| `just tb {ver} run {idx}`     | 处理单个 TB 文件                             |
| `just ec {ver} x run {idx}`   | 处理单个 EC X光机文件                        |
| `just ec {ver} src run {idx}` | 处理单个 EC 放射源文件                       |
| `just tbfit {ver}`            | TB 二维面拟合                                |
| `just ecfit {ver}`            | E-C 关系拟合                                 |
| `just new-payload {ver}`      | 为新载荷生成配置和目录骨架                   |

## 仓库结构

```
.
├── README.md
├── AGENTS.md                       # 面向 AI 助手的工程约定与接入方法论
├── docs/                           # 详细文档（清单见 docs/README.md）
│   ├── README.md                   # 文档目录说明与清单
│   ├── data.md                     # 数据准备与目录约定、接入步骤骨架
│   ├── results.md                  # 结果产物与 QA 指标
│   └── 12B_13B/                    # 12B/13B 类载荷的数据说明（每类载荷一个目录）
│       └── data.md                 # 点位对照表、各数据文件的特殊情况
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

- [X] 12B 简易结果
- [ ] 12B 完整结果（TB: 备份目录全 54 点位温偏面拟合，适用偏压 ≥27.5V 内残差 <3.5%；EC: 与其他载荷相同的 K 边拆段二次拟合（EC_low/EC_high 输出 schema 一致），锚点为 Am241/Na22/Cs137/Co60（双高斯拟合 1332 keV）+ X光机 20-100 kV（低管压点的能量按管压赋值、有已知系统偏差，仅作参考；见 docs/12B_13B/data.md））
- [ ] 13B 完整结果
- [ ] 10B, 11B 塑闪结果

## Contributor

- wqd
- 丘智勇
- Lee
- yhl
- xsq(Haruka Kujo)
