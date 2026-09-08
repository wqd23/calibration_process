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
just check {ver}                     # 验证部署（新配置/manifest 层）
```

**跑数据（显式 workflow）：**
```bash
just all {ver}                       # 处理该版本全部 TB + EC 数据（单拟合 + 全局拟合）
just fit-one {ver} tb {id}           # 单 measurement 单拟合
just global {ver} tb                 # TB 二维面拟合
just global {ver} ec                 # E-C 关系拟合
just discover {ver} {branch}         # 扫描目录 -> manifest 草稿
```

新版本由**显式 workflow** 驱动：配置在 `src/calibration_process/configs/{ver}/`（YAML + strict schema），每一步都调用同一个未改动的科学内核，因此结果与历史实现一致。

## 标定的基础流程

GRID 载荷标定分成 **TB**（温度-偏压）与 **EC**（能量-道址）两条线，处理链路都是：

```text
[Raw Measurement Bundle]        # 一个 measurement = science + hk + aux + metadata（manifest 记录）
        │  discover → 人工确认 manifest
        ▼
[单文件处理]  read → TV/temp-bias 校正 → 能谱 → 峰拟合   # file_lib / util_lib（科学内核）
        ▼
[SingleFitResult]  peak / center / sigma / resolution / errors
        ▼
[TBPoint / ECPoint]  build_tb_points / build_ec_points（含 channel 过滤、能量映射）
        ▼
[版本特有 correction（若有）]   例如 12B 的 bias≥27.5V、11B 的 lmfit
        ▼
[全局拟合]  TB：峰位~（温度,偏压）二维面；EC：按 K 边拆两段二次拟合 + 分辨率拟合
        ▼
[最终产物]  JSON / NumPy / 图   # 与历史格式兼容
```

- **TB** 单谱峰拟合**不做**校正（corr = 1），温压依赖交给二维面建模。
- **EC** 单谱在生成阶段做**TV 校正**：以 TB 二维面为参考（默认 25°C/28.5V），
  把每个 EC 文件修正到该参考点（`get_spectrum` 里对 amp 乘因子）。
- EC 再分成**放射源**与**X 光机**两条分支（背景/选点不同，见 [workflows.md](docs/workflows.md)）。
- 每个版本的实际选择规则与特殊处理见 [docs/workflows.md](docs/workflows.md)。

## 增加一个新载荷

最小流程（新配置 + 新数据，**无需改科学内核**）：

```bash
calib scaffold {ver} --data-dir /path/to/data   # 生成 configs/{ver}/*.yaml + 目录 + 软链 raw_data
# 1) 编辑 configs/{ver}/payload.yaml：reader / bin_width / adc_max / 分支路径 / channel_count
# 2) 新增 reader（若包格式不同）：lib_reader/src/lib_reader/reader{ver}/，在 lib_reader/__init__.py 注册
calib discover {ver} tb                        # 扫描 -> manifest 草稿（tb / ec_source / ec_xray 各一次）
# 3) 人工确认 manifest（去重、剔除坏点、设 use / channels.use）
# 4) 填 configs/{ver}/fit_range_tb.yaml / fit_range_ec_source.yaml / fit_range_ec_xray.yaml
calib check {ver}                               # 校验配置/manifest 层与数据链接
calib all {ver}                                 # 全流程：单拟合 + TB/EC 全局拟合
```

- 若该版本的流程与已有版本**完全相同**，在 `workflows/versions/` 里复用对应文件的
  `enumerate_measurements`，并去 `workflows/registry.py` 注册即可（版本选择只发生一次）。
- 若流程**真的不同**（不同的选点/背景/分辨率/模型），在 `workflows/versions/v{ver}.py`
  里显式写出差异，并把差异点写进 [docs/workflows.md](docs/workflows.md)。
- 配置都是 **strict schema**：未知字段/非法枚举/重复 id/指向缺失文件 都会在
  `calib check` 或 workflow 开始前失败，报错会带上 version/branch/measurement/字段。
- 参考：不同读者 / 包格式判定、HK 截取、找峰、Co60 双高斯等工程细节见
  [docs/data.md](docs/data.md)；每个载荷现有 workflow 见 [docs/workflows.md](docs/workflows.md)。

## 文档

详细文档在 [docs/](docs/) 目录（清单与写作约定见 [docs/README.md](docs/README.md)）：

| 文档 | 内容 |
|------|------|
| [docs/workflows.md](docs/workflows.md) | 每个载荷各自的显式 workflow（选点、reader、背景轮转、分辨率、特殊处理） |
| [docs/workflow_matrix.md](docs/workflow_matrix.md) | 历史 workflow 完整审计矩阵（每个版本实际做了什么） |
| [docs/deploy.md](docs/deploy.md) | 部署指南：新机器上装环境、挂数据、`calib check` 验证 |
| [docs/data.md](docs/data.md) | 数据准备与目录约定：新配置/manifest 格式、新载荷接入流程与方法论 |
| [docs/results.md](docs/results.md) | 结果产物与 QA 指标说明 |
| [docs/12B_13B/data.md](docs/12B_13B/data.md) | 12B/13B 类载荷的数据说明（点位对照表、各文件的特殊情况） |
| [docs/payloads_data.md](docs/payloads_data.md) | 其它 7 个载荷的数据说明（简化版：通用结构 + 各版本选点/排除） |

每类载荷的数据说明单独放一个目录（如 `docs/12B_13B/`），新增载荷时仿照添加。

## 命令速查

| 命令 | 说明 |
|------|------|
| `just` | 查看所有命令 |
| `just init {ver} {path}` | 初始化版本：软链数据 + 建目录 |
| `just check {ver}` | 校验新配置/manifest 层与数据链接（`--fix` 建缺失目录） |
| `just all {ver}` | 处理该版本全部数据（单拟合 + TB/EC 全局拟合） |
| `just fit-one {ver} {branch} {id}` | 处理单个 measurement 单拟合 |
| `just fit {ver} {branch}` | 处理单个分支单拟合 |
| `just global {ver} {branch}` | TB/EC 全局拟合 |
| `just discover {ver} {branch}` | 扫描目录生成 manifest 草稿 |
| `just list {ver} {branch}` | 列出已确认 measurement |
| `just config {ver} {branch} {id}` | 查看单个 measurement 的 resolved 配置 |
| `just new-payload {ver}` | 为新载荷生成 YAML 配置/目录骨架 |
| `just compare {ver}` | 新流程 vs 冻结 legacy oracle 差异对比 |

## 仓库结构

```
.
├── README.md
├── docs/                           # 详细文档（清单见 docs/README.md）
│   ├── README.md                   # 文档目录说明与清单
│   ├── workflows.md                # 每个载荷各自的显式 workflow
│   ├── workflow_matrix.md          # 历史 workflow 审计矩阵
│   ├── deploy.md                   # 部署指南
│   ├── data.md                     # 数据准备与目录约定、接入方法论
│   ├── results.md                  # 结果产物与 QA 指标
│   ├── payloads_data.md            # 其它 7 个载荷的数据说明（简化版）
│   └── 12B_13B/                    # 12B/13B 类载荷的数据说明（详细版）
│       └── data.md                 # 点位对照表、各数据文件的特殊情况
├── data/{ver}/                     # 各版本数据目录
│   ├── raw_data -> /path/to/data   # 原始数据软链
│   ├── single_process/             # 拟合产物（+ 历史 JSON oracle）
│   ├── tb_logs/                    # TB 面拟合产物
│   └── ec_logs/                    # E-C 拟合产物
├── lib_reader/                     # 各版本数据读取库（workspace 子包）
├── lib_plot/                       # 绘图库（workspace 子包）
├── src/calibration_process/
│   ├── cli.py                      # calib 命令入口
│   ├── pipeline.py                 # 高层编排（discover/list/fit/global/all）
│   ├── config_schema.py            # strict Pydantic schema（payload/analysis/fit_range/manifest）
│   ├── manifest.py                 # manifest 发现与加载（运行时不再扫目录）
│   ├── runtime.py                  # 解析后的运行时配置（resolved context）
│   ├── deploy.py                   # 部署校验 + 新载荷脚手架（calib check/scaffold）
│   ├── products.py                 # typed 中间产物（SingleFitResult/TBPoint/ECPoint）
│   ├── file_lib.py                 # 单文件读取与拟合（protected kernel）
│   ├── util_lib.py                 # 工具函数（拟合、缓存、阈值，protected kernel）
│   ├── configs/{ver}/              # 每版本 YAML 配置 + manifest
│   └── workflows/
│       ├── common.py               # 复用内核的 stage（single fit/points/global）
│       ├── registry.py             # 版本 -> workflow 模块（选择只发生一次）
│       └── versions/v{ver}.py      # 每版本显式 workflow（文件选择规则）
├── pyproject.toml
├── uv.lock
└── justfile
```

## TODO

- [x] 12B 简易结果
- [x] 12B 完整结果（TB: 备份目录全 54 点位温偏面拟合，适用偏压 ≥27.5V 内残差 <3.5%；EC: 与其他载荷相同的 K 边拆段二次拟合（EC_low/EC_high 输出 schema 一致），锚点为 Am241/Na22/Cs137/Co60（双高斯拟合 1332 keV）+ X光机 20-100 kV（低管压点的能量按管压赋值、有已知系统偏差，仅作参考；见 docs/12B_13B/data.md））
- [ ] 13B 完整结果
- [ ] 10B, 11B 塑闪结果

## Contributor

- wqd
- 丘智勇
- Lee
- yhl
- xsq(Haruka Kujo)
