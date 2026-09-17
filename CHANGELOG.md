# Changelog

本项目从当前提交起维护 Changelog。条目**不自动生成**，只在需要时人工添加（例如重构、新版本接入、行为变更等）。历史提交不再补记。

版本号遵循 [Semantic Versioning](https://semver.org/)，并见 `pyproject.toml`。

## [0.4.0] - 2026-09-18

### Added
- **接入 14B / 13B / 15B 三个载荷，均跑到 TB 单谱 + 二维全局拟合**（新 reader、
  `configs/{ver}/`、`workflows/versions/v{ver}.py`）：
  - 14B：新增 `reader14B`（584B `grid1x_ft_packet` + 178B `hk_packet`，XML 为
    `yunyao_packet.xml`；该 XML 里的 224B HK 与实际数据不符，实际 HK 复用 N1 的
    178B 定义），数据 `14B15B/03`。
  - 13B：复用 12B 的 grid1x reader，数据为 12B、13B 目录下的设备 **073**。
  - 15B：复用 14B reader，数据为 `14B15B/02`。
- `tests/golden/GRIDN1`、`tests/golden/12B` 之外的自洽基准思路不变；本轮新增的是
  载荷接入与初步 TB 结果，未改动单谱/TB 数学模型。

> **13B、15B 的 TB 结果未经人工核验**，只是流程接入与初步结果，核验前不要用于
> 正式标定。已知的排除项与残差：
> - 13B：`m20C_280` 与 `m20C_283` 逐字节重复（280 真实数据丢失）、`m20C_285` 的 HK
>   平台实为 28.7V/−8°C 的重测，两者 `use:false`；26.5V 行贴近阈值沿、单谱 redchi
>   很大，二维拟合用 `bias_min_filter: 26.75` 排除。保留点相对残差 max < 5%。
> - 15B：`30C_265`（最高温+最低偏压，四通道单谱 QA fail）`use:false`；保留点相对
>   残差 max < 5%。

### Changed
- **N1 / 12B / 14B 报出的 `bias` 不再扣 499Ω 串联电阻压降**（`bias = vMon`）：跨版本
  检查（`docs/intermediate_data.md` §6.1）显示这几型的 HK 监测值就等于 PID 设定值，
  扣 I·R 会让结果系统性低于文件名；依据待硬件确认。为此 L2 处理缓存 schema 5→6，
  N1 reader golden 按新约定重新冻结，12B 与 legacy oracle 的逐字节比对暂时豁免。
- `reader12` / `reader14B` 的 L1/L2 缓存版本改为按原始路径推导，使同一个 reader 能
  服务 12B/13B、14B/15B 两组数据（各自缓存到 `data/{ver}/`）。
- 13B 的 `fit_range_tb.yaml` 以 12B 窗口中心为先验、在 13B 直方图上重新找峰生成；
  15B 沿用 14B 的找峰规则。

## [0.3.0] - 2026-09-14

### Changed
- **流水线式改造完成**：单谱拟合（L3）→ 构造点（L4）→ 全局拟合（L5）走显式 workflow
  （`workflows/versions/v{ver}.py` 每版本一个）+ 人工确认 manifest + `configs/{ver}/*.yaml`
  （strict Pydantic schema）；reader 统一到注册表并产出 L1 忠实帧 / L2 处理结果缓存；
  `calib all {ver} --until L1..L5` 可在任一层停下。原有科学数学（单谱/TB/EC 模型、reader
  数据结构）不变，本次是把接入与编排显式化、可回归。
- GRIDN1 按流水线组织为 `configs/GRIDN1/{GAGG,CLYC,EC,Neutron}` 与 `data/GRIDN1/*` 四个子版本。
- 12B、GRIDN1/GAGG、GRIDN1/CLYC 去掉 `bias_min_filter`、`skip_qa_fail` 等默认过滤，改为
  先全量单谱拟合、再由 QA 报告逐通道显式排除（写入 manifest 的 `channels.use:false`）。
- GRIDN1/EC：ch1/ch2（GAGG）保持 Gd K 边分段二次，ch0/ch3（CLYC）改为单条线性；分辨率
  改用非负约束拟合 `polyfit_nn`（裸 `polyfit` 在点少时会给出负系数，曲线在高能端断开）。

### Added
- **GRIDN1 初步标定结果**（自洽基准，`tests/golden/GRIDN1/`）：
  - TB：GAGG（Am241，ch1/ch2）与 CLYC（Na22，ch0–3）温度-偏压二维面，相对残差 max < 5%。
  - EC：源 4 点（Cs137/Co60）+ X 光 13 个能量，四通道 E-C 相对偏差 max ≤ 4.4%。
  - 中子束流 TB 快照（固定偏压、简并，仅作自洽回归基准）。
- GRIDN1 低能 X 光（15–35 keV，ft 包）接入 `ec.x_path_low` + `xray_reader_low`；本轮不参与
  E-C 拟合（峰贴近阈值、跨通道本底留大负凹陷），原因见 `docs/N1/data.md`。
- `scripts/fit_qa.py`：只读的单谱拟合 QA 报告（分类、redchi、边界命中、残差等）。
- `grid_common` 共享数值包；`scripts/gen_golden.py` 支持为有意变更的版本重建自洽 golden。

### Fixed
- 12B 救回 10 个原本全 null 的低偏压 TB 点；`30C_265`（文件异常）与 `20C_265` ch2 显式排除。
- EC 绘图支持各通道能量点集不同（GRIDN1 ch1/ch2 有 X 光、ch0/ch3 仅源锚点），GRIDN1 也能出图。
- `_center_fit` 对"点数 = 阶数+1"的段退化为插值拟合，避免 `global_ec` 崩溃。

### Removed
- GRIDN1/EC 的 Th228（583/2614 keV）锚点：该源测量太弱，拟合面积显著性 `a/a_err` 仅
  1.5–2.9，2614 峰实为噪声涨落，583 也仅为连续谱上的弱肩。

## [0.2.0] - 2026-09-09

### Changed
- 重构了整个项目框架与 workflow：显式 workflow（`workflows/versions/v{ver}.py` 每版本一个）+ 人工确认的 manifest + `configs/{ver}/*.yaml`（YAML + strict Pydantic schema）+ 统一 `calib` CLI。
- 将隐藏在版本继承、配置文件、目录扫描中的处理流程，改为显式、分层、可复用、可测试的 scientific workflow。
- **原有的科学数据处理流程与处理代码保持不变**：单谱拟合、TB/EC 数学模型、reader、正式绘图、中间/最终产物格式均与迁移前一致；本次为纯架构重构，结果经 `legacy vs new` 全量回归验证一致。
