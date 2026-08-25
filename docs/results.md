# 结果产物与 QA 指标

修订日期：2026-08-25

## 产物清单

| 产物 | 路径 | 格式 | 说明 |
|------|------|------|------|
| TB 单谱拟合 | `data/{ver}/single_process/TB_fit_result/*.pickle` | pickle | 每个数据文件一个，含 4 通道拟合结果 |
| EC 单谱拟合 | `data/{ver}/single_process/EC_fit_result/*.pickle` | pickle | 每个数据文件一个，含 4 通道拟合结果 |
| TB 面拟合 | `data/{ver}/tb_logs/*_temp_bias_fit.json` | JSON | 4 通道的温度-偏压二维面拟合参数 |
| E-C 系数 | `data/{ver}/ec_logs/*_ec_coef_*.json` | JSON | 每通道一组能量-道址关系系数 |
| E-C 数据 | `data/{ver}/ec_logs/*_ec_data_*.npy` | numpy | 拟合用的原始数据点 |
| 拟合图 | `data/{ver}/single_process/single_fit_fig/*.png` | PNG | 每个数据文件每个通道一张 |

## pickle 内容

单谱拟合 pickle 包含以下字段：

- `file`：原始数据文件路径
- `fit_result`：4 元素列表，每个元素为该通道的拟合结果 dict（详见下文），`null` 表示未配置拟合区间
- `spectrum` / `x`：能谱数据和道址
- `tel`：遥测数据（温度、偏压等）
- `config`：本次拟合使用的配置

### fit_result dict 字段

| 字段 | 类型 | 说明 |
|------|------|------|
| `a`, `b`, `c` | float | 峰面积、峰位、峰宽（sigma） |
| `a_err`, `b_err`, `c_err` | float | 对应参数的误差 |
| `rate`, `rate_err` | float | 计数率及其误差 |
| `resolution`, `resolution_err` | float | 能量分辨率（FWHM/峰位）及其误差 |
| `bkg` | dict | 本底拟合信息和函数 |
| `redchi` | float | 约化卡方（χ²/ndf），拟合质量的核心指标 |
| `ndf` | int | 自由度（数据点数 - 拟合参数数） |
| `success` | bool | 拟合是否收敛 |
| `boundary_hit` | list | 顶到边界的参数列表，如 `["peak_sigma@max"]` |
| `qa_flag` | str | 质量判定：`"ok"` / `"warn"` / `"fail"` |

## QA 指标

### redchi（约化卡方）

拟合残差相对于误差棒的归一化度量。理想值为 1；大于 1 说明拟合残差大于误差棒预期；远大于 1 说明拟合质量差或误差棒被低估。

TB 和 EC 的 redchi 分布不同（EC 数据物理更复杂，redchi 天然偏高），阈值分别设置。

### success

lmfit 拟合器是否报告收敛。`False` 通常意味着拟合发散或参数无法确定。

### boundary_hit

参数值是否顶到了人为设定的边界（如峰位顶到拟合区间端点、峰宽顶到上限）。顶到边界说明拟合结果不可靠，需要检查拟合区间配置或数据质量。

### qa_flag

三级判定，基于 redchi 与阈值的比较：

| 标签 | 含义 |
|------|------|
| `ok` | redchi ≤ warn 阈值，拟合质量正常 |
| `warn` | warn < redchi ≤ fail 阈值，建议检查 |
| `fail` | redchi > fail 阈值或拟合失败，需要人工审查 |

阈值配置在 `data/{ver}/single_process/qa_thresholds.json`。默认值：
- TB：warn=2.0, fail=5.0
- EC：warn=5.0, fail=50.0

### TB 二维面拟合

`tb_logs/*_temp_bias_fit.json` 包含 5 个物理参数（G0, k, V0, b, c）及其误差，以及 `redchi` 和 `ndf`。该 redchi 反映的是温度-偏压二维面上峰位的整体拟合质量，受系统误差主导，数值通常在 40–150 范围内（远大于单谱 redchi），属于正常现象。
