# 中间数据与流程可定制性

修订日期：2026-09-08

本文回答四个问题：处理流程里有哪些中间数据（格式、怎么读）、能不能被本/其它项目
当输入做额外分析、能不能在本项目里定制一条"与众不同"的 pipeline、以及 pickle
是否依赖当前 Python 环境（是否需要完全相同的类）。

## 1. 中间数据清单（格式 / 位置 / 怎么读）

| 数据 | 位置 | 格式 | 怎么读 |
|------|------|------|--------|
| 原始数据解析缓存 | `.cache/` | cachier 自带缓存（dill blob / 逐文件 dill，reader09/10/12 为 `separate_files`，reader04/03B/05B/07 为共享大 blob） | **不要直接读文件**；用 reader 函数 `single_readXX(path)`，缓存命中即返回解析结果 |
| 单谱拟合 pickle | `single_process/{TB,EC}_fit_result/*.pickle` | **dill** pickle 的 dict | `dill.load(path)` 或 `util_lib.pickle_load(path)` |
| 单谱拟合图 | `single_process/single_fit_fig/*.png` | PNG | `matplotlib.image.imread` |
| TB 二维面拟合 | `tb_logs/*_temp_bias_fit.json` | JSON（每通道 G0/k/V0/b/c+err+redchi+ndf） | `json.load` |
| EC 系数 | `ec_logs/*_ec_coef_sci_ch{i}.json` | JSON（EC_low/EC_high/resolution_* 及误差） | `json.load` |
| EC 数据点 | `ec_logs/*_ec_data_ch{i}.npy` | NumPy `[[energy...],[center...]]` shape (2, n) | `numpy.load` |
| EC / resolution 图 | `ec_logs/*.png` | PNG | `matplotlib.image.imread` |
| 配置 / manifest | `configs/{ver}/*.yaml` | YAML（strict） | `yaml.safe_load` + Pydantic schema |
| QA 阈值 | `single_process/qa_thresholds.json` | JSON | `json.load` |

**单谱 pickle 的内容**（最"完整"的中间产物）：

```python
import dill
d = dill.load("single_process/EC_fit_result/15keV.pickle")
# d = {
#   'file': 原始路径,
#   'fit_result': [ {a,b,c,a_err,b_err,c_err,rate,rate_err,resolution,resolution_err,
#                    bkg,redchi,ndf,success,boundary_hit,qa_flag} × 4 或 null ],
#   'spectrum': [4 个 ndarray], 'x': [4 个 ndarray],
#   'tel': {tempSipm/bias/timestamp/... 每通道数组},
#   'config': {'read_config','spectrum_config'(含 corr 闭包),'fit_config'}
# }
```
其中 `spectrum`/`x`/`tel` 是 NumPy 数组，`fit_result` 才是你分析最常用的。

## 2. 能不能把它当输入做额外分析？

**能，推荐程度按"跨环境可移植性"排：**

1. **`ec_logs/*.npy`、各种 `*.json` → 最稳**：纯 NumPy/JSON，**无任何 Python 类依赖**，
   拷贝到任何项目、任何语言都能读。适合跨项目复用（如 E-C 曲线、分辨率、TB 参数）。
2. **`*.pickle`（单谱结果）→ 本项目内 / 装了这个包的环境**：内容最全（谱、遥测、
   逐通道拟合、配置），但要 **dill** 且能 import 本包（见第 4 节）。适合在本仓库或
   装了本包的同事机器上做二次谱分析。
3. **原始数据 + reader → 最自由**：用 `lib_reader` 的 `single_readXX(path)` 直接拿
   原始 `(sci, tel)`，想怎么处理都行，不依赖任何中间文件。

示例见 [scripts/extra_analysis_example.py](../scripts/extra_analysis_example.py) 的
`demo_read_intermediate()`：分别读 pickle / npy / json 并取值。

## 3. 能不能定制一条与众不同的 pipeline？

**能。** 管线是**可独立调用的函数**，不是黑盒：
- `calib` 的每个 stage 都在 `workflows/common.py`（`run_single_fit`、
  `build_tb_points`、`build_ec_points`、`global_tb`、`global_ec`）与
  `pipeline.py`（`fit_branch`/`global_*`/`all`）里，可被脚本直接 import。
- 中间边界是 **typed 数据**（`SingleFitResult`/`TBPoint`/`ECPoint`，见
  `products.py`），所以你可以在"第一步单拟合 → 第二步造点"之间**插入任何加工**，
  再把结果传给下游的全局拟合。

例如：`只做前两步，中间加一道加工，再走后面的步骤`：

```python
from calibration_process import pipeline
from calibration_process import manifest as man
from calibration_process.workflows import common as stages

rt = pipeline.load_rt("09")
manifest = man.load_manifest(pipeline._manifest_path("09", "tb"))

# ① 单拟合（可只跑前 N 个 measurement）
for m in man.filtered_measurements(manifest)[:3]:
    stages.run_single_fit(rt, "tb", m, pipeline.output_root("09"))

# ② 载入所有单拟合并构建 TB 点
items = [(m, stages.load_single_fp_from_store(rt, "tb", m, pipeline.output_root("09")))
         for m in man.filtered_measurements(manifest)]
per_channel = stages.build_tb_points(rt, items)

# ③ 自定义加工：例如按温度过滤 / 改某通道 / 加偏移
custom = [[p for p in pts if p.temperature < 35.0] for pts in per_channel]

# ④ 交给现有全局拟合
stages.global_tb(rt, custom, pipeline.output_root("09") / "tb_logs_custom")
```

完整可运行示例见 `scripts/extra_analysis_example.py`（已跑通：读中间产物 + 插入
"按温度过滤"这一步 + 用过滤后的点做 TB 全局拟合）。

> 也可以直接用 `calib fit-one {ver} tb {id}` 处理单个 measurement，逐点检查后再
> 手动串起来。粒度控制（fit-one / fit / global / all）见 justfile。

## 4. Pickle 依赖当前 Python 环境吗？需要完全相同的类吗？

**依赖，但"依赖"的具体含义是：**

- 它是个 **dill pickle**（`util_lib.pickle_save` 用 `import dill as pickle`）。读取要
  用 **`dill.load`**（`util_lib.pickle_load`）；纯标准库 `pickle` 可能无法还原里面的
  闭包/数据类。
- dill 通过 **模块路径** 还原对象（例如
  `calibration_process.file_lib.Read_config`、
  `calibration_process.runtime` 里的 corr 闭包、`calibration_process.products`
  的 dataclass）。所以**本包必须能 import**。换到没有这个包的环境，dill 会报
  `ModuleNotFoundError`（我们迁移时删了 `operation.py`，旧 oracle pickle 就因此
  无法加载，是同一个原因）。
- **不要求"内存里完全相同的实例"**：字段值按值存进 `__dict__` 还原；类本身取自
  当前包。所以只要**模块路径不变、类名一致**，即使类字段有变化，实例的旧字段值
  通常也能还原（dataclass 实例的 `__dict__` 保留）。

**结论**：
- 想跨环境/跨项目重放，**优先 .npy / .json**（零类依赖）。
- 只有在"本仓库 / 装了本包且用了 dill"的环境里才可靠地读 pickle。
- 若只想做谱级二次分析且不想绑死在这套类上，稳妥做法是**存成 npz/CSV** 或在读
  pickle 后把需要字段另存为 npy/json；本项目若要加这类"可移植中间产物"，可再加
  一个 `calib export` 子命令（目前未做）。

## 相关
- 结构：`README.md` 仓库结构
- 每版本选点：`workflows.md` / `payloads_data.md`
- 产物含义：`results.md`
- 示例脚本：`scripts/extra_analysis_example.py`
