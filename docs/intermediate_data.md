# 中间数据与流程可定制性

修订日期：2026-09-08

本文回答四个问题：处理流程里有哪些中间数据（格式、怎么读）、能不能被本/其它项目
当输入做额外分析、能不能在本项目里定制一条"与众不同"的 pipeline、以及 pickle
是否依赖当前 Python 环境（是否需要完全相同的类）。

## 1. 中间数据清单（格式 / 位置 / 怎么读）

| 数据 | 位置 | 格式 | 怎么读 |
|------|------|------|--------|
| L1 解析缓存（全版本） | `data/{ver}/l1_cache/` | polars parquet（zstd）+ `meta.json` + `index.json` | 见本文「5. L1 parquet 缓存」；亦可脱离本项目用纯 polars 读 parquet |
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

## 5. L1 parquet 缓存（全版本）

修订日期：2026-09-10

为了替代 dill，本仓库给 reader 结果加了一层 L1 parquet 缓存，统一放在
`data/{ver}/l1_cache/`。按缓存边界分两种：

- **帧缓存**（`schema_ver=1`，10B/11B/12B）：缓存 `readSci`/`readHK` 的**原始解析
  结果**（每事件一行），命中后**重跑** `single_readXX` 的既有后处理（amp、遥测换算、
  稳定段截取、按通道分组），对科学结果透明。
- **最终输出缓存**（`schema_ver=2`，03B/04/05B/07/09）：缓存整个 `(sci, tel)` 的
  **最终输出**，命中即返回（不重跑）。参差的四通道量按"四通道叠加"压成规整 parquet：
  每个通道键把 4 个通道数组**首尾相连**成一列、并记录各通道长度（`sublens`），读回时
  再按边界切回 4 个通道；共享数组与空条目分别存表/存 `meta.json`。

一次性灌满全部版本用 `python scripts/warm_l1_cache.py [ver ...]`（只做读出、不拟合）。

### 布局与缓存键

```
data/{ver}/l1_cache/
  index.json                  # 原始文件 -> 缓存映射（每条记录一个原始文件）
  cache/{key}/
    events.parquet            # readSci 的 raw 解析结果（每事件一行）
    tel.parquet               # readHK  的 raw 解析结果（每 HK 记录一行）
    meta.json                 # schema_ver / reader / ver / kind / dtypes / rows
```

缓存键为 `sha256(f"{ver}|{reader}|{kind}|{raw_path}|{parse_kwargs 的排序 json}")[:16]`，
`parse_kwargs` 只含影响解析的参数（如 `readSci` 的 `mode`；`readHK` 为空），
**不含** `bin_width`/`adc_max`/`fit_range` 等后续拟合参数。`overwrite_cache=True`
会强制重解析并覆盖缓存。kind 为 `sci`/`tel`，二者键不同、各占一个 `cache/{key}/`
目录。

### 为什么不建议直接读缓存文件

- parquet 里存的是**原始整数/bool 数组**（`amp` 不缓存，由后处理重算），所以要读
  出能用于分析的"峰"仍要走 `single_readXX`。
- `index.json` 只是"哪个原始文件对应哪个 parquet"的发现表；`meta.json` 里记录各列
  的 numpy dtype，用于读取时若有宽度提升则按原 dtype 回退。

### 跨项目读取（纯 polars，无本包 import）

```python
import json, polars as pl
idx = json.load(open("data/12B/l1_cache/index.json"))
rec = next(r for r in idx if r["kind"] == "sci")
df = pl.read_parquet(rec["events"])          # 每事件一行
amp = df["data_max"] - df["data_base"].to_numpy() / 4.0
```

### 「四通道叠加」缓存布局（03B/04/05B/07/09）

```
data/{ver}/l1_cache/cache/{key}/
  meta.json               # schema_ver=2, reader/ver/kind, sections{sci,tel}
  sci__chan__0.parquet    # 通道键（同 sublens 的合一张表）：每键一列=concat(4通道) + __channel__
  sci__flat__0.parquet    # 共享 1-D 数组（如 effectiveCount/missingCount）
  tel__chan__0.parquet
  tel__flat__0.parquet
```

`meta.json` 的 `sections.<name>.files` 记录每张表是 `channel`（含 `sublens`/`keys`）还是
`flat`（含 `shape`/`keys`），`empties` 记录空条目（空 list / None）。读回外部工具时按
`sublens` 用 `np.split` 切回四通道即可；`__channel__` 列给出每行归属的通道号。

> 缓存键与 `schema_ver` 详见 `lib_reader/src/lib_reader/l1_cache.py`；命中校验同时检查
> `schema_ver` 与 `kind`，任何格式/语义变化都必须 bump 常量以作废旧缓存。

## 6. 读取层统一（packet_parser / frame_io）

修订日期：2026-09-11

各载荷的原始数据在**字节层**都是"一串按包成帧的数据"，只是存储方式有三类：
直接二进制（10B/11B/12B）、UDP 包裹的二进制（03B/05B，内层才是科学帧）、
空格分隔的十六进制文本（04/07/09，本质是字节流的 hexdump）。据此把读取的**帧层**
统一：

- `lib_reader/packet_parser.py` + `parity_check.py`：XML 驱动的包解析器
  （`parse_grid_data_new`），已合并原来 reader10/11/12 的三份副本；每个载荷只保留
  自己的 `grid_packet.xml`，调用时显式传 `xml_file`。
- `lib_reader/frame_io.py`：字节来源层——`load_binary`（二进制）、`load_hex_text`
  （hex 文本→字节）。03B/05B 的 UDP 分块与 HK/timeline 解码在
  `lib_reader/reader05/readout.py`（见下）。

各载荷的帧定义与适配层（均已切到统一解析器）：

| 载荷族 | 帧 XML | 适配层（帧表→`(sci, tel)`） | 说明 |
|--------|--------|------------------------------|------|
| 10B/11B/12B | 各 `reader{10,11,12}/grid_packet.xml` | `readerXX/read.py` 的后处理 | 原始帧缓存（schema_ver=1） |
| 04/07/09 | `reader07/grid_packet.xml` | `reader07/frame_adapter.py` | 三版合并为一个参数化 reader；`hex_sci_packet`（主事件 + 43 子事件，步进 11B）、`hex_tel_packet`（7×70B） |
| 03B/05B | `reader05/grid_packet.xml` | `reader05/frame_adapter.py` + `reader05/readout.py` | `sci_wf_packet`（waveform）/`sci_ft_packet`（feature 20×24B）；HK/timeline 解码与 UTC 拟合在自包含的 `readout.py` |

要点与坑：

- **04 与 07 的差异是两个科学常量**：`internal_resistance`（04=2.1、07=1.1）与 iMon
  除数（04 除以 2.0、07 除以 1.0）。合并时做成版本参数 `_PARAMS`，不"顺手统一"。
- **C 组科学包不是纯 XML 直出**：主事件在包头，43 个子事件在 repeat 区（步进 11B），
  适配层用 stable-argsort 按通道分组以复现 legacy 的"逐包、先主事件后子事件"顺序。
- **B 组 UDP 分块读取会造成帧丢失**：legacy 以 `maxUdpReadout=10000` 个 UDP 包为一批
  读取，跨批边界的科学帧被丢弃。适配层复用 `readout.extractSciRawData` 按同样的批次
  切分，**忠实复现**这一行为（否则 03B 输出会与 legacy/oracle 不一致）。
- **`cutFileRef` 逐文件时间截断**：部分 03B X 光机文件带按文件名指定的 time cut，
  适配层在末尾复现该截断。
- HK/timeline 解码（`extractHKData*`/`extractTimelineData*`/`getUTC`）与 `findPackPos`/
  `doFitLin` 已从 legacy 单体逐字抽到 `reader05/readout.py`（约 300 行），B 族不再依赖
  大文件。

验证：
- `tests/test_frame_io.py`：`load_binary` / `load_hex_text` 字节来源。
- `tests/test_reader_golden.py`：B/C 各路径的截断真实样本冻结为 legacy `(sci,tel)`
  输出，统一 reader 重跑逐字段（dtype + 值）对齐（见第 7 节）。

迁移收尾已删除：C 组单体 `gridBasicFunctions02.py`（reader04/07）、B 组单体
`gridBasicFunctions.py`/`gridParametersCommon.py`/`gridProcessFunctions*/`、
`extractSciEvents`/`dataReadout`/`fitBaseline` 等编排函数、`*_legacy` 包装与全部
`@cachier`（含依赖）。B 族现在只依赖 `reader05/readout.py` 的小函数；`util_lib` 用到的
`getSpectrum`/`gehrelsErr`/`residualTempbias2D` 迁到 `reader05/fit_utils.py`。

## 7. Reader golden（自包含回归）

`tests/test_reader_golden.py` + `tests/golden/reader/<sample>/`：每个样本提交一小段
截断真实原始文件（03B/05B 连同 HK/TimeLine/config 兄弟文件），并把 **legacy reader**
的 `(sci, tel)` 输出冻结成 `expected.npz` + `structure.json`。测试用**当前统一 reader**
在样本上重跑并逐字段对齐冻结值，零 `raw_data` 依赖、全新 clone 可跑。覆盖：C 组 hex
解码（07/04）、B 组 waveform noUdp（05B normal）、大小端 HK（05B xray）、waveform UDP
（03B src）、feature UDP 与逐文件 time cut（03B xray）。重新生成：
`python scripts/gen_reader_golden.py`。

## 相关
- 结构：`README.md` 仓库结构
- 每版本选点：`workflows.md` / `payloads_data.md`
- 产物含义：`results.md`
- 示例脚本：`scripts/extra_analysis_example.py`
