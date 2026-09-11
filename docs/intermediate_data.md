# 中间数据与流程可定制性

修订日期：2026-09-08

本文回答四个问题：处理流程里有哪些中间数据（格式、怎么读）、能不能被本/其它项目
当输入做额外分析、能不能在本项目里定制一条"与众不同"的 pipeline、以及 pickle
是否依赖当前 Python 环境（是否需要完全相同的类）。

## 1. 中间数据清单（格式 / 位置 / 怎么读）

| 数据 | 位置 | 格式 | 怎么读 |
|------|------|------|--------|
| L1 忠实帧缓存（全版本） | `data/{ver}/l1/<key>/<kind>.parquet` + `<kind>.meta.json` + `index.json` | polars parquet（zstd） | 见本文「5. L1/L2 缓存」；亦可脱离本项目用纯 polars 读 parquet |
| L2 处理结果缓存（全版本） | `data/{ver}/l2/<key>/`（若干 parquet）+ `meta.json` + `index.json` | polars parquet（zstd，四通道叠加） | 见本文「5. L1/L2 缓存」 |
| 单谱拟合参数（可移植） | `single_process/{TB,EC}_fit_result/*.fit.json` | JSON（`fit_result` 4 通道参数） | `json.load` |
| 单谱能谱（可移植） | `single_process/{TB,EC}_fit_result/*.spectrum.parquet` | polars parquet（长表：channel/bin/x/spectrum/spectrum_err） | `polars.read_parquet` |
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

## 5. L1/L2 缓存（全版本）

修订日期：2026-09-12

读取分成两层，各写一层 parquet 缓存，都放在 `data/{ver}/` 下：

- **L1 忠实帧**（`data/{ver}/l1/<key>/<kind>.parquet`，`schema_ver=4`）：把原始帧里
  每个字段原样落成"天然单元"的一行——科学数据一个粒子一行，HK/时间线一次采样
  一行。包含帧头元数据、数据字段（如 `data_sum`）、CRC 数值，以及解析时算出的
  `crc_check`。处理（CRC 过滤、run 切分、单位换算、按通道分组）不属于 L1。
- **L2 处理结果**（`data/{ver}/l2/<key>/`，`schema_ver=5`）：各版本 reader 的
  `(sci, tel)` 最终输出。参差的四通道量按"四通道叠加"压成规整 parquet：每个键把
  4 个通道数组首尾相连成一列并记录各通道长度（`sublens`），读回时按边界切回；
  共享的 1-D/2-D 数组与空条目分别存表/存 `meta.json`。每个键还记录**逐通道
  dtype**（空通道默认 float64，不能与整型通道混成一个 dtype）。

读取入口：`single_readXX`（`lib_reader`）内部就是"取 L1 → 处理成 L2"；命中 L2 时
直接返回、不再碰 L1。`lib_reader.read_frames(path, ver, kind)` 只做 L1 解码，供只
想拿原始帧的项目使用。

一次性灌满全部版本用 `python scripts/warm_l1_cache.py [ver ...]`（只做读出、不拟合）。

### 布局与缓存键

```
data/{ver}/l1/<key>/
  sci.parquet  hk.parquet  tl.parquet     # 按 kind 各一张（存在哪些由载荷决定）
  <kind>.meta.json                         # schema_ver/reader/ver/kind/dtypes/shapes/rows
  （sci.meta.json 里 03B UDP 文件另有 legacy_drop_frame_idx）
data/{ver}/l1/index.json                   # 原始文件 -> key 映射
data/{ver}/l2/<key>/
  sci__chan__0.parquet  sci__flat__0.parquet  tel__chan__0.parquet ...
  meta.json                                # sections{files(sublens/shape/keys/dtypes),empties}
data/{ver}/l2/index.json
```

缓存键为 `sha256(f"{ver}|{reader}|{raw_path}|{parse_kwargs 的排序 json}")[:16]`，
`parse_kwargs` 只含影响解析的参数（如 `readSci` 的 `mode`、B 的 `feature_mode`），
**不含** `bin_width`/`adc_max`/`fit_range` 等后续拟合参数。`overwrite_cache=True`
会强制重算并覆盖缓存。

### 为什么不建议直接读 L1 文件

- L1 parquet 里是**逐字段的原始整数/bool 数组**（`amp` 由 L2 用 `data_max-data_base`
  重算），要拿能用于分析的"峰"仍要走 L2（`single_readXX`）。
- `meta.json` 记录各列 numpy dtype 与 2-D 形状，读回时按原值还原；`index.json` 只是
  "哪个原始文件对应哪个 key"的发现表。

### 跨项目读取（纯 polars，无本包 import）

```python
import json, polars as pl
idx = json.load(open("data/12B/l1/index.json"))
rec = next(r for r in idx if r["kind"] == "sci")
df = pl.read_parquet("data/12B/l1/" + rec["key"] + "/sci.parquet")   # 每粒子一行
amp = df["data_max"] - df["data_base"].to_numpy() / 4.0
```

`index.json` 的每条记录带相对路径字段 `file`，可直接使用。

### 03B 的 legacy 丢帧

旧实现用 `maxUdpReadout=10000` 分批读 UDP，跨批边界的科学帧会被丢弃。新实现 L1 按
完整格式解码全部帧，并把旧实现会丢的帧号写进 `sci.meta.json` 的
`legacy_drop_frame_idx`；L2 据此丢帧，保证输出与 legacy 逐字节一致。

> 缓存键与 `schema_ver` 详见 `lib_reader/src/lib_reader/l1_cache.py`；任何格式/语义
> 变化都必须 bump 常量以作废旧缓存。

## 5b. L3 拟合参数 / L4 输入

修订日期：2026-09-12

- **L3**：`run_single_fit` 除原有 dill pickle 外，新增可移植的
  `<stem>.fit.json`（4 通道拟合参数）与 `<stem>.spectrum.parquet`（长表
  `channel/bin/x/spectrum/spectrum_err`）。pickle 保留给旧消费方与冻结 oracle。
- **L4**：`load_single_fp_from_store` 从 `fit.json` 读 `fit_result`，TB 遥测
  （`tempSipm`/`bias`）改从 reader 的 **L2 处理结果**取（按同一 read spec 重解析，
  命中 L2 即返回），不再依赖 pickle 里的 `tel`。


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
| 10B/11B/12B | 各 `reader{10,11,12}/grid_packet.xml` | `readerXX/read.py` | L1 帧缓存 + L2 处理输出 |
| 04/07/09 | `reader07/grid_packet.xml` | `reader07/frame_adapter.py` | 三版合并为一个参数化 reader；`hex_sci_packet`（主事件 + 43 子事件，步进 11B）、`hex_tel_packet`（7×70B） |
| 03B/05B | `reader05/grid_packet.xml` | `reader05/frame_adapter.py` + `reader05/readout.py` | `sci_wf_packet`（waveform）/`sci_ft_packet`（feature 20×24B）；HK/timeline 解码与 UTC 拟合在自包含的 `readout.py` |

要点与坑：

- **04 与 07 的差异是两个科学常量**：`internal_resistance`（04=2.1、07=1.1）与 iMon
  除数（04 除以 2.0、07 除以 1.0）。这些常量放在 `configs/{ver}/reader.yaml`，由
  `file_lib` 读出后传给 reader（直接调用 reader 时用内建默认值），不"顺手统一"。
- **C 组科学包不是纯 XML 直出**：主事件在包头，43 个子事件在 repeat 区（步进 11B），
  适配层用 stable-argsort 按通道分组以复现 legacy 的"逐包、先主事件后子事件"顺序。
- **B 组 UDP 分块读取会造成帧丢失**：legacy 以 `maxUdpReadout=10000` 个 UDP 包为一批
  读取，跨批边界的科学帧被丢弃。新实现 L1 解码全部帧、把丢失帧号记进
  `legacy_drop_frame_idx`，由 L2 丢弃以复现 legacy 行为（见第 5 节）。
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
`getSpectrum`/`gehrelsErr`/`residualTempbias2D` 及 `resolutionFunction`/`headtime` 已迁到
中立包 `grid_common`（见第 8 节）。

版本级读取参数（`engine` + handler 注册表 + 常量）放在
`configs/{ver}/reader.yaml`，由 `calib check` 做 strict 校验；`file_lib.__read` 经
`lib_reader.READERS` 注册表统一分发，不再有 `if ending == ...` 的硬编码。

## 7. Reader golden（自包含回归）

`tests/test_reader_golden.py` + `tests/golden/reader/<sample>/`：每个样本提交一小段
截断真实原始文件（03B/05B 连同 HK/TimeLine/config 兄弟文件），并把 **legacy reader**
的 `(sci, tel)` 输出冻结成 `expected.npz` + `structure.json`。测试用**当前统一 reader**
在样本上重跑并逐字段对齐冻结值，零 `raw_data` 依赖、全新 clone 可跑。覆盖：C 组 hex
解码（07/04）、B 组 waveform noUdp（05B normal）、大小端 HK（05B xray）、waveform UDP
（03B src）、feature UDP 与逐文件 time cut（03B xray）。重新生成：
`python scripts/gen_reader_golden.py`。

## 8. 分层流水线与共享包

修订日期：2026-09-12

流水线现在是**声明式步骤**，可用 `--until` 停在任意层：

```
L1 忠实帧（reader 解码）        data/{ver}/l1/
L2 处理输出（reader 处理）      data/{ver}/l2/
L3 单谱拟合                    single_process/*/*.{fit.json,spectrum.parquet,pickle}
L4 构造点（TB/EC）             内存中（由 L3+L2 得到）
L5 全局拟合                    tb_logs/、ec_logs/
```

- `calib all {ver} --until L2` 只跑读出（建 L1/L2 缓存），不拟合；
  `--until L3` 跑完单谱拟合即停；`--until L4` 构造点后停；默认 `L5` 跑完整流程。
- L1/L2 由同一次 reader 调用一起产生；L4 的 TB 遥测来自 L2，`fit_result` 来自 L3。

共享数值（谱直方图 `getSpectrum`、`gehrelsErr`、residual/response 函数、
`resolutionFunction`、`headtime`）集中在中立包 **`grid_common`**，供
`calibration_process` 与 `lib_plot` 共用。这样 `lib_plot` 不再反向 import
`calibration_process`，workspace 内没有包循环；`lib_reader` 也不再持有拟合工具。

## 相关
- 结构：`README.md` 仓库结构
- 每版本选点：`workflows.md` / `payloads_data.md`
- 产物含义：`results.md`
- 示例脚本：`scripts/extra_analysis_example.py`
