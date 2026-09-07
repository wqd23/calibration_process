# Historical Workflow Audit (workflow_matrix)

修订日期：2026-09-08

本文档是架构迁移前的历史 workflow 审计文档（**不参与 runtime**）。它描述每一个历史版本实际上做了什么，而不是"看起来应该做什么"。无法从代码确定的行为标注为 `UNKNOWN`。

审计来源：

- `src/calibration_process/process.py`（CLI 入口与 `OPERATION_SPEC`）
- `src/calibration_process/cmd.py`（`VersionProcessOp` / `VersionProcessOp10B` / `VersionProcessOp12B`）
- `src/calibration_process/operation.py`（各版本 TB / EC Operation 类）
- `src/calibration_process/file_lib.py`（单文件读取 + 峰拟合）
- `src/calibration_process/util_lib.py`（工具函数）
- `src/calibration_process/config.json`（版本配置）
- `lib_reader/`（各版本 reader）

---

## 0. 共享的科学内核（Protected Kernel）

所有版本共用以下科学函数（本次迁移不改）：

- 单文件处理：`file_lib.File_operation_05b`，其 `get_spectrum()` / `peak_fit()` 是唯一的单文件科学路径。
- 峰拟合：`util.peak_fit`（Gaussian + 可选背景）或 `fitting.peak_fit`（支持 peak_form / compton，仅 10B 之后的 `File_operation_10b` 使用）。
- TB 二维面拟合：`util.temp_bias_fit_curvefit`（默认）或 `util.temp_bias_lmfit`（11B）。
- EC 能量标定：`Operation.center_fit`（polyfit deg=2）+ `Operation.resolution_fit`（每版本不同）。
- 最终序列化：`util.pickle_save` / `util.json_save` / `np.save`。

---

## 1. 入口与调度

`process.py` 中 `OPERATION_SPEC` 决定每个版本的 `(wrapper, tb_op, ec_op, extra_kwargs)`：

| Version | wrapper | fp_method | suffix |
|---|---|---|---|
| 03B | `VersionProcessOp` | "03" | dat |
| 04 | `VersionProcessOp` | "04" | txt |
| 05B | `VersionProcessOp` | `None` | dat |
| 07 | `VersionProcessOp` | "07" | txt |
| 10B | `VersionProcessOp10B` | "10" | dat |
| 11B | `VersionProcessOp10B` | "11" | dat |
| 09 | `VersionProcessOp` | "09" | txt |
| 12B | `VersionProcessOp12B` | "12" | dat |

**关键调度事实（重要）：**

- `VersionProcessOp.__init__` 中，`tb` 分支的 `file_list_op` **不传 fp_method**，`ec['x']` 传 fp_method，`ec['src']` 不传。
- `operation.process(op, file, fp_method=None)` 内：
  ```python
  fp_method = __get_fp05B if fp_method is None else __get_fp03B
  ```
  因此：
  - **TB 单文件**、**EC source 单文件** 永远走 `__get_fp05B`（单 4 通道文件）。
  - **EC X-ray 单文件**：fp_method 非 None 的版本（03B/04/07/09/10B/11B/12B）走 `__get_fp03B`（每个通道一个独立文件，4 文件重建）；fp_method=None 的版本（05B）走 `__get_fp05B`（单 4 通道文件）。

这个差异是历史 orchestration 的核心差异之一，不能统一。

`fp_method` 实际取值与 `file_lib.__read` 的 `ending` 是**两套独立的字符串**，不要混淆：
- `fp_method`：决定使用 `__get_fp05B` 还是 `__get_fp03B`。
- `Read_config.ending`：决定 `single_readXX` 哪个 reader。

---

## 2. 矩阵总览

> 表内"reader"列为 `Read_config.ending` 实际选中的 reader；"TV corr 位置"为 temp/bias 校正发生层级；"特殊处理"记录非通用逻辑。

| Version | Raw inputs | Reader (ending) | TB flow | EC source flow | EC X-ray flow | TV correction position | Special selection | Special handling |
|---|---|---|---|---|---|---|---|---|
| 03B | sci `.dat` + scienceConfig `.json`（config_file） | `03b`（normal/src: `03b-src`） | 峰拟合→点收集→`temp_bias_fit_curvefit` | `03b-src` 读法 | 每能量 4 通道文件，`__get_fp03B` 重建 | **EC spectrum 层**（corr 在 `get_spectrum`） | TB 排除 `baseline`/`CI`/`50C`；X 用 `_rundata`+非 CI；src 硬编码 3 个 + 硬编码 bkg | `resolution_lmfit`；`energy_split_low=49, high=51`；EC src bkg 硬编码映射 |
| 04 | sci `.txt` | `04` | 峰拟合→点收集→`temp_bias_fit_curvefit` | `04` 读法 | 每能量 4 通道文件，`__get_fp03B` | **EC spectrum 层** | TB 全 `.txt`；X 排除 `_18p0_`；src 硬编码 4 个 + 硬编码 bkg | `resolution_ExprFit`；split low=49 high=55 |
| 05B | sci `.dat`（normal/src）、sci+config+time_cut（xray） | `normal` / `xray` | 峰拟合→点收集→`temp_bias_fit_curvefit` | `normal` 读法，带 bkg | xray 单 4 通道文件，`__get_fp05B`，bkg=同文件 rotated time_cut | **EC spectrum 层** | TB 排除 `50C`；X 排除 `XM_22`，`_observe.dat`；src 排除 `bkg` | 唯一 xray 走 fp05B；xray 用 `config_file` + `time_cut.json` + `bkg_time_cut`（rotated）；split low=49 high=52 |
| 07 | sci `.txt` | `07` | 峰拟合→点收集→`temp_bias_fit_curvefit` | `07` 读法 | 每能量 4 通道文件，`__get_fp03B` | **EC spectrum 层** | TB 全 `.txt`；X 排除 `40p0/15p0/12p0/99p9/90p1`；src 动态按 `bkg{src}` | `resolution_polyfit`；split low=49 high=55 |
| 10B | sci `observe.dat` + `hk.dat`（read10 自行 replace observe→hk） | `10b` | 峰拟合→点收集→`temp_bias_fit_curvefit` | `10b` 读法 | 每能量 4 通道文件，`__get_fp03B` | **EC spectrum 层** | TB 排除 `50C_265`；X 排除 `hk`；src 硬编码 4 个 + 硬编码 bkg | **EC 只有 3 通道（CHN_NUM=3）**，后补丁 ch3=ch0；`resolution_polyfit`；默认 split |
| 11B | sci `observe.dat` + `hk.dat`（getHK 按 stem 前缀配对） | `11b` (mode=wf) | 峰拟合→点收集→**`temp_bias_lmfit`** | `11b` 读法 | 每能量 4 通道文件，`__get_fp03B` | **EC spectrum 层** | TB 特殊 glob 4 子目录 + 多个显式 remove + `_50_Cs_2` 过滤；X 排除 `hk`，并 remove 能量 `20`；src 硬编码 4 个（Ba133 无 bkg） | TB 用 `temp_bias_lmfit`；EC 仅 3 通道补丁；split low=49 high=55 |
| 09 | sci `.txt` | `09`（= reader07 的 cache 别名） | 峰拟合→点收集→`temp_bias_fit_curvefit` | `09` 读法 | 每能量 4 通道文件，`__get_fp03B` | **EC spectrum 层** | TB 全 `.txt` 但**显式 remove 6 个文件**；X 排除 `65keV_`；src 硬编码 3 个 + bkg | `resolution_polyfit`；`fit_range.get` / `bkg_form.get` 带默认值兜底；split low=49 high=55 |
| 12B | sci `.dat` + 显式 `hk`（tb_file_map.json 或 (kv,ch) 配对） | `12b` (mode=ft) | 峰拟合→点收集→`temp_bias_fit_curvefit`（覆盖初始 p0）；load 时 bias<27.25 处丢弃 | `12b` 读法 | 每能量 4 通道文件，`__get_fp03B` | **EC spectrum 层** | TB 用 `tb_file_map.json`，EXCLUDE (-20,275)/(-20,285)，`sci_half`/`hk_bias` 特殊 | EC global 走 `VersionProcessOp12B.ecfit`；X 过滤 `old`，hk 完整性+fit_range 完整性；src 共享 bkg `0611env.dat`；split low=49 high=55 |

---

## 3. 逐版本详细说明

### 3.1 输入文件组成

- **03B/05B（旧包格式，binary/struct）**：TB 读单个 rundata `.dat`（03B 还需 scienceConfig `.json` 作为 config_file；05B 不需要）；EC xray（05B）额外需要 `config.json`（x_config）+ `time_cut.json`。
- **04/07/09（hex/ASCII `.txt`）**：单文件；`07`/`09` 不需配对文件；`04` 单文件。
- **10B/11B/12B（新包格式）**：science 文件（`observe`/`*.dat`）+ HK 文件。配对方**因版本而异**：
  - 10B：`observe` → 替换为 `hk`。
  - 11B：按 stem 前缀（`getHK` 在目录里找 `{idx}*hk*`）。
  - 12B：TB 用 `tb_file_map.json` 显式指定 `hk_file`；X-ray 用 `(kv, ch)` glob；可传 `hk_path`。

### 3.2 TB 单文件流程（所有版本一致，仅 reader/bin/range 不同）

```text
process(op.tb, file, fp_method=None)
→ __get_fp05B → File_operation_05b(path, read_config, bkg_read_config="", spectrum_config, fit_config)
    ├─ read_out()   → __read(read_config)  → reader → (sci, tel)
    ├─ get_spectrum() → corr = [lambda t,b:1]*4 (无 TV 校正); count histogram; 无 bkg
    ├─ peak_fit()   → util.peak_fit(..., bkg_form=fit_config.bkg_form="lin")
    └─ save()       → pickle TB_fit_result/{stem}.pickle
```

- TB 单文件峰型固定 Gaussian，背景固定 `"lin"`（`Fit_config` 默认 `bkg_form="lin"`，TB 从不为每个文件覆盖）。
- **TB 无 TV/temp-bias 校正**（`corr = [1]*4`）。TB 的 temp/bias 依赖由全局二维拟合建模，而不是在单文件层校正。

### 3.3 TB 全局拟合（`temp_bias_fit`）

```text
load_data(): 遍历 TB_fit_result/{stem}.pickle
  → fit["b"], fit["b_err"]（peak center）
  → tel 中 len==4 的字段取 ch i → temp/bias 的 mean/std
  → data = [center, center_err, temp_mean, temp_std, bias_mean, bias_std]
temp_bias_fit(data_all):
  → util.temp_bias_fit_curvefit(center, center_err, temp, bias, p0, maxfev)
     （11B 用 util.temp_bias_lmfit）
  → 每通道出 temp_bias_fit_{ich}.png + 写 tb_logs/{ts}_temp_bias_fit.json
```

- 默认 `TB_FIT_P0=None, TB_FIT_MAXFEV=10000`。
- **12B 覆盖**：`TB_FIT_P0=[-0.02,0.07,24.4,-35.0,-1000.0]`，`TB_FIT_MAXFEV=100000`；`load_data()` 额外 `[data[data[:,4]>=27.25] for ...]`（仅保留 bias≥27.5V 行参与拟合）。

### 3.4 EC 单文件流程

```text
EC src:  process(op.ec.src, file, fp_method=None) → __get_fp05B → File_operation_05b
EC xray: process(op.ec.x,  file, fp_method=<ver>)  → __get_fp03B → 每通道一个 File_operation_05b + __dict_4ch_reconstruct
```

- `File_operation_05b`（EC branch）的 `get_spectrum()`：**这里有 TV/temp-bias 校正**。

### 3.5 TV correction（temp/bias 校正）真实位置

**处于 spectrum（单文件）层**：`file_lib.get_spectrum()` 内部，在 `count` histogram 之前对 amp 乘以因子：

```python
corr = [c_f(np.mean(temp), np.mean(bias)) for c_f, temp, bias in zip(spectrum_config.corr, tel["tempSipm"], tel["bias"])]
amp  = [factor * a for factor, a in zip(corr, sci["amp"])]
```

其中 `spectrum_config.corr` 在 `EC_operation_*.__init__` 中构造：

```python
ref_temp, ref_bias = 25, 28.5
ref_func = [lambda t,b: tempbias2DFunction(t,b, c["G0"],c["k"],c["V0"],c["b"],c["c"]) for c in tb_result]
self.corr = [lambda t,b: f(ref_temp, ref_bias)/f(t,b) for f in ref_func]
```

即：以 `tb_result_path` 提供的 TB 全局拟合为参考（T=25°C, bias=28.5V），把每个 EC 单文件校正到该参考点。`tb_result_path` 在 config 里指向 `single_process/{ts}_temp_bias_fit.json`（**注意：TB 全局输出写在 `tb_logs`，而 EC 读取的是 `single_process` 下的一个独立参考文件**——二者数值在 09 上完全相等，但路径是解耦的，属于历史耦合点，需记录）。

- TB 分支：无校正（corr=1）。
- 05B 与其它版本差异在 bkg：05B xray 的 bkg 是**同一文件**用 rotated `bkg_time_cut`（`{k:[v[1],v[2],v[3],v[0]]}`）再读一次；其它版本的 xray bkg 是**相邻通道文件**（`read_config[1:4]+read_config[0]` 或 `[ch1,ch2,ch0,ch0]`）。

### 3.6 EC 全局拟合（`ec_fit`）

三处不同实现：

1. `VersionProcessOp.ecfit` + `EC_operation_*.ec_fit`（03B/04/05B/07/09）：`util.get_fit_dict(save_path, energy, suffix)` 扫描 `EC_fit_result` 目录，用三种正则（`^(\d+)p(\d+)` / `[xX][mM]\S+observe` / `src\S+`）+ fallback 映射 pickle→能量。
2. `VersionProcessOp10B.ecfit`（10B/11B）：按 `energy.keys()` 的 `observe` 关键字直接分组读 pickle。
3. `VersionProcessOp12B.ecfit`（12B）：按 `k.isdigit()` 区分 xray(`kV`) 与 src。

`ec_fit` 主体（`EC_operation_05B.ec_fit`）：
```text
result = src_result + x_result;  energy = src_energy + x_energy
按 energy 排序
q_low = energy < split_low;  q_high = energy >= split_high
每通道:
  center_fit: np.polyfit(center, energy, deg=2, w=1/center_err)   # EC_low / EC_high
  resolution_fit: 分版本 polyfit / ExprFit / lmfit               # resolution_low/high
写 ec_logs/{ts}_ec_coef_sci_ch{i}.json  + ec_data_ch{i}.npy  + plot.ec_plot
```

- **10B/11B**：`CHN_NUM=3`，仅打包 ch0/1/2，然后 `center=[c0,c1,c2,c0]`、`result=[r0,r1,r2,r0]` 等补丁成 4 通道（`ec_plot` 需要 4 通道）。这是历史行为，需保留。

### 3.7 version 特殊 selection 汇总（`files`）

| Version | TB files | EC x files | EC src files |
|---|---|---|---|
| 03B | `rundata` in & not `baseline` & not `CI` & not `50C` | `_rundata` in & not `CI` → `x_list=set(f.split("_")[1])` | 硬编码 `src_Na22...` / `src_Am241...` / `src_Cs137...`（3 个）+ 硬编码 bkg |
| 04 | `.txt` | `_ch` in & not `_18p0_` → `set(f.split("_")[3])` | 硬编码 4 个 src + 4 个 bkg |
| 05B | `rundata` in & not `50C` | `_observe.dat` in & not `XM_22` | `rundata` in & not `bkg` |
| 07 | `.txt` | `_ch` in & not incl `40p0/15p0/12p0/99p9/90p1` → `set(f.split("_")[2])` | `src` in & not `_bk_` & not `bkg` |
| 09 | `.txt` **且 remove 6 个指定文件** | `_ch` in & not `65keV_` → `set(f.split("_")[0])` | 硬编码 3 个 src + 3 个 bkg |
| 10B | `observe` in & not `50C_265` | `_ch` in & not `hk` → `set(f.split("_")[2])` | 硬编码 4 个 src + 4 个 bkg（`073/085/089/077_observe_*.dat`） |
| 11B | glob 4 子目录 + 显式 remove + 过滤 `_50_Cs_2` | `_ch` in & not `hk` → `set(f.split("_")[2])` & remove `20` | 硬编码 4 个 src + 4 个 bkg（Ba133 bkg 为 `""`） |
| 12B | `tb_file_map.json`，EXCLUDE `(-20,275)/(-20,285)` | `.dat` & `_ch` & not `old` → `set(f.split("_")[1])`，再过滤 hk 完整 + fit_range 完整 | 硬编码 4 个 src，共享 bkg `0611env.dat` |

### 3.8 reader/ending 对应关系（`file_lib.__read`）

| ending | reader | 说明 |
|---|---|---|
| `normal` | `single_read05b_normal(path)` | 05B TB / src |
| `xray` | `single_read05b_xray(path,config_file,time_cut)` | 05B xray |
| `03b` | `single_read03b(path,config_file)` | 03B TB / xray |
| `03b-src` | `src_read03b(path,config_file)` | 03B src |
| `04` | `single_read04(path)` | 04 全部 |
| `07` | `single_read07(path)` | 07 全部 |
| `09` | `single_read09(path)` | 09 全部（= `reader07.read` 的独立 cache 别名） |
| `10b` | `single_read10(path)` | 10B 全部 |
| `11b` | `single_read11(path, mode=kwarg["mode"])`（默认 wf） | 11B 全部 |
| `12b` | `single_read12(path, **kwarg)`（默认 ft） | 12B 全部 |

---

## 4. 已知的历史耦合点 / 疑似问题（记录，不修）

1. **EC 参考耦合**：TB 全局输出写到 `tb_logs/{ts}_temp_bias_fit.json`，但 EC 的 `corr` 读取 `single_process/{ts}_temp_bias_fit.json`（config 里 `tb_result_path`）。二者路径解耦；数值一致（09 上已核对，4/4 通道 MATCH）。`UNKNOWN`：这是历史人工约定的产物，无法从代码判断是否应合并。
2. **`09` 的 TB 文件顺序**：`self.files = [splitext==.txt]` 后再 `remove` 6 个文件，依赖 `os.listdir` 顺序；`remove` 的 6 个文件在未来目录内容变化时可能 `ValueError`。
3. **`09` 的 `fit_range.get` / `bkg_form.get` 兜底**：EC xray_config 用 `fit_range.get(energy_name, [[None,None]]*4)` 与 `bkg_form.get(energy_name,"lin")`，其余版本用直接下标。这属于版本差异兜底，历史行为。
4. **`04` / `09` 的 `src_list` 硬编码为 .txt，`03B` 为 .dat**：各版本源极其不同，无法统一，除非使用 manifest。
5. **`10B/11B` EC 隐藏 3 通道**：`CHN_NUM=3` + `center/result/src_result/x_result` 补丁为 4 通道。下游 `ec_plot` 需要 4 通道，但科学上可能只有 3 通道有效。`UNKNOWN`：真实物理意图无法从代码判断。
6. **`process()` 中 `fp_method` 与 `ending` 两套命名**：易混淆，是历史积累产物。

---

## 5. 审计结论（Gate A）

- 所有版本共享**单一**单文件科学内核（`File_operation_05b`），差异全部集中在：
  1. reader 选择（ending）；
  2. 文件 selection 与 order（扫描规则 / 硬编码列表）；
  3. EC X-ray 的 4 通道重建 vs 单 4 通道文件（fp05B vs fp03B）；
  4. EC 的「bkg 来自相邻通道文件 vs 同一文件 rotated time_cut」；
  5. TB 全局拟合函数（curvefit vs lmfit；p0/maxfev 覆盖；12B 的 bias 过滤）；
  6. EC 全局拟合实现（3 种 ecfit）+ 通道数（4 vs 3+补丁）；
  7. `resolution_fit` 每版本不同（polyfit / ExprFit / lmfit）。

这正是本次重构要"显式化"的差异集合。
