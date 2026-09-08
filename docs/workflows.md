# 各载荷的显式 workflow

修订日期：2026-09-08

本文描述当前代码里**每个载荷版本**实际执行的 workflow。它们都由
`src/calibration_process/workflows/versions/v{ver}.py` 显式定义（文件选择规则、
reader、背景轮转、分辨率拟合方法、等等），并复用同一个未改动的科学内核。
若某一版本与其它版本有差异，都能从该版本自身的 workflow 文件里直接读出来。

> 完整的历史审计矩阵见 [workflow_matrix.md](workflow_matrix.md)；
> 科学内核（单谱拟合 / TB 二维面 / EC 分段拟合 / 绘图）见 [data.md](data.md) 与 [results.md](results.md)。

概括几个跨版本的概念：

- **TB**：温度-偏压标定。对每个（温度, 偏压）点位做单谱峰拟合，得到峰位，再对
  4 个通道各自做「峰位 ~ (温度, 偏压)」的二维面拟合。
- **EC-source**：放射源能量-道址标定。用已标定好的 TB 面把每个放射源数据修正到
  参考温压（默认 25°C / 28.5V），拟合全能峰，得到 能量↔道址 点。
- **EC-xray**：X 光机能量-道址标定。同上，对每个管压点拟合。
- **EC 全局**：把 E-C 点按 Gd K 吸收边（~50 keV）拆成低/高两段，各做二次拟合，
  同时做分辨率拟合。
- **TV 校正位置**：温度/偏压校正都发生在**单谱生成阶段**（`get_spectrum` 里对
  amp 乘因子），不发生在点/全局层。

---

## 03B

- Reader：`03b`；EC 放射源用 `03b-src`。
- TB：目录里 `rundata` 且在文件名中**不含** `baseline`/`CI`/`50C`，bin_width 6、
  adc_max 16384；二维面用共享 curvefit。
- EC-source：3 个硬编码源（Na22/Am241/Cs137）+ 硬编码本底；resolution 用 `lmfit`。
- EC-xray：`_rundata` 且不含 `CI`，能量名取文件名第 2 段（`split("_")[1]`）；
  **circle** 背景轮转（ch i 用 ch (i+1)），resolution `lmfit`；拆分 49/51 keV。

## 04

- Reader：`04`（hex/ASCII `.txt`）。
- TB：目录里所有 `.txt`，bin_width 6、adc_max 65535；二维面共享 curvefit。
- EC-source：4 个硬编码源（Co60/Na22/Cs137/Am241）+ 硬编码本底；resolution `ExprFit`。
- EC-xray：`_ch` 且不含 `_18p0_`，能量名取文件名第 4 段（`split("_")[3]`）；
  **circle** 背景轮转，resolution `ExprFit`；拆分 49/55 keV。

## 05B

- Reader：TB/源用 `normal`，X 光机用 `xray`（需要 `x_config` 与 `time_cut`）。
- TB：`rundata` 且不含 `50C`，bin_width 6、adc_max 16384；二维面共享 curvefit。
- EC-source：`rundata` 且不含 `bkg`，本底按源名与 `src_bkg` 配对；resolution 共享
  polyfit。
- EC-xray：**每个管压一个 4 通道文件**（不是每通道一个文件），文件名含
  `_observe.dat` 且不含 `XM_22`；**单文件**路径，背景是同文件按**循环移位的时间窗**
  （`bkg_time_cut`）再读一次；resolution polyfit；拆分 49/52 keV。
  - 这是唯一走"单 4 通道文件"的版本（`xray_single_file = true`）。

## 07

- Reader：`07`（hex/ASCII `.txt`）。
- TB：目录里所有 `.txt`，bin_width 6、adc_max 65535；二维面共享 curvefit。
- EC-source：文件名含 `src`，不含 `_bk_`/`bkg`；本底按 `bkg{源}` 动态配对；
  resolution polyfit。
- EC-xray：`_ch`，且不含 `40p0/15p0/12p0/99p9/90p1`，能量名取文件名第 3 段
  （`split("_")[2]`）；**circle** 背景轮转，resolution polyfit；拆分 49/55 keV。

## 09

- Reader：`09`（= reader07 的独立缓存别名，hex/ASCII `.txt`）。
- TB：目录里所有 `.txt`，**显式剔除 6 个文件**，bin_width 6、adc_max 65535；
  二维面共享 curvefit。
- EC-source：3 个硬编码源（Na22/Co60/Cs137）+ 硬编码本底；resolution polyfit。
- EC-xray：`_ch` 且不含 `65keV`，能量名取文件名第 1 段（`split("_")[0]`）；
  **circle** 背景轮转，resolution polyfit；拆分 49/55 keV。
  - 注意 09 的 `fit_range`/`bkg_form` 用 `.get(...)` 带默认值兜底（历史行为）。

## 10B

- Reader：`10b`（新包格式 `grid1x_wf_packet`，science+HK 自动配对）。
- TB：文件名含 `observe` 且不含 `50C_265`，bin_width 6、adc_max 16384；
  二维面共享 curvefit。
- EC-source：4 个硬编码源 + 硬编码本底；resolution polyfit。
- EC-xray：`_ch` 且不含 `hk`，能量名取文件名第 3 段（`split("_")[2]`）；
  **fixed** 背景轮转 `[ch1,ch2,ch0,ch0]`；resolution polyfit；拆分 49/55 keV。
- **EC 物理 3 通道**（`channel_count=3`）：单拟合仍读出 4 通道（fp03B 重建），
  但 EC 全局只取 ch0/1/2，绘图层把 ch3 复制成 ch0。

## 11B

- Reader：`11b`（`grid1x_wf_packet`，HK 按 stem 配对）。
- TB：**glob 4 个温度偏压子目录**的 `*_observe*.dat`，并**显式剔除若干文件**、
  过滤 `_50_Cs_2`；bin_width 6、adc_max 16384；二维面用 **lmfit**。
- EC-source：4 个硬编码源（Ba133 **无背景**）；resolution polyfit。
- EC-xray：`_ch` 且不含 `hk`，能量名取第 3 段，**剔除能量 `20`**；**fixed** 背景
  轮转；resolution polyfit；拆分 49/55 keV；3 通道 + 绘图补齐。
- 备注：TB 的 fit_range 按**文件 stem** 键控，而 `just list` 展示带扩展名的
  basename，manifest 用 `fit_key` 元数据衔接（详见 workflow_matrix.md）。

## 12B

- Reader：`12b`（`grid1x_ft_packet`，HK 按 stem 或 (kV,ch) 配对）。
- TB：文件集来自 **`tb_file_map.json`**（备份目录映射），**剔除
  (-20,275)/(-20,285)** 两个点位；两个点共用同一 observe 文件，用 `sci_half`
  按事件序 / `hk_bias` 按目标偏压切半；bin_width 6、adc_max 16384；二维面用
  **定制初值** `[-0.02,0.07,24.4,-35.0,-1000.0]`、maxfev 100000，且**仅用
  bias ≥ 27.5 V** 的点参与拟合（`load_data` 过滤）。
- EC-source：4 个硬编码源，**共用**环境本底 `0611env.dat`；resolution polyfit。
- EC-xray：`_ch` 且不含 `old`，能量名取第 2 段；要求 **HK 配对完整** 且
  **4 通道 fit_range 都完整** 才保留；**fixed** 背景轮转；resolution polyfit；
  拆分 49/55 keV。
- 逐通道 **None** fit range 表示该通道在该点不拟合（历史行为，保留）。

---

## 如何读一个版本的 workflow

1. `configs/{ver}/payload.yaml`：reader、bin_width、adc_max、channel_count、
   X 光机过滤/背景轮转、能量分界、温度参考、EC 的 TB 参考路径等**版本级**参数。
2. `configs/{ver}/analysis.yaml`：背景/峰型默认与逐点 override。
3. `configs/{ver}/fit_range_*.yaml`：逐 measurement 每通道拟合区间。
4. `configs/{ver}/*_manifest.yaml`：**人工确认**的 measurement 列表（science/hk/aux/
   metadata/use/channels）。
5. `workflows/versions/v{ver}.py`：**文件选择规则**（哪个文件进哪些分支）与版本特有
   逻辑；`enumerate_measurements` 复现了历史选择。
6. `workflows/common.py`：被所有版本共用的 stage（单拟合 / 点构建 / 全局拟合），
   直接调用未改动的科学内核。
