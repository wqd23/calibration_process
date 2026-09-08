# 各载荷数据说明（除 12B/13B，简化版）

修订日期：2026-09-08

本文是 **12B/13B 之外**各载荷的数据说明。

- 12B/13B 由于其数据最特殊（同一 observe 文件被两个点位共用、大量
  补测/重测/逐字节复制、HK 带爬升与回温段、点位需文件↔条件对照表），单独建了
  [12B_13B/data.md](12B_13B/data.md)。
- 其余 7 个载荷的数据结构**高度一致**：都是「TB 目录 + X光机目录 + 放射源目录」
  三块，每个**科学文件对应一个实验条件**，靠文件名过滤选点，目录里混有大量
  **辅助文件**（HK / TimeLine / SciConfig / ivscan / 其它未用实验），由选择规则
  跳过。因此本文用同一套模板概括，再给每个版本一个简表。**未在工程里核实的
  数据级"异常"不写，标 UNKNOWN 或不写。**

> 选择规则由 `workflows/versions/v{ver}.py` 的 `enumerate_measurements` 决定，
> 与本文一一对应；配置在 `configs/{ver}/*.yaml`。

## 通用模板

```text
data/{ver}/raw_data/
├── <TB 目录>      # 每个 (温度,偏压) 点位一个科学文件；同目录有 *_HK / *_TimeLine / *_SciConfig / ivscan 等辅助文件
├── <X光机目录>    # 每个管压一个科学文件；同目录有 *_HK / *_TimeLine / ivscan 等辅助文件
└── <放射源目录>   # 每个源一个科学文件 + 本底文件；同目录有 IV / CI_on / ivscan 等未用实验
```

- 科学文件是真正被读取的；辅助文件（HK/TimeLine/SciConfig/ivscan/CI/IV）由
  选择规则排除，**不进 manifest**。
- 每个版本在自己的 `v{ver}.py` 里用**文件名关键字**过滤出科学文件，并把需要的
  本底（源）用硬编码列表或按名字动态配对。
- 说明里的"排除"指**不会出现在该分支**的点/文件；若已在 manifest 里被标记
  `use:false`（如 10B 的 90 keV），表示该点整体不参与。

---

## 03B

- 目录：TB `20210501_tempbias_03B`、X光机 `20210429_Xray_03B`、源 `20210504_source_03B`。
- reader：`03b`（EC 源用 `03b-src`）；`.dat`。
- TB：文件名含 `rundata` 且**不含** `baseline`/`CI`/`50C`；bin_width 6、adc_max 16384；
  二维面共享 curvefit。56 个点位（8 个通道拟合区间为 None，即该通道不拟合）。
- X光机：`jly_{能量}p{小数}_ch{通道}_30s_rundata*.dat`；文件名含 `_rundata` 且
  不含 `CI`；能量名取第 2 段（`split("_")[1]`）；**circle** 背景轮转；resolution
  `lmfit`；拆分 49/51 keV。8 个能量点。
- 放射源：3 个硬编码源（Na22/Am241/Cs137）+ 硬编码本底；目录里另有
  `VIresult*.dat`/`Vbrresult*.dat` 等未用文件。resolution `lmfit`。

## 04

- 目录：TB `20210501_tempbias_Am241_GRID04`、X光机 `20210429_Xray_GRID04/data_GRID`、
  源 `20210504_source_GRID04`。
- reader：`04`；hex/ASCII `.txt`。
- TB：所有 `.txt`；bin_width 6、adc_max 65535；二维面共享 curvefit。54 个点位。
- X光机：`210429..._COM12_jly_{能量}p{小数}_ch{通道}_30s-Data.txt`；`_ch` 且不含
  `_18p0_`；能量名取第 4 段（`split("_")[3]`）；**circle** 背景轮转；resolution
  `ExprFit`；拆分 49/55 keV。7 个能量点。
- 放射源：4 个硬编码源（Co60/Na22/Cs137/Am241）+ 硬编码本底；目录里另有
  `COM6_src_IV-Data.txt`/`CI_on...` 等未用文件。resolution `ExprFit`。

## 05B

- 目录：TB `温度偏压实验`、X光机 `X光机实验-天格_reset`、源 `放射源实验`。
- reader：TB/源用 `normal`，X 光机用 `xray`（需 `x_config` + `time_cut`）；`.dat`。
- TB：`rundata` 且不含 `50C`；bin_width 6、adc_max 16384；二维面共享 curvefit。
  56 个点位（14 个通道 None）。
- X光机：**每个管压一个 4 通道文件**（不是每通道一个文件），文件名
  `XM_{管压}_{序号}_observe.dat`，同目录 `_observeHK`/`_observeTimeLine`；
  含 `_observe.dat` 且不含 `XM_22`；**单文件**路径 + **循环移位时间窗**背景
  （`bkg_time_cut`）再读一次；resolution polyfit；拆分 49/52 keV。16 个点。
- 放射源：`src_{源}_{距离}_rundata*_dat`（含 `rundata` 且不含 `bkg`）；本底按
  `src_bkg` 与源名配对。目录含 HK/SciConfig/TimeLine 等辅助文件。resolution polyfit。

## 07

- 目录：TB `北师大正样温度偏压标定数据-20211124`、X光机 `X光机标定实验/正样`、
  源 `北师大豁免源标定数据-20211127`。
- reader：`07`；hex/ASCII `.txt`。
- TB：所有 `.txt`；bin_width 6、adc_max 65535；二维面共享 curvefit。56 个点位
  （1 个通道 None）。
- X光机：`{时间}_{代码}_bnu01_{能量}p{小数}_2m_ch{通道}_COM3-Data.txt`；`_ch` 且
  不含 `40p0/15p0/12p0/99p9/90p1`；能量名取第 3 段（`split("_")[2]`）；**circle**
  背景轮转；resolution polyfit；拆分 49/55 keV。12 个能量点。
- 放射源：文件名含 `src` 且不含 `_bk_`/`bkg`；本底按 `bkg{源}` **动态**配对
  （如 `src_Cs137...` ↔ `src_bkgCs...`）；另有 `src_bkgCo...`。resolution polyfit。

## 09

- 目录：TB `temp_bias`、X光机 `Xray`、源 `src`。
- reader：`09`（= reader07 的独立缓存别名）；hex/ASCII `.txt`。
- TB：所有 `.txt`，**显式剔除 6 个文件**（见 `v09.TB_EXCLUDE`，含 `0825_0C_290`、
  几个 `0826_30C_265`、一个带空格的 `(1)` 文件）；bin_width 6、adc_max 65535；
  二维面共享 curvefit。48 个点位。
- X光机：`{能量}keV_ch{通道}_0x{编号}_2m.txt`；`_ch` 且不含 `65keV`；能量名取
  第 1 段（`split("_")[0]`）；**circle** 背景轮转；resolution polyfit；拆分
  49/55 keV。13 个能量点。
- 放射源：3 个硬编码源（Na22/Co60/Cs137，能量 511/1332/662 keV）+ 硬编码本底
  （Na22_bkg、以及 Co60/Cs137 **共用** `0827_10C_285_bkg_20m_0x00CF.txt`）。
  resolution polyfit。
- 注意：09 的 `fit_range`/`bkg_form` 用 `.get(...)` 带默认值兜底（历史行为）。

## 10B

- 目录：TB `tb_data`、X光机 `x_data`、源 `src_data`。
- reader：`10b`（新包格式 `grid1x_wf_packet`，science+HK 自动配对）；`.dat`。
- TB：文件名含 `observe` 且不含 `50C_265`；bin_width 6、adc_max 16384；二维面共享
  curvefit。63 个点位（14 个通道 None）。
- X光机：`{序号}_observe_{能量}_ch{通道}.dat`，同目录 `{序号}_hk_...`/ivscan；
  `_ch` 且不含 `hk`；能量名取第 3 段（`split("_")[2]`）；**fixed** 背景轮转
  `[ch1,ch2,ch0,ch0]`；resolution polyfit；拆分 49/55 keV。21 个点。
  **单一 90 keV 点被标记 `use:false`**：该点 ch3 的 `gaus` 背景拟合失败，旧版
  绘图对失败结果取 `["a"]` 抛 `KeyError`，因此**旧版本就处理不了这个点**。新
  版本沿用同一个科学内核与绘图，行为一致，仍无法产生该点；作为"历史不可处理点"
  显式排除（`use:false`），不做额外处理。
- 放射源：4 个硬编码源（Cs137/Na22/Am241/Co60）+ 硬编码本底；目录含
  `{序号}_observe_CI_off/on.dat`、`ivscan` 等未用文件。
- **EC 物理 3 通道**（`channel_count=3`）：单拟合仍读 4 通道（fp03B 重建），
  但 EC 全局只取 ch0/1/2，绘图把 ch3 复制成 ch0。

## 11B

- 目录：TB `tempbias`（其下 4 个温度偏压子目录）、X光机 `GRID-11B_x_data`、
  源 `src/能量道址实验——1.15/数据/`。注意 TB 目录下有 `__MACOSX` 目录（跳过）。
- reader：`11b`（`grid1x_wf_packet`，HK 按 stem 配对）；`.dat`。
- TB：**glob 4 个子目录**的 `*_observe*.dat`，排除路径含 `on`/`off` 的文件，
  **显式剔除若干文件**（`v11B.TB_REMOVE`，含文件名带尾部空格、重测等），并过滤
  `_50_Cs_2`；bin_width 6、adc_max 16384；二维面用 **lmfit**。54 个点位
  （3 个通道 None；fit_range 按文件 **stem** 键控）。
- X光机：`{序号}_observe_{能量}_ch{通道}.dat`；`_ch` 且不含 `hk`；能量名取第 3 段，
  **剔除能量 `20`**；**fixed** 背景轮转；resolution polyfit；拆分 49/55 keV。
  21 个点。
- 放射源：4 个硬编码源（Ba133 **无背景**、Cs137、Eu152、Co60），Ba133 的
  `aux` 为空；其余有本底文件。目录含 `ivscan`/`chargeon` 等未用文件。
- **EC 物理 3 通道**（同 10B，绘图层补齐）。

---

## 对照小结

- 除 12B/13B 外，其余载荷都是**一科学文件 = 一条件**，靠文件名关键字选点；
  它们的"数据复杂"都体现在**目录里混有大量辅助/未用文件**（HK/TimeLine/
  SciConfig/ivscan/CI/IV）与**少量排除点/None 通道/3 通道**，而非 12B 那种
  "文件↔条件需查表"。
- 如果你发现某个版本的数据确实需要逐文件核对的映射表（类似 12B），再把它单独
  拆成一个子目录详述；目前**没有**。
