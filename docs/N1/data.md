# N1（GRIDN1）数据说明

修订日期：2026-09-12

本文记录 GRIDN1 载荷标定数据的来源、文件组织与各文件特殊情况。GRIDN1 是
中子/伽马双探测器载荷（GAGG 晶体通道测伽马、CLYC 晶体通道测中子），与 GRID B
系列同族但包格式和数据组织不同。数据处理方法论见 [../data.md](../data.md)；
本载荷的早期处理流水线在另一个仓库 `gridN_cali`（分层 L0/L1/L2 + Snakemake），
本仓库的包定义 XML 与人工拟合区间来自它，但**直接从原始二进制文件处理，不依赖
它的中间产物**。

## 目录与版本划分

GRIDN1 统一收在 `configs/GRIDN1/` 与 `data/GRIDN1/` 下，按组成部分分目录，
版本名即路径：

| 子版本 | 分支 | 数据软链 | 说明 |
|--------|------|----------|------|
| `GRIDN1/GAGG` | `tb` | `raw_data -> 260303温度偏压/GAGG` | Am241，ft 包，ch1/2 |
| `GRIDN1/CLYC` | `tb` | `raw_data -> 260303温度偏压/CLYC` | Na22，wf 512 包，ch0–3 |
| `GRIDN1/EC` | `ec_source`,`ec_xray` | `ec_src -> 260326放射源/src`、`ec_xray -> 260202计量院标定` | 能量标定 |
| `GRIDN1/Neutron` | `tb`,`neutron` | `raw_data -> 260322中子束流` | 中子束流（256 点 wf） |

原始数据在 `/home/wqd/cali_data/GRIDN1/data/`。命令示例：
`calib check GRIDN1/GAGG`、`calib all GRIDN1/CLYC --until L3`、
`calib global GRIDN1/EC ec`。

## 数据来源与组织

### 温度偏压（TB）——`260303温度偏压/`

- `GAGG/`：**Am241 源**（59.5 keV），结果只针对 GAGG 通道 ch1/ch2。负温度
  （m20C/m10C）是每个偏压一个文件（`m20C-265-125.event.dat` = −20°C、26.5V），
  非负温度（0C/10C/20C/30C）是**偏压扫描文件**（`0C-265-290-139.event.dat`，
  一个文件依次测 8 个偏压点）。
- `CLYC/`：**Na22 源**（511 keV），四个通道全部参与拟合。全部 6 个温度都是
  偏压扫描文件，部分扫描末尾有补测偏压段（如 m20C 多 277/279/281 三个回测点）。

排除文件：名字含 `To`（温度转换段）、`CI`（电流扫描）、`test`（测试轮）。
`260302温度偏压测试/` 是更早一轮，未使用。

### 中子束流——`260322中子束流/`

固定偏压 **28.5V 的温度扫描**，命名为 `{温度档}-28.5-{idx}.event.dat`
（如 `0degC-28.5-191`、`0-75C-28.5-193`）。文件覆盖一段温度区间，实测 SiPM
温度约 26–27°C（HK 的 `sipm_temp`），偏压恒定。

### EC——源与 X 光机

- 源：`260326放射源/src`（Cs137/Co60/Th228 等 ft 文件，每个含四通道）。
- X 光：`260202计量院标定`，逐通道准直文件 `{E}-ch{N}-{idx}.event.dat`
  （E = 管压 keV），能量 40/45/47/55/60/65/70/75/80/90。

## 数据格式（与 12B 的关键差别）

- GAGG 是 584 字节 `grid1x_ft_packet`（38 事件/包、每事件 14 字节，比 12B 多一个
  `data_ccm`），解析用 `grid_packet.xml`、`multi_evt=38, multi_step=14`。
- CLYC 是 1080 字节 `grid1x_wf_packet`（512 点波形，`grid_packet_512wave.xml`），
  用 `parse_grid_data_iter` 分块解析并滤掉 CRC 错误事件（reader 模式 `wf`）。
- 中子束流是 568 字节 `grid1x_wf_packet`（256 点波形），reader 模式 `wf256`。
- 三种包的 `data_max/data_base/data_sum` 都在，幅度均按
  `amp = data_max − data_base/4` 计算（注意 `/4`，与旧 A 族的 `/1` 不同）。
- HK 是 178 字节 `hk_packet`（定义在 `yingtian_packet.xml`；与 12B 的 187 字节
  `grid1x_hk_packet` 头魔数相同但字段布局不同）。
- **N1 科学包带真实 utc**（12B 地面数据 utc 全 0），所以 sci 与 hk 按最近 utc
  对齐，扫描文件按"utc 间隔大于 4 秒的安全间隙"切成偏压段，用 `seg_bias` 按
  实测偏压选段（0.1V 容差）。
- 包定义 XML 复制自 `gridN_cali/reader/`，解析引擎与 reader11/12 相同。

## TB：两个数据集分开拟合

两轮用不同放射源（GAGG 轮 Am241 59.5 keV，CLYC 轮 Na22 511 keV），同一通道在
两轮拟合的是不同能量的峰，峰位 ADC 差约 9–10 倍，不能放进同一个二维面。新架构
用两个子版本分别做二维面拟合，再用 `scripts/merge_tb.py` 按通道归并成标准 4 元
`temp_bias_fit.json`：

- `GRIDN1/GAGG` 拟合 ch1/ch2（`tb.channels=[1,2]`）；
- `GRIDN1/CLYC` 拟合 ch0–3（`tb.channels=[0,1,2,3]`）；
- `GRIDN1/EC` 的 `payload.ec.tb_ref_path` 指向合并文件。

通道归属：ch0/ch3 取 CLYC、ch1/ch2 取 GAGG（即 legacy 的 `CHANNEL_DS`）。逐通道
初值 `tb_fit_p0_by_channel` 取自 `gridN_cali` 的 L2 结果。

Na22 的 1274 keV 峰评估过、放弃：这批数据采集幅度在 ~8700 ADC 截断（各通道均有
溢出尖峰），1274 峰位（511×2.494）在 28.3V 以上全部越界；ch0/ch3 谱被噪声主导
也看不到 1274 峰。

## EC：单条二次

N1 的 E-C 不做 Gd K 边拆段，用**单条二次**（`ec_form: "quadratic"`）：ch1/ch2
（GAGG）由放射源高能线与 X 光机低中能点共同拟合；ch0/ch3（CLYC）没有 Gd K 边，
不套用 B 方案的 K 边假设。X 光逐通道准直，背景用旋转 `[1,2,0,0]`
（`xray_bkg_rotation: fixed`）；源文件不减本底。X 光各通道来自不同文件，速率按
**逐通道时间跨度**计算（`rate_span: channel`），否则本底扣减会出负值。

## 中子：TB 简并，仅作自洽基准

中子数据固定 28.5V，没有偏压自由度，且谱形非高斯，5 参数二维面拟合是**简并**的
（`redchi=inf`、`ndf=0`）；`0degC-28.5-191` 还没有时间对齐的 HK。因此
`GRIDN1/Neutron` 的 TB 结果只作为**自洽回归基准**（`tests/golden/GRIDN1/Neutron/`），
不是物理标定；中子的粒子甄别（PSD）另见版本模块 `vN1.selection` 的占位钩子。

## 拟合区间来源与坏点

`GRIDN1/{GAGG,CLYC}/fit_range_tb.yaml` 由 `gridN_cali/data/TB_{GAGG,CLYC}/L1_cfg/*_fit_cfg.yaml`
转出（人工区间，`can_fit: false` 的通道置 null），共 107 个点位（48 GAGG + 59 CLYC，
含 11 个 CLYC 补测段）。`GRIDN1/EC` 的区间来自 `ec_fit_range.json`，源/X 光分别
存于 `fit_range_ec_source.yaml` / `fit_range_ec_xray.yaml`。

已知坏点（数据本身问题）：

- **CLYC m20C 扫描的 265/270/275 段 ch0/ch3**：低增益通道在最低偏压下峰无法区分
  （ch0 峰位非单调），已在 `fit_range_tb.yaml` 置 null。
- GAGG 轮（Am241）只针对 ch1/ch2，ch0/ch3 区间多为 null，二维面也只做 ch1/ch2。

## 处理结果（2026-09-12，新架构）

- TB：两数据集共 107 点位；`merge_tb.py` 合并后与 `gridN_cali` 参考
  `20260826185532_temp_bias_fit.json` 对比，**GAGG ch1/ch2 系数逐位一致（0.00%）**，
  **CLYC ch0/ch3 在 ±1.6% 内**（k/V0≈0；二维面在标定范围内差异 0.0000%，属简并
  参数化差异）。点数与 `gridN_cali` 表一致（GAGG 41/43，CLYC 53/58/58/55）。
- EC：源 6 点 + X 光 10 点，四通道单二次系数有限；`gridN_cali`/`feat/gridN1` 未提交
  N1 EC 系数，按物理合理性验收（ch1 662 keV 处约 737 keV）。
- 中子：20 点自洽 TB 快照（见上）。
