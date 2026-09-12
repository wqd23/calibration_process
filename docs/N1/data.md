# N1（GRIDN1）数据说明

修订日期：2026-09-13

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
- X 光机（分两部分，方向不同）：
  - `260202计量院标定`：40/45/47/49/51/53/55/60/65/70/75/80/90 keV，ch0–3 齐全，
    是 `{E}-ch{N}-{idx}.event.dat` 的 **wf 512** 包；接在 `ec.x_path`。
  - `260129计量院标定`：15/20/25/30/35 keV，**只有 ch1/ch2/ch3（无 ch0）**，
    而且同名前缀是 **ft** 包（`grid1x_ft_packet`，与 260202 的 wf 不同）；接在
    `ec.x_path_low`，并用 `ec.xray_reader_low: n1` 指定 ft reader。
  - 文件命名为 `{E}-ch{N}-{idx}.event.dat`（E = 管压 keV）。`_start`、`.cut`、
    `CI`、无能量前缀（如 `087.event.dat`）的文件都排除。

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

## EC：GAGG 分段二次，CLYC 单条二次

N1 的 E-C 按通道形式不同（`ec.ec_form`）：

- **ch1/ch2（GAGG）**：按 Gd K 边（50.2 keV）把能量轴拆成两段，各做一条二次，
  拆分点为 `energy_split_low=49.0`、`energy_split_high=55.0`；落在死区
  `[49,55)` 的 X 光点（49/51/53）不参与拟合，但仍做单谱拟合留档。
- **ch0/ch3（CLYC）**：没有 Gd K 边，用**单条二次**（`ec_form: "quadratic"`）。

X 光逐通道准直，背景用旋转 `[1,2,0,0]`（`xray_bkg_rotation: fixed`）；源文件不减
本底。X 光各通道来自不同文件，速率按**逐通道时间跨度**计算（`rate_span: channel`），
否则本底扣减会出负值。

低能 X 光（260129，15–35 keV，ft 包）已接入 manifest 与 `energy_map`，但**本轮不
参与 E-C 拟合**：它的峰落在阈值附近，跨通道本底扣减在阈值处留下一个很大的负凹陷，
ch2/ch3 的谱被噪声主导，单谱拟合中心不稳定（redchi 常在 10 以上、σ 撞边界）。
按"宁缺毋滥"，这些点在 `fit_range_ec_xray.yaml` 里置 null，低段只保留 40/45/47 三个
锚点（三点正好定死一条二次）。

放射源区间这次重开了窗：`fit_range_ec_source.yaml` 里 ch1/ch2 的旧窗口是按更高的
增益设的，与实际峰位不符（662 keV 峰落在窗口下沿、拟合中心撞边界）；现按实测谱
重开 Cs137（662）、Co60（1173/1332）、Th228（583）的 ch1/ch2 窗口，ch0/ch3 不变。

## 中子：TB 简并，仅作自洽基准

中子数据固定 28.5V，没有偏压自由度，且谱形非高斯，5 参数二维面拟合是**简并**的
（`redchi=inf`、`ndf=0`）；`0degC-28.5-191` 还没有时间对齐的 HK。因此
`GRIDN1/Neutron` 的 TB 结果只作为**自洽回归基准**（`tests/golden/GRIDN1/Neutron/`），
不是物理标定；中子的粒子甄别（PSD）另见版本模块 `vN1.selection` 的占位钩子。

## 拟合区间来源与坏点

`GRIDN1/{GAGG,CLYC}/fit_range_tb.yaml` 由 `gridN_cali/data/TB_{GAGG,CLYC}/L1_cfg/*_fit_cfg.yaml`
转出（人工区间，`can_fit: false` 的通道置 null），共 107 个点位（48 GAGG + 59 CLYC，
含 11 个 CLYC 补测段）。`GRIDN1/EC` 的区间来自 `ec_fit_range.json`，源/X 光分别
存于 `fit_range_ec_source.yaml` / `fit_range_ec_xray.yaml`（X 光含新增的 49/51/53
和 15–35；后者置 null，见上）。

已知坏点（数据本身问题）：

- **CLYC m20C 扫描的 265/270/275 段 ch0/ch3**：低增益通道在最低偏压下峰无法区分
  （ch0 峰位非单调），已在 `fit_range_tb.yaml` 置 null。
- GAGG 轮（Am241）只针对 ch1/ch2，ch0/ch3 区间多为 null，二维面也只做 ch1/ch2。

## 处理结果（2026-09-13，新架构）

- TB：两数据集共 107 点位。2026-09-13 关闭 `skip_qa_fail`（原先会默认丢掉
  `qa_flag=fail` 的点），改为全部点先做单谱拟合、QA 报告出来后再逐通道人工排除。
  本次排除（写在 manifest 的 `channels.{ch}.use:false`）：GAGG `0C/10C/20C_265` 的
  ch1/ch2、CLYC `10C_265` 的 ch1/ch2、`0C_270` 的 ch0、`m10C_265` 的 ch0。
  排除后二维面相对残差 max < 5%：GAGG ch1/ch2 为 4.1%/2.8%，CLYC ch0–3 为
  4.6%/4.3%/4.0%/3.6%。
- EC：源 6 点 + X 光 13 个能量（40/45/47/49/51/53/55/60/65/70/75/80/90）进入单谱
  拟合，其中 49/51/53 落在死区、不参与 E-C 拟合。ch1/ch2 分段二次、ch0/ch3 单条
  二次。四通道 E-C 相对偏差 max < 5%（ch0 2.9%、ch1 2.6%、ch2 4.4%、ch3 3.8%）。
  ch1/ch2 的低段只有 40/45/47 三个点，正好定死一条二次；`_center_fit` 对这种
  "点数 = 阶数+1" 的情形退化为不带协方差的插值拟合（误差棒记 0，系数不受影响）。
  N1 没有 legacy EC 系数可对，按物理合理性验收，并作为自洽 golden
  （`tests/golden/GRIDN1/EC/`）。
- 中子：20 点自洽 TB 快照（见上）。
