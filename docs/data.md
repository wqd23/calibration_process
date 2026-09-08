# 数据准备与目录约定

修订日期：2026-08-25

术语约定：**TB** 指温度-偏压标定（temperature-bias，测探测器响应对温度和
工作偏压的依赖）；**EC** 指能量-道址标定（energy-channel，建立 ADC 道址
与光子能量的对应关系）；**HK** 指工程遥测参数（housekeeping，含温度、
偏压等慢监测量）；**ADC** 指模数转换后的道址，即未做能量修正的原始幅度单位。

## 目录结构

配置已迁移为 **YAML + manifest**（在 `src/calibration_process/configs/{ver}/`），
新流程只读这些 YAML；`data/{ver}/single_process/` 下仍保留历史 JSON：
`fit_range.json` / `bkg_form.json` / `ec_energy.json` / `time_cut.json` /
`tb_file_map.json` / `{ts}_temp_bias_fit.json` 作为**迁移前的 oracle / 审计**（不再被新流程读，
但 `qa_thresholds.json` **仍是 QA 阈值的来源**，新流程 `common.py` 通过
`util.load_qa_thresholds` 读取）。

```
data/{ver}/                       # 数据目录（数据 + 产物）
├── raw_data -> /path/to/data     # 软链，指向原始标定数据
├── single_process/               # 单谱拟合产物 + 历史 JSON（oracle）
│   ├── fit_range.json ...        # （历史 oracle，新流程不再读）
│   ├── TB_fit_result/            # TB 单谱拟合产物（pickle）
│   ├── EC_fit_result/            # EC 单谱拟合产物（pickle）
│   └── single_fit_fig/           # 拟合图（png）
├── tb_logs/                      # TB 二维面拟合产物
└── ec_logs/                      # E-C 关系拟合产物

src/calibration_process/configs/{ver}/   # 新配置（唯一事实来源）
├── payload.yaml                 # 版本级科学/reader 参数
├── analysis.yaml                # 人类可编辑：背景/峰型默认与逐点 override
├── fit_range_tb.yaml / fit_range_ec_source.yaml / fit_range_ec_xray.yaml
├── tb_manifest.yaml / ec_source_manifest.yaml / ec_xray_manifest.yaml
```

`just init {ver} {path}` 负责**软链 raw_data + 创建输出目录**；`just check {ver}`
（`calib check`）校验新配置/manifest 层与数据链接。

## 配置文件

### payload.yaml

版本级参数：`reader`、`bin_width`、`adc_max`、`channel_count`（10B/11B 为 3）、
TB 的二维面初值/过滤（`tb_fit_p0`/`tb_fit_maxfev`/`bias_min_filter`/`tb_file_map`）、
EC 的各分支路径（`x_path`/`src_path`）、能量分界（`energy_split_low/high`）、
温度参考（`ref_temp`/`ref_bias`）、EC 的温度偏压参考（`tb_ref_path`）、
能量映射（`energy_map`）、X 光机的过滤/背景轮转（`xray_drop_*`、`xray_bkg_rotation`、
`xray_name_index`）、以及 05B 的单文件 X 光机标志（`xray_single_file` + `time_cut`）。

示例（09）：
```yaml
tb: {reader: "09", bin_width: 6, adc_max: 65535.0, science_dir: raw_data/temp_bias}
ec:
  reader: "09"
  energy_split_low: 49.0
  energy_split_high: 55.0
  xray_bkg_rotation: circle
  energy_map: {"15keV": 15.0, "0827_10C_285_Co60_20m_0x010F.txt": 1332.0}
```

### analysis.yaml

每个分支一份 `default_bkg`（默认背景）与 `background_overrides`（逐 measurement 的
背景/峰型覆盖）。可选值：`"lin"`/`"quad"`/`"exp"`/`"gaus"`/`null`。

### fit_range_*.yaml

逐 measurement 的每通道拟合区间。每个 measurement 是 4 元素列表（ch0–ch3），
每元素是 `[low, high]`（ADC 道址范围），`null` 表示该通道不拟合（12B 用到）。

```yaml
measurements:
  "0825_0C_265_4m_0x00C5.txt":
    - [235.2, 775.2]     # ch0
    - [235.2, 775.2]     # ch1
    - [235.2, 775.2]     # ch2
    - [235.2, 775.2]     # ch3
```

### *_manifest.yaml

**人工确认**的 measurement 列表。每个条目至少含 `id`、`branch`、`science_files`、
可选 `hk_files`/`aux_files`/`metadata`，以及 `use`（整条是否启用）和
`channels: {chN: {use: false}}`（单通道剔除）。

```yaml
measurements:
  - id: 0825_0C_265_4m_0x00C5.txt
    branch: tb
    science_files: [raw_data/temp_bias/0825_0C_265_4m_0x00C5.txt]
    use: true
```

> 这些 YAML 都是 **strict schema**：未知字段、非法枚举、重复 id、指向缺失文件
> 都会在 workflow 开始前报错（报错带 version/branch/measurement/字段）。

## 新版本接入流程

```bash
calib scaffold {ver} --data-dir /path/to/data   # 生成 configs/{ver}/*.yaml + 目录 + 软链 raw_data
# 1) 编辑 configs/{ver}/payload.yaml：reader / bin_width / adc_max / 分支路径 / channel_count
# 2) 若包格式不同，新增 reader：lib_reader/src/lib_reader/reader{ver}/，并在 lib_reader/__init__.py 注册
calib discover {ver} tb                        # 扫描 -> manifest 草稿（tb / ec_source / ec_xray 各一次）
# 3) 人工确认 manifest（去重、剔除坏点、设 use / channels.use）
# 4) 填 configs/{ver}/fit_range_tb.yaml / fit_range_ec_source.yaml / fit_range_ec_xray.yaml
calib check {ver}                               # 校验
calib all {ver}                                 # 全流程：单拟合 + TB/EC 全局拟合
```

> 若该版本流程与已有版本完全相同，可复用对应 `workflows/versions/v{ver}.py` 的
> `enumerate_measurements`，并在 `workflows/registry.py` 注册；若流程不同，在
> `workflows/versions/v{ver}.py` 显式写出差异（并更新 docs/workflows.md）。
> 每个载荷各自的 workflow 见 [workflows.md](workflows.md)；各载荷的**数据/选点**说明见
> [payloads_data.md](payloads_data.md)（除 12B/13B，见 [12B_13B/data.md](12B_13B/data.md)）。

### 接入细节与经验（以 12B 为例）

**1. 先判定数据包格式，再写 reader。**
用 `xxd 文件 | head` 看文件头字节，对照解析代码 `grid_packet.xml` 里各包的
`head` 属性确定包类型。注意同一系列不同代载荷的包定义可能不同：

- 12B/13B 的 `.dat` 是 `grid1x_ft_packet`（特征包，528 字节，41 事件/包，
  每事件 12 字节：timestamp 4B + data_max 2B + data_base 2B + data_sum 4B），
  而 11B 的 `.dat` 是 `grid1x_wf_packet`（波形包）；两者的 ft 包字段布局也不同
  （11B 的 timestamp 是 8 字节、无 data_sum）
- 12B/13B 的 `.hk` 是 187 字节的 `grid1x_hk_packet`，11B 是 139 字节的
  `hk_grid1x_packet`，sipm 字段偏移不同（82 vs 111）

最稳妥的做法是把**随数据附带的解析代码**（如数据目录里的
`python_Grid11B解析代码/`）中的 `grid_packet.xml` 整体复制进新 reader，
再从最近代的 reader（reader11）复制 `parse_grid_data.py` / `parity_check.py`
（仓库版已把 import 改成相对导入，数据目录版是绝对导入），read 逻辑仿照
`reader{ver}/read.py` 编写。

**2. 新 reader 必查四件事。**

- 科学包类型与 `multi_evt`/`multi_step`（12B 用 ft + 41/12）
- HK 包 tag 与长度
- `.dat` ↔ `.hk` 配对规则：12B 的 TB/放射源是同名 `.hk`；X 光机是
  `{idx}_{kV}_ch{n}`，dat 与 hk 的 idx 不一致，要按 (kV, ch) 配对，
  配对不上的能量点整体跳过
- 地面数据 `utc` 全 0，不能用 11B 的 utc 时间cut；HK 文件包含偏压爬升段，
  备份的 hk 还可能在同一偏压下带很长的回温段（如 12B 的 ecu_1_088.hk 在
  29V 下从 −16.4°C 回温到 +29.5°C），所以要按（偏压, 温度）二维稳定段
  截取：reader12 对逐记录通道均值做 0.1V × 1°C 分箱，取记录数最多的箱；
  一个 hk 文件覆盖多个偏压点位时（如 12B 的 093 覆盖 27.0/27.5V 两段），
  用 `hk_bias` 参数显式指定目标偏压；一个 observe 文件覆盖两个点位时
  （前半/后半），用 `sci_half` 按事件序切分
- 不要用"全记录平均"读 HK：爬升段会把平均偏压拉低几 V。12B 曾据此误判
  放射源轮偏压 20.4V 低于击穿电压，实际各轮都升到了 28.5V

**3. fit_range 用直方图找峰生成，不要手猜。**

- 源峰可能贴近阈值沿（低增益时 Am241 峰在 ~150-500 ADC），要把直方图
  范围压小、bin 放小（4 ADC）才看得见；检测时先找阈值沿，再找沿后谷值，
  峰取谷值之后平滑谱的最大值，prominence（峰/谷）不足的通道置
  null 跳过，宁缺毋滥
- 拟合窗口要窄：以峰位为中心 ±2.8σ（σ 由峰右半高宽估计），下沿夹在
  谷值之上，把阈值沿排除在窗口外；窗口太宽会把右侧连续谱的下降尾巴
  包进来，线性本底拟合会把峰位拉偏
- X 光机数据增益低，光子信号埋在噪声包附近，要用同 kV 点其他通道的文件
  互扣本底（`payload.ec.xray_bkg_rotation` = `fixed [1,2,0,0]` / `circle` 就是干这个的）才能
  露出峰；注意扣除残余会在峰低能侧产生假的负凹陷，拟合窗口下沿要避开
- 低 kV 点（12B 的 ≤65 kV）的峰包不是干净的全能峰：管压越低，连续谱
  端点离阈值越近，谱形被阈值截断，且端点以上有堆积尾巴（12B 的 40 kV
  数据计数延伸到 ~130 keV 表观能量，真实光子不可能超过管压）。这类点
  按"能量=管压"赋值会系统偏高（12B 实测 +8%~+35%，越低越严重），不能
  用作 EC 锚点。判断一个峰包能不能用：看它随管压的移动是否落在已有
  E-C 的延长线上、端点以上有没有尾巴。12B 各点位的具体状态见
  [12B_13B/data.md](12B_13B/data.md)
- Co60 有两个全能峰（1173、1332 keV），相距仅 ~3σ，拟合错峰会把整个
  高能端带歪（12B 曾把 1173 当 1332，导致 662 keV 残差恒为 +2.3%）。
  判峰方法：双峰位置比应等于能量比（1173/1332 = 0.88）；用 Na22 的
  1275 keV 峰（不参与拟合）做独立验证。拟合用双高斯：窗口覆盖两个峰，
  `bkg_form` 填 `"gaus"`（主峰高斯 + 高斯"本底"吸收另一个峰，配置层面
  即可实现，无需改代码）；主峰初值取自数据右半，所以窗口要让最右边的
  峰占据右半
- TB 文件先做 md5 查重：12B 有两个 dat 是前一点位的逐字节复制
  （真实数据丢失），要在 operation 的文件列表里剔除

**3b. TB 数据质量筛查：同一点位复测一致性。**

接入后先列出每个文件的实测 (温度, 偏压, 峰位) 表。12B 的实际温度与文件名
档位差异很大（SiPM 自热，m0C 档实测 +15.5°C），必须以 HK 实测值为准。
要重点排除两类坏文件，它们会把整个二维面拉歪（12B 实测剔除两个文件后
残差从 ±20% 降到 <2%）：

- HK 里有多个偏压平台段（科学数据采在前面一段，平台均值却是后面一段，
  表现为同一 (T,B) 的两个文件峰位差 30%）
- 运行中温度漂移（12B 正常文件逐条 HK 温度 std ~0.03°C，漂移文件
  ~0.9°C）

二维面拟合不收敛时先检查初值：共享的 `temp_bias_fit_curvefit` 默认初值
是为其他版本调的，12B 用自己的初值（`payload.yaml` 的 `tb.tb_fit_p0` /
`tb.tb_fit_maxfev`，`global_tb` 传给 `temp_bias_fit_curvefit` 的 `p0`/`maxfev`）。
验收标准：相对残差 max < 5%。

**3c. 低增益点位的峰不要轻易放弃：用模型反推重试。**

自动找峰（prominence 阈值）会放弃峰与阈值沿部分重叠的点位。但这些峰
肉眼可见、且物理上必然存在。流程：先用高置信度点位拟合出二维响应，再对每个
缺测通道用模型在实测 (T, B) 下反推预期峰位，开窗口 [谷值, 预期+2.5σ] 做
高斯+线性拟合，验收条件：拟合中心与预期偏差 ≤12%、redchi<5、
center_err<5%。重试成功的点写回 fit_range 后重跑 `calib global {ver} tb`，迭代到没有新增
为止（12B 由此多救回 60 个通道，含全部 27.0V 行）。

注意 12B 的 27.0V 行拟合中心对窗口不敏感（窗口扫描中心变化 <1 ADC），
但对二维面有 +3~+6% 的系统偏高——实测增益随过偏压的变化比模型的
Vov² 形式平缓，是模型形式本身的偏差，不是拟合错误。因此二维拟合限定
在偏压 ≥27.5V 范围（`payload.yaml` 的 `bias_min_filter`，`global_tb` 在点过滤时使用），27.0V 行的单谱
拟合结果仍保留在 pickle 中可参考。EC 运行偏压 28.5V，在适用范围内。

**4. EC 的两个易踩的坑。**

- 能量映射（`payload.ec.energy_map`，原 `ec_energy.json`）的键命名影响 EC 全局
  对 X 光机/放射源条目的区分：10B/11B 用文件名里是否含 `observe`，12B 的放射源
  文件名不含 observe，改用"键是否纯数字（kV 值）"区分。新架构里这个"哪个文件是
  source/xray"由 manifest 的 `branch` 决定（`ec_source` vs `ec_xray`），不再依赖
  文件名启发式。
- EC 拟合按 Gd K 吸收边（50.2 keV）拆成低/高两段二次拟合，拆分阈值
  （`energy_split_low`/`energy_split_high`）要按实际点位设置：12B 用
  49/55 keV，落在死区里的点（如 51/53 kV）不参与拟合。若某载荷所有点
  都在拆分线同一侧，空的那组会让拟合崩溃，需调整阈值
- TB 二维面拟合只在 TB 数据覆盖的温压范围内可靠；接入后先检查各 EC 点
  的修正因子是否合理（数量级 1），不合理多半不是模型外推问题，而是
  EC 轮 HK 读数被爬升段污染（见第 2 点），先修 HK 截取再下结论。
  12B 最终用全温度段 TB 拟合结果做修正，修正因子 ~1.1

**5. 跑完看 QA。**
`single_process/{TB,EC}_fit_result/*.pickle` 里每个通道有 `qa_flag`
（ok/warn/fail）和 `redchi`，汇总确认通过情况；拟合失败或跳过的点要在
结果说明里如实记录原因。
