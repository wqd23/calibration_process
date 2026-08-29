# AGENTS.md / 工程约定与接入方法论

本文件面向在本仓库工作的 AI 助手（以及需要接入新载荷的人）。这里只写
「为了让工作不出错而必须遵守」的规则、命令与踩坑经验；面向人阅读的背景
知识在 [docs/](docs/) 的对应文档里。代码注释与 commit message 一律用英文。

## 常用命令

环境用 [uv](https://docs.astral.sh/uv/) 管理、命令用 [just](https://github.com/casey/just) 运行：

```bash
just                # 列出所有命令 + 可用版本（读 config.json 的键）
uv sync             # 安装依赖到 .venv
just init {ver} {path}    # 首次接入：软链 raw_data + 建输出目录
just check {ver} [--fix]  # 验证部署完整性（--fix 自动建缺失输出目录）
just tb {ver} run {idx}   # 处理单个 TB 文件（idx: 0-3 或 all）
just ec {ver} x run {idx} # 处理单个 EC X 光机文件
just ec {ver} src run {idx}# 处理单个 EC 放射源文件
just all {ver}            # 处理该版本全部 TB + EC
just tbfit {ver}          # TB 二维面拟合
just ecfit {ver}          # E-C 关系拟合
just new-payload {ver}    # 为新载荷生成 config 条目 + 目录 + reader 骨架
```

lint / 测试：

```bash
ruff check src lib_reader lib_plot    # 代码风格
pytest tests/                           # 单元测试（tests/ 目前可为空）
coverage run --source src/ -m pytest tests/ && coverage report -m
```

## 代码约定

- 不加注释（除非必要且用英文）；docstring 简短。
- `calibration_process`/`lib_reader`/`lib_plot` 是 uv workspace 子包
  （见 `pyproject.toml` 的 `[tool.uv.workspace]`），包内部用相对导入。
- 环境由 uv 管理，CLI 入口统一为 `python -m calibration_process`（`process.py`
  惰性加载，`cmd.py` 做命令行包装）。
- 配置集中在 `src/calibration_process/config.json`，路径以 `data/{ver}/raw_data/...`
  相对引用；`_defaults` 是 tb/ec 的模板段。
- 数据软链到 `data/{ver}/raw_data`（不复制），产物输出到
  `data/{ver}/single_process/`、`tb_logs/`、`ec_logs/`。

## 接入新载荷的方法论与坑

这是最容易出错的部分。接一个新版本时按下面的顺序走，每一条都是踩过的坑。

### 1. 先判定数据包格式，再写 reader

用 `xxd 文件 | head` 看文件头字节，对照随数据附带的解析代码
`grid_packet.xml` 里各包的 `head` 属性确定包类型。同一系列不同代载荷的
包定义可能不同：

- 12B/13B 的 `.dat` 是 `grid1x_ft_packet`（特征包，528 字节，41 事件/包，
  每事件 12 字节：timestamp 4B + data_max 2B + data_base 2B + data_sum 4B），
  而 11B 的 `.dat` 是 `grid1x_wf_packet`（波形包）；两者 ft 包字段布局也不同
  （11B 的 timestamp 是 8 字节、无 data_sum）。
- 12B/13B 的 `.hk` 是 187 字节的 `grid1x_hk_packet`，11B 是 139 字节的
  `hk_grid1x_packet`，sipm 字段偏移不同（82 vs 111）。

最稳妥的做法：把数据目录里随附的 `grid_packet.xml` 整体复制进新 reader，
再从最近代的 reader（如 reader11）复制 `parse_grid_data.py` / `parity_check.py`
（仓库版已改相对导入，数据目录版是绝对导入），read 逻辑仿照已有的
`reader{ver}/read.py` 编写。

### 2. 新 reader 必查四件事

- 科学包类型与 `multi_evt`/`multi_step`（12B 用 ft + 41/12）。
- HK 包 tag 与长度（不同代不同，见上）。
- `.dat` ↔ `.hk` 配对规则：TB/放射源通常是同名 `.hk`；X 光机形如
  `{idx}_{kV}_ch{n}`，dat 与 hk 的 idx 可能不一致，要按 (kV, ch) 配对，
  配对不上的能量点整体跳过。
- 地面数据 `utc` 全 0，不能沿用 11B 的 utc 时间 cut。
  HK 文件含偏压爬升段，备份的 hk 还可能在同偏压下带很长的回温段（如 12B 的
  `ecu_1_088.hk` 在 29V 下从 −16.4°C 回温到 +29.5°C）。要按（偏压, 温度）
  二维稳定段截取：`reader12` 对逐记录通道均值做 0.1V × 1°C 分箱，取记录数
  最多的箱。一个 hk 文件覆盖多个偏压点位时（如 `093` 覆盖 27.0/27.5V）用
  `hk_bias` 显式指定目标偏压；一个 observe 文件覆盖两个点位时用 `sci_half`
  按事件序切分。

> 不要用「全记录平均」读 HK：爬升段会把平均偏压拉低几 V。12B 曾据此误判
> 放射源轮偏压 20.4V 低于击穿电压，实际各轮都升到了 28.5V。

### 3. fit_range 用直方图找峰生成，不要手猜

- 源峰可能贴近阈值沿（低增益时 Am241 峰在 ~150-500 ADC）：把直方图范围压小、
  bin 放小（4 ADC）才看得见。检测时先找阈值沿，再找沿后谷值，峰取谷值之后
  平滑谱的最大值；prominence（峰/谷）不足的通道置 null 跳过，宁缺毋滥。
- 拟合窗口要窄：以峰位为中心 ±2.8σ（σ 由峰右半高宽估计），下沿夹在谷值之上，
  把阈值沿排除在窗口外。窗口太宽会把右侧连续谱的下降尾巴包进来，线性本底拟合
  会把峰位拉偏。
- X 光机数据增益低，光子信号埋在噪声包附近，要用同 kV 点其他通道的文件互扣
  本底（`xray_config` 里的 bkg 轮转 `[1,2,0,0]`）才能露出峰；注意扣除残余会在
  峰低能侧产生假的负凹陷，拟合窗口下沿要避开。
- 低 kV 点（12B 的 ≤65 kV）的峰包不是干净的全能峰：管压越低连续谱端点离阈值
  越近、谱形被阈值截断，端点以上还有堆积尾巴（40 kV 数据计数延伸到 ~130 keV
  表观能量）。这类点按「能量=管压」赋值会系统偏高（12B 实测 +8%~+35%，越低越
  严重），不能用作 EC 锚点。判断能否用：看它随管压的移动是否落在已有 E-C 的
  延长线上、端点以上有没有尾巴。
- Co60 有两个全能峰（1173、1332 keV），相距仅 ~3σ，拟合错峰会把高能端带歪
  （12B 曾把 1173 当 1332，导致 662 keV 残差恒为 +2.3%）。判峰方法：双峰位置比
  应等于能量比（1173/1332 = 0.88）；用 Na22 的 1275 keV 峰（不参与拟合）做独立
  验证。拟合用双高斯：窗口覆盖两峰，`bkg_form` 填 `"gaus"`（主峰高斯 + 高斯
  "本底"吸收另一峰，配置层面即可实现）；主峰初值取自数据右半，所以窗口要让
  最右边的峰占据右半。
- TB 文件先做 md5 查重：12B 有两个 dat 是前一点位的逐字节复制（真实数据丢失），
  要在 operation 的文件列表里剔除。

### 4. TB 数据质量筛查：同一点位复测一致性

接入后先列出每个文件的实测 (温度, 偏压, 峰位) 表。12B 的实际温度与文件名档位
差异很大（SiPM 自热，m0C 档实测 +15.5°C），必须以 HK 实测值为准。重点排除两类
坏文件（12B 剔除两个后残差从 ±20% 降到 <2%）：

- HK 里有多个偏压平台段（科学数据采在前面一段、平台均值却是后面一段，表现为
  同一 (T,B) 的两个文件峰位差 30%）。
- 运行中温度漂移（正常文件逐条 HK 温度 std ~0.03°C，漂移文件 ~0.9°C）。

二维面拟合不收敛时先检查初值：共享的 `temp_bias_fit_curvefit` 默认初值是为其他
版本调的，新版本用自己的初值（如 `TB_operation_12B` 的类属性
`TB_FIT_P0`/`TB_FIT_MAXFEV`，基类 `temp_bias_fit` 会把它们传给
`temp_bias_fit_curvefit` 的 `p0`/`maxfev`）。验收标准：相对残差 max < 5%。

### 5. 低增益点的峰不要轻易放弃：用模型反推重试

自动找峰（prominence 阈值）会放弃与阈值沿部分重叠的点位，但这些峰肉眼可见、
物理上必然存在。流程：先用高置信度点位拟合出二维响应，再对每个缺测通道用模型
在实测 (T, B) 下反推预期峰位，开窗口 [谷值, 预期+2.5σ] 做高斯+线性拟合，验收
条件：拟合中心与预期偏差 ≤12%、redchi<5、center_err<5%。重试成功的点写回
fit_range 后重跑 tbfit，迭代到没有新增为止（12B 由此多救回 60 个通道）。

注意：12B 的 27.0V 行拟合中心对窗口不敏感（窗口扫描中心变化 <1 ADC），但对
二维面有 +3~+6% 的系统偏高——实测增益随过偏压的变化比模型的 Vov² 形式平缓，
是模型形式本身的偏差。因此二维拟合需限定在偏压 ≥27.5V（`load_data` 过滤），
27.0V 行的单谱结果仍保留在 pickle 中可参考。EC 运行偏压 28.5V，在适用范围内。

### 6. EC 的两个坑

- `ec_energy.json` 的键命名影响 `cmd.py` 里 ecfit 对 X 光机/放射源条目的区分：
  10B/11B 用文件名是否含 `observe` 区分，12B 的放射源文件名不含 observe，改用
  「键是否纯数字（kV 值）」区分（`VersionProcessOp12B`）。
- EC 拟合按 Gd K 吸收边（50.2 keV）拆成低/高两段二次拟合，拆分阈值
  （`energy_split_low`/`energy_split_high`）按实际点位设置：12B 用 49/55 keV，
  落在死区里的点（如 51/53 kV）不参与拟合。若某载荷所有点都在拆分线同一侧，
  空的那组会让拟合崩溃，需调整阈值。
- TB 二维面拟合只在 TB 数据覆盖的温压范围内可靠；接入后先检查各 EC 点的修正
  因子是否合理（数量级 1）。不合理多半不是模型外推问题，而是 EC 轮 HK 读数被
  爬升段污染（见第 2 点），先修 HK 截取再下结论（12B 最终用全温度段 TB 拟合
  结果做修正，修正因子 ~1.1）。

### 7. 跑完看 QA

`single_process/{TB,EC}_fit_result/*.pickle` 里每个通道有 `qa_flag`
（ok/warn/fail）和 `redchi`，汇总确认通过情况；拟合失败或跳过的点要在结果说明里
如实记录原因。各载荷数据的状态与特殊情况（哪些文件有坑）记录在
`docs/{载荷族}/data.md` 里（如 `docs/12B_13B/data.md`）。
