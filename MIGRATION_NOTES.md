# Migration Notes — historical bugs recorded, not fixed

修订日期：2026-09-08

本次架构迁移遵循"记录，不修复"原则：新 workflow 首先复现 legacy 行为。下列问题均为历史实现中存在、但在本次迁移中**刻意保留**的行为。修复它们需要单独的科学验证，不在本次纯架构迁移范围内。

每个条目给出：现象、触发条件、影响、以及新架构为何保留。

---

## B1. EC 温度偏压参考（corr）与 TB 输出路径解耦

- **现象**：TB 全局拟合把结果写到 `data/<ver>/tb_logs/{ts}_temp_bias_fit.json`，而 EC 单文件的 temp/bias 校正因子（`corr`）读取的是 `config.json` 里 `ec.tb_result_path` 指向的 `data/<ver>/single_process/{ts}_temp_bias_fit.json`（一个**独立的**参考文件）。两个路径不同，内容在 09 上逐通道完全相等（已核对 4/4 MATCH）。
- **影响**：TB 重新拟合后不会自动更新 EC 的参考；EC 只认 `single_process` 下那个固定文件。迁移后 EC 的 `corr` 仍从 `payload.yaml` 的 `ec.tb_ref_path` 读取。
- **保留原因**：路径解耦是历史人工约定，`UNKNOWN` 其真实意图；强行合并会改变 EC 结果。

## B2. 09 的 TB 文件枚举依赖 `os.listdir` 顺序 + 硬编码 remove

- **现象**：`TB_operation_09.__init__` 用 `[f for f in os.listdir(path) if splitext(f)[1]=='.txt']` 再 `.remove()` 6 个文件。`os.listdir` 顺序是文件系统相关；若目录内容变化导致某个待 remove 文件不存在，`list.remove` 会抛 `ValueError`。
- **影响**：结果顺序不稳定（本次在同一机器上稳定）。
- **保留原因**：新 `v09.py` 用同样的 `os.listdir` + 相同 `TB_EXCLUDE` 列表复现。目录顺序以 manifest 固定后，科学不受影响。

## B3. `plot.ec_plot` 中 `center[i][xpoint]` 长度假设（潜在越界）

- **现象**：`ec_plot` 内 `center[i][xpoint]`，其中 `xpoint` 由 `x_energy`（长度 = x 点数）导出，而 `center[i]` 长度 = x 点数 + src 点数。只有 `src` 能量全部 **大于** x 能量时（排序后 src 在数组末尾）才会恰好对齐。对 09 成立（x 15–60 keV，src 511/662/1332 keV）。
- **影响**：若某个版本 src 能量不全部高于 x，散点图会错位（拟合曲线不受影响）。
- **保留原因**：这是正式 plotting function（Protected Kernel），且属于显示层；新 `global_ec` 以与 legacy 相同的参数形状调用 `plot.ec_plot`，保持行为一致。

## B4. 09 的 EC 三个源中 Co60 与 Cs137 共用同一背景文件

- **现象**：`EC_operation_09.src_bkg` 中 Co60 与 Cs137 都指向 `0827_10C_285_bkg_20m_0x00CF.txt`。这意味着 Co60 与 Cs137 使用同一个背景谱。
- **影响**：是否物理上正确无法从代码确认。
- **保留原因**：`UNKNOWN`。新 manifest 的 `aux_files` 按 legacy 硬编码映射复现。

## B5. 09 的 EC fit_range / bkg_form 默认值兜底

- **现象**：`EC_operation_09.xray_config` 使用 `fit_range.get(energy_name, [[None,None]]*4)` 与 `bkg_form.get(energy_name, "lin")`，其它版本用直接下标。
- **影响**：缺失配置不会报错，而是静默用默认值（可能导致无拟合/无背景），掩盖配置错误。新 `v09.py` 在 fit_range 文件中显式列出全部 13 个能量，缺失时由 Pydantic strict 校验报错（比 legacy 更早失败）。
- **保留原因**：legacy 行为如此；新架构用严格 schema 暴露缺失（这是架构改进，不是科学改动）。

## B6. EC X-ray 使用相邻通道作为背景

- **现象**：对 03B/04/07/09/10B/11B/12B，通道 `ch i` 的背景读的是通道 `(i+1)%4` 的同一个 channel 数据（`read_config[1:4]+[read_config[0]]` 或 `[ch1,ch2,ch0,ch0]`）；仅 05B 使用"同文件 rotated time_cut"。
- **影响**：把相邻通道当背景是历史设计，物理合理性待评估。
- **保留原因**：这是对各版本最重要的 orchestration 差异之一，已在 `_rotate_bkg`（`circle`/`fixed`）中按版本保留。

## B7. `fp_method` 与 `Read_config.ending` 两套命名

- **现象**：`process.py` 的 `fp_method`（05B 用 None，其余用 "03/04/07/09/10/11/12"）决定使用单 4 通道文件路径还是 4 通道重建；`Read_config.ending` 决定 reader。二者是独立的字符串语义。
- **影响**：易混淆，是历史积累。
- **保留原因**：新架构把这一差异建模为 `payload.ec.xray_bkg_rotation` + 分支重建，并显式记录在 `workflow_matrix.md`。

## B8. 10B/11B EC 仅 3 通道 + 补丁成 4 通道

- **现象**：`EC_operation_10B/11B.ec_fit` 用 `CHN_NUM=3`，随后 `center=[c0,c1,c2,c0]`、`result=[r0,r1,r2,r0]`、`src_result/x_result` 同样把 ch3 补成 ch0。
- **影响**：下游 `plot.ec_plot` 需要 4 通道，但 ch3 实际上是 ch0 的拷贝。物理意图 `UNKNOWN`。
- **保留原因**：这是历史行为。新架构通过 `payload.ec.channel_count` + `global_ec` 的补齐逻辑保留（09 为 4 通道，不使用补齐）。

## B9. 12B 的逐通道 `None` fit range

- **现象**：12B 的 `fit_range.json` 中，某些 measurement 的某些通道 range 为 `None`（该通道不参与拟合）。`File_operation_05b.peak_fit` 对 `None` 通道产出 `fit_result=None`，下游 `load_data` / `ec_fit` 跳过。
- **影响**：这些通道在 TB/EC 中不产生点；点构建时被跳过。
- **保留原因**：这是历史行为。新 `FitRangeSet` 允许逐通道 `null`，`build_tb_points` / `build_ec_points` 跳过 `None` 拟合。

## B10. 12B EC-Xray 完整性过滤（HK 配对 + fit range 完整）

- **现象**：12B `EC_operation_12B` 对 xray 能量做两级过滤：(a) `__x_hk_complete`（4 通道文件都能配对到 HK，否则丢弃整点）；(b) `all(r is not None for r in fit_range[e])`（4 通道 fit range 完整才保留）。
- **影响**：能量点要么 4 通道齐全、要么整点丢弃；缺失的通道不会被半途使用。
- **保留原因**：新 `v12B.py` 以 `xray_require_hk` / `xray_require_fit_range` 复现。

## B12. 05B 的 X 光机单 4 通道文件 + 时间窗（time_cut）

- **现象**：05B 的 EC X-ray 每个管压在**同一个文件**里带 4 通道（不是每通道一个文件）；背景是**同一个文件**按通道循环移位的时间窗（`bkg_time_cut = {k:[v[1],v[2],v[3],v[0]]}`）再读一次；reader 为 `xray` 并需要 `x_config`。能量分界为 49/52 keV。
- **影响**：这是与其它版本（03B/04/07/09/10B/11B/12B 的每通道文件）完全不同的 X 光机路径。
- **保留原因**：新 `v05B.py` 以 `xray_single_file` + `time_cut` + `xray_config_file` 复现；`common.single_run_spec` 对 `xray_single_file` 走单文件 fp05B 路径。
- **备注**：`xray_config_file` 保留 config.json 的**完整** `data/<ver>/...` 字符串（而非 manifest 相对路径），以与 legacy 的 `config_file` 字段逐字符一致。

## B11. 12B TB 的定制拟合参数与 bias 过滤

- **现象**：`TB_operation_12B` 覆盖 `TB_FIT_P0=[-0.02,0.07,24.4,-35.0,-1000.0]`、`TB_FIT_MAXFEV=100000`，且 `load_data` 把参与 2D 拟合的点限制在 `bias>=27.25 V`。文档注释说明这是为了适配 12B 增益响应。
- **影响**：这是科学方法参数，不是 bug；对 12B 是刻意的 version override。
- **保留原因**：新 `payload.yaml` 的 `tb.tb_fit_p0 / tb_fit_maxfev / bias_min_filter` 原样携带，`global_tb` 读取它们。
