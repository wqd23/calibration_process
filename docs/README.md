# docs/ — 技术文档

修订日期：2026-09-16

本目录存放面向人阅读的技术文档：背景知识讲解与专题技术备忘，读者默认是有物理或工程背景、但不熟悉本系统细节的人。项目本身的用法、配置与状态记录仍以根目录 README.md 为准，docs 里的文档被 README 相应章节引用，不重复其内容。

写作约定：开头标注修订日期；先讲背景与思路、再给细节；术语首次出现时给一句话解释；关键结论给出真实可查的出处链接；代码注释与提交信息仍按仓库约定用英文。

## 文档清单

| 编号 | 文档 | 内容 |
|------|------|------|
| D2 | [data.md](data.md) | 数据准备与目录约定：新配置（`configs/{ver}/*.yaml`（含 `reader.yaml`）+ manifest）的格式、新载荷接入流程、以及方法论经验（找峰、HK 截取、质量筛查、拟合技巧） |
| D3 | [results.md](results.md) | 结果产物清单与质量判定指标（pickle 字段、`fit.json`/`spectrum.parquet`、redchi、qa_flag）的含义 |
| D4 | [workflows.md](workflows.md) | 每个载荷各自的显式 workflow：文件选择、reader、背景轮转、分辨率方法、特殊处理，一步可从 `v{ver}.py` 读出 |
| D5 | [workflow_matrix.md](workflow_matrix.md) | 历史 workflow 完整审计矩阵（重构前每个版本实际做了什么，逐版本表格） |
| D6 | [12B_13B/data.md](12B_13B/data.md) | 12B/13B 类载荷的数据说明：温度偏压点位与文件对照表、X 光机/放射源各数据文件的状态与特殊情况 |
| D7 | [payloads_data.md](payloads_data.md) | 其它 7 个载荷的数据说明（简化版）：通用结构 + 每个版本的目录/reader/命名/选择与排除/已知坑 |
| D8 | [intermediate_data.md](intermediate_data.md) | 中间数据分层（L1 忠实帧 / L2 处理输出 / L3 `fit.json`+`spectrum.parquet`、pickle）、怎么读、能否跨项目分析、如何定制 pipeline；第 5/5b 节为 L1/L2 缓存与 L3/L4，第 6/7 节为读取层统一与 reader golden，第 8 节为分层流水线与 `grid_common` |
| D9 | [N1/data.md](N1/data.md) | N1（GRIDN1）的数据说明：GAGG/CLYC 双数据集、两种包格式（ft/wf 512）、扫描文件分段、坏点清单 |
| D10 | [test_fixtures.md](test_fixtures.md) | 测试数据与 legacy oracle：三类测试数据（git 内 golden、raw_data、仓库外 oracle）、oracle 的 8 个版本与存放/备份位置、软链与 `CALIB_ORACLE_DIR` 用法、raw_data 路径表、维护约定 |

**数据说明的约定**：数据特殊的载荷（如 12B/13B）单独放一个子目录详述
（数据从哪来、每个文件对应什么条件、哪些文件有坑：截断/复制/补测等）。
数据结构高度一致的载荷（本次其余 7 个：一科学文件=一条件，靠文件名选点）
用 [payloads_data.md](payloads_data.md) 的简化模板统一说明即可，无需各建一个
复杂子目录。新增载荷时先判断属于哪种，再仿照建。

