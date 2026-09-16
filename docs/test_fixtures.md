# 测试数据与 legacy oracle

修订日期：2026-09-16

本文说明跑完整测试时需要的三类数据：**git 内自带的 golden**、**本机的原始数据
（raw_data）**、以及**仓库之外的 legacy oracle**。重点是最后一类：它是什么、
放在哪、怎么配置、为什么不能进 git。回归策略的背景见根目录
[`AGENTS.md`](../AGENTS.md) 的「回归策略」一节。

## 一、三类测试数据

| 类别 | 在哪 | 谁需要 |
|------|------|--------|
| golden（小体积冻结基准） | **在 git 里**：`tests/golden/`（约 2.6 MB） | `test_golden.py`、`test_reader_golden.py`、`test_l1_frames.py` 等，clone 即可跑 |
| raw_data（各载荷原始数据） | 本机 `/home/wqd/cali_data/`（不在 git） | reader / L1 / L2 / 端到端测试；缺数据时部分测试跳过、部分直接失败 |
| legacy oracle（旧实现的完整产物） | 本机固定目录 `/home/wqd/cali_data/calib_test_fixtures/oracle`（不在 git） | `tests/test_pipeline_run.py` 的逐字段回归对照 |

`tests/golden/` 是自包含的，因此**只跑单元与 golden 测试时不需要任何外部数据**。
要跑完整测试（尤其是 `test_pipeline_run.py`），才需要 raw_data 和 oracle。

## 二、legacy oracle 是什么

oracle 是迁移到显式 workflow **之前**的旧实现，在 8 个载荷版本上跑出的完整产物
快照，冻结下来作为新实现的回归基准：

- 版本：`03B`、`04`、`05B`、`07`、`09`、`10B`、`11B`、`12B`
- 每个版本包含：`TB_fit_result/`、`EC_fit_result/`（单谱拟合 pickle）、
  `tb_logs/`、`ec_logs/`（全局拟合的 json/npy/图）、`single_fit_fig/`（图）；
  其中 12B 还带 L3 可移植产物 `*.fit.json` 与 `*.spectrum.parquet`
- 体积：原始目录约 394 MB；`tar.zst` 压缩备份约 141 MB
- `tests/test_pipeline_run.py` 目前只回归 **09 和 12B** 两个版本
  （`VERSIONS` 列表），其余版本的 oracle 作为存量基准保留，暂未被测试引用

**它不可再生**：旧实现（`operation.py` 等）已经删除，无法用现在的代码重新生成
这份产物。因此它必须被当作只读冻结数据长期保存，并且要有备份。

**GRIDN1 没有 legacy oracle**：GRIDN1 是迁移之后才接入的载荷，没有对应的旧产物；
它的回归基准是 git 内的 `tests/golden/GRIDN1/`。

**它不放进 git**：一是体积大（算上 LFS/历史会拖慢 clone），二是它对标定数据
有敏感性，不适合放到公开仓库。仓库里只保留"如何取得和使用"的说明与脚本。

## 三、存放与备份

固定目录（仓库外，跨 worktree 共用）：

```
/home/wqd/cali_data/calib_test_fixtures/
├── oracle/                          # oracle 实体（8 个版本，约 394 MB）
├── oracle_legacy.tar.zst            # 压缩备份（约 141 MB）
├── oracle_legacy.tar.zst.sha256     # 备份的校验和
└── repo_all_refs_20260916.bundle    # 一次性迁移备份：仓库全部本地分支与 tag
```

`repo_all_refs_*.bundle` 是仓库迁移时留下的全 refs 快照，可用
`git clone <bundle> <dir>` 恢复任意本地分支；确认不再需要后可自行删除。

每个 worktree 里用软链指向实体：

```bash
bash scripts/setup_test_fixtures.sh              # 默认给当前目录建 .oracle 软链
bash scripts/setup_test_fixtures.sh /path/to/worktree
```

脚本会：检查固定目录 → 缺失时校验 sha256 并从 `oracle_legacy.tar.zst` 解压
（也支持用 `CALIB_ORACLE_URL` 下载）→ 建立 `.oracle` 软链。路径可用环境变量覆盖：

| 变量 | 默认值 |
|------|--------|
| `CALIB_FIXTURES_DIR` | `/home/wqd/cali_data/calib_test_fixtures` |
| `CALIB_ORACLE_ARCHIVE` | `$CALIB_FIXTURES_DIR/oracle_legacy.tar.zst` |
| `CALIB_ORACLE_URL` | 空（仅当本地没有备份时才需要） |

`.oracle` 已被 `.gitignore` 忽略（目录和软链都忽略），所以软链不会出现在
`git status` 里。

## 四、在测试里使用

### 方式一：worktree 下的软链（默认）

`tests/test_pipeline_run.py` 默认从仓库根目录的 `.oracle/<ver>` 查找 oracle。
执行上面的脚本建立软链即可。

### 方式二：`CALIB_ORACLE_DIR` 环境变量

不想建软链时，可以直接把 oracle 根目录告诉测试：

```bash
CALIB_ORACLE_DIR=/home/wqd/cali_data/calib_test_fixtures/oracle \
    pytest tests/test_pipeline_run.py
```

两种方式等价；软链适合长期开发，环境变量适合临时或不想改动工作区的场景。

### 从备份恢复

```bash
cd /home/wqd/cali_data/calib_test_fixtures
sha256sum -c oracle_legacy.tar.zst.sha256
tar --zstd -xf oracle_legacy.tar.zst        # 解出 oracle/
```

## 五、raw_data 的挂载位置

完整测试还需要各版本的原始数据。用 `just init {ver} {path}` 建软链，
当前机器上的路径如下：

| 版本 | raw_data 路径 |
|------|---------------|
| 03B | `/home/wqd/cali_data/03B` |
| 04 | `/home/wqd/cali_data/04` |
| 05B | `/home/wqd/cali_data/05B` |
| 07 | `/home/wqd/cali_data/07` |
| 09 | `/home/wqd/cali_data/09` |
| 10B | `/home/wqd/cali_data/10B` |
| 11B | `/home/wqd/cali_data/11B` |
| 12B | `/home/wqd/cali_data/12B、13B/data` |
| GRIDN1/GAGG | `/home/wqd/cali_data/GRIDN1/data/260303温度偏压/GAGG` |
| GRIDN1/CLYC | `/home/wqd/cali_data/GRIDN1/data/260303温度偏压/CLYC` |
| GRIDN1/Neutron | `/home/wqd/cali_data/GRIDN1/data/260322中子束流` |
| GRIDN1/EC | 用 `ec_src`/`ec_xray`/`ec_xray_low` 三个软链，见 [N1/data.md](N1/data.md) |

缺少 raw_data 时：`test_frame_io`、`test_l1_cache`、`test_l1_processed`
等会跳过（有 skip guard）；`test_versions_enumerate_matches_manifest`
等会直接失败，因为它要扫描各版本的原始目录来核对 manifest。

## 六、维护约定

- oracle 是**只读冻结基准**：不要用 `just oracle {ver}` 的输出覆盖它。
  那个命令是"把当前产物快照成 oracle"，与这份 legacy oracle 是两回事。
- 至少保留两份：固定目录里的实体 + `oracle_legacy.tar.zst` 压缩备份。
- 校验和随备份一起保存，恢复前先 `sha256sum -c`。
- 如果确需更新 oracle 内容（原则上不应发生），要作为独立决策记录在
  CHANGELOG，并同步更新 `tests/golden/` 中受影响的冻结值。