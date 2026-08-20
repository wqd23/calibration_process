# 部署指南

本指南面向首次在新机器上部署本项目的组内成员。

## 前置环境

本项目使用 [uv](https://docs.astral.sh/uv/) 管理 Python 环境和依赖，使用 [just](https://github.com/casey/just) 作为命令运行器。Python 版本由 uv 自动管理（`.python-version` 中指定），不需要手动安装。

**安装 uv：**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**安装 just：**
```bash
# macOS
brew install just
# Linux (cargo)
cargo install just
# 或参考 https://github.com/casey/just#installation
```

## 安装步骤

```bash
git clone <repo-url>
cd calibration_process
uv sync                          # 安装依赖到 .venv
just init {ver} {data-path}      # 软链数据 + 创建输出目录
just check {ver}                 # 验证部署是否完整
```

其中 `{ver}` 为载荷版本号（如 `03B`），`{data-path}` 为该版本原始标定数据的绝对路径。

`just`（无参数）可查看所有可用命令和版本列表。

## 验证部署

`just check {ver}` 会逐项检查：
- **raw_data 软链**：是否存在、目标是否可达
- **输入路径**：数据目录、fit_range / bkg_form / ec_energy 配置文件是否存在
- **输出目录**：TB_fit_result / EC_fit_result / single_fit_fig / tb_logs / ec_logs 是否存在（缺失时加 `--fix` 自动创建）

全部 OK 时输出 `[ver] READY`，exit code 0。有缺失时输出 `[ver] NOT READY`，exit code 1。

## 已知问题

**cachier 缓存锁**：`lib_reader` 使用 cachier 缓存读取结果到 `.cache/` 目录。03B/04/05B/07 使用共享缓存数据库，多进程并发读写同一数据库会触发 portalocker 锁冲突。当前解决方案：这 4 个版本在批跑时内部串行（已由脚手架处理），不同版本之间可以并行。

如遇锁冲突（报错 `AlreadyLocked`），删除 `.cache/.lib_reader*` 文件后重跑即可。
