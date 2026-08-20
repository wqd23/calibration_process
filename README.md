# GRID 标定数据处理

本仓库用于 GRID 载荷标定数据的处理、拟合与质量检验。

## 快速上手

**前置依赖：** [uv](https://docs.astral.sh/uv/)（Python 环境管理）和 [just](https://github.com/casey/just)（命令运行器）。
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh   # 安装 uv
cargo install just                                  # 安装 just（或 brew install just）
```

**部署：**
```bash
git clone <repo-url> && cd calibration_process
uv sync                              # 安装依赖
just init {ver} {data-path}          # 软链数据 + 创建目录（ver 用 just 查看）
just check {ver}                     # 验证部署
```

**跑数据：**
```bash
just all {ver}                       # 处理全部 TB + EC 数据
just tbfit {ver}                     # TB 二维面拟合
just ecfit {ver}                     # E-C 关系拟合
```

详细文档见 [docs/](docs/) 目录。

## 命令速查

| 命令 | 说明 |
|------|------|
| `just` | 查看所有命令和可用版本 |
| `just init {ver} {path}` | 初始化版本：软链数据 + 建目录 |
| `just check {ver}` | 验证部署完整性（`--fix` 自动创建缺失目录） |
| `just all {ver}` | 处理该版本全部数据 |
| `just tb {ver} run {idx}` | 处理单个 TB 文件 |
| `just ec {ver} x run {idx}` | 处理单个 EC X光机文件 |
| `just ec {ver} src run {idx}` | 处理单个 EC 放射源文件 |
| `just tbfit {ver}` | TB 二维面拟合 |
| `just ecfit {ver}` | E-C 关系拟合 |
| `just new-payload {ver}` | 为新载荷生成配置和目录骨架 |

## 仓库结构

```
.
├── README.md
├── docs/                           # 详细文档
│   ├── deploy.md                   # 部署指南
│   ├── data.md                     # 数据准备与目录约定
│   └── results.md                  # 结果产物与 QA 指标
├── data/{ver}/                     # 各版本数据目录
│   ├── raw_data -> /path/to/data   # 原始数据软链
│   ├── single_process/             # 配置文件和拟合产物
│   ├── tb_logs/                    # TB 面拟合产物
│   └── ec_logs/                    # E-C 拟合产物
├── lib_reader/                     # 各版本数据读取库（workspace 子包）
├── lib_plot/                       # 绘图库（workspace 子包）
├── src/calibration_process/
│   ├── process.py                  # CLI 入口（惰性加载）
│   ├── cmd.py                      # 命令行包装
│   ├── operation.py                # TB/EC 操作类
│   ├── file_lib.py                 # 单文件读取与拟合
│   ├── fitting.py                  # 通用峰型拟合
│   ├── util_lib.py                 # 工具函数（拟合、缓存、阈值）
│   ├── check.py                    # 部署校验
│   ├── scaffold.py                 # 新版本脚手架
│   └── config.json                 # 各版本配置（模板化）
├── pyproject.toml
├── uv.lock
└── justfile
```

## TODO

- [x] 12B 简易结果
- [ ] 12B, 13B 完整结果
- [ ] 10B, 11B 塑闪结果

## Contributor

- wqd
- 丘智勇
- Lee
- yhl
- xsq(Haruka Kujo)
