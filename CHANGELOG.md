# Changelog

本项目从当前提交起维护 Changelog。条目**不自动生成**，只在需要时人工添加（例如重构、新版本接入、行为变更等）。历史提交不再补记。

版本号遵循 [Semantic Versioning](https://semver.org/)，并见 `pyproject.toml`。

## [0.2.0] - 2026-09-09

### Changed
- 重构了整个项目框架与 workflow：显式 workflow（`workflows/versions/v{ver}.py` 每版本一个）+ 人工确认的 manifest + `configs/{ver}/*.yaml`（YAML + strict Pydantic schema）+ 统一 `calib` CLI。
- 将隐藏在版本继承、配置文件、目录扫描中的处理流程，改为显式、分层、可复用、可测试的 scientific workflow。
- **原有的科学数据处理流程与处理代码保持不变**：单谱拟合、TB/EC 数学模型、reader、正式绘图、中间/最终产物格式均与迁移前一致；本次为纯架构重构，结果经 `legacy vs new` 全量回归验证一致。
