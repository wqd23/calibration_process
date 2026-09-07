# Migration Checklist

修订日期：2026-09-08

| Item | Status | Evidence |
|---|---|---|
| legacy oracle frozen | ✅ | `.oracle/09` (TB 48 + EC 16 pickles, tb/ec logs, figures) |
| workflow matrix complete | ✅ | `docs/workflow_matrix.md` |
| YAML schema implemented | ✅ | `config_schema.py` (strict/Pydantic) |
| old config migration verified | ✅ | `migration.py` + `test_migration_writes_valid_yaml` |
| manifest generated and verified | ✅ | `configs/09/*_manifest.yaml` (TB 48 / src 3 / xray 13) |
| preview works without fit range | ⬜ | preview stage intentionally not built in this phase (see note) |
| single fit regression passes | ✅ | `test_regression_09.py` + pipeline run (Level 1) |
| TB points regression passes | ✅ | `test_tb_points_match_legacy_load_data` (Level 2) |
| EC points regression passes | ✅ | `test_ec_global_outputs_match_legacy` (Level 2/3) |
| final TB output passes | ✅ | tb_logs JSON exact match |
| final EC output passes | ✅ | ec_logs JSON + npy exact match |
| single-fit figures checked | ✅ | pixel-identical |
| TB figures visually checked | ✅ | pixel-identical |
| EC figures visually checked | ✅ | `~0.03%` scatter-point-only diff |
| representative version full-run passes | ✅ | `test_full_pipeline_matches_oracle` |
| all historical versions full-run pass | ⬜ | beyond Gate B (requires migrating other versions) |
| legacy retired | ⬜ | beyond Gate D |

---

**Note on preview**: the plan's `preview` stage (single-file preview that works even
without a fit range, plus log-log/auto-peak) is a separate feature that this
phase deliberately does **not** build. The single-fit path (which is preview's
sibling but operates on confirmed ranges) is the focus of Gate B. This will not
affect the scientific-equivalence proof.
