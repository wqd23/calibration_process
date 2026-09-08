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
| single fit regression passes | ✅ | `test_regression.py` config-level (09/04/07/12B); full-run 09 + 12B |
| TB points regression passes | ✅ | full-run 09 + 12B; config-level 04/07 |
| EC points regression passes | ✅ | full-run 09 + 12B; config-level 04/07 |
| final TB output passes | ✅ | tb_logs JSON exact match (09 + 12B) |
| final EC output passes | ✅ | ec_logs JSON + npy exact match (09 + 12B) |
| single-fit figures checked | ✅ | pixel-identical |
| TB figures visually checked | ✅ | pixel-identical |
| EC figures visually checked | ✅ | `~0.03%` scatter-point-only diff |
| representative version full-run passes | ✅ | `test_all_versions_match_oracle` (09 + 12B) |
| all historical versions config-level pass | ✅ | `test_config_resolution_matches_legacy` all 8 versions |
| all historical versions full-run pass | ✅ (except 10B EC) | 09/12B/03B/04/05B/07/11B PASS; 10B single+TB-global PASS, EC global blocked by 90 keV legacy defect |
| legacy retired | ⬜ | beyond Gate D (needs consumer confirmation) |

---

**Note on preview**: the plan's `preview` stage (single-file preview that works even
without a fit range, plus log-log/auto-peak) is a separate feature that this
phase deliberately does **not** build. The single-fit path (which is preview's
sibling but operates on confirmed ranges) is the focus of Gate B. This will not
affect the scientific-equivalence proof.
