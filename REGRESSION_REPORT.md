# Regression Report — legacy vs new (version 09)

修订日期：2026-09-08

Comparison mode: `legacy` (frozen `.oracle/09`) vs `new` (`.new/09`, produced by the
explicit workflow). All numeric comparisons are **exact** (no tolerance was
needed for any scientific field).

## Matrix

| Version | TB | EC-src | EC-xray | Final outputs | Figures | Notes |
|---|---|---|---|---|---|---|
| 09 | PASS | PASS | PASS | PASS | PASS | no tolerance exception |

## Detail

### Level 1 — single fit (per measurement / channel)
- TB: 48 pickles, 4 channels each — peak `b`, sigma `c`, resolution, errors all exact.
- EC: 16 pickles (13 x-ray + 3 src) — exact.
- single-fit figure set: 64 files, names + dimensions match; pixel-identical.
- `test_config_resolution_matches_legacy[tb|ec_source|ec_xray]` proves the resolved
  `Read_config / Spectrum_config / Fit_config` are identical to legacy for every
  measurement (the strongest equivalence argument).

### Level 2 — intermediate points
- TB points (per channel: center, center_err, temp, temp_err, bias, bias_err) match
  legacy `load_data()` exactly.
- EC global fit reads the same points; since the global fit is deterministic, the
  exact final output implies identical point content.

### Level 3 — final products
- `tb_logs/{ts}_temp_bias_fit.json` — exact.
- `ec_logs/{ts}_ec_coef_sci_ch{0..3}.json` — exact.
- `ec_logs/{ts}_ec_data_ch{0..3}.npy` — exact (np.array_equal).
- File sets + counts identical.

### Figures
- TB figures (`temp_bias_fit_*`): pixel-identical (max diff 0.0).
- single-fit figures: pixel-identical.
- EC figures (`ec_fit_*`, `resolution_fit_*`): max pixel diff 1.0 at
  `mean=0.0002–0.0023` (over 0.03%–0.34% of pixels) — scatter marker ordering
  only; fit curves identical. Accepted under the plan's "PNG need not be
  byte-identical; visually equivalent" rule.

## Tolerance exceptions

None. No scientific field required a tolerance; every float/np comparison was
exact.

## Unfixed historical bugs

See `MIGRATION_NOTES.md` (B1–B8). All preserved by the new workflow.

## How to reproduce locally

```bash
# needs the full raw data + a frozen legacy oracle
uv sync
pytest tests/ -q                # 32 tests (incl. full pipeline run)
coverage report                # new-code coverage 94%
```
