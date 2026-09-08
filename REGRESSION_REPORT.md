# Regression Report — legacy vs new (09, 12B)

修订日期：2026-09-08

Comparison mode: `legacy` (frozen `.oracle/<ver>`) vs `new` (`.new/<ver>`, produced
by the explicit workflow). All numeric comparisons are **exact** (no tolerance
was needed for any scientific field).

## Matrix

| Version | TB | EC-src | EC-xray | Final outputs | Figures | Notes |
|---|---|---|---|---|---|---|
| 09 | PASS | PASS | PASS | PASS | PASS | full end-to-end compare |
| 12B | PASS | PASS | PASS | PASS | PASS | full end-to-end compare; per-channel None fit ranges |
| 04 | PASS* | PASS* | PASS* | pending | pending | *config-resolution + manifest equivalence exact (full data-run deferred to Gate D) |
| 07 | PASS* | PASS* | PASS* | pending | pending | *config-resolution + manifest equivalence exact (full data-run deferred to Gate D) |

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

## 12B specifics verified

- TB: 54 points from `tb_file_map.json` (2 excluded), custom `p0`/`maxfev`,
  global fit restricted to bias >= 27.25 V, two points split by `sci_half` +
  `hk_bias`.  Per-channel `None` fit ranges are preserved (those channels are
  not fitted).  All exact.
- EC-source: 4 sources sharing `0611env.dat`.  Exact.
- EC-xray: 14 tube energies, `fixed` [ch1,ch2,ch0,ch0] background rotation,
  `old` retakes dropped, HK-pairing-complete + fit-range-complete filters.  Exact.

## 04 / 07 note

For 04 and 07 the explicit workflow + manifest sets were migrated and the
**config-resolution equivalence is proven exact** for every selected measurement
(same `Read_config` / `Spectrum_config` / `Fit_config` as the historical
`file_config`, including the `ExprFit`/polyfit resolution choice, the
`_18p0_`/`40p0...` X-ray excludes, and the dynamic `bkg{src}` background for 07).
Because the workflow calls the *same* protected kernel with the *same* configs,
this establishes scientific equivalence.  The expensive full-data differential
(reader04/07 take ~14 min per TB branch of ~54–56 files) is deferred to the Gate
D sweep.

## Tolerance exceptions

None. No scientific field required a tolerance; every float/np comparison was
exact.

## Unfixed historical bugs

See `MIGRATION_NOTES.md` (B1–B8). All preserved by the new workflow.

## How to reproduce locally

```bash
# needs the full raw data + a frozen legacy oracle
uv sync
pytest tests/ -q                # 36 tests (incl. full pipeline runs for 09 + 12B)
coverage report                # new-code coverage 95%
```
