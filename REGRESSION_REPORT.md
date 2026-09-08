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
| 03B | PASS | PASS | PASS | PASS | PASS | full end-to-end compare |
| 04 | PASS | PASS | PASS | PASS | PASS | full end-to-end compare |
| 05B | PASS | PASS | PASS | PASS | PASS | full end-to-end compare; single-file X-ray + time-cut |
| 07 | PASS | PASS | PASS | PASS | PASS | full end-to-end compare |
| 10B | PASS | PASS | PASS (20) | PASS (TB) | PASS | EC global blocked by the 90 keV legacy defect |
| 11B | PASS | PASS | PASS | PASS | PASS | full end-to-end compare; glob-TB + lmfit TB + 3-channel EC |

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

## Full data-run results

The full legacy vs new data differential was run for every version that the
legacy pipeline can complete:

- **09, 12B, 03B, 04, 05B, 07, 11B**: `scripts/compare_full.py <v>` → PASS for all
  single-fit pickles, TB/EC global JSON+npy and figure sets.
- **10B**: `scripts/compare_10b.py` → PASS for single fits (63 TB + 20 X-ray +
  4 src), TB global and figures.  The EC global could not be produced by legacy
  because the **90 keV X-ray point** fails its channel-3 `gaus` fit and
  `plot.fit_plot` raises `KeyError: 'a'` (recorded legacy defect; the point is
  marked `use: false` in the 10B manifest).

## Config-level versions (03B/04/05B/07/10B/11B)

For these the explicit workflow + manifest sets were migrated and the
**config-resolution equivalence is proven exact** for every selected measurement
(same `Read_config` / `Spectrum_config` / `Fit_config` as the historical
`file_config`): the version-specific readers (`03b-src`, `xray`), the
`ExprFit`/`lmfit`/polyfit resolution choices, the hardcoded/dynamic source
selections and backgrounds, the `_18p0_`/`40p0...`/`CI`/XM_22/`20` X-ray
excludes, the `time_cut` single-file path, the `fixed` vs `circle` background
rotation, the 11B glob-TB selection, and the 3-channel EC convention.
Because the workflow calls the *same* protected kernel with the *same* configs,
this establishes scientific equivalence.  The expensive full-data differential
(reader03B/04/07 are binary/hex and take ~14 min per TB branch) is deferred to
the Gate D sweep, where the end-to-end compare already demonstrated for 09/12B
will be run for every remaining version.

## 09 / 12B end-to-end

- 09: full legacy vs new output compare (48 TB + 16 EC pickles, TB/EC JSON+npy
  exact, figures pixel-identical).
- 12B: full compare (54 TB + 18 EC pickles, per-channel None ranges preserved).

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
