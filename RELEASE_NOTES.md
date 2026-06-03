# Release Notes

## v1.0.3-review - 2026-06-03

This layout-only update refreshes the manuscript PDF after the v1.0.2-review
archive. Catalogue files, numerical products, selection-function products, and
headline science quantities are unchanged from v1.0.2-review.

- Removed a forced page break before the single-column inner-Galaxy table block
  so the Bailer-Jones overlap paragraph is no longer stranded on an almost
  blank page.
- Kept the Sgr A* candidate table and its cumulative-probability figure together
  before the robustness section, preserving table/figure order while reducing
  visible whitespace.
- Regenerated `main.pdf` and the descriptive manuscript copy
  `Humble_2026_Gaia_DR3_slow_vgrf_catalogue_manuscript.pdf`.

## v1.0.2-review - 2026-06-01

This review-stage update packages the referee-panel provenance reconciliation
and reproduce-cycle audit after the v1.0.1 barred-orbit correction. Headline
membership counts and the default catalogue definition are unchanged.

### Reproducibility product addition

- Added per-source GaiaUnlimited DR3-RVS parent-count diagnostics to
  `analysis_products/expanded_selection_function.fits`: `sf_parent_count` and
  `sf_prior_dominated_n_lt10`. Cells absent from the GaiaUnlimited
  `dr3-rvs-nk.h5` grid are encoded as parent count zero, matching the
  `p=0.5` prior-mean fill used by `DR3RVSSelectionFunction.query(...,
  fill_nan=True)`.
- Regenerated `expanded_selection_function_summary.*` and
  `expanded_selection_function_lown_audit.csv` from those per-source columns.
  The regenerated low-parent-count fractions reproduce the manuscript table
  values to rounding; no manuscript science numbers changed.

### Reproduce cycle

Operator-directed provenance pass for the point-estimate orbit products.

- Regenerated the Tier~A+B+C point-estimate static and barred orbit catalogue
  with `trajsize=40001` and parabola-interpolated cylindrical pericentres at
  `dR=0` sign changes. The tier counts remain unchanged
  (`N_A+B=517`, `N_A+B+C=1835`) and the bridger counts remain unchanged
  (static = 3, barred Hunter/Sormani default = 22).
- The 40,001-point static readout gives median `R_peri = 0.1123 kpc` and
  `R_peri < 100 pc = 824`. The barred default gives
  `R_peri < 100 pc = 1543`, while its median pericentre remains
  readout-density sensitive (`trajsize=20001` gives 20.3 pc; `trajsize=40001`
  gives 14.1 pc). For this reason, manuscript text and tables treat barred
  median pericentres as resolution-limited diagnostics and emphasize the more
  stable threshold counts and bridger counts.
- Regenerated the matched point-estimate sensitivity products at the same
  40,001-point/interpolated readout: static halo mass variants, barred
  Hunter/Sormani pattern speeds 24/28/33/37.5/41, bar-angle sweep,
  halo-flattening sweep, and the Portail+2017 M2M bar cross-check.
- Re-seeded the expanded Sgr A* approacher audit from the regenerated orbit
  catalogue. The strict negative result remains unchanged:
  `P(r_sph,min < 10 pc) > 0.5` count = 0. The softer
  `P(r_sph,min < 100 pc) > 0.5` diagnostic is now 6 for the 10 seeded
  candidate/potential rows. The Sgr A* Monte Carlo engine itself already uses
  `trajsize=40001` and quadratic interpolation of `r_sph^2(t)`.
- Audited the full orbit-Monte-Carlo layer (`phase14x`): it remains the
  separate static-potential, 5,000-realisation-per-star product with
  `trajsize=1001`, giving the headline posterior-median values
  `median R_peri = 151 pc`, `median R_apo = 7.45 kpc`, and
  `median e = 0.949`. No sampling-sensitive deep-tail MC statistic from this
  product is quoted as a headline, so the reproduce cycle documents this
  sampling rather than rerunning the full MC layer.

### Referee-panel Round 1 provenance reconciliation

Provenance reconciliation of the barred orbit summary (panel decisions D2R / 9c648515).

- Historical note: the barred values in `catalogues/expanded_orbit_summary.json`
  were stale relative
  to the corrected v1.0.1 barred product: the sidecar JSON reported the pre-correction
  barred figures (median R_peri = 0.0556 kpc; R_peri < 100 pc = 1305; bridgers = 32),
  whereas the shipped `catalogues/catalogue_expanded_orbits_tierABC.fits`, the manuscript
  tables (tab_reach, tab_sensitivity), and the v1.0.1 corrections below all use the
  corrected barred default (median R_peri = 0.0676 kpc; R_peri < 100 pc = 1226;
  bridgers = 22). The JSON has been regenerated from the shipped FITS so the sidecar now
  agrees with the catalogue and manuscript.
- Superseded trajectory-size provenance: the subsequent reproduce-cycle addendum
  above reconciles the point-estimate generators and regenerated products at
  `trajsize=40001` with interpolated pericentres.

## v1.0.1-review - 2026-05-13

This review-stage update corrects the barred Hunter+2024 orbit-integration
convention used in the barred sensitivity products.

### Changed

- Barred integrations now use the fixed present-day Hunter+2024 full-bar
  potential with `Omega` supplied to `agama.orbit`, avoiding double
  application of the rotating-frame convention.
- Regenerated the Tier A+B+C point-estimate orbit catalogue.
- Regenerated the expanded potential-sensitivity products.
- Regenerated the expanded Sgr A* approacher Monte Carlo refinement.
- Updated manuscript tables and prose to match the regenerated products.
- Rebuilt `main.pdf` from the updated manuscript source.

### Updated Numerical Values

- Default barred bridgers: 22.
- Bar-pattern-speed bridger sensitivity: 6--24 (across $\Omega_p \in \{24, 28, 33, 37.5, 41\}$ km/s/kpc).
- Barred `R_peri < 100 pc`: 1,226.
- Barred `min r_sph < 10 pc`: 0.
- Sgr A* strict `P(r_sph,min < 10 pc) > 0.5`: 0.
- Sgr A* soft `P(r_sph,min < 100 pc) > 0.5`: 3.

### Unchanged Interpretation

The default manuscript interpretation remains anchored to the static
Hunter+2024 axisymmetric potential. The barred products remain sensitivity
tests for inner-Galaxy geometry rather than the default catalogue definition.

### Zenodo Note

If this package is deposited to Zenodo as an update, publish it as a new
version of the existing record and update any version-specific DOI fields
after Zenodo mints the new version DOI. Use the concept DOI when citing the
evolving catalogue across versions.
