# Release Notes

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
- Bar-pattern-speed bridger sensitivity: 9--24.
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
