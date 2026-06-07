# Phase 14 Validation Products

This directory contains final review-stage validation and sensitivity
products for the public catalogue release. The `phase14` name is a pipeline
provenance label, not an astronomy standard.

For archive users, the important products are:

- `expanded_orbit_mc/`: full orbit Monte Carlo products for the Tier A+B+C
  catalogue.
- `expanded_potential_sensitivity_*`: halo-mass and barred-potential
  point-estimate sensitivity reruns.
- `expanded_sgrA_refinement_*`: Monte Carlo refinement of the closest
  Sgr A* approacher seeds.
- `expanded_selection_function.*`: Gaia DR3 RVS selection-function values,
  per-source parent-count diagnostics, and summaries.
- `injection_recovery/wp7_probability_calibration_*`: truth-labelled
  GeDR3mock probability-score diagnostic products.
- `external_rv/`: APOGEE DR17 and GALAH DR3 radial-velocity consistency
  audit for the sparse footprint overlap with Tier A+B+C.
- `expanded_chemistry_summary.json` and related spectroscopy files:
  chemistry cross-match/context products.
- `cluster_search_*`: local nearest-neighbour/compact-pair diagnostics.
- `filion2025_comparison.*`: external low-azimuthal-velocity comparison
  against the Filion et al. (2025) APOGEE/Gaia sample.

The manuscript uses these files as validation and robustness products; the
headline catalogue products are in `catalogues/`, with column definitions in
`COLUMNS.md`.
