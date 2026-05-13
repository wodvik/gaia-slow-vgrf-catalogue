# Analysis Products

This directory contains final review-stage validation and sensitivity
products for the public catalogue release.

For archive users, the important products are:

- `expanded_orbit_mc/`: full orbit Monte Carlo products for the Tier A+B+C
  catalogue.
- `expanded_potential_sensitivity_*`: halo-mass and barred-potential
  point-estimate sensitivity reruns.
- `expanded_sgrA_refinement_*`: Monte Carlo refinement of the closest
  Sgr A* approacher seeds.
- `expanded_selection_function.*`: Gaia DR3 RVS selection-function values
  and summaries.
- `expanded_chemistry_summary.json` and related spectroscopy files:
  chemistry cross-match/context products.
- `cluster_search_*`: local nearest-neighbour/compact-pair diagnostics.

The manuscript uses these files as validation and robustness products; the
headline catalogue products are in `catalogues/`, with column definitions in
`COLUMNS.md`.
