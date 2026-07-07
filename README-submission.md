# Submission Reproducibility Runbook

This runbook describes the first public review-stage release for the
probabilistic Gaia DR3 slow-Vgrf catalogue. It is intended for referees
and archive curators who want to understand which products can be
reproduced from the reviewer bundle alone and which require external raw
catalogues.

## Release Identity

- Release tag: `v1.0.8-review`
- Manuscript: `main.tex`, compiled to `main.pdf`
- Reviewer bundle: `gaia_slow_vgrf_catalogue_v1.0.8_review.zip`
- Zenodo concept DOI: `10.5281/zenodo.20116134`
- Primary catalogue: Tier A+B, `P(Vgrf < 25 km/s) > 0.84`, 541 stars
- Broader orbit-summary catalogue: Tier A+B+C, `P > 0.50`, 1,952 stars
- Propagated candidate pool: 20,829 stars
- Broad parent buffer scanned from local Gaia DR3 mirror: 5,867,654 unique source IDs

## Environment

Create the software environment from:

```bash
conda env create -f environment.yml
conda activate gaia2026-v2
```

AGAMA was run from WSL for the orbit calculations. The full
orbit Monte Carlo used `OMP_NUM_THREADS=32`.

AGAMA can require platform-specific compilation. If the pip install in
`environment.yml` is not sufficient on a new machine, install AGAMA using
the upstream instructions for that platform and then rerun the orbit
steps from the same conda environment.

## Included Reproducibility Files

- `config.yml`: central numerical configuration and release paths.
- `environment.yml`: conda/pip dependency specification.
- `ADQL_selection.sql`: Gaia Archive query equivalent to the Gaia DR3
  6D parent-source selection before local buffering and downstream
  processing.
- `COLUMNS.md`: data dictionary for the public catalogue products.
- `LICENSE.md`: data and code license statement.
- `scripts/`: review-stage pipeline scripts.
- `catalogues/`: FITS/CSV catalogue products.
- `potentials/`: AGAMA potential `.ini` files used by the shipped orbit
  rerun scripts.
- `phase14/expanded_orbit_mc/`: full orbit Monte Carlo products.
- `phase14/injection_recovery/`: compact GeDR3mock injection-recovery
  summaries and probability-score diagnostics; the public rerun scripts are
  `scripts/wp7_injection_recovery.py` and
  `scripts/phase14ak_probability_calibration.py`.
- `phase14/external_rv/`: APOGEE DR17 and GALAH DR3 radial-velocity
  consistency audit products; `scripts/phase14al_external_rv_check.py`
  rebuilds them by querying public VizieR TAP tables.
- `phase14/README.md`: explanation of the final validation/sensitivity
  product directory and its pipeline-stage naming.
- `RELEASE_NOTES.md`: review-stage version history and corrected barred
  sensitivity values.

## Inspect From The Bundle Alone

The bundle contains the released catalogue products, manuscript source,
figures, tables, orbit-MC summaries, chemistry cross-match products, and
validation summaries. From the bundle alone, a reviewer can:

1. Inspect and join the FITS/CSV catalogues by `source_id`.
2. Verify the primary and broader tier counts.
3. Recompute manuscript summary statistics from the released tables.
4. Recompile the manuscript PDF from `main.tex`.
5. Inspect the full 5,000-realisation orbit-MC output for all 1,952
   Tier A+B+C stars and the 10,000-realisation convergence subset.
6. Inspect the truth-labelled mock score diagnostic and the sparse external
   APOGEE/GALAH radial-velocity cross-check.

## Reuse From The Bundle Alone

The public FITS/CSV/MRT products are ready for downstream science without a
Gaia Archive rescan. The following checks and product-level operations are
repo-relative:

```bash
cd gaia_slow_vgrf_catalogue_v1.0.8_review
make validate-release                      # full pre-submission validation
python tests/smoke_regression.py --bundle-root .
python scripts/make_mrt_tables.py --bundle-root .
```

`make validate-release` (equivalently `python scripts/validate_release.py`)
chains three integrity checks: the primary-count smoke test, the tier-count
provenance guard (`scripts/check_provenance.py`, which also scans the
manuscript and tables for stale counts), and a SHA-256 verification of the
released figures, tables, catalogues, MRT products, and bundled potential
configuration files against
`release_checksums.sha256` (regenerate with `make checksums`). A pinned pip
dependency lock is provided in `requirements-lock.txt` (conda-forge and
source-compiled dependencies such as AGAMA are specified in the open
`environment.yml` and noted in the lock), and the omitted large
parent-buffer intermediate is
documented with row counts, velocity partitions, and reconstruction
instructions in `private_inputs/parent_buffer_manifest.json`.

The smoke harness pins the final v1.0.8-review counts
(`master=20,829`, `Tier A=289`, `Tier A+B=541`,
`Tier A+B+C=1,952`, corrected point-estimate `<25 km/s=2,755`).

## Full Rerun Boundary And Private Inputs

The released bundle supports inspection, manuscript recompilation, checksum
verification, MRT regeneration, and candidate-to-release reruns from the
shipped products once AGAMA is installed. The full raw Gaia-to-release rebuild
has a narrower boundary: it requires external public catalogues and one large
local intermediate product:

- ESA Gaia Archive / Gaia DR3 `gaia_source`
- Bailer-Jones et al. (2021) EDR3 distance catalogue
- APOGEE DR17 and GALAH DR3 spectroscopic catalogues
- Gaia XP/GSP-Phot parameter products used for context
- The broad local parent-buffer CSV, omitted from the reviewer zip for
  size: `D:/GAIA/parent_scan/gaia_parent_buffer_vgrf200_full.csv`

The parent-buffer CSV is reproducible from the local Gaia DR3 source
mirror with the Phase 0 scanner. It is not needed to inspect or use the
released catalogue products. The downstream 20,829-row
`private_inputs/expanded_candidates_mc_tiered.csv` is shipped, so orbit,
action, threshold-table, and candidate-to-release validation scripts do not
require a Gaia mirror.

Portable paths for raw private inputs and shipped rerun inputs are listed in
`config.yml` under `input.*` and default to `private_inputs/...`. The exact production
locations are retained only under `input.production_private_inputs`:

| Production path | Portable placeholder | Purpose |
|---|---|---|
| `D:/GAIA/csv` (`/mnt/d/GAIA/csv` under WSL) | `private_inputs/gaia_dr3_source_mirror` | Local Gaia DR3 source mirror for the parent scan. |
| `D:/GAIA/parent_scan/gaia_parent_buffer_vgrf200_full.csv` | `private_inputs/gaia_parent_buffer_vgrf200_full.csv` | Non-redistributed broad parent buffer. |
| `D:/GAIA/parent_scan/expanded_candidates_mc_tiered.csv` | `private_inputs/expanded_candidates_mc_tiered.csv` | Materialized MC-tiered expanded candidate table for full reruns. |

## Full-Rerun Command Ledger

The parent-buffer scan was run against the local Gaia DR3 CSV
mirror at `D:/GAIA/csv`, writing chunk-by-chunk to:

```text
D:/GAIA/parent_scan/
```

On another machine, either place equivalent files under `private_inputs/` or
pass explicit `--input-csv` / `--out-dir` arguments. No released FITS/MRT
catalogue requires those paths for normal use.

The final full orbit Monte Carlo was run from WSL with:

```bash
cd gaia_slow_vgrf_catalogue_v1.0.8_review
OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python scripts/phase14x_expanded_mc_orbits.py \
  --input-csv private_inputs/expanded_candidates_mc_tiered.csv \
  --out-dir phase14/expanded_orbit_mc \
  --n-samp 5000 \
  --chunk-size 25 \
  --n-convergence 10000 \
  --convergence-stars 100 \
  --convergence-chunk-size 10
```

The output files are:

- `phase14/expanded_orbit_mc/expanded_catalogue_mc_orbits.fits`
- `phase14/expanded_orbit_mc/expanded_catalogue_mc_orbits.csv`
- `phase14/expanded_orbit_mc/expanded_catalogue_mc_orbits_summary.json`
- `phase14/expanded_orbit_mc/expanded_catalogue_mc_orbits_convergence_10000.csv`

The corrected barred-potential sensitivity products were regenerated with:

```bash
python scripts/phase0g_expanded_orbits.py \
  --input-csv private_inputs/expanded_candidates_mc_tiered.csv
python scripts/phase14v_expanded_potential_sensitivity.py
python scripts/phase14q_expanded_sgrA_refinement.py \
  --n-samp 5000 --chunk-size 500 --n-per-potential 5
```

These commands update the point-estimate orbit catalogue, potential
sensitivity tables, and Sgr A* approacher refinement products used by the
current manuscript.

The bundled AGAMA potential definitions are in `potentials/`; the public orbit
scripts prefer those files and fall back to the author's private production
work directory only when run inside the original development tree.

The GeDR3mock injection-recovery and probability-score diagnostics are shipped
as:

```bash
python scripts/wp7_injection_recovery.py
python scripts/phase14ak_probability_calibration.py
```

It can reuse cached mock products if present under
`phase14/injection_recovery/`; otherwise it queries the public GeDR3mock TAP
service and rebuilds the recovery table/figure.

The external radial-velocity audit can be rebuilt with:

```bash
python scripts/phase14al_external_rv_check.py
```

This step queries public APOGEE DR17 and GALAH DR3 tables through VizieR TAP.

## Known Limitations

The reviewer bundle is a reproducibility and inspection package, not a full
redistribution of the Gaia DR3 source table. Private input paths are documented
explicitly in `config.yml`; public script defaults are bundle-relative where
the required input can be supplied locally.

The final Vc=229 km/s catalogue counts are adopted in this bundle. The
version DOI remains pending until Zenodo mints the new version-specific DOI;
the concept DOI above should be used until that release freeze is complete.

The catalogue has been evaluated with GaiaUnlimited's
Castro-Ginard et al. (2023) Gaia DR3 RVS selection function. The per-star
table is `phase14/expanded_selection_function.fits`, with summaries in
`phase14/expanded_selection_function_summary.csv` and
`phase14/expanded_selection_function_summary.json`. Inverse-selection
weighted sums are retained as contextual observability diagnostics, not
as volume-complete primary population counts.

The GeDR3mock probability-score table is a truth-labelled mock diagnostic, not
an empirical calibration of the real Gaia DR3 population prior. The APOGEE and
GALAH radial-velocity audit is intentionally conservative but sparse and
footprint-selected, so it is used as a robustness cross-check rather than as a
whole-catalogue replacement for Gaia DR3 RVS.
