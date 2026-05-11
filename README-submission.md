# Submission Reproducibility Runbook

This runbook describes the expanded review-stage release for the
probabilistic Gaia DR3 slow-Vgrf catalogue. It is intended for referees
and archive curators who want to understand which products can be
reproduced from the reviewer bundle alone and which require external raw
catalogues.

## Release Identity

- Release tag: `expanded-fullmc-20260510`
- Manuscript: `main.tex`, compiled to `main.pdf`
- Reviewer bundle: `_review_bundle_expanded_20260510_fullmc.zip`
- Headline catalogue: Tier A+B, `P(Vgrf < 25 km/s) > 0.84`, 517 stars
- Broader orbit-summary catalogue: Tier A+B+C, `P > 0.50`, 1,835 stars
- Expanded propagated candidate pool: 20,829 stars
- Broad parent buffer scanned from local Gaia DR3 mirror: 5,867,654 unique source IDs

## Environment

Create the software environment from:

```bash
conda env create -f environment.yml
conda activate gaia2026-v2
```

AGAMA was run from WSL for the orbit calculations. The full expanded
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
- `phase14/expanded_orbit_mc/`: full expanded orbit Monte Carlo products.

## What Can Be Reproduced From The Bundle Alone

The bundle contains the released catalogue products, manuscript source,
figures, tables, orbit-MC summaries, chemistry cross-match products, and
validation summaries. From the bundle alone, a reviewer can:

1. Inspect and join the FITS/CSV catalogues by `source_id`.
2. Verify the headline and broader tier counts.
3. Recompute manuscript summary statistics from the released tables.
4. Recompile the manuscript PDF from `main.tex`.
5. Inspect the full 5,000-realisation orbit-MC output for all 1,835
   Tier A+B+C stars and the 10,000-realisation convergence subset.

## What Requires External Data

The full raw-to-release rebuild requires external public catalogues and
one large local intermediate product:

- ESA Gaia Archive / Gaia DR3 `gaia_source`
- Bailer-Jones et al. (2021) EDR3 distance catalogue
- APOGEE DR17 and GALAH DR3 spectroscopic catalogues
- Gaia XP/GSP-Phot parameter products used for context
- The broad local parent-buffer CSV, omitted from the reviewer zip for
  size: `D:/GAIA/parent_scan/gaia_parent_buffer_vgrf200_full.csv`

The parent-buffer CSV is reproducible from the local Gaia DR3 source
mirror with the Phase 0 scanner. It is not needed to inspect or use the
released catalogue products.

## Production Command Ledger

The expanded parent-buffer scan was run against the local Gaia DR3 CSV
mirror at `D:/GAIA/csv`, writing chunk-by-chunk to:

```text
D:/GAIA/parent_scan/
```

On another machine, remap the local paths in `config.yml` and command
arguments before rerunning scripts. The production defaults document this
release environment: Windows paths such as `D:/GAIA/...` correspond to
WSL paths such as `/mnt/d/GAIA/...`, and
`C:/Users/humbl/GAIA2026` corresponds to `/mnt/c/Users/humbl/GAIA2026`.
No released FITS/MRT catalogue requires those paths for normal use.

The final full expanded orbit Monte Carlo was run from WSL with:

```bash
cd /mnt/c/Users/humbl/GAIA2026
OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 release/v2/scripts/phase14x_expanded_mc_orbits.py \
  --out-dir /mnt/c/Users/humbl/GAIA2026/release/v2/phase14/expanded_orbit_mc \
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

## Known Limitations

The reviewer bundle is a reproducibility and inspection package, not a
full redistribution of the Gaia DR3 source table. Some scripts retain
local default paths documenting the exact production environment; the
inputs and outputs are also listed explicitly in `config.yml` so they can
be remapped on another machine.

The expanded catalogue has been evaluated with GaiaUnlimited's
Castro-Ginard et al. (2023) Gaia DR3 RVS selection function. The per-star
table is `phase14/expanded_selection_function.fits`, with summaries in
`phase14/expanded_selection_function_summary.csv` and
`phase14/expanded_selection_function_summary.json`. Inverse-selection
weighted sums are retained as contextual observability diagnostics, not
as volume-complete headline population counts.
