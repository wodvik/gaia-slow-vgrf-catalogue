# Gaia DR3 Slow-Vgrf Catalogue

This repository contains the review-stage catalogue, scripts, and derived
products for the manuscript:

**A Gaia DR3 Catalogue of Stars with Very Low Galactic Rest-Frame Speeds**

The headline sample is a probabilistic Gaia DR3 catalogue of stars with very
low Galactic rest-frame speed, defined by `P(Vgrf < 25 km/s) > 0.84`.  The
repository is currently intended for private review.  A DOI-backed Zenodo
release will be created when the manuscript and data products are frozen.

## Contents

- `paper/manuscript_draft.pdf` - current review draft of the manuscript.
- `release/v2/ADQL_selection.sql` - equivalent Gaia Archive query for the
  parent Gaia DR3 6D source selection.
- `release/data/` - compact upstream CSV inputs needed by the retained
  scripts, including the enriched slow-star table, matched control-band
  inputs, and APOGEE/GALAH cross-match tables.
- `release/v2/config.yml` - numerical settings used by the pipeline.
- `release/v2/environment.yml` and top-level `environment.yml` - software
  environment.
- `release/v2/phase1/` - catalogue construction and velocity-threshold
  membership products.
- `release/v2/phase2/` - Gaia DR3 RVS selection-function products.
- `release/v2/phase3_agama/` - adopted Galactic potential configuration files
  and checks.
- `release/v2/phase4/` - point-estimate orbit catalogue.
- `release/v2/phase5/` - matched-control, chemistry, mock, and population
  context products.
- `release/v2/phase6/` - Monte Carlo orbit and robustness products.
- `release/v2/phase12_diagnostics/` and `release/v2/phase14/` - focused
  diagnostic and sensitivity products used by the manuscript.
- `release/v2/scripts/` - Python scripts used to generate the catalogue and
  diagnostics.
- `release/v2/mrt/` - AAS-style ASCII machine-readable versions of the Tier A
  and Tier A+B headline catalogues, plus the manuscript sample-table source.
- `docs/COLUMNS.md` - column dictionary for the principal FITS products.
- `docs/STAR_COUNTS.md` - plain-English ledger of the catalogue counts.
- `docs/DATA_PRODUCTS.md` - short map from manuscript claims to files.

## Headline Files

The most useful starting points are:

- `release/v2/phase1/catalogue_tierAB.fits` - 334 headline Tier A+B stars.
- `release/v2/phase1/catalogue_tierA.fits` - 214 high-confidence Tier A stars.
- `release/v2/phase1/catalogue_v2.fits` - 2,859-source master candidate table.
- `release/v2/phase4/catalogue_v2_orbits.fits` - point-estimate orbit outputs.
- `release/v2/phase6/catalogue_mc_orbits.fits` - Monte Carlo orbit posterior
  summaries for Tier A+B+C.
- `release/v2/phase5/control_orbits.fits` - matched-control orbit catalogue.
- `release/v2/phase5/chemistry_v2.fits` - APOGEE/GALAH survey-match rows.
- `release/v2/phase5/populations_v2.fits` - chemistry-context labels.
- `release/v2/mrt/catalogue_tierAB_mrt.txt` - ASCII machine-readable version
  of the 334-star Tier A+B headline catalogue.

The chemistry products are survey-match tables, not strictly one row per Gaia
source: a small number of stars have both APOGEE and GALAH rows.  The catalogue
products used for the headline sample and orbit summaries are one row per Gaia
`source_id`.

## Quick Read Example

```python
from astropy.table import Table

cat = Table.read("release/v2/phase1/catalogue_tierAB.fits")
print(len(cat))
print(cat["source_id", "P_vgrf_below_25", "vgrf_default", "tier"][:5])
```

## Environment

Create the analysis environment with:

```bash
conda env create -f environment.yml
conda activate gaia2026-v2
```

Some scripts query public external services or require optional astronomy
packages such as `agama`, `gaiaunlimited`, and `astroquery`.  The released FITS,
CSV, and JSON products are included so readers can inspect the manuscript-level
outputs without re-running the full pipeline.

## Source Data

The source catalogues are public survey products:

- Gaia DR3 through the ESA Gaia Archive.
- Bailer-Jones et al. (2021) photogeometric distances through VizieR.
- APOGEE DR17 through SDSS.
- GALAH DR3 through the GALAH Survey.
- Gaia XP / GSP-Phot parameter products as described in the manuscript.

The raw Gaia DR3 archive mirror used during development is not included.

## Citation

This private review package should not be cited yet.  When the manuscript is
ready for submission, a tagged GitHub release and Zenodo archive will be created
and the DOI will replace this note.
