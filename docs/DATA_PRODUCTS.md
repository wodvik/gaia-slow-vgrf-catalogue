# Data Products

This repository is a private review package for the manuscript.  The final
public release will be archived on Zenodo with a persistent DOI.

## Catalogue Construction

- `release/v2/ADQL_selection.sql` - equivalent Gaia Archive parent-source
  selection.
- `release/v2/phase1/catalogue_v2.fits` - master 2,859-source candidate table.
- `release/v2/phase1/catalogue_tierAB.fits` - 334-star headline sample with
  `P(Vgrf < 25 km/s) > 0.84`.
- `release/v2/phase1/catalogue_tierA.fits` - 214-star high-confidence core with
  `P(Vgrf < 25 km/s) > 0.95`.

## Selection Function

- `release/v2/phase2/catalogue_v2_sf.fits` - per-source Gaia DR3 RVS
  selection-function values and inverse weights.
- `release/v2/phase2/sf_healpix_nside64.fits` - sky-pixel selection-function
  summary.

## Orbit Products

- `release/v2/phase4/catalogue_v2_orbits.fits` - point-estimate orbit outputs
  in the adopted Galactic potential.
- `release/v2/phase6/catalogue_mc_orbits.fits` - Monte Carlo orbit posterior
  summaries.
- `release/v2/phase6/orbits_static_dt0p1.fits` and
  `release/v2/phase6/orbits_barred_dt0p1.fits` - fine-sampled orbit diagnostics.
- `release/v2/phase3_agama/*.ini` - potential configuration files.

## Controls and Context

- `release/v2/phase5/control_orbits.fits` - matched-control orbit catalogue.
- `release/v2/phase5/geometric_mock.fits` - geometric-bias mock used for
  interpreting the eccentricity versus velocity trend.
- `release/v2/phase5/chemistry_v2.fits` - APOGEE/GALAH survey-match chemistry
  rows.
- `release/v2/phase5/populations_v2.fits` - subset-qualified chemodynamic
  context labels.

## Robustness and Diagnostics

- `release/v2/phase6/*` - solar-parameter, halo-mass, bar-pattern-speed,
  distance-method, and Monte Carlo robustness products.
- `release/v2/phase12_diagnostics/*` - retrograde and pericentre-coordinate
  diagnostics.
- `release/v2/phase14/*` - focused sensitivity products for distance coupling,
  control reweighting, central-approach interpolation, and clustering checks.

