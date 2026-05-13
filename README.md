# Gaia DR3 Slow-Vgrf Catalogue

This is the first public review-stage release for the probabilistic Gaia DR3
slow-Galactic-rest-frame-speed catalogue.

Zenodo DOI: [10.5281/zenodo.20116135](https://doi.org/10.5281/zenodo.20116135)

Release version: `v1.0.1-review` (see `RELEASE_NOTES.md`).

## Key Counts

- Propagated candidate pool: 20,829 stars
- Tier A (`P(Vgrf < 25 km/s) > 0.95`): 276 stars
- Tier A+B headline catalogue (`P > 0.84`): 517 stars
- Tier A+B+C orbit-summary catalogue (`P > 0.50`): 1,835 stars
- Corrected point-estimate `Vgrf < 25 km/s`: 2,591 stars

## Primary Catalogue Products

Some product filenames retain expanded as an internal build label; this is the first public release.

- `catalogues/catalogue_expanded_master.fits`
- `catalogues/catalogue_expanded_tierA.fits`
- `catalogues/catalogue_expanded_tierAB.fits`
- `catalogues/catalogue_expanded_tierABC.fits`
- `catalogues/catalogue_expanded_orbits_tierABC.fits`
- `phase14/expanded_orbit_mc/expanded_catalogue_mc_orbits.fits`
- `phase14/expanded_selection_function.fits`
- `mrt/catalogue_tierA_mrt.txt`
- `mrt/catalogue_tierAB_mrt.txt`

`COLUMNS.md` defines the public FITS/CSV/MRT columns. Gaia `source_id`
values should be treated as integer identifiers and never converted to
floating point.

## Manuscript Files

- `main.tex` and `main.pdf` are the compile-ready manuscript source and
  default LaTeX output names.
- `Humble_2026_Gaia_DR3_slow_vgrf_catalogue_manuscript.pdf` is a
  descriptive copy of the same PDF for archive browsing.

## Verification

The catalogue source IDs were independently checked against
the ESA Gaia Archive `gaiadr3.gaia_source` table by `source_id`.

- Tier A+B: 517/517 recovered, 0 discrepancies
- Tier A+B+C: 1,835/1,835 recovered, 0 discrepancies
- Master catalogue: 20,829/20,829 recovered, 0 discrepancies

The verification report and script are in `online_verification/`.

## Reproducibility

`phase14/` contains final review-stage validation and sensitivity
products. The name is a pipeline provenance label; see `phase14/README.md`.

`README-submission.md` describes which products can be reproduced from
the release package alone and which require external catalogues or the
large local Gaia parent-buffer scan. `environment.yml` records the
Python environment used for the public products; AGAMA orbit calculations
may require platform-specific installation.

## License

See `LICENSE.md`. Data and documentation are released under CC BY 4.0;
code is released under the MIT License.
