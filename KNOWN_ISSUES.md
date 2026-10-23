# Known Issues in v1.2.0-review

Last updated: 2026-10-23. These issues were found in a post-release audit of
`v1.2.0-review`. The catalogue data products themselves (the tier assignments,
membership probabilities and per-star orbit products) are not affected. The
problems are stale numbers in the bundled manuscript, two mislabelled columns,
and one documentation error. All of them will be corrected in the next
release. Every corrected value below was either recomputed from the files
shipped in this release (the catalogue products and
`catalogues/retier_summary.json`) or taken from its adopted-tier tables.

## 1. Which catalogue files are adopted

`README.md` lists the `catalogues/catalogue_expanded_*` files under "Primary
Catalogue Products". Those files carry the **superseded forward-score tiers**
(289 / 541 / 1,952). The adopted catalogue is in:

- `catalogues/catalogue_retier_master.fits`
- `catalogues/catalogue_retier_tierA.fits`
- `catalogues/catalogue_retier_tierAB.fits` (primary sample, 276 stars)
- `catalogues/catalogue_retier_tierABC.fits` (621 stars)
- `catalogues/catalogue_retier_orbits_tierABC.fits`

In these files `tier` is the adopted population-prior tier. `tier_forward`
and both probability columns are included as well (see `COLUMNS.md`).

## 2. `tier` in two phase14 products is the forward definition

In `phase14/expanded_orbit_mc/expanded_catalogue_mc_orbits.fits` and
`phase14/expanded_selection_function.fits`, the `tier` column (and the
accompanying probability) is the **forward-score** tier. That contradicts
`COLUMNS.md`, which says the column is copied from the master catalogue. To
get the adopted tiers, join on `source_id` to
`catalogues/catalogue_retier_master.fits`. The orbit Monte Carlo file covers
the 1,952 forward Tier A+B+C stars, which include all 621 adopted Tier A+B+C
stars.

## 3. Forward-era numbers that survive in the bundled manuscript

The abstract and the membership sections use the adopted tiers. The
following passages still carry values from the superseded 1,952-star
forward-score sample:

| Location in `main.tex` | Printed (forward-era) | Adopted Tier A+B+C (N = 621) |
|---|---|---|
| Table `tab:orbit_summary` (all rows; caption says N = 621) | median e 0.949; median R_peri 154 pc (p16-p84 89-274) | median e 0.965 (0.939-0.981); median R_peri 153.5 pc (83-255); median R_apo 8.24 kpc (7.11-9.62), max 14.58 kpc; median z_max 3.97 kpc |
| Conclusions, orbit summary sentence | e 0.949, R_peri 154 pc, R_apo 7.35 kpc | e 0.965, R_peri 153 pc, R_apo 8.24 kpc (as in the abstract) |
| Sky distribution (Section 3) | 39.7 +/- 1.1 % toward the Galactic-centre direction (l in [330, 30] deg); 6.2 % toward the anti-centre (l in [150, 210] deg) | 14.2 % and 7.1 %; a uniform distribution in longitude gives 16.7 % each. The Galactic-centre concentration is a property of the forward-score sample and **does not hold** for the adopted catalogue. |
| Robustness appendix, RV-uncertainty defensive cut | "171 of the 420 Gold stars" | The Gold subset of the adopted catalogue has 226 stars (220 after the `ipd_frac_multi_peak` cut). The 171/420 figures are forward-era. |
| Table `tab:uncertainty_ledger`, rows "Observed-control ratio" (7.4 and 25.6) | observed-to-null excess ratios | These ratios were **withdrawn** in v1.2.0 and should be disregarded. The same applies to the Interpretation-section reference to the slow sample's "count relative to smooth velocity-field continuations of the controls". |
| Inner-Galaxy close-approach appendix, Table `tab:sgrA_candidates` | candidate close approachers | None of the listed candidates is in the adopted Tier A+B+C. They belong to the forward-score sample only. |

## 4. RVS selection-function inputs

The released selection-function values (`phase14/expanded_selection_function*`)
were obtained by querying the GaiaUnlimited DR3 RVS selection function with
G_RVS and BP-RP colour. That model is defined on G and G-G_RP. The
completeness figures derived from it should therefore be treated as
approximate until the corrected computation in the next release. The effect
has not yet been quantified.

## 5. Archive versions

The Zenodo concept DOI
([10.5281/zenodo.20116134](https://doi.org/10.5281/zenodo.20116134)) currently
resolves to `v1.0.8-review`, which uses the forward-score tiers. This GitHub
repository at `v1.2.0-review` is the most recent release.
