# Known Issues in v1.2.0-review

First posted 2026-10-23; revised 2026-10-24. These issues were found in a
post-release audit of `v1.2.0-review`. Item 1 changes which tier definition
users should adopt. The other items are documentation and manuscript errors.
All of them will be corrected in the next release.

## 1. The population-prior tiers in this release are biased low: use the forward tiers

`v1.2.0-review` redefined the membership tiers on a "population-prior"
probability (`P_vgrf_below_25_popprior`; `tier` in the
`catalogues/catalogue_retier_*` files; 173 / 276 / 621 stars). A post-release
re-analysis, confirmed by a second, independent analysis, found that this
correction is **biased low by roughly a factor of two**.

The cause is a double count. The per-star Monte Carlo Vgrf histograms that
the correction uses as likelihoods are in fact posteriors under a prior that
is flat in 3-D velocity, which is proportional to v² in speed. The
phase-space volume is therefore counted twice. In mock catalogues built on
the real stars' measurement errors, the shipped correction under-recovered
the slow population by 44-72 % for every population shape tested.

The forward Monte Carlo score is calibrated to within a few per cent when
the phase-space density is flat at low speed, and the data show exactly
that. A corrected reconstruction gives an expected slow count within about
4 % of the forward value. These mock and corrected results come from the
re-analysis and will be published with the next release.

**Until the next release, use the forward-score tiers.** They ship in this
release:

- `catalogues/catalogue_expanded_tierA.fits`: Tier A, 289 stars
  (`P(Vgrf < 25 km/s) > 0.95`)
- `catalogues/catalogue_expanded_tierAB.fits`: Tier A+B, 541 stars
  (`P > 0.84`, primary sample)
- `catalogues/catalogue_expanded_tierABC.fits`: Tier A+B+C, 1,952 stars
  (`P > 0.50`)
- `catalogues/catalogue_expanded_orbits_tierABC.fits`: orbit products for
  all 1,952
- `catalogues/catalogue_expanded_master.fits`: the 20,829-star candidate pool

Equivalently, select on `tier_forward` or `P_vgrf_below_25_forward` in
`catalogues/catalogue_retier_master.fits`. Score-implied purities for the
forward tiers are 98.7 % / 94.3 % / 72.7 %.

Do not use `tier` or `P_vgrf_below_25_popprior` from the `catalogue_retier_*`
files. The population-prior claims in the bundled manuscript are withdrawn.
These include:

- the abstract's 276 / 621-star samples with 95 % / 79 % purity;
- the statement that the correction lowers a forward-defined sample's purity
  from 73 % to 42 %;
- the subsection "The Population Prior Halves the Expected Slow Count".

## 2. The bundled manuscript mixes the two tier definitions

The abstract, the Sect. 2.3 purity paragraph, Sect. 3.4 and the tables that
quote N = 276 or N = 621 describe the population-prior tiers of item 1. Some
passages instead retain values for the forward-score tiers:

- Table `tab:orbit_summary`: its caption says N = 621, but its values (median
  e 0.949, median R_peri 154 pc) are for the forward 1,952-star Tier A+B+C.
- The Conclusions orbit-summary sentence: e 0.949, R_peri 154 pc,
  R_apo 7.35 kpc.
- The sky distribution: 39.7 % toward the Galactic-centre direction.
- The RV-uncertainty defensive cut: "171 of the 420 Gold stars".
- The inner-Galaxy close-approach candidates (`tab:sgrA_candidates`), which
  are members of the forward Tier A+B+C.

The abstract's orbit summary (e 0.965, R_peri 153 pc, R_apo 8.24 kpc) is for
the population-prior Tier A+B+C. The forward-tier values are the ones listed
above. Note that the forward-tier values are subject to the other corrections
planned for the next release; see item 5 for an example.

## 3. `tier` in two phase14 products

`phase14/expanded_orbit_mc/expanded_catalogue_mc_orbits.fits` and
`phase14/expanded_selection_function.fits` carry the **forward-score** tier,
which is the recommended definition (item 1). `COLUMNS.md` says these columns
are copied from the master catalogue, whose `tier` is the population-prior
definition; that statement is wrong. The orbit Monte Carlo file covers all
1,952 forward Tier A+B+C stars.

## 4. Withdrawn excess ratios

Table `tab:uncertainty_ledger` still lists "Observed-control ratio" rows
(7.4 and 25.6). These observed-to-null excess ratios were withdrawn in
v1.2.0 and should be disregarded. The same applies to the
Interpretation-section reference to the slow sample's "count relative to
smooth velocity-field continuations of the controls".

## 5. RVS selection-function inputs

The released selection-function values (`phase14/expanded_selection_function*`)
were obtained by querying the GaiaUnlimited DR3 RVS selection function with
G_RVS and BP-RP colour. That model is defined on G and G-G_RP. The
completeness figures derived from it should be treated as approximate until
the corrected computation in the next release. The effect has not yet been
quantified.

## 6. Archive versions

The Zenodo concept DOI
([10.5281/zenodo.20116134](https://doi.org/10.5281/zenodo.20116134)) currently
resolves to `v1.0.8-review`, which uses the forward-score tiers recommended in
item 1. This repository at `v1.2.0-review` additionally ships the
population-prior columns described in item 1.
