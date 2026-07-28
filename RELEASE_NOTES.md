# Release Notes

## v1.2.0-review - 2026-07-27

**Catalogue redefinition.** Tier membership is now defined on a
population-prior-corrected probability rather than the forward Monte Carlo
score. Adopted counts change from Tier A / A+B / A+B+C = 289 / 541 / 1,952 to
**173 / 276 / 621**. Both probabilities and both tier labels ship, so either
definition is reproducible from the released files.

### Why the tiers changed

The forward score answers "what fraction of this star's error realisations
fall below 25 km/s". It carries no information about the parent population.
Because the candidate pool rises steeply with `Vgrf` (188 sources below
10 km/s against 7,635 between 40 and 50), far more stars are available to
scatter *down* into the slow window than up out of it - the classical
Eddington bias. Scored correctly, the previously released Tier A+B and
Tier A+B+C samples have purities of 80% and 42%, not the 94% and 73% implied
by the forward scores. Those purity figures were wrong, which is what
motivated the redefinition rather than a change of presentation.

The correction is **perfectly nested**: no star enters any tier it was not
already in (Spearman 0.997, Kendall 0.963), so the ranking of the catalogue is
unchanged and only the thresholds move. Equivalently, the forward cuts of
0.95 / 0.84 / 0.50 correspond to population-prior cuts of 0.987 / 0.957 /
0.813. On the new definition the tier purities are 99.0% / 95.4% / 78.6%.

### New analysis

- Latent `Vgrf`-distribution reconstruction by iterative deconvolution
  (Richardson 1972; Lucy 1974; equivalently EM, Dempster et al. 1977) with
  smoothed-EM regularisation (Silverman et al. 1990), following the recent
  application of this approach by Banik et al. (2026) to the oldest-star age
  extremum. Expected `N(Vgrf<25)` falls from 2,556 to 1,234 (+110/-146).
  Scripts `phase16b`-`phase16e`.
- The reconstructed latent distribution passes smoothly through the 25 km/s
  boundary (curvature of `ln P` below 0.02 per bin across 20-30 km/s against a
  maximum of 0.08 elsewhere), confirming the catalogue is a tail selection
  rather than a distinct population.
- Quality-cut orthogonality audit (`phase16a`): pass rate of each cut as a
  function of `Vgrf` across the full 20,829-source pool, demonstrating the cuts
  are not circular. New table `tab_cut_orthogonality`.
- New quality cut `ipd_frac_multi_peak <= 1`, targeting marginally resolved
  companions between the RUWE and non-single-star regimes.

### Removed

- The observed/predicted "excess" ratios (25.6x / 7.4x / 11.0x), together with
  `tab_null_models` and `fig_ndf_expectation`. The ratio had no defensible
  parent normalisation: the matched-control library is band-balanced by
  construction, so its sampling rate relative to the Gaia 6D catalogue varies
  by more than two orders of magnitude across velocity bands, and the
  prediction was anchored to 10,078 control stars while the numerator counted
  observed slow stars from the complete parent scan. The matched-control
  comparison is retained as a distributional (shape) test, where no
  normalisation is required. Both files remain in the bundle but are no longer
  cited by the manuscript.

### Regenerated

- 11 figures on the new tiers; `tab_selection_function`,
  `tab_covariance_stress`; all four MRT tables (173 / 276 / 621 / 621 rows).
- Selection-function diagnostics re-aggregated from the shipped per-source
  table (no GaiaUnlimited re-query needed, the retier being nested):
  Tier A+B `f_n<10` = 81.3%, `f_sum,n<10` = 87.5%.
- Covariance stress restated: **276 of 276** Tier A+B members remain above
  P = 0.84 under maximal copula coupling (median |dP| = 0.0010), stronger than
  the forward-era 527 of 541 because the new thresholds sit further from the
  boundary.
- Orbit summaries re-aggregated over the nested subsets; no orbit was
  re-integrated. Monte Carlo posterior medians: eccentricity 0.965 (was
  0.949), `R_peri` 153 pc (was 154), `R_apo` 8.24 kpc (was 7.35).

### Fixed

- `phase16f` initially subset the orbit table's rows without relabelling its
  `tier` column, leaving forward labels inside the retiered row set. Totals
  were correct but per-tier splits, legends and marker styles were not. Caught
  by inspecting a rendered figure; an assertion now guards it.
- Population-prior probabilities are clipped to [0, 1]; the EM normalisation
  had left 22 values an epsilon above unity.
- `phase14t` wrote absolute local paths into a released JSON; now
  bundle-relative.
- `check_provenance` validated only one tier definition; it now checks both and
  asserts the retier is a strict nested subset of the forward tiers.

### Reporting

- Chemistry is reported on both samples: the 51-star population-prior subset
  and the 125-star forward-defined subset. Literature-region proportions are
  stable across the size change (Splash/GSE/Aurora/disc/unclassified =
  22/9/10/1/6 against 51/23/22/2/19).

## v1.1.0-review - 2026-07-17 (draft)

Manuscript-clarification revision begun from the complete
v1.0.9-review bundle.  This entry will be expanded as the revision
proceeds.

- Figure 6 now states explicitly that a plotted pericentre is the
  minimum Galactocentric-radius point on a continuous trajectory, not a
  physical reversal of the star's velocity at the Galactic centre.
- The manuscript now defines the orbit integrations as fully
  three-dimensional while distinguishing the cylindrical turning-point
  summaries $R_{\rm peri}=\min\sqrt{x^2+y^2}$ and
  $R_{\rm apo}=\max\sqrt{x^2+y^2}$ from the spherical Galactic-centre
  closest approach $r_{\rm sph,min}=\min\sqrt{x^2+y^2+z^2}$.
- The abstract, introductory interpretation, orbit methods, primary
  results, figure and table captions, limitations, conclusions, and
  Sgr A* appendix now use this terminology consistently.  Small
  $R_{\rm peri}$ is described as Galactic spin-axis reach rather than,
  by itself, a close passage to Sgr A*.
- The Sgr A* candidate calculation and all numerical values are
  unchanged: its seeds and Monte Carlo probabilities already use the
  three-dimensional spherical closest approach.  Auto-generated table
  sources were updated so future regeneration preserves the clarified
  terminology.

## v1.0.9-review - 2026-07-11

Manuscript-format release. No data, catalogue, script, or numerical
changes of any kind; primary counts are unchanged (Tier A / A+B /
A+B+C = 289 / 541 / 1,952).

- The bundled manuscript source is now typeset with the Astronomy &
  Astrophysics document class (`aa.cls` v9.4) instead of `aastex631`.
  `main.tex` compiles standalone within this bundle (30 pages, 0
  errors) and uses the A&A structured abstract
  (Context/Aims/Methods/Results/Conclusions).
- Release metadata is now journal-agnostic: the release description no
  longer names a target journal or a manuscript tracking number.
  Historical entries below are retained verbatim.
- `aastex631.cls` removed from the bundle; `aa.cls`, `lineno.sty`, and
  `linenoaa.sty` added.

## v1.0.8-review - 2026-07-06

Manuscript-typesetting fix release, cut on the day the paper was
submitted to The Astronomical Journal (manuscript AAS78818). No data,
catalogue, script, or numerical changes of any kind.

- The aastex631 `acknowledgments` environment mis-measures its height
  when the `linenumbers` class option is active, overprinting the
  `\facilities`/`\software` block onto the GALAH acknowledgement
  paragraph. The environment is replaced with a plain
  `\section*{Acknowledgments}`; pages 26-27 verified clean by render.
- The `linenumbers` class option (required by AAS at submission) is now
  enabled in the shipped manuscript source; the compiled manuscript is
  37 pages.
- The bundled `main.pdf` and descriptive manuscript copy now match the
  PDF submitted to the journal.

## v1.0.7-review - 2026-07-04

Pre-submission touch-up pass following a second independent numerical
audit of the v1.0.6 bundle (14 check groups recomputed from the shipped
products; all headline and load-bearing numbers reproduced, including the
previously pending Anderson-Darling p<0.001 re-run, which confirms at
A^2 = 531 from the released eccentricity arrays). Primary catalogue
counts are unchanged (Tier A / A+B / A+B+C = 289 / 541 / 1,952); no
science, figure, or catalogue-file changes.

**Manuscript fixes** (`main.tex`, `tables/`):
- 8-Gyr convergence check: stale sample-median deltas replaced with the
  values recomputed from `phase14/wp5_8gyr_convergence_per_star.csv`
  (+5.5 pc R_peri, +0.017 kpc R_apo, +0.039 kpc z_max; were +1.6 pc /
  +0.010 / +0.060).
- Bar-induced chaos: "one of three reach candidates" corrected to "none
  of the three" (shipped lambda*T_orb = 0.63 / 0.49 / 0.00, none
  order-unity; `wp5_chaos_per_star.csv` and `barred_ftlyap_norm` agree).
- tab_gold_dynamical: Tier A+B median z_max 3.92 -> 3.93 (recomputed
  3.9345).
- tab_potential: barred e>0.95 cells Omega_p=24: 93.7 -> 93.9 (93.85)
  and Omega_p=28: 93.0 -> 92.7 (92.67), matching
  `expanded_potential_sensitivity_summary.csv`.
- Anderson-Darling vs single-ellipsoid mock: sentence now discloses that
  the mock eccentricity realisations are not archived in the release
  bundle (the control-band comparison is fully re-runnable).
- NSS wording: "predominantly SB1" scoped to the full Tier A+B+C
  extension (Tier A+B alone splits 4 SB1 / 4 Orbital /
  1 AstroSpectroSB1).

**Packaging / metadata**:
- COLUMNS.md: new "Matched-Control Reweighting Products" section pointing
  reproducers to `phase14/expanded_control_weights.csv`
  (`expanded_kw_weight`) as the operative kernel weight -- the legacy
  `kw_weight` column inside `control_orbits.fits` does not reproduce the
  manuscript's kernel-weighted medians.
- Version bumped to 1.0.7-review (CITATION.cff, config.yml, READMEs);
  checksums and review zip regenerated; `main.pdf` recompiled.

**Presentation pass** (same day, from a fresh-eyes referee read; no
science or catalogue-file changes):
- Every figure and table is now called out in the text (8 figures and 8
  tables previously had no in-text reference); the one-row Bailer-Jones
  overlap table is folded into a sentence (table count 33 -> 32).
- The accidental full-width prose block in the inner-reach discussion
  (\onecolumngrid around five paragraphs) is restored to two-column
  typesetting; `tab_reach` converted to a regular single-column float.
- Mandatory survey acknowledgements added (Gaia/ESA DPAC, SDSS-IV/APOGEE
  DR17, GALAH DR3 incl. AAT program list), plus `\facilities`,
  `\software`, `\correspondingauthor`, and `\shorttitle`/`\shortauthors`
  running heads; the aastex631 acknowledgments line-number quirk is
  suppressed with `\nolinenumbers`.
- Gehrels (1986) added to the bibliography (cited three times but
  missing); the reference list is now a single alphabetized sequence
  (93 entries); Wright & Binney (2026) entry de-cluttered.
- Conclusions wording: "calibrated probabilistic sample" ->
  "probability-scored sample", consistent with the Section 2.3
  calibration disclaimer.
- `fig_ndf_expectation` regenerated: in-plot annotation now reads
  N_exp(<25)=76.4 (was a stale 76.5 disagreeing with text and table);
  the shipped `phase14ab` script now renders its `N_EXP_25` constant.
- `fig_chaos_vs_rperi` regenerated from the shipped
  `wp5_chaos_per_star.csv`: the stale in-plot "Table 15" legend is
  replaced by table-number-free labels and the internal "WP-5" tag is
  dropped from the title.
- Assorted wording fixes: circular sparsity sentence tightened;
  "low-speed ball" -> "low-speed region"; dangling Section-4.5 forward
  pointers now target Appendix A (Table solar_variants) with a new
  Solar-parameter summary paragraph; Gold dynamical pericentre row
  converted from kpc to pc.
- The compiled manuscript is now 38 pages (survey acknowledgements add
  roughly one page).

**Smooth-DF eccentricity mock shipped** (same day): the second
Anderson-Darling comparison (Tier A+B+C eccentricities vs the
single-ellipsoid smooth-continuation mock) is now fully re-runnable from
bundle products. `phase14/control_orbits.fits` gains the Cartesian
velocity columns `vx_kms`/`vy_kms`/`vz_kms`; the new
`scripts/phase14am_smooth_df_ecc_mock.py` regenerates the mock
realisation (`phase14/smooth_df_eccentricity_mock.csv`, seed 20260523)
and writes `phase14/smooth_df_ecc_mock_summary.json` with both AD
statistics (slow vs 25-50 control: A^2 = 531.1; slow vs mock:
A^2 = 300.2; both asymptotic p < 0.001, matching the manuscript). The
manuscript parenthetical now points to the shipped realisation instead
of disclosing its absence.

## v1.0.6-review - 2026-06-28

Pre-submission correctness, consistency, and packaging pass over the v1.0.5
bundle following an independent end-to-end audit. Primary catalogue counts
are unchanged (Tier A / A+B / A+B+C = 289 / 541 / 1,952); no scientific
result or figure changes.

**Manuscript fixes** (`main.tex`):
- Abstract labels the headline triple (ecc 0.949, R_peri 154 pc, R_apo
  7.35 kpc) as Monte Carlo posterior medians, distinguishing it from the
  point-estimate median R_peri (115.5 pc) used in the body.
- Pericentre readout made internally consistent: the point-estimate methods
  text now states the 40,001-point parabola-interpolated readout (was
  "2,001 stored samples"); Gold-subset pericentres corrected to 114.8 pc
  (Gold) / 112.5 pc (full Tier A+B) (were 121 / 120 pc); the face-on caption
  Tier A+B+C pericentre corrected to 115.5 pc (was 112 pc).
- Eccentricity Anderson-Darling comparisons now state rejection at
  asymptotic p<0.001 without the specific A^2 statistics (192.27 / 36.97),
  which did not reproduce from the released eccentricity arrays (pending
  re-run).
- Observability-kernel 200-260 km/s control pericentre corrected to
  2.25 kpc (was 2.20; now matches the table and shipped summary).
- Gold sigma_RV<2 subset median Vgrf corrected to 16.4 km/s vs 15.8 km/s
  (full Gold) (were 16.57 / 15.72).
- Tier D defined by P<=0.50 rather than by a "remaining" subtraction.
- Reference author lists corrected: Belokurov et al. (2018; Deason, not
  Sherwin) and Sormani et al. (2022; Gerhard, Portail, Vasiliev, Clarke).

**Packaging / metadata**:
- requirements-lock.txt: removed private editable installs; documented the
  conda/compiled scientific dependencies (agama, gaiaunlimited, healpy,
  hdbscan, dustmaps, diptest).
- CITATION.cff: added ORCID and e-mail; version bumped to 1.0.6-review.
- NSS non-single-star flag merged: `nss_two_body` + `nss_solution_type` added to
  the master/tier/orbit catalogues, the four MRT tables, COLUMNS.md, and
  `make_mrt_tables.py` (9/541 Tier A+B, 23/1,952 Tier A+B+C; SB1=16/Orbital=6/
  AstroSpectroSB1=1), so the manuscript's "marked in the released catalogue"
  statement is backed by data.
- Potential-sensitivity sweep extended (`phase14/expanded_potential_sensitivity_*`,
  now 18 variants + `frac_ecc_gt_0p95`): added Omega_p=24/28 (N_bridge 8/7); a
  genuine bar-angle sweep that corrects the tab_potential bar-angle N_bridge
  range from 17-22 to 18-28; and a genuine halo-flattening q_z=0.80-1.05 sweep
  (median R_peri 114-118 pc), replacing the prior placeholder row whose numbers
  mirrored the halo-mass bracket. tab_potential and appendix prose updated.

**Final v1.0.6 packaging**:
- The `nss_two_body` cross-match flag is present in the catalogues and
  COLUMNS.md, so the "marked in the released catalogue" NSS statement is backed
  by data.
- The extended potential-sensitivity helper has been folded into
  `scripts/phase14v_expanded_potential_sensitivity.py`; the shipped script now
  regenerates all 18 variants, including the Hunter+2024 halo-flattening
  q_z sweep.
- The bundle directory and archive are packaged as
  `gaia_slow_vgrf_catalogue_v1.0.6_review`, with regenerated SHA-256 checksums
  and review zip.
- `main.pdf` has been recompiled from the corrected `main.tex`.

## v1.0.5-review - 2026-06-07

Response to two referee reports -- the external deep-review (15 numbered
issues) and the data-bundle review (8 issues). The manuscript builds cleanly;
`make validate-release` passes (smoke + provenance guard + SHA-256
checksums). Primary catalogue counts are unchanged
(Tier A / A+B / A+B+C = 289 / 541 / 1,952).

**Analysis additions** (new shipped scripts, run against the released products):
- Matched-control covariate balance + entropy-balancing weighted inference
  (`phase14ac`); the monotonic pericentre trend survives unweighted, kernel,
  propensity-IPW, and entropy-balanced schemes.
- Stronger smooth-DF nulls (`phase14ad`): the slow-tail excess is
  null-dependent (Gaussian mixture 7.4x, KDE 11x, single ellipsoid 25.6x)
  with bootstrap intervals; the single ellipsoid is demoted to a heuristic.
- Distance-covariance stress test expanded from 200 to the full 3,342-source
  tier-sensitive band with a confusion matrix (`phase14ae`); primary Tier A+B
  is 97.4% stable under maximal parallax-distance coupling.
- Radial-velocity / unresolved-multiplicity audit (`phase14af`): Gaia DR3
  non-single-star cross-match (1.2-1.7% flagged) + per-tier RV-quality table.
- External radial-velocity audit (`phase14al`): sparse APOGEE DR17/GALAH DR3
  overlap shows median Gaia-minus-external offsets of -0.10 and +0.16 km/s
  with robust scatters of 1.10 and 1.72 km/s.
- Truth-labelled probability-score diagnostic (`phase14ak`): GeDR3mock bins
  confirm monotonic score behaviour while explicitly avoiding an empirical
  real-population calibration claim.
- Uncertainty ledger of primary numbers (`phase14ag`).
- McMillan (2017) comparison potential (`phase14ah`): reproduces the static
  Hunter+2024 orbit summary, confirming the radial/compact result is not
  specific to the Hunter model.
- Filion et al. (2025) external comparison (`phase14aj`): resolves the 69
  APOGEE/Gaia low-azimuthal-velocity stars to Gaia DR3, applies the present
  low-\vgrf{} tiering and RVS-quality audit, and reports that none enters
  Tier A+B+C.

**Provenance & reproducibility**: pre-submission guard
(`scripts/check_provenance.py`), `make validate-release`, pinned
`requirements-lock.txt`, parent-buffer manifest, bundled AGAMA potential `.ini`
files, and SHA-256 release checksums.

**Manuscript**: catalogue-first framing tightened; orbit claims split into
robust (high-e, compact pericentre) vs potential-sensitive (central reach,
barred long-apocentre inner-reach candidates); resonance overlay recast as qualitative; turning-point
framing in the introduction; citations updated (Wright & Binney 2026 and
Zhang et al. 2024 to published; added Hunt et al. 2016, Babusiaux et al.
2023, McMillan 2017, Palicio et al. 2023; removed an unused reference).
Exploratory robustness material (sensitivity batteries, Sgr A* approachers,
resonances, DBSCAN, distance-covariance confusion matrix, and the Filion et
al. external comparison) moved to an appendix; probability calibration,
external-RV, null-model, and compact-pericentre language made more explicitly
diagnostic/modest; appendix tables converted to floats for clean typesetting.

**Data-bundle review response** (second referee, 8 issues):
- Chemistry columns merged into every catalogue FITS and `COLUMNS.md`
  (`mh_gspphot` Gaia GSP-Phot [M/H]; `feh_spec`/`feh_spec_err`/`alpha_spec`
  spectroscopic abundances from the APOGEE DR17 / GALAH DR3 cross-match;
  `chem_survey`; `chem_population`), via
  `scripts/phase14u2_merge_chemistry_into_catalogues.py`. Tier A+B+C coverage
  is 1,279 GSP-Phot metallicities and a 117-star spectroscopic alpha subset
  (median [Fe/H] = -0.81). Photometric and spectroscopic metallicities are
  kept in separate columns. [Issue 2]
- Machine-readable tables converted to AAS/CDS byte-by-byte MRT format
  (astropy-generated by `scripts/make_mrt_tables.py`, fully regenerable). Tier
  MRTs now carry the chemistry columns, and a new
  `mrt/catalogue_orbits_tierABC_mrt.txt` adds orbit summaries + actions; sparse
  fields render as blank fixed-width columns with CDS `?` null markers.
- Full-sample convergent action audit
  (`scripts/phase14ai_full_action_audit.py`): orbit-time-averaged actions
  (`J_*_timeavg`) and a per-star `action_reliability_flag` now cover all
  1,952 Tier A+B+C stars instead of a 100-star control. The azimuthal action
  is exact and J_R reliable (median 8%), but J_z is unreliable for these
  plunging orbits (median 62%; 1,164 stars `sampled_poor`); actions are
  presented as caveated diagnostics. [Issue 1]
- Nominal membership purity stated in the abstract and tier definitions:
  Tier A 98.7%, Tier A+B 94.3%, Tier A+B+C 72.7%. [Issue 4]
- Barred-frame energy is not conserved (relative drift median 2%, 99th
  percentile 12%) but the Jacobi integral is conserved to ~5e-7; disclosed in
  the orbit-model section. [Issue 6]
- 25 km/s threshold sensitivity quantified: 1,420 / 2,755 / 4,594 stars at
  point-estimate Vgrf < 20 / 25 / 30 km/s. [minor]
- Absolute filesystem paths relativised to `<GAIA_ROOT>` in the generated
  catalogue/orbit summaries. [Issue 8]

**Zenodo**: version bumped to v1.0.5-review (`config.yml`, `README.md`,
`README-submission.md`, `CITATION.cff`); mint a new Zenodo version at upload.
The concept DOI `10.5281/zenodo.20116134` continues to resolve to the latest
version.

## v1.0.4-review - 2026-06-05

This review-stage release adopts the final `Vc=229.0 km/s` Galactocentric
frame convention throughout the bundle, replacing the earlier 232 km/s
working convention. Catalogue products, manuscript text, tables, smoke
checks, and release documentation have been aligned to the regenerated
v1.0.4 products.

### Numerical Refresh

| Quantity | 232 km/s working draft | v1.0.4 adopted 229 km/s |
|---|---:|---:|
| Tier A / A+B / A+B+C | 276 / 517 / 1,835 | 289 / 541 / 1,952 |
| Gold obs. subset / point-estimate `<25 km/s` | 397 / 2,591 | 420 / 2,755 |
| MC median `e`, `R_peri`, `R_apo` | 0.949 / 151 pc / 7.45 kpc | 0.949 / 154 pc / 7.35 kpc |
| Smooth-DF excess | 24.0x | 25.5x |
| Tier A+B low-parent-count / `V_eff` factor | 83.8% / 2.13x | 83.2% / 2.09x |
| Chemistry alpha subset / median `[Fe/H]` | 113 / -0.80 | 117 / -0.81 |

### Manuscript And Product Updates

- Updated the Galactocentric transformation text to use
  `Vcirc = 229.0 km/s`, consistent with the Eilers et al. (2019)
  rotation-curve value and the axisymmetrised Hunter et al. (2024)
  potential used for orbit integrations.
- Added the present-day radial-phase sanity check:
  `scripts/phase14y_radial_phase_mc.py`,
  `phase14/radial_phase/radial_phase_mc_summary.json`, and
  `tables/v15/tab_radial_phase.tex`.
- Added a private diagnostic zero-gravity/rectilinear baseline:
  `scripts/phase15a_rectilinear_baseline.py`,
  `phase14/rectilinear_baseline/rectilinear_baseline_summary.json`, and
  `tables/v15/tab_rectilinear_baseline.tex`. This diagnostic is shipped
  for inspection but is not inserted as a manuscript result.
- Regenerated the Tier A, Tier A+B, and Tier A+B+C MRT products from the
  final v1.0.4 catalogue counts.
- Regenerated stale figure/provenance products that still carried the earlier
  232 km/s counts, including `fig01`, `fig02`, `fig15`,
  `fig_ndf_expectation`, `expanded_context_figures_summary.json`, and the
  static-energy cache used by the E--`L_z` figure; removed superseded
  `wp2_veff.*` and `wp2_potential_sensitivity_*` duplicates in favour of the
  current `expanded_*` products.
- Pinned the public smoke-regression count checks to the v1.0.4 values.
- Tightened the terminal availability-note layout before the references
  so the Software Availability note no longer opens with excessive
  whitespace.
- Regenerated `main.pdf` and the descriptive manuscript copy
  `Humble_2026_Gaia_DR3_slow_vgrf_catalogue_manuscript.pdf`.

## v1.0.3-review - 2026-06-03

This layout-only update refreshes the manuscript PDF after the v1.0.2-review
archive. Catalogue files, numerical products, selection-function products, and
primary science quantities are unchanged from v1.0.2-review.

- Removed a forced page break before the single-column inner-Galaxy table block
  so the Bailer-Jones overlap paragraph is no longer stranded on an almost
  blank page.
- Kept the Sgr A* candidate table and its cumulative-probability figure together
  before the robustness section, preserving table/figure order while reducing
  visible whitespace.
- Regenerated `main.pdf` and the descriptive manuscript copy
  `Humble_2026_Gaia_DR3_slow_vgrf_catalogue_manuscript.pdf`.

## v1.0.2-review - 2026-06-01

This review-stage update packages the referee-panel provenance reconciliation
and reproduce-cycle audit after the v1.0.1 barred-orbit correction. In that
earlier 232 km/s working convention, primary membership counts and the
default catalogue definition were unchanged relative to v1.0.1.

### Reproducibility product addition

- Added per-source GaiaUnlimited DR3-RVS parent-count diagnostics to
  `phase14/expanded_selection_function.fits`: `sf_parent_count` and
  `sf_prior_dominated_n_lt10`. Cells absent from the GaiaUnlimited
  `dr3-rvs-nk.h5` grid are encoded as parent count zero, matching the
  `p=0.5` prior-mean fill used by `DR3RVSSelectionFunction.query(...,
  fill_nan=True)`.
- Regenerated `expanded_selection_function_summary.*` and
  `expanded_selection_function_lown_audit.csv` from those per-source columns.
  The regenerated low-parent-count fractions reproduce the manuscript table
  values to rounding; no manuscript science numbers changed.

### Reproduce cycle

Operator-directed provenance pass for the point-estimate orbit products.

- Regenerated the Tier~A+B+C point-estimate static and barred orbit catalogue
  with `trajsize=40001` and parabola-interpolated cylindrical pericentres at
  `dR=0` sign changes. In that earlier 232 km/s reproduce cycle, the tier
  counts remained unchanged
  (`N_A+B=517`, `N_A+B+C=1835`) and the long-apocentre inner-reach counts remain unchanged
  (static = 3, barred Hunter/Sormani default = 22).
- The 40,001-point static readout gives median `R_peri = 0.1123 kpc` and
  `R_peri < 100 pc = 824`. The barred default gives
  `R_peri < 100 pc = 1543`, while its median pericentre remains
  readout-density sensitive (`trajsize=20001` gives 20.3 pc; `trajsize=40001`
  gives 14.1 pc). For this reason, manuscript text and tables treat barred
  median pericentres as resolution-limited diagnostics and emphasize the more
  stable threshold counts and long-apocentre inner-reach counts.
- Regenerated the matched point-estimate sensitivity products at the same
  40,001-point/interpolated readout: static halo mass variants, barred
  Hunter/Sormani pattern speeds 24/28/33/37.5/41, bar-angle sweep,
  halo-flattening sweep, and the Portail+2017 M2M bar cross-check.
- Re-seeded the expanded Sgr A* approacher audit from the regenerated orbit
  catalogue. The strict negative result remains unchanged:
  `P(r_sph,min < 10 pc) > 0.5` count = 0. The softer
  `P(r_sph,min < 100 pc) > 0.5` diagnostic is now 6 for the 10 seeded
  candidate/potential rows. The Sgr A* Monte Carlo engine itself already uses
  `trajsize=40001` and quadratic interpolation of `r_sph^2(t)`.
- Audited the full orbit-Monte-Carlo layer (`phase14x`): it remains the
  separate static-potential, 5,000-realisation-per-star product with
  `trajsize=1001`, giving the primary posterior-median values
  `median R_peri = 151 pc`, `median R_apo = 7.45 kpc`, and
  `median e = 0.949`. No sampling-sensitive deep-tail MC statistic from this
  product is quoted as a primary claim, so the reproduce cycle documents this
  sampling rather than rerunning the full MC layer.

### Referee-panel Round 1 provenance reconciliation

Provenance reconciliation of the barred orbit summary (panel decisions D2R / 9c648515).

- Historical note: the barred values in `catalogues/expanded_orbit_summary.json`
  were stale relative
  to the corrected v1.0.1 barred product: the sidecar JSON reported the pre-correction
  barred figures (median R_peri = 0.0556 kpc; R_peri < 100 pc = 1305; long-apocentre reach candidates = 32),
  whereas the shipped `catalogues/catalogue_expanded_orbits_tierABC.fits`, the manuscript
  tables (tab_reach, tab_sensitivity), and the v1.0.1 corrections below all use the
  corrected barred default (median R_peri = 0.0676 kpc; R_peri < 100 pc = 1226;
  long-apocentre reach candidates = 22). The JSON has been regenerated from the shipped FITS so the sidecar now
  agrees with the catalogue and manuscript.
- Superseded trajectory-size provenance: the subsequent reproduce-cycle addendum
  above reconciles the point-estimate generators and regenerated products at
  `trajsize=40001` with interpolated pericentres.

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

- Default barred long-apocentre reach candidates: 22.
- Bar-pattern-speed long-apocentre reach sensitivity: 6--24 (across $\Omega_p \in \{24, 28, 33, 37.5, 41\}$ km/s/kpc).
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
