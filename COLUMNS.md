# Catalogue Column Dictionary

This file defines the columns used in the first public review-stage FITS and
CSV catalogue products. Gaia `source_id` values are stored as integer
identifiers and should not be converted to floating point.

## Master And Tier Subset Catalogues

These columns appear in:

- `catalogues/catalogue_expanded_master.fits`
- `catalogues/catalogue_expanded_tierA.fits`
- `catalogues/catalogue_expanded_tierAB.fits`
- `catalogues/catalogue_expanded_tierABC.fits`

| Column | Unit | Description |
|---|---:|---|
| `source_id` | none | Gaia DR3 source identifier. |
| `ra` | deg | ICRS right ascension from Gaia DR3. |
| `dec` | deg | ICRS declination from Gaia DR3. |
| `l` | deg | Galactic longitude from Gaia DR3. |
| `b` | deg | Galactic latitude from Gaia DR3. |
| `parallax` | mas | Gaia DR3 parallax before zero-point correction. |
| `parallax_error` | mas | Gaia DR3 parallax uncertainty. |
| `parallax_over_error` | none | Gaia DR3 parallax signal-to-noise ratio. |
| `zpcorr_value_uas` | microarcsec | Lindegren et al. (2021) parallax zero-point correction applied to the source. |
| `zpcorr_valid` | boolean | Whether the source lies inside the configured Lindegren et al. (2021) validity window. |
| `parallax_zpcorr` | mas | Parallax after applying the zero-point correction. |
| `dist_pc` | pc | Adopted distance, usually Bailer-Jones et al. (2021) photogeometric distance. |
| `dist_lo_pc` | pc | Lower credible-distance bound used for asymmetric distance draws. |
| `dist_hi_pc` | pc | Upper credible-distance bound used for asymmetric distance draws. |
| `dist_source` | none | Distance provenance label. |
| `pmra` | mas/yr | Gaia DR3 proper motion in right ascension, including cos(dec). |
| `pmra_error` | mas/yr | Uncertainty on `pmra`. |
| `pmdec` | mas/yr | Gaia DR3 proper motion in declination. |
| `pmdec_error` | mas/yr | Uncertainty on `pmdec`. |
| `radial_velocity` | km/s | Gaia DR3 line-of-sight radial velocity. |
| `radial_velocity_error` | km/s | Gaia DR3 radial-velocity uncertainty. |
| `rv_quality` | none | Internal Gaia DR3 RVS-quality class: `ok`, `marginal`, or `poor`, using the same thresholds as the manuscript quality checks. |
| `rvs_quality_ok` | boolean | True when `rv_quality == "ok"`; this is the RVS-quality gate used in the Gold subset. |
| `rv_chisq_pvalue` | none | Gaia DR3 radial-velocity chi-square p-value used in the RVS-quality diagnostic. |
| `rv_expected_sig_to_noise` | none | Gaia DR3 expected radial-velocity signal-to-noise ratio used in the RVS-quality diagnostic. |
| `rv_nb_transits` | none | Number of Gaia RVS transits used for the radial velocity. |
| `rv_amplitude_robust` | km/s | Gaia DR3 robust radial-velocity amplitude statistic. |
| `rv_template_teff` | K | Effective temperature of the Gaia RVS template. |
| `phot_g_mean_mag` | mag | Gaia DR3 mean G-band magnitude. |
| `bp_rp` | mag | Gaia DR3 BP-RP colour. |
| `grvs_mag` | mag | Gaia DR3 RVS magnitude. |
| `ruwe` | none | Gaia DR3 renormalised unit weight error. |
| `legacy_v_total_grf` | km/s | Preliminary Galactic-rest-frame speed used only for legacy preselection provenance. |
| `vgrf_zpcorr_inv` | km/s | Point-estimate Galactic-rest-frame speed using zero-point-corrected inverse-parallax distance. |
| `vgrf_bj_or_inv` | km/s | Point-estimate Galactic-rest-frame speed using Bailer-Jones distance where available, otherwise inverse parallax. |
| `vgrf_default` | km/s | Final adopted point-estimate Galactic-rest-frame speed. |
| `P_vgrf_below_25` | dimensionless | Monte Carlo probability that `Vgrf < 25 km/s`; this is a probability and has no velocity unit. |
| `mc_realisations` | none | Number of velocity-threshold Monte Carlo realisations used for this source. |
| `tier` | none | Probability tier: A (`P>0.95`), B (`0.84<P<=0.95`), C (`0.50<P<=0.84`), D (point-estimate below 25 km/s but `P<=0.50`), or X. |
| `source_in_old_preselection` | boolean | Whether the source was present in the preliminary 2,859-source development preselection. |
| `parent_scan_file` | none | Gaia mirror CSV file from which the source was recovered during the parent-buffer scan. |

## Point-Estimate Orbit Catalogue

These columns appear in:

- `catalogues/catalogue_expanded_orbits_tierABC.fits`
- `catalogues/catalogue_expanded_orbits_tierABC.csv`

This product contains one row per Tier A+B+C star. Orbit quantities are
point-estimate summaries in the adopted Hunter+2024 static potential and
the barred Hunter/Sormani sensitivity potential; they are model outputs,
not observables.

| Column | Unit | Description |
|---|---:|---|
| `source_id` | none | Gaia DR3 source identifier. |
| `tier` | none | Probability tier copied from the master catalogue. |
| `P_vgrf_below_25` | dimensionless | Monte Carlo threshold-membership probability. |
| `rv_quality` | none | Internal Gaia DR3 RVS-quality class copied from the master catalogue. |
| `rvs_quality_ok` | boolean | True when `rv_quality == "ok"`. |
| `source_in_old_preselection` | boolean | Whether the source was present in the preliminary 2,859-source development preselection. |
| `vgrf_default_exact` | km/s | Final adopted point-estimate Galactic-rest-frame speed from the catalogue pipeline. |
| `vgrf_default_orbit_ic` | km/s | Galactic-rest-frame speed recomputed from the orbit-integration initial condition. |
| `dist_pc_final_screen` | pc | Adopted distance used for the orbit initial condition. |
| `dist_source_final_screen` | none | Distance provenance label used for the orbit initial condition. |
| `x_kpc`, `y_kpc`, `z_kpc` | kpc | Galactocentric Cartesian position used as the orbit initial condition. |
| `vx_kms`, `vy_kms`, `vz_kms` | km/s | Galactocentric Cartesian velocity used as the orbit initial condition. |
| `static_R_peri_kpc` | kpc | Minimum cylindrical Galactocentric radius in the static potential integration. |
| `static_R_apo_kpc` | kpc | Maximum cylindrical Galactocentric radius in the static potential integration. |
| `static_z_max_kpc` | kpc | Maximum absolute Galactocentric height in the static potential integration. |
| `static_min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius in the static potential integration. |
| `static_ecc` | dimensionless | Eccentricity proxy `(R_apo - R_peri)/(R_apo + R_peri)` from the static integration. |
| `static_n_peri` | none | Number of sampled cylindrical pericentre passages in the static integration. |
| `static_E_drift_rel` | dimensionless | Relative endpoint energy drift in the static integration. |
| `static_E_range_rel` | dimensionless | Relative sampled energy range in the static integration. |
| `static_EJ_drift_rel` | dimensionless | Jacobi-energy drift diagnostic; retained for schema symmetry and zero for the non-rotating static case. |
| `static_EJ_range_rel` | dimensionless | Jacobi-energy range diagnostic; retained for schema symmetry and zero for the non-rotating static case. |
| `static_Lz_kpc_kms` | kpc km/s | Mean Galactocentric angular momentum about the z-axis in the static integration. |
| `static_star_idx` | none | Row index used internally by the static orbit integration. |
| `barred_R_peri_kpc` | kpc | Minimum cylindrical Galactocentric radius in the barred-potential sensitivity integration. |
| `barred_R_apo_kpc` | kpc | Maximum cylindrical Galactocentric radius in the barred-potential sensitivity integration. |
| `barred_z_max_kpc` | kpc | Maximum absolute Galactocentric height in the barred-potential sensitivity integration. |
| `barred_min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius in the barred-potential sensitivity integration. |
| `barred_ecc` | dimensionless | Eccentricity proxy from the barred-potential sensitivity integration. |
| `barred_n_peri` | none | Number of sampled cylindrical pericentre passages in the barred-potential sensitivity integration. |
| `barred_E_drift_rel` | dimensionless | Relative endpoint energy drift diagnostic in the rotating barred integration. |
| `barred_E_range_rel` | dimensionless | Relative sampled energy range diagnostic in the rotating barred integration. |
| `barred_EJ_drift_rel` | dimensionless | Relative endpoint Jacobi-energy drift in the rotating barred integration. |
| `barred_EJ_range_rel` | dimensionless | Relative sampled Jacobi-energy range in the rotating barred integration. |
| `barred_Lz_kpc_kms` | kpc km/s | Mean Galactocentric angular momentum about the z-axis in the barred integration. |
| `barred_star_idx` | none | Row index used internally by the barred orbit integration. |
| `J_R`, `J_z`, `J_phi` | kpc km/s | AGAMA Staeckel-fudge radial, vertical, and azimuthal actions in the static potential. |
| `Omega_R`, `Omega_z`, `Omega_phi` | km/s/kpc | AGAMA radial, vertical, and azimuthal orbital frequencies in the static potential. |
| `res_ratio_OmegaR_over_dPhi` | dimensionless | Resonance diagnostic `Omega_R/(Omega_phi - Omega_p)` evaluated relative to the default bar pattern speed. |

## Orbit Monte Carlo Catalogue

These columns appear in:

- `analysis_products/expanded_orbit_mc/expanded_catalogue_mc_orbits.fits`
- `analysis_products/expanded_orbit_mc/expanded_catalogue_mc_orbits.csv`

| Column | Unit | Description |
|---|---:|---|
| `source_id` | none | Gaia DR3 source identifier. |
| `tier` | none | Probability tier copied from the master catalogue. |
| `P_vgrf_below_25` | dimensionless | Monte Carlo threshold-membership probability. |
| `vgrf_default_exact` | km/s | Final adopted point-estimate Galactic-rest-frame speed used for the orbit-MC product. |
| `mc_realisations` | none | Number of velocity-threshold Monte Carlo realisations for the source. |
| `orbit_mc_realisations` | none | Number of orbit Monte Carlo realisations; 5,000 for every row in this product. |
| `R_peri_kpc_p16`, `R_peri_kpc_p50`, `R_peri_kpc_p84` | kpc | 16th, 50th, and 84th percentiles of cylindrical pericentre in the adopted static Hunter+2024 potential. |
| `R_apo_kpc_p16`, `R_apo_kpc_p50`, `R_apo_kpc_p84` | kpc | 16th, 50th, and 84th percentiles of cylindrical apocentre. |
| `z_max_kpc_p16`, `z_max_kpc_p50`, `z_max_kpc_p84` | kpc | 16th, 50th, and 84th percentiles of maximum absolute Galactocentric height. |
| `ecc_p16`, `ecc_p50`, `ecc_p84` | dimensionless | 16th, 50th, and 84th percentiles of orbital eccentricity. |
| `min_r_sph_kpc_p16`, `min_r_sph_kpc_p50`, `min_r_sph_kpc_p84` | kpc | 16th, 50th, and 84th percentiles of minimum spherical Galactocentric radius. |
| `J_R_p16`, `J_R_p50`, `J_R_p84` | kpc km/s | 16th, 50th, and 84th percentiles of radial action from AGAMA. |
| `J_z_p16`, `J_z_p50`, `J_z_p84` | kpc km/s | 16th, 50th, and 84th percentiles of vertical action from AGAMA. |
| `J_phi_p16`, `J_phi_p50`, `J_phi_p84` | kpc km/s | 16th, 50th, and 84th percentiles of azimuthal action/angular momentum proxy from AGAMA; sign distinguishes prograde and retrograde motion in the adopted convention. |

## Selection-Function Catalogue

These columns appear in:

- `analysis_products/expanded_selection_function.fits`

| Column | Unit | Description |
|---|---:|---|
| `source_id` | none | Gaia DR3 source identifier. |
| `tier` | none | Probability tier copied from the master catalogue. |
| `l` | deg | Galactic longitude copied from the master catalogue for positional context. |
| `b` | deg | Galactic latitude copied from the master catalogue for positional context. |
| `grvs_mag` | mag | Gaia RVS magnitude supplied as the GaiaUnlimited `g` coordinate. |
| `bp_rp` | mag | Gaia BP-RP colour supplied as the GaiaUnlimited colour coordinate; missing colours are filled with the finite catalogue median for evaluation. |
| `sf_value` | dimensionless | Castro-Ginard et al. (2023) Gaia DR3 RVS selection-function value returned by GaiaUnlimited. |
| `sf_weight` | dimensionless | Contextual inverse-selection weight, computed as `1/max(sf_value, 0.02)` for finite values. |
| `sf_invalid` | boolean | True when GaiaUnlimited returned a non-finite selection-function value for the source. |
| `P_vgrf_below_25` | dimensionless | Monte Carlo probability that `Vgrf < 25 km/s`, copied from the master catalogue. |
