# Column Dictionary

This file summarizes the principal machine-readable products in this repository. FITS files are readable with `astropy.table.Table.read`.

## `release/v2/phase1/catalogue_tierAB.fits`

Headline Tier A+B catalogue. Rows: 334. Columns: 26.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `ra` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `dec` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax_error` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `zpcorr_uas` | microarcsec | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `zpcorr_valid` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax_zpcorr` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `dist_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_lo_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_hi_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_source` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmra` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmdec` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmra_error` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmdec_error` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `radial_velocity` | km/s | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `radial_velocity_error` | km/s | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `rv_quality` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `vgrf_default` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_grav22` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_lsr6` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_rb20` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `P_vgrf_below_25` | km/s | Monte Carlo probability that Vgrf is below 25 km/s. |
| `mc_realisations` | count or identifier | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |

## `release/v2/phase1/catalogue_tierA.fits`

Tier A high-confidence subset. Rows: 214. Columns: 26.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `ra` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `dec` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax_error` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `zpcorr_uas` | microarcsec | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `zpcorr_valid` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax_zpcorr` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `dist_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_lo_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_hi_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_source` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmra` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmdec` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmra_error` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmdec_error` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `radial_velocity` | km/s | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `radial_velocity_error` | km/s | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `rv_quality` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `vgrf_default` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_grav22` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_lsr6` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_rb20` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `P_vgrf_below_25` | km/s | Monte Carlo probability that Vgrf is below 25 km/s. |
| `mc_realisations` | count or identifier | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |

## `release/v2/phase1/catalogue_v2.fits`

Master candidate catalogue. Rows: 2859. Columns: 26.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `ra` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `dec` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax_error` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `zpcorr_uas` | microarcsec | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `zpcorr_valid` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `parallax_zpcorr` | mas | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `dist_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_lo_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_hi_pc` | pc | Adopted Bailer-Jones photogeometric distance posterior summary. |
| `dist_source` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmra` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmdec` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmra_error` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `pmdec_error` | mas/yr | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `radial_velocity` | km/s | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `radial_velocity_error` | km/s | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `rv_quality` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `vgrf_default` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_grav22` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_lsr6` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `vgrf_rb20` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `P_vgrf_below_25` | km/s | Monte Carlo probability that Vgrf is below 25 km/s. |
| `mc_realisations` | count or identifier | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |

## `release/v2/phase2/catalogue_v2_sf.fits`

Selection-function table. Rows: 2859. Columns: 8.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `l` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `b` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `grvs_mag` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `sf_value` | dimensionless | Gaia DR3 RVS selection-function value. |
| `sf_weight` | dimensionless | Inverse selection-function weight. |
| `sf_invalid` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |

## `release/v2/phase4/catalogue_v2_orbits.fits`

Point-estimate orbit catalogue. Rows: 2859. Columns: 44.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |
| `P_vgrf_below_25` | km/s | Monte Carlo probability that Vgrf is below 25 km/s. |
| `vgrf_default` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `rv_quality` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `x_kpc` | kpc | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `y_kpc` | kpc | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `z_kpc` | kpc | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `vx_kms` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `vy_kms` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `vz_kms` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `static_R_peri_kpc` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `static_R_apo_kpc` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `static_z_max_kpc` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `static_min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `static_ecc` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `static_n_peri` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `static_E_mean` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `static_E_drift_rel` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `static_E_range_rel` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `static_EJ_drift_rel` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `static_EJ_range_rel` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `static_Lz_kpc_kms` | kpc km/s | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `static_star_idx` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `barred_R_peri_kpc` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `barred_R_apo_kpc` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `barred_z_max_kpc` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `barred_min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `barred_ecc` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `barred_n_peri` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `barred_E_mean` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `barred_E_drift_rel` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `barred_E_range_rel` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `barred_EJ_drift_rel` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `barred_EJ_range_rel` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `barred_Lz_kpc_kms` | kpc km/s | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `barred_star_idx` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `J_R` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_z` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_phi` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `Omega_R` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `Omega_z` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `Omega_phi` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `res_ratio_OmegaR_over_dPhi` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |

## `release/v2/phase5/control_orbits.fits`

Matched-control orbit catalogue. Rows: 10078. Columns: 20.

| Column | Unit | Description |
|---|---|---|
| `band` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `ra` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `dec` | deg | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `G` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `dist_pc_v1` | pc | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `b_deg` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `vgrf_v1` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `R_peri_kpc` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `R_apo_kpc` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `z_max_kpc` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `ecc` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `J_R` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_z` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_phi` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `Omega_R` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `Omega_phi` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `res_ratio` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `kw_weight` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |

## `release/v2/phase5/chemistry_v2.fits`

APOGEE/GALAH chemistry matches. Rows: 141. Columns: 14.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |
| `P_vgrf_below_25` | km/s | Monte Carlo probability that Vgrf is below 25 km/s. |
| `vgrf_default` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `rv_quality` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `Teff` | K | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `logg` | dex | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `FeH` | dex | Metallicity estimate or selected metallicity proxy. |
| `e_FeH` | dex | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `MgFe` | dex | Abundance ratio or alpha-element proxy from external spectroscopy. |
| `AlFe` | dex | Abundance ratio or alpha-element proxy from external spectroscopy. |
| `MnFe` | dex | Abundance ratio or alpha-element proxy from external spectroscopy. |
| `alpha_proxy` | dex | Abundance ratio or alpha-element proxy from external spectroscopy. |
| `survey` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |

## `release/v2/phase5/populations_v2.fits`

Chemistry and population labels. Rows: 2865. Columns: 13.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |
| `P_vgrf_below_25` | km/s | Monte Carlo probability that Vgrf is below 25 km/s. |
| `vgrf_default` | km/s | Galactic rest-frame speed under the named solar-parameter convention. |
| `MH_xp` | dex | Metallicity estimate or selected metallicity proxy. |
| `FeH` | dex | Metallicity estimate or selected metallicity proxy. |
| `MgFe` | dex | Abundance ratio or alpha-element proxy from external spectroscopy. |
| `alpha_proxy` | dex | Abundance ratio or alpha-element proxy from external spectroscopy. |
| `survey` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `best_FeH` | dex | Metallicity estimate or selected metallicity proxy. |
| `best_alpha` | dex | Abundance ratio or alpha-element proxy from external spectroscopy. |
| `chem_source` | category | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `population` | category | Literature-region label used for subset-qualified chemodynamic context. |

## `release/v2/phase6/catalogue_mc_orbits.fits`

Monte Carlo orbit posterior summaries. Rows: 632. Columns: 27.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |
| `n_realisations` | count or identifier | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `R_peri_kpc_p16` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `R_peri_kpc_p50` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `R_peri_kpc_p84` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `R_apo_kpc_p16` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `R_apo_kpc_p50` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `R_apo_kpc_p84` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `z_max_kpc_p16` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `z_max_kpc_p50` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `z_max_kpc_p84` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `ecc_p16` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `ecc_p50` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `ecc_p84` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `min_r_sph_kpc_p16` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `min_r_sph_kpc_p50` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `min_r_sph_kpc_p84` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `J_R_p16` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_R_p50` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_R_p84` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_z_p16` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_z_p50` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_z_p84` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_phi_p16` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_phi_p50` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `J_phi_p84` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |

## `release/v2/phase6/bar_speed_subsample.fits`

Bar-pattern-speed sensitivity subsample. Rows: 999. Columns: 10.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |
| `Omega_p` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `speed_label` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `R_peri_kpc` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `R_apo_kpc` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `z_max_kpc` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `ecc` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `EJ_drift_rel` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |

## `release/v2/phase6/orbits_barred_dt0p1.fits`

Fine-sampled barred orbit diagnostics. Rows: 2859. Columns: 15.

| Column | Unit | Description |
|---|---|---|
| `star_idx` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |
| `fine_R_peri_kpc` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `fine_R_apo_kpc` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `fine_z_max_kpc` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `fine_min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `fine_ecc` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `fine_EJ_drift_rel` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |
| `fine_EJ_range_rel` | kpc km/s | Action or angular-momentum-like quantity from agama orbit/action calculation. |
| `old_R_peri_kpc` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `old_R_apo_kpc` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `old_z_max_kpc` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `old_ecc` | dimensionless | Orbital eccentricity in the adopted Galactic potential. |
| `old_EJ_drift` | see description | Column carried by the pipeline product; see manuscript methods and scripts for derivation. |

## `release/v2/phase6/orbits_static_dt0p1.fits`

Fine-sampled static orbit diagnostics. Rows: 632. Columns: 8.

| Column | Unit | Description |
|---|---|---|
| `source_id` | count or identifier | Gaia DR3 source identifier. |
| `tier` | category | Velocity-threshold membership tier: A, B, C, D, or X as defined in the manuscript. |
| `fine_R_peri_kpc` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `fine_R_apo_kpc` | kpc | Orbital apocentre in the adopted Galactic potential. |
| `fine_z_max_kpc` | kpc | Maximum vertical excursion in the adopted Galactic potential. |
| `fine_min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius reached in the integration. |
| `old_R_peri_kpc` | kpc | Orbital pericentre in the adopted Galactic potential. |
| `old_min_r_sph_kpc` | kpc | Minimum spherical Galactocentric radius reached in the integration. |

