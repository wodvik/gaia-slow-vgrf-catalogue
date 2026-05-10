-- Gaia DR3 6D parent-source selection used by the slow-Vgrf pipeline.
--
-- The local pipeline was run from a local mirror of gaia_source.  This
-- query is the equivalent Gaia Archive selection for reproducing the
-- parent source table before downstream zero-point correction,
-- distance joining, Vgrf calculation, and probabilistic tiering.

SELECT
  source_id,
  ra,
  dec,
  parallax,
  parallax_error,
  pmra,
  pmra_error,
  pmdec,
  pmdec_error,
  ra_dec_corr,
  ra_parallax_corr,
  ra_pmra_corr,
  ra_pmdec_corr,
  dec_parallax_corr,
  dec_pmra_corr,
  dec_pmdec_corr,
  parallax_pmra_corr,
  parallax_pmdec_corr,
  pmra_pmdec_corr,
  radial_velocity,
  radial_velocity_error,
  rv_method_used,
  rv_expected_sig_to_noise,
  rv_renormalised_gof,
  rv_chisq_pvalue,
  grvs_mag,
  phot_g_mean_mag,
  phot_bp_mean_mag,
  phot_rp_mean_mag,
  bp_rp,
  ruwe,
  duplicated_source,
  l,
  b
FROM gaiadr3.gaia_source
WHERE radial_velocity IS NOT NULL
  AND parallax > 0
  AND parallax_over_error > 5
  AND pmra IS NOT NULL
  AND pmdec IS NOT NULL;
