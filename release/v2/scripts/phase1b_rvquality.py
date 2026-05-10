"""Phase 1B' — RVS quality vetting (Katz+2023 § 5.4 thresholds).

The slow_stars catalogue already carries the rv_* columns from the Gaia
DR3 main table; this step just classifies each star into rv_quality
∈ {ok, marginal, poor} and writes a small companion catalogue.

Katz+2023 (DR3 RVS paper) flags:
  - rv_chisq_pvalue < 0.01            -> poor (constant-RV chi2 model fit)
  - rv_amplitude_robust > 10 km/s and rv_nb_transits < 10 -> marginal
  - rv_expected_sig_to_noise < 2     -> poor
  - rv_template_teff outside 3500–7000 K and grvs_mag > 12 -> marginal
Otherwise: ok.

Outputs
-------
release/v2/phase1/catalogue_rvquality.fits
release/v2/phase1/gate1Bp_rvquality.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT_DIR = REPO / "release/v2/phase1"
SRC_CSV = REPO / CONFIG["input"]["source_csv"]


def classify(df: pd.DataFrame) -> np.ndarray:
    n = len(df)
    quality = np.full(n, "ok", dtype=object)

    pval = df["rv_chisq_pvalue"].to_numpy(dtype=float)
    snr_rv = df["rv_expected_sig_to_noise"].to_numpy(dtype=float)
    amp = df["rv_amplitude_robust"].to_numpy(dtype=float)
    n_trans = df["rv_nb_transits"].to_numpy(dtype=float)
    teff = df["rv_template_teff"].to_numpy(dtype=float)
    grvs = df["grvs_mag"].to_numpy(dtype=float)

    poor = (
        (np.isfinite(pval) & (pval < 0.01))
        | (np.isfinite(snr_rv) & (snr_rv < 2.0))
    )
    marginal = (
        (np.isfinite(amp) & np.isfinite(n_trans)
         & (amp > 10.0) & (n_trans < 10))
        | (np.isfinite(teff) & np.isfinite(grvs)
           & ((teff < 3500) | (teff > 7000)) & (grvs > 12.0))
    )
    quality[marginal] = "marginal"
    quality[poor] = "poor"
    return quality


def main() -> int:
    df = pd.read_csv(SRC_CSV)
    print(f"[1B'] {len(df)} rows")
    quality = classify(df)

    out = Table()
    out["source_id"] = df["source_id"].to_numpy()
    out["rv_quality"] = quality.astype(str)
    out["rv_chisq_pvalue"] = df["rv_chisq_pvalue"].to_numpy()
    out["rv_expected_sig_to_noise"] = df["rv_expected_sig_to_noise"].to_numpy()
    out["rv_amplitude_robust"] = df["rv_amplitude_robust"].to_numpy()
    out["rv_nb_transits"] = df["rv_nb_transits"].to_numpy()
    out["rv_template_teff"] = df["rv_template_teff"].to_numpy()
    out["grvs_mag"] = df["grvs_mag"].to_numpy()
    out["radial_velocity"] = df["radial_velocity"].to_numpy()
    out["radial_velocity_error"] = df["radial_velocity_error"].to_numpy()

    out_path = OUT_DIR / "catalogue_rvquality.fits"
    out.write(out_path, format="fits", overwrite=True)
    print(f"[1B'] wrote {out_path.relative_to(REPO)}")

    # Sanity: median Vgrf of full sample vs poor-excluded.
    vgrf_all = df["V_grf"].to_numpy(dtype=float)
    poor_mask = quality == "poor"
    vgrf_clean = vgrf_all[~poor_mask]
    summary = {
        "n_total": int(len(df)),
        "n_ok": int((quality == "ok").sum()),
        "n_marginal": int((quality == "marginal").sum()),
        "n_poor": int(poor_mask.sum()),
        "median_vgrf_all": float(np.nanmedian(vgrf_all)),
        "median_vgrf_no_poor": float(np.nanmedian(vgrf_clean)),
    }
    (OUT_DIR / "gate1Bp_rvquality.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
