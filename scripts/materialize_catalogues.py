"""
Materialize the expanded MC-tiered candidate table into public catalogue products.

This is the bridge from the parent-buffer rebuild back into the paper/release
tree. It writes a compact FITS master catalogue plus Tier A and Tier A+B subset
FITS files.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table


REPO = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = Path("D:/GAIA/parent_scan/expanded_candidates_mc_tiered.csv")
DEFAULT_OUT = REPO / "catalogues"

CORE_COLUMNS = [
    "source_id",
    "ra",
    "dec",
    "l",
    "b",
    "parallax",
    "parallax_error",
    "parallax_over_error",
    "zpcorr_value_uas",
    "zpcorr_valid",
    "parallax_zpcorr",
    "dist_pc_final_screen",
    "dist_lo_pc_final_screen",
    "dist_hi_pc_final_screen",
    "dist_source_final_screen",
    "pmra",
    "pmra_error",
    "pmdec",
    "pmdec_error",
    "radial_velocity",
    "radial_velocity_error",
    "rv_quality",
    "rvs_quality_ok",
    "rv_chisq_pvalue",
    "rv_expected_sig_to_noise",
    "rv_nb_transits",
    "rv_amplitude_robust",
    "rv_template_teff",
    "phot_g_mean_mag",
    "bp_rp",
    "grvs_mag",
    "ruwe",
    "legacy_v_total_grf",
    "vgrf_zpcorr_inv",
    "vgrf_bj_or_inv",
    "vgrf_default_exact",
    "P_vgrf_below_25",
    "mc_realisations",
    "tier",
    "source_in_old_preselection",
    "parent_scan_file",
]


def classify_rvs_quality(df: pd.DataFrame) -> np.ndarray:
    """Apply the Gaia DR3 RVS-quality classes used by the catalogue build."""
    n = len(df)
    quality = np.full(n, "ok", dtype=object)

    pval = pd.to_numeric(df["rv_chisq_pvalue"], errors="coerce").to_numpy(dtype=float)
    snr_rv = pd.to_numeric(df["rv_expected_sig_to_noise"], errors="coerce").to_numpy(dtype=float)
    amp = pd.to_numeric(df["rv_amplitude_robust"], errors="coerce").to_numpy(dtype=float)
    n_trans = pd.to_numeric(df["rv_nb_transits"], errors="coerce").to_numpy(dtype=float)
    teff = pd.to_numeric(df["rv_template_teff"], errors="coerce").to_numpy(dtype=float)
    grvs = pd.to_numeric(df["grvs_mag"], errors="coerce").to_numpy(dtype=float)

    poor = (
        (np.isfinite(pval) & (pval < 0.01))
        | (np.isfinite(snr_rv) & (snr_rv < 2.0))
    )
    marginal = (
        (np.isfinite(amp) & np.isfinite(n_trans) & (amp > 10.0) & (n_trans < 10))
        | (
            np.isfinite(teff)
            & np.isfinite(grvs)
            & ((teff < 3500.0) | (teff > 7000.0))
            & (grvs > 12.0)
        )
    )
    quality[marginal] = "marginal"
    quality[poor] = "poor"
    return quality


def to_table(df: pd.DataFrame) -> Table:
    cols = [c for c in CORE_COLUMNS if c in df.columns]
    out = df[cols].copy()
    out = out.rename(
        columns={
            "dist_pc_final_screen": "dist_pc",
            "dist_lo_pc_final_screen": "dist_lo_pc",
            "dist_hi_pc_final_screen": "dist_hi_pc",
            "dist_source_final_screen": "dist_source",
            "vgrf_default_exact": "vgrf_default",
        }
    )
    return Table.from_pandas(out)


def summary(df: pd.DataFrame) -> dict:
    old = df["source_in_old_preselection"].astype(bool)
    tier = df["tier"].astype(str)
    tiers = ["A", "B", "C", "D", "X"]
    return {
        "n_processed": int(len(df)),
        "n_unique_source_id": int(df["source_id"].nunique()),
        "n_old_preselection": int(old.sum()),
        "n_outside_old": int((~old).sum()),
        "tier_counts": {t: int((tier == t).sum()) for t in tiers},
        "tier_counts_old_preselection": {t: int(((tier == t) & old).sum()) for t in tiers},
        "tier_counts_outside_old": {t: int(((tier == t) & ~old).sum()) for t in tiers},
        "headline_tier_A_plus_B": int(tier.isin(["A", "B"]).sum()),
        "orbit_summary_tier_A_plus_B_plus_C": int(tier.isin(["A", "B", "C"]).sum()),
        "point_estimate_vgrf_lt25": int((df["vgrf_default_exact"] < 25.0).sum()),
        "point_estimate_vgrf_lt25_outside_old": int(((df["vgrf_default_exact"] < 25.0) & ~old).sum()),
        "legacy_v_total_grf_lt25": int((df["legacy_v_total_grf"] < 25.0).sum()),
        "mc_realisations_counts": {
            str(int(k)): int(v) for k, v in df["mc_realisations"].value_counts().sort_index().items()
        },
        "rv_quality_counts": {
            str(k): int(v) for k, v in df["rv_quality"].value_counts().sort_index().items()
        },
        "P_vgrf_below_25_range": [
            float(df["P_vgrf_below_25"].min()),
            float(df["P_vgrf_below_25"].max()),
        ],
    }


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input_csv)
    df["source_id"] = df["source_id"].astype("int64")
    df["tier"] = df["tier"].astype(str)
    df["rv_quality"] = classify_rvs_quality(df)
    df["rvs_quality_ok"] = df["rv_quality"] == "ok"
    df = df.sort_values(["tier", "P_vgrf_below_25", "vgrf_default_exact", "source_id"],
                        ascending=[True, False, True, True]).reset_index(drop=True)

    master = to_table(df)
    master_path = out_dir / "catalogue_expanded_master.fits"
    master.write(master_path, overwrite=True)

    tier_a = to_table(df[df["tier"] == "A"].copy())
    tier_a_path = out_dir / "catalogue_expanded_tierA.fits"
    tier_a.write(tier_a_path, overwrite=True)

    tier_ab = to_table(df[df["tier"].isin(["A", "B"])].copy())
    tier_ab_path = out_dir / "catalogue_expanded_tierAB.fits"
    tier_ab.write(tier_ab_path, overwrite=True)

    tier_abc = to_table(df[df["tier"].isin(["A", "B", "C"])].copy())
    tier_abc_path = out_dir / "catalogue_expanded_tierABC.fits"
    tier_abc.write(tier_abc_path, overwrite=True)

    s = summary(df)
    s.update(
        {
            "input_csv": str(args.input_csv),
            "master_fits": str(master_path),
            "tierA_fits": str(tier_a_path),
            "tierAB_fits": str(tier_ab_path),
            "tierABC_fits": str(tier_abc_path),
        }
    )
    (out_dir / "expanded_catalogue_summary.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))
    return s


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", default=str(DEFAULT_INPUT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
