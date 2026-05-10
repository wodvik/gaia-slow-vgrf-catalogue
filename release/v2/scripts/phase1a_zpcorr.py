"""Phase 1A — Lindegren+2021 zero-point correction per star.

Reads release/data/slow_stars_enriched_orbits.csv, applies zpt.get_zpt
to every source, and writes catalogue_zpcorr.fits + a histogram of the
correction values to release/v2/phase1/.

L21 validity window (from config.yml zpt_validity):
  6 < G < 21
  1.1 < nu_eff < 1.9 (5p)  OR  1.24 < pseudocolour < 1.72 (6p)
  astrometric_params_solved in {31, 95}

Stars outside the validity window get zpcorr_valid=False and zpcorr_value
left as the value zpt returns (with _warnings=False to silence per-star
spam) -- i.e. we record what L21 says, but flag it as extrapolation.
We do NOT extrapolate them into the science sample without the flag.

Outputs:
  release/v2/phase1/catalogue_zpcorr.fits
  release/v2/phase1/gate1A_zpcorr_hist.png
  release/v2/phase1/gate1A_zpcorr.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from astropy.table import Table

# --- repo paths ---
REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
SRC_CSV = REPO / CONFIG["input"]["source_csv"]
OUT_DIR = REPO / "release/v2/phase1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ZPV = CONFIG["zpt_validity"]


def main() -> int:
    print(f"[1A] reading {SRC_CSV.relative_to(REPO)}")
    df = pd.read_csv(SRC_CSV)
    n = len(df)
    print(f"[1A] {n} rows")
    assert n == CONFIG["input"]["n_stars_expected"], (
        f"row count {n} != expected {CONFIG['input']['n_stars_expected']}"
    )

    g = df["phot_g_mean_mag"].to_numpy(dtype=float)
    nu = df["nu_eff_used_in_astrometry"].to_numpy(dtype=float)
    psc = df["pseudocolour"].to_numpy(dtype=float)
    ecl = df["ecl_lat"].to_numpy(dtype=float)
    asp = df["astrometric_params_solved"].to_numpy(dtype=int)
    plx = df["parallax"].to_numpy(dtype=float)
    plx_err = df["parallax_error"].to_numpy(dtype=float)

    # zpt expects nu_eff for 5p (asp=31), pseudocolour for 6p (asp=95).
    # The function accepts NaN for the unused channel; we leave the source
    # column values in place so it can branch internally.
    from zero_point import zpt
    zpt.load_tables()

    # Validity flag (per L21 paper §2 + zpt docs).
    is_5p = asp == 31
    is_6p = asp == 95
    asp_ok = is_5p | is_6p
    g_ok = (g > ZPV["G_min"]) & (g < ZPV["G_max"])
    nu_ok = (nu > ZPV["nu_eff_min"]) & (nu < ZPV["nu_eff_max"])
    ps_ok = (psc > ZPV["pseudocolour_min"]) & (psc < ZPV["pseudocolour_max"])
    colour_ok = np.where(is_5p, nu_ok, np.where(is_6p, ps_ok, False))
    valid = asp_ok & g_ok & colour_ok & np.isfinite(plx)

    # Compute correction in mas for every row that has a valid asp
    # (zpt requires asp in {31,95}); for asp not in that set, set 0.
    zpcorr_mas = np.full(n, np.nan, dtype=float)
    if asp_ok.any():
        # zpt complains on NaN; clip for the call but record validity flag.
        g_in = np.where(np.isfinite(g), g, 12.0)
        nu_in = np.where(np.isfinite(nu), nu, 1.5)
        ps_in = np.where(np.isfinite(psc), psc, 1.5)
        ecl_in = np.where(np.isfinite(ecl), ecl, 0.0)
        asp_in = np.where(asp_ok, asp, 31)
        # _warnings=False to silence per-star validity warnings.
        zpcorr_mas = zpt.get_zpt(g_in, nu_in, ps_in, ecl_in, asp_in,
                                 _warnings=False)
        zpcorr_mas = np.where(asp_ok, zpcorr_mas, np.nan)

    # L21 returns the bias *to subtract*: parallax_zpcorr = parallax - zp.
    parallax_zpcorr = plx - zpcorr_mas
    # Where the correction is invalid, fall back to raw parallax (and flag).
    parallax_zpcorr = np.where(valid, parallax_zpcorr, plx)
    parallax_zpcorr_error = plx_err.copy()  # L21 correction is ~deterministic

    out = Table()
    out["source_id"] = df["source_id"].to_numpy()
    out["parallax"] = plx
    out["parallax_error"] = plx_err
    out["zpcorr_value_uas"] = zpcorr_mas * 1000.0
    out["zpcorr_valid"] = valid
    out["parallax_zpcorr"] = parallax_zpcorr
    out["parallax_zpcorr_error"] = parallax_zpcorr_error
    out["astrometric_params_solved"] = asp
    out["phot_g_mean_mag"] = g
    out["nu_eff_used_in_astrometry"] = nu
    out["pseudocolour"] = psc
    out["ecl_lat"] = ecl

    out_path = OUT_DIR / "catalogue_zpcorr.fits"
    out.write(out_path, format="fits", overwrite=True)
    print(f"[1A] wrote {out_path.relative_to(REPO)} ({len(out)} rows)")

    # ---- Gate 1A diagnostics ----
    zp_valid = out["zpcorr_value_uas"][valid]
    pct = np.percentile(zp_valid, [16, 50, 84])
    frac_invalid = float((~valid).sum()) / n
    summary = {
        "n_total": int(n),
        "n_valid": int(valid.sum()),
        "frac_zpcorr_invalid": frac_invalid,
        "zpcorr_uas_p16": float(pct[0]),
        "zpcorr_uas_p50": float(pct[1]),
        "zpcorr_uas_p84": float(pct[2]),
        "n_5p": int(is_5p.sum()),
        "n_6p": int(is_6p.sum()),
        "n_other_asp": int((~asp_ok).sum()),
        "expected_p50_uas_range": [-35, -25],
        "anomalous": not (-50 < pct[1] < -10),
    }
    (OUT_DIR / "gate1A_zpcorr.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(zp_valid, bins=80, color="steelblue", edgecolor="white")
    for q, label in zip(pct, ["p16", "median", "p84"]):
        ax.axvline(q, color="black", lw=1, ls="--",
                   label=f"{label} = {q:.1f} μas")
    ax.set_xlabel("L21 parallax zero-point (μas)")
    ax.set_ylabel("count")
    ax.set_title(f"Phase 1A — L21 ZP for {valid.sum()}/{n} valid stars")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "gate1A_zpcorr_hist.png", dpi=140)
    plt.close(fig)
    print(f"[1A] wrote gate1A histogram")

    if summary["anomalous"]:
        print("[1A] ANOMALY: median ZP outside [-50, -10] μas; STOP.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
