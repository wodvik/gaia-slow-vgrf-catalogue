"""Phase 1D — adaptive Monte Carlo tiering of P(Vgrf < 25 km/s).

For every star:
  pass 1 -- 500 realisations from
            (parallax, pmra, pmdec) ~ MVN(mu, Sigma_3x3)
            radial_velocity ~ N(rv, rv_err)
            distance ~ split-normal from BJ photogeo posterior
            Vgrf computed via astropy Galactocentric (default variant).
  pass 2 (refine) -- 5,000 realisations for stars whose pass-1 P
                     falls in (0.30, 0.70).
  pass 3 (ultra)  -- 10,000 realisations for stars whose point-estimate
                     |Vgrf - 25| < 2 km/s.

Outputs
-------
release/v2/phase1/catalogue_v2.fits  -- full master table
release/v2/phase1/gate1D_tiering.json -- counts + 100-star convergence
release/v2/phase1/gate1D_convergence_500_vs_5000.png

Tiers (default solar params, Phase 0 lock-in):
  A: P > 0.95   |   B: P > 0.84  [HEADLINE]   |   C: P > 0.50
  D: point-est < 25 but P <= 0.50
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT_DIR = REPO / "release/v2/phase1"
SRC_CSV = REPO / CONFIG["input"]["source_csv"]
DIST_PATH = OUT_DIR / "catalogue_dist.fits"
ZP_PATH = OUT_DIR / "catalogue_zpcorr.fits"
RVQ_PATH = OUT_DIR / "catalogue_rvquality.fits"
VGRF_PATH = OUT_DIR / "catalogue_vgrf.fits"
BJ_CACHE = REPO / CONFIG["compute"]["cache_dir"] / "bailer_jones_2021.fits"

VGRF_CUT = float(CONFIG["vgrf_cutoff_kms"])
TIERS = CONFIG["tiers"]
MC = CONFIG["mc"]
DEFAULT = CONFIG["solar_variants"]["default"]


def galcen_frame() -> coord.Galactocentric:
    return coord.Galactocentric(
        galcen_distance=DEFAULT["R0_kpc"] * u.kpc,
        z_sun=DEFAULT["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            DEFAULT["U_kms"] * u.km/u.s,
            (DEFAULT["Vc_kms"] + DEFAULT["V_kms"]) * u.km/u.s,
            DEFAULT["W_kms"] * u.km/u.s,
        ),
    )


def split_normal_sample(med: np.ndarray, lo: np.ndarray, hi: np.ndarray,
                        n_samples: int, rng: np.random.Generator) -> np.ndarray:
    """Sample a split-normal posterior. Shape: (n_stars, n_samples)."""
    n_stars = len(med)
    sig_lo = np.maximum(med - lo, 1e-3)  # pc; floor to avoid degenerate
    sig_hi = np.maximum(hi - med, 1e-3)
    u_ = rng.standard_normal((n_stars, n_samples))
    sig = np.where(u_ < 0,
                   sig_lo[:, None],
                   sig_hi[:, None])
    return med[:, None] + u_ * sig


def covariance_3x3(parallax_err, pmra_err, pmdec_err,
                   plx_pmra_corr, plx_pmdec_corr, pmra_pmdec_corr):
    """Build (n_stars, 3, 3) covariance for (parallax, pmra, pmdec)."""
    n = len(parallax_err)
    sig = np.stack([parallax_err, pmra_err, pmdec_err], axis=-1)  # (n,3)
    rho = np.zeros((n, 3, 3))
    rho[:, 0, 0] = 1.0
    rho[:, 1, 1] = 1.0
    rho[:, 2, 2] = 1.0
    rho[:, 0, 1] = rho[:, 1, 0] = plx_pmra_corr
    rho[:, 0, 2] = rho[:, 2, 0] = plx_pmdec_corr
    rho[:, 1, 2] = rho[:, 2, 1] = pmra_pmdec_corr
    sig_outer = sig[:, :, None] * sig[:, None, :]
    return rho * sig_outer


def cholesky_psd(C):
    """Cholesky for a stack of (3,3) PSD matrices, regularising if needed."""
    n = C.shape[0]
    out = np.zeros_like(C)
    fail = 0
    for i in range(n):
        try:
            out[i] = np.linalg.cholesky(C[i])
        except np.linalg.LinAlgError:
            # add tiny diagonal jitter
            out[i] = np.linalg.cholesky(C[i] + 1e-9 * np.eye(3))
            fail += 1
    if fail:
        print(f"[1D] cholesky regularised on {fail} stars")
    return out


def mc_pass(df_sub: pd.DataFrame, n_samples: int, rng: np.random.Generator,
            galcen: coord.Galactocentric) -> tuple[np.ndarray, np.ndarray]:
    """Run one MC pass on a subset of stars.

    Returns (P, vgrf_samples_med) where P is shape (n_stars,) and
    vgrf_samples_med is the median of the per-star Vgrf samples.
    """
    n = len(df_sub)

    # 1) Sample (parallax_zpcorr, pmra, pmdec) jointly.
    mu = np.stack([
        df_sub["parallax_zpcorr"].to_numpy(dtype=float),
        df_sub["pmra"].to_numpy(dtype=float),
        df_sub["pmdec"].to_numpy(dtype=float),
    ], axis=-1)  # (n,3)
    cov = covariance_3x3(
        df_sub["parallax_error"].to_numpy(dtype=float),
        df_sub["pmra_error"].to_numpy(dtype=float),
        df_sub["pmdec_error"].to_numpy(dtype=float),
        df_sub["parallax_pmra_corr"].to_numpy(dtype=float),
        df_sub["parallax_pmdec_corr"].to_numpy(dtype=float),
        df_sub["pmra_pmdec_corr"].to_numpy(dtype=float),
    )
    L = cholesky_psd(cov)  # (n,3,3)
    z = rng.standard_normal((n, n_samples, 3))   # (n, S, 3)
    samp = mu[:, None, :] + np.einsum("nij,nsj->nsi", L, z)  # (n,S,3)
    plx_s = samp[:, :, 0]
    pmra_s = samp[:, :, 1]
    pmdec_s = samp[:, :, 2]

    # 2) RV independent.
    rv_mu = df_sub["radial_velocity"].to_numpy(dtype=float)
    rv_err = df_sub["radial_velocity_error"].to_numpy(dtype=float)
    rv_s = rv_mu[:, None] + rng.standard_normal((n, n_samples)) * rv_err[:, None]

    # 3) Distance from BJ split-normal (or fallback gaussian on 1/plx).
    dist_med = df_sub["dist_pc"].to_numpy(dtype=float)
    dist_lo = df_sub["dist_lo_pc"].to_numpy(dtype=float)
    dist_hi = df_sub["dist_hi_pc"].to_numpy(dtype=float)
    dist_s = split_normal_sample(dist_med, dist_lo, dist_hi, n_samples, rng)
    dist_s = np.maximum(dist_s, 1.0)  # floor at 1 pc

    # 4) Tile RA/Dec (no astrometric error broadcast; ra/dec_error is μas).
    ra = np.broadcast_to(df_sub["ra"].to_numpy(dtype=float)[:, None],
                         (n, n_samples)).reshape(-1)
    dec = np.broadcast_to(df_sub["dec"].to_numpy(dtype=float)[:, None],
                          (n, n_samples)).reshape(-1)

    icrs = coord.SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        distance=dist_s.reshape(-1) * u.pc,
        pm_ra_cosdec=pmra_s.reshape(-1) * u.mas/u.yr,
        pm_dec=pmdec_s.reshape(-1) * u.mas/u.yr,
        radial_velocity=rv_s.reshape(-1) * u.km/u.s,
        frame="icrs",
    )
    g = icrs.transform_to(galcen)
    vx = g.v_x.to_value(u.km/u.s)
    vy = g.v_y.to_value(u.km/u.s)
    vz = g.v_z.to_value(u.km/u.s)
    vgrf = np.sqrt(vx*vx + vy*vy + vz*vz).reshape(n, n_samples)

    P = np.mean(vgrf < VGRF_CUT, axis=1)
    vgrf_med = np.median(vgrf, axis=1)
    return P, vgrf_med


def main() -> int:
    t0 = time.time()
    src = pd.read_csv(SRC_CSV)
    zp = Table.read(ZP_PATH).to_pandas()[
        ["source_id", "parallax_zpcorr", "parallax_zpcorr_error",
         "zpcorr_value_uas", "zpcorr_valid"]
    ]
    dist = Table.read(DIST_PATH).to_pandas()[
        ["source_id", "dist_pc", "dist_lo_pc", "dist_hi_pc", "dist_source"]
    ]
    rvq = Table.read(RVQ_PATH).to_pandas()[
        ["source_id", "rv_quality"]
    ]
    vgrf = Table.read(VGRF_PATH).to_pandas()

    df = (src.merge(zp, on="source_id")
              .merge(dist, on="source_id")
              .merge(rvq, on="source_id")
              .merge(vgrf[["source_id", "vgrf_default", "vgrf_grav22",
                            "vgrf_lsr6", "vgrf_rb20"]],
                     on="source_id"))
    print(f"[1D] {len(df)} stars after merges")

    rng = np.random.default_rng(MC["random_seed"])
    galcen = galcen_frame()

    # ----- pass 1: base 500 -----
    n0 = MC["base_realisations"]
    print(f"[1D] pass 1 ({n0} realisations) on all {len(df)} stars")
    t = time.time()
    P, vgrf_med = mc_pass(df, n0, rng, galcen)
    n_real = np.full(len(df), n0, dtype=int)
    print(f"[1D] pass 1 done in {time.time()-t:.1f}s")

    # ----- pass 2: refine in transition band -----
    refine_mask = (P > MC["refine_lo"]) & (P < MC["refine_hi"])
    n_refine = int(refine_mask.sum())
    print(f"[1D] pass 2 ({MC['refine_realisations']} realisations) "
          f"on {n_refine} transition-band stars")
    if n_refine > 0:
        t = time.time()
        P_ref, _ = mc_pass(df.loc[refine_mask].reset_index(drop=True),
                           MC["refine_realisations"], rng, galcen)
        P[refine_mask.to_numpy() if hasattr(refine_mask, "to_numpy") else refine_mask] = P_ref
        n_real[refine_mask.to_numpy() if hasattr(refine_mask, "to_numpy") else refine_mask] = MC["refine_realisations"]
        print(f"[1D] pass 2 done in {time.time()-t:.1f}s")

    # ----- pass 3: ultra near the cut -----
    near_cut_mask = np.abs(df["vgrf_default"].to_numpy() - VGRF_CUT) < MC["ultra_window_kms"]
    n_ultra = int(near_cut_mask.sum())
    print(f"[1D] pass 3 ({MC['ultra_realisations']} realisations) "
          f"on {n_ultra} near-cut stars")
    if n_ultra > 0:
        t = time.time()
        P_ult, _ = mc_pass(df.loc[near_cut_mask].reset_index(drop=True),
                           MC["ultra_realisations"], rng, galcen)
        P[near_cut_mask] = P_ult
        n_real[near_cut_mask] = MC["ultra_realisations"]
        print(f"[1D] pass 3 done in {time.time()-t:.1f}s")

    # ----- assign tiers -----
    tier = np.full(len(df), "X", dtype=object)
    point_below = df["vgrf_default"].to_numpy() < VGRF_CUT
    tier[point_below] = "D"
    tier[(P > TIERS["C"])] = "C"
    tier[(P > TIERS["B"])] = "B"
    tier[(P > TIERS["A"])] = "A"
    tier[(~point_below) & (P <= TIERS["C"])] = "X"  # not in any tier

    # ----- assemble master catalogue -----
    out = Table()
    out["source_id"] = df["source_id"].to_numpy()
    out["ra"] = df["ra"].to_numpy()
    out["dec"] = df["dec"].to_numpy()
    out["parallax"] = df["parallax"].to_numpy()
    out["parallax_error"] = df["parallax_error"].to_numpy()
    out["zpcorr_uas"] = df["zpcorr_value_uas"].to_numpy()
    out["zpcorr_valid"] = df["zpcorr_valid"].to_numpy()
    out["parallax_zpcorr"] = df["parallax_zpcorr"].to_numpy()
    out["dist_pc"] = df["dist_pc"].to_numpy()
    out["dist_lo_pc"] = df["dist_lo_pc"].to_numpy()
    out["dist_hi_pc"] = df["dist_hi_pc"].to_numpy()
    out["dist_source"] = df["dist_source"].to_numpy().astype(str)
    out["pmra"] = df["pmra"].to_numpy()
    out["pmdec"] = df["pmdec"].to_numpy()
    out["pmra_error"] = df["pmra_error"].to_numpy()
    out["pmdec_error"] = df["pmdec_error"].to_numpy()
    out["radial_velocity"] = df["radial_velocity"].to_numpy()
    out["radial_velocity_error"] = df["radial_velocity_error"].to_numpy()
    out["rv_quality"] = df["rv_quality"].to_numpy().astype(str)
    out["vgrf_default"] = df["vgrf_default"].to_numpy()
    out["vgrf_grav22"] = df["vgrf_grav22"].to_numpy()
    out["vgrf_lsr6"] = df["vgrf_lsr6"].to_numpy()
    out["vgrf_rb20"] = df["vgrf_rb20"].to_numpy()
    out["P_vgrf_below_25"] = P
    out["mc_realisations"] = n_real
    out["tier"] = tier.astype(str)

    out_path = OUT_DIR / "catalogue_v2.fits"
    out.write(out_path, format="fits", overwrite=True)
    print(f"[1D] wrote {out_path.relative_to(REPO)}")

    # ----- gate 1D -----
    counts = {t_: int((tier == t_).sum()) for t_ in ("A", "B", "C", "D", "X")}
    summary = {
        "n_processed": int(len(df)),
        "tier_counts": counts,
        "headline_tier_B": counts["B"] + counts["A"],   # B includes A
        "tier_C_total":   counts["C"] + counts["B"] + counts["A"],
        "v1_published_count": 1248,
        "v1_csv_count": 2859,
        "vgrf_cut_kms": VGRF_CUT,
        "mc_base": n0,
        "mc_refined_count": n_refine,
        "mc_ultra_count": n_ultra,
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    # ----- 100-star convergence test (500 vs 5000) -----
    sample_idx = rng.choice(len(df), size=100, replace=False)
    df100 = df.iloc[sample_idx].reset_index(drop=True)
    P_500a, _ = mc_pass(df100, 500, np.random.default_rng(MC["random_seed"]+1), galcen)
    P_5000, _ = mc_pass(df100, 5000, np.random.default_rng(MC["random_seed"]+2), galcen)
    dP = P_5000 - P_500a
    summary["convergence_100_star"] = {
        "max_abs_dP": float(np.max(np.abs(dP))),
        "median_abs_dP": float(np.median(np.abs(dP))),
        "p84_abs_dP": float(np.percentile(np.abs(dP), 84)),
    }

    (OUT_DIR / "gate1D_tiering.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(P_500a, P_5000, s=18, alpha=0.7)
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("P(Vgrf<25)  -- 500 realisations")
    ax.set_ylabel("P(Vgrf<25)  -- 5000 realisations")
    ax.set_title(f"convergence: max|ΔP| = {summary['convergence_100_star']['max_abs_dP']:.3f}")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "gate1D_convergence_500_vs_5000.png", dpi=140)
    plt.close(fig)

    if counts["B"] + counts["A"] < 200:
        print("[1D] ANOMALY: Tier B count < 200; STOP.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
