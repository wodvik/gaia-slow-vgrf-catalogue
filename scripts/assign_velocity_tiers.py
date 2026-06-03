"""
Expanded-candidate point-estimate and MC tiering pass.

This script takes the Phase 0D expanded candidate table, recomputes exact
Astropy Galactocentric point-estimate Vgrf values, and assigns
P(Vgrf < 25 km/s) tiers using the same adaptive schedule as Phase 1D:

  - 500 realisations for all expanded candidates
  - 5000 realisations for base-pass P in (0.30, 0.70)
  - 10000 realisations for point-estimate |Vgrf - 25| < 2 km/s

Outputs are written to D:/GAIA/parent_scan by default so the large rebuild
artifacts stay off the project drive.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u


REPO = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO / "data" / "slow_stars_expanded_candidates_vgrf50.csv"
DEFAULT_OUT_DIR = Path("D:/GAIA/parent_scan")

VGRF_CUT = 25.0
TIERS = {"A": 0.95, "B": 0.84, "C": 0.50}
MC = {
    "base_realisations": 500,
    "refine_realisations": 5000,
    "ultra_realisations": 10000,
    "refine_lo": 0.30,
    "refine_hi": 0.70,
    "ultra_window_kms": 2.0,
    "random_seed": 20260502,
}
DEFAULT_SOLAR = {
    "R0_kpc": 8.178,
    "z_sun_pc": 25.0,
    "Vc_kms": 232.0,
    "U_kms": 11.1,
    "V_kms": 12.24,
    "W_kms": 7.25,
}


def galcen_frame() -> coord.Galactocentric:
    return coord.Galactocentric(
        galcen_distance=DEFAULT_SOLAR["R0_kpc"] * u.kpc,
        z_sun=DEFAULT_SOLAR["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            DEFAULT_SOLAR["U_kms"] * u.km / u.s,
            (DEFAULT_SOLAR["Vc_kms"] + DEFAULT_SOLAR["V_kms"]) * u.km / u.s,
            DEFAULT_SOLAR["W_kms"] * u.km / u.s,
        ),
    )


def dist_columns(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    bj_ok = df["rpgeo"].notna().to_numpy()
    med = df["distance_inv_zpcorr_pc"].to_numpy(dtype=float).copy()
    lo = med.copy()
    hi = med.copy()
    source = np.full(len(df), "inv_parallax_zpcorr", dtype=object)

    med[bj_ok] = df.loc[bj_ok, "rpgeo"].to_numpy(dtype=float)
    lo[bj_ok] = df.loc[bj_ok, "rpgeo_lo"].to_numpy(dtype=float)
    hi[bj_ok] = df.loc[bj_ok, "rpgeo_hi"].to_numpy(dtype=float)
    source[bj_ok] = "bailer_jones_2021_photogeo"

    fallback = ~bj_ok
    if fallback.any():
        plx = df.loc[fallback, "parallax_zpcorr"].to_numpy(dtype=float)
        plx_err = df.loc[fallback, "parallax_error"].to_numpy(dtype=float)
        inv_pc = np.where(plx > 0, 1000.0 / plx, np.nan)
        sigma_pc = np.abs(inv_pc * (plx_err / plx))
        med[fallback] = inv_pc
        lo[fallback] = inv_pc - sigma_pc
        hi[fallback] = inv_pc + sigma_pc

    lo = np.where(np.isfinite(lo), lo, med)
    hi = np.where(np.isfinite(hi), hi, med)
    return med, lo, hi, source


def point_vgrf(df: pd.DataFrame, dist_pc: np.ndarray, batch_size: int) -> np.ndarray:
    frame = galcen_frame()
    out = np.full(len(df), np.nan, dtype=float)
    for start in range(0, len(df), batch_size):
        stop = min(start + batch_size, len(df))
        sub = df.iloc[start:stop]
        d = dist_pc[start:stop]
        finite = (
            np.isfinite(d)
            & np.isfinite(sub["ra"].to_numpy(dtype=float))
            & np.isfinite(sub["dec"].to_numpy(dtype=float))
            & np.isfinite(sub["pmra"].to_numpy(dtype=float))
            & np.isfinite(sub["pmdec"].to_numpy(dtype=float))
            & np.isfinite(sub["radial_velocity"].to_numpy(dtype=float))
        )
        if not finite.any():
            continue
        ss = sub.loc[finite]
        dd = d[finite]
        icrs = coord.SkyCoord(
            ra=ss["ra"].to_numpy(dtype=float) * u.deg,
            dec=ss["dec"].to_numpy(dtype=float) * u.deg,
            distance=dd * u.pc,
            pm_ra_cosdec=ss["pmra"].to_numpy(dtype=float) * u.mas / u.yr,
            pm_dec=ss["pmdec"].to_numpy(dtype=float) * u.mas / u.yr,
            radial_velocity=ss["radial_velocity"].to_numpy(dtype=float) * u.km / u.s,
            frame="icrs",
        )
        g = icrs.transform_to(frame)
        vx = g.v_x.to_value(u.km / u.s)
        vy = g.v_y.to_value(u.km / u.s)
        vz = g.v_z.to_value(u.km / u.s)
        local_idx = np.flatnonzero(finite)
        out[start + local_idx] = np.sqrt(vx * vx + vy * vy + vz * vz)
    return out


def split_normal_sample(
    med: np.ndarray, lo: np.ndarray, hi: np.ndarray, n_samples: int, rng: np.random.Generator
) -> np.ndarray:
    sig_lo = np.maximum(med - lo, 1e-3)
    sig_hi = np.maximum(hi - med, 1e-3)
    u0 = rng.standard_normal((len(med), n_samples))
    sig = np.where(u0 < 0, sig_lo[:, None], sig_hi[:, None])
    return med[:, None] + u0 * sig


def covariance_3x3(df: pd.DataFrame) -> np.ndarray:
    parallax_err = np.nan_to_num(df["parallax_error"].to_numpy(dtype=float), nan=0.0)
    pmra_err = np.nan_to_num(df["pmra_error"].to_numpy(dtype=float), nan=0.0)
    pmdec_err = np.nan_to_num(df["pmdec_error"].to_numpy(dtype=float), nan=0.0)
    sig = np.stack([parallax_err, pmra_err, pmdec_err], axis=-1)

    plx_pmra = np.nan_to_num(df["parallax_pmra_corr"].to_numpy(dtype=float), nan=0.0)
    plx_pmdec = np.nan_to_num(df["parallax_pmdec_corr"].to_numpy(dtype=float), nan=0.0)
    pmra_pmdec = np.nan_to_num(df["pmra_pmdec_corr"].to_numpy(dtype=float), nan=0.0)
    for arr in (plx_pmra, plx_pmdec, pmra_pmdec):
        np.clip(arr, -0.999, 0.999, out=arr)

    n = len(df)
    rho = np.zeros((n, 3, 3))
    rho[:, 0, 0] = 1.0
    rho[:, 1, 1] = 1.0
    rho[:, 2, 2] = 1.0
    rho[:, 0, 1] = rho[:, 1, 0] = plx_pmra
    rho[:, 0, 2] = rho[:, 2, 0] = plx_pmdec
    rho[:, 1, 2] = rho[:, 2, 1] = pmra_pmdec
    return rho * (sig[:, :, None] * sig[:, None, :])


def cholesky_psd(cov: np.ndarray) -> np.ndarray:
    out = np.zeros_like(cov)
    for i in range(cov.shape[0]):
        try:
            out[i] = np.linalg.cholesky(cov[i])
        except np.linalg.LinAlgError:
            out[i] = np.linalg.cholesky(cov[i] + 1e-9 * np.eye(3))
    return out


def mc_pass(
    df: pd.DataFrame,
    dist_med: np.ndarray,
    dist_lo: np.ndarray,
    dist_hi: np.ndarray,
    n_samples: int,
    rng: np.random.Generator,
    star_batch: int,
) -> np.ndarray:
    frame = galcen_frame()
    P = np.full(len(df), np.nan, dtype=float)
    for start in range(0, len(df), star_batch):
        stop = min(start + star_batch, len(df))
        sub = df.iloc[start:stop].reset_index(drop=True)
        n = len(sub)

        mu = np.stack(
            [
                sub["parallax_zpcorr"].to_numpy(dtype=float),
                sub["pmra"].to_numpy(dtype=float),
                sub["pmdec"].to_numpy(dtype=float),
            ],
            axis=-1,
        )
        cov = covariance_3x3(sub)
        L = cholesky_psd(cov)
        z = rng.standard_normal((n, n_samples, 3))
        samp = mu[:, None, :] + np.einsum("nij,nsj->nsi", L, z)
        pmra_s = samp[:, :, 1]
        pmdec_s = samp[:, :, 2]

        rv_mu = sub["radial_velocity"].to_numpy(dtype=float)
        rv_err = np.nan_to_num(sub["radial_velocity_error"].to_numpy(dtype=float), nan=0.0)
        rv_err = np.maximum(rv_err, 0.0)
        rv_s = rv_mu[:, None] + rng.standard_normal((n, n_samples)) * rv_err[:, None]

        dist_s = split_normal_sample(
            dist_med[start:stop], dist_lo[start:stop], dist_hi[start:stop], n_samples, rng
        )
        dist_s = np.maximum(dist_s, 1.0)

        ra = np.broadcast_to(sub["ra"].to_numpy(dtype=float)[:, None], (n, n_samples)).reshape(-1)
        dec = np.broadcast_to(sub["dec"].to_numpy(dtype=float)[:, None], (n, n_samples)).reshape(-1)
        icrs = coord.SkyCoord(
            ra=ra * u.deg,
            dec=dec * u.deg,
            distance=dist_s.reshape(-1) * u.pc,
            pm_ra_cosdec=pmra_s.reshape(-1) * u.mas / u.yr,
            pm_dec=pmdec_s.reshape(-1) * u.mas / u.yr,
            radial_velocity=rv_s.reshape(-1) * u.km / u.s,
            frame="icrs",
        )
        g = icrs.transform_to(frame)
        vx = g.v_x.to_value(u.km / u.s)
        vy = g.v_y.to_value(u.km / u.s)
        vz = g.v_z.to_value(u.km / u.s)
        vgrf = np.sqrt(vx * vx + vy * vy + vz * vz).reshape(n, n_samples)
        P[start:stop] = np.mean(vgrf < VGRF_CUT, axis=1)

        done = stop
        print(
            f"MC pass S={n_samples}: {done:,}/{len(df):,} stars complete",
            flush=True,
        )
    return P


def assign_tiers(P: np.ndarray, point: np.ndarray) -> np.ndarray:
    tier = np.full(len(P), "X", dtype=object)
    point_below = point < VGRF_CUT
    tier[point_below] = "D"
    tier[P > TIERS["C"]] = "C"
    tier[P > TIERS["B"]] = "B"
    tier[P > TIERS["A"]] = "A"
    tier[(~point_below) & (P <= TIERS["C"])] = "X"
    return tier


def summarize(df: pd.DataFrame, tier: np.ndarray, P: np.ndarray, n_real: np.ndarray) -> dict[str, Any]:
    old = df["source_in_old_preselection"].astype(bool).to_numpy()
    counts = {name: int((tier == name).sum()) for name in ["A", "B", "C", "D", "X"]}
    counts_old = {name: int(((tier == name) & old).sum()) for name in ["A", "B", "C", "D", "X"]}
    counts_new = {name: int(((tier == name) & ~old).sum()) for name in ["A", "B", "C", "D", "X"]}
    return {
        "n_processed": int(len(df)),
        "n_old_preselection": int(old.sum()),
        "n_outside_old": int((~old).sum()),
        "tier_counts": counts,
        "tier_counts_old_preselection": counts_old,
        "tier_counts_outside_old": counts_new,
        "headline_tier_A_plus_B": int(((tier == "A") | (tier == "B")).sum()),
        "orbit_summary_tier_A_plus_B_plus_C": int(((tier == "A") | (tier == "B") | (tier == "C")).sum()),
        "point_estimate_vgrf_lt25": int((df["vgrf_default_exact"] < VGRF_CUT).sum()),
        "point_estimate_vgrf_lt25_outside_old": int(((df["vgrf_default_exact"] < VGRF_CUT) & ~old).sum()),
        "mc_realisations_counts": {str(k): int((n_real == k).sum()) for k in sorted(set(n_real.tolist()))},
        "P_vgrf_below_25_min": float(np.nanmin(P)),
        "P_vgrf_below_25_max": float(np.nanmax(P)),
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    t0 = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    df["source_id"] = df["source_id"].astype("int64")
    df["source_in_old_preselection"] = df["source_in_old_preselection"].astype(bool)

    dist_med, dist_lo, dist_hi, dist_source = dist_columns(df)
    print(f"Loaded {len(df):,} expanded candidates")

    point = point_vgrf(df, dist_med, args.point_batch)
    df["dist_pc_final_screen"] = dist_med
    df["dist_lo_pc_final_screen"] = dist_lo
    df["dist_hi_pc_final_screen"] = dist_hi
    df["dist_source_final_screen"] = dist_source
    df["vgrf_default_exact"] = point
    point_path = out_dir / "expanded_candidates_point_vgrf.csv"
    df.to_csv(point_path, index=False)
    print(f"Wrote point estimates: {point_path}")
    print(
        f"Point-estimate Vgrf<{VGRF_CUT}: {int((point < VGRF_CUT).sum())} "
        f"({int(((point < VGRF_CUT) & ~df['source_in_old_preselection'].to_numpy()).sum())} outside old)"
    )

    rng = np.random.default_rng(MC["random_seed"])
    P = mc_pass(
        df,
        dist_med,
        dist_lo,
        dist_hi,
        MC["base_realisations"],
        rng,
        args.base_star_batch,
    )
    n_real = np.full(len(df), MC["base_realisations"], dtype=int)
    np.savez(out_dir / "expanded_mc_after_base.npz", P=P, n_real=n_real)

    refine_mask = (P > MC["refine_lo"]) & (P < MC["refine_hi"])
    print(f"Refine pass: {int(refine_mask.sum())} transition-band stars")
    if refine_mask.any():
        idx = np.flatnonzero(refine_mask)
        P_ref = mc_pass(
            df.iloc[idx].reset_index(drop=True),
            dist_med[idx],
            dist_lo[idx],
            dist_hi[idx],
            MC["refine_realisations"],
            rng,
            args.refine_star_batch,
        )
        P[idx] = P_ref
        n_real[idx] = MC["refine_realisations"]
        np.savez(out_dir / "expanded_mc_after_refine.npz", P=P, n_real=n_real)

    near_cut_mask = np.abs(point - VGRF_CUT) < MC["ultra_window_kms"]
    print(f"Ultra pass: {int(near_cut_mask.sum())} near-cut point-estimate stars")
    if near_cut_mask.any():
        idx = np.flatnonzero(near_cut_mask)
        P_ult = mc_pass(
            df.iloc[idx].reset_index(drop=True),
            dist_med[idx],
            dist_lo[idx],
            dist_hi[idx],
            MC["ultra_realisations"],
            rng,
            args.ultra_star_batch,
        )
        P[idx] = P_ult
        n_real[idx] = MC["ultra_realisations"]
        np.savez(out_dir / "expanded_mc_after_ultra.npz", P=P, n_real=n_real)

    tier = assign_tiers(P, point)
    df["P_vgrf_below_25"] = P
    df["mc_realisations"] = n_real
    df["tier"] = tier
    out_csv = out_dir / "expanded_candidates_mc_tiered.csv"
    df.to_csv(out_csv, index=False)

    summary = summarize(df, tier, P, n_real)
    summary["elapsed_seconds"] = round(time.time() - t0, 1)
    summary["input_csv"] = str(args.input_csv)
    summary["output_csv"] = str(out_csv)
    summary["point_csv"] = str(point_path)
    summary["mc_schedule"] = MC
    summary_path = out_dir / "expanded_candidates_mc_tiering_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", default=str(DEFAULT_INPUT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--point-batch", type=int, default=5000)
    parser.add_argument("--base-star-batch", type=int, default=150)
    parser.add_argument("--refine-star-batch", type=int, default=30)
    parser.add_argument("--ultra-star-batch", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
