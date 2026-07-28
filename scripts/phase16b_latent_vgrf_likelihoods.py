"""Phase 16B -- per-star Vgrf likelihood grids.

Phase 0E stores only the scalar membership probability P(Vgrf < 25 km/s) for
each candidate. The latent-distribution reconstruction in Phase 16C needs the
full shape of each star's Vgrf posterior, not just its mass below the
threshold. This script re-runs the same Monte Carlo forward model as Phase 0E
and records, for every candidate, a histogram of its Vgrf realisations on a
fixed grid.

The forward model is identical to Phase 0E: draws from the Gaia
parallax--proper-motion covariance submatrix with the Lindegren et al. (2021)
zero-point correction, a split-normal draw on the Bailer-Jones et al. (2021)
photogeometric distance (inverse-parallax fallback), a Gaussian draw on the
DR3 radial velocity, and an exact Astropy Galactocentric transformation.

Outputs: <out-dir>/vgrf_likelihood_grid.npz
  grid_edges  (NBIN+1,)  bin edges in km/s
  L           (NSTAR, NBIN) row-normalised over the FULL grid; rows sum to
              the fraction of realisations landing inside the grid
  source_id   (NSTAR,)
  n_samples   (NSTAR,)
  vgrf_mean, vgrf_std, p_below_25  (NSTAR,)
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import pandas as pd

BUNDLE = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = BUNDLE / "private_inputs" / "expanded_candidates_mc_tiered.csv"
DEFAULT_OUT_DIR = BUNDLE / "phase14" / "latent_deconvolution"

VGRF_CUT = 25.0
GRID_MAX = 120.0
GRID_STEP = 1.0

# Heavier than the Phase 0E tiering schedule: a histogram needs more draws per
# star than a single tail probability does.
MC = {
    "base_realisations": 2000,
    "refine_realisations": 20000,
    "refine_lo": 5.0,    # |point Vgrf - 25| below this -> refine
    "random_seed": 20260727,
}
DEFAULT_SOLAR = {
    "R0_kpc": 8.178,
    "z_sun_pc": 25.0,
    "Vc_kms": 229.0,
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


def dist_columns(df: pd.DataFrame):
    bj_ok = df["rpgeo"].notna().to_numpy()
    med = df["distance_inv_zpcorr_pc"].to_numpy(dtype=float).copy()
    lo, hi = med.copy(), med.copy()
    med[bj_ok] = df.loc[bj_ok, "rpgeo"].to_numpy(dtype=float)
    lo[bj_ok] = df.loc[bj_ok, "rpgeo_lo"].to_numpy(dtype=float)
    hi[bj_ok] = df.loc[bj_ok, "rpgeo_hi"].to_numpy(dtype=float)
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
    return med, lo, hi


def split_normal_sample(med, lo, hi, n_samples, rng):
    sig_lo = np.maximum(med - lo, 1e-3)
    sig_hi = np.maximum(hi - med, 1e-3)
    u0 = rng.standard_normal((len(med), n_samples))
    sig = np.where(u0 < 0, sig_lo[:, None], sig_hi[:, None])
    return med[:, None] + u0 * sig


def covariance_3x3(df: pd.DataFrame) -> np.ndarray:
    sig = np.stack([
        np.nan_to_num(df["parallax_error"].to_numpy(dtype=float), nan=0.0),
        np.nan_to_num(df["pmra_error"].to_numpy(dtype=float), nan=0.0),
        np.nan_to_num(df["pmdec_error"].to_numpy(dtype=float), nan=0.0),
    ], axis=-1)
    c01 = np.clip(np.nan_to_num(df["parallax_pmra_corr"].to_numpy(dtype=float), nan=0.0), -0.999, 0.999)
    c02 = np.clip(np.nan_to_num(df["parallax_pmdec_corr"].to_numpy(dtype=float), nan=0.0), -0.999, 0.999)
    c12 = np.clip(np.nan_to_num(df["pmra_pmdec_corr"].to_numpy(dtype=float), nan=0.0), -0.999, 0.999)
    n = len(df)
    rho = np.zeros((n, 3, 3))
    rho[:, 0, 0] = rho[:, 1, 1] = rho[:, 2, 2] = 1.0
    rho[:, 0, 1] = rho[:, 1, 0] = c01
    rho[:, 0, 2] = rho[:, 2, 0] = c02
    rho[:, 1, 2] = rho[:, 2, 1] = c12
    return rho * (sig[:, :, None] * sig[:, None, :])


def cholesky_psd(cov: np.ndarray) -> np.ndarray:
    out = np.zeros_like(cov)
    for i in range(cov.shape[0]):
        try:
            out[i] = np.linalg.cholesky(cov[i])
        except np.linalg.LinAlgError:
            out[i] = np.linalg.cholesky(cov[i] + 1e-9 * np.eye(3))
    return out


def mc_histograms(df, dist_med, dist_lo, dist_hi, n_samples, rng, edges, star_batch):
    """Return (hist counts (n,NBIN), mean, std, p_below_25) for one MC schedule."""
    frame = galcen_frame()
    nbin = len(edges) - 1
    H = np.zeros((len(df), nbin))
    mean = np.full(len(df), np.nan)
    std = np.full(len(df), np.nan)
    pbel = np.full(len(df), np.nan)

    for start in range(0, len(df), star_batch):
        stop = min(start + star_batch, len(df))
        sub = df.iloc[start:stop].reset_index(drop=True)
        n = len(sub)
        mu = np.stack([
            sub["parallax_zpcorr"].to_numpy(dtype=float),
            sub["pmra"].to_numpy(dtype=float),
            sub["pmdec"].to_numpy(dtype=float),
        ], axis=-1)
        L = cholesky_psd(covariance_3x3(sub))
        z = rng.standard_normal((n, n_samples, 3))
        samp = mu[:, None, :] + np.einsum("nij,nsj->nsi", L, z)

        rv_mu = sub["radial_velocity"].to_numpy(dtype=float)
        rv_err = np.maximum(np.nan_to_num(sub["radial_velocity_error"].to_numpy(dtype=float), nan=0.0), 0.0)
        rv_s = rv_mu[:, None] + rng.standard_normal((n, n_samples)) * rv_err[:, None]

        dist_s = np.maximum(split_normal_sample(
            dist_med[start:stop], dist_lo[start:stop], dist_hi[start:stop], n_samples, rng), 1.0)

        ra = np.broadcast_to(sub["ra"].to_numpy(dtype=float)[:, None], (n, n_samples)).reshape(-1)
        dec = np.broadcast_to(sub["dec"].to_numpy(dtype=float)[:, None], (n, n_samples)).reshape(-1)
        icrs = coord.SkyCoord(
            ra=ra * u.deg, dec=dec * u.deg, distance=dist_s.reshape(-1) * u.pc,
            pm_ra_cosdec=samp[:, :, 1].reshape(-1) * u.mas / u.yr,
            pm_dec=samp[:, :, 2].reshape(-1) * u.mas / u.yr,
            radial_velocity=rv_s.reshape(-1) * u.km / u.s, frame="icrs")
        g = icrs.transform_to(frame)
        vgrf = np.sqrt(g.v_x.to_value(u.km / u.s) ** 2
                       + g.v_y.to_value(u.km / u.s) ** 2
                       + g.v_z.to_value(u.km / u.s) ** 2).reshape(n, n_samples)

        for k in range(n):
            H[start + k] = np.histogram(vgrf[k], bins=edges)[0]
        mean[start:stop] = vgrf.mean(axis=1)
        std[start:stop] = vgrf.std(axis=1)
        pbel[start:stop] = np.mean(vgrf < VGRF_CUT, axis=1)
        print(f"  S={n_samples}: {stop:,}/{len(df):,}", flush=True)
    return H, mean, std, pbel


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", default=str(DEFAULT_INPUT))
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    ap.add_argument("--star-batch", type=int, default=400)
    args = ap.parse_args()

    t0 = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(MC["random_seed"])

    df = pd.read_csv(args.input_csv, low_memory=False)
    ok = (df["radial_velocity"].notna() & df["pmra"].notna()
          & df["pmdec"].notna() & df["vgrf_default_exact"].notna())
    df = df.loc[ok].reset_index(drop=True)
    print(f"Loaded {len(df):,} candidates with complete 6D input")

    edges = np.arange(0.0, GRID_MAX + GRID_STEP, GRID_STEP)
    dist_med, dist_lo, dist_hi = dist_columns(df)

    point = df["vgrf_default_exact"].to_numpy(dtype=float)
    refine = np.abs(point - VGRF_CUT) < MC["refine_lo"]
    print(f"Refined (|Vgrf-25| < {MC['refine_lo']}): {int(refine.sum()):,}")

    H = np.zeros((len(df), len(edges) - 1))
    mean = np.full(len(df), np.nan)
    std = np.full(len(df), np.nan)
    pbel = np.full(len(df), np.nan)
    nsamp = np.zeros(len(df), dtype=int)

    for mask, S in [(~refine, MC["base_realisations"]), (refine, MC["refine_realisations"])]:
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            continue
        sub = df.iloc[idx].reset_index(drop=True)
        h, m, s, p = mc_histograms(sub, dist_med[idx], dist_lo[idx], dist_hi[idx],
                                   S, rng, edges, args.star_batch)
        H[idx], mean[idx], std[idx], pbel[idx] = h, m, s, p
        nsamp[idx] = S

    # Row-normalise by the number of draws so each row sums to the fraction of
    # realisations landing inside the grid (i.e. mass above GRID_MAX is lost,
    # and that loss is recorded rather than hidden by renormalisation).
    L = H / nsamp[:, None]
    inside = L.sum(axis=1)

    np.savez_compressed(
        out_dir / "vgrf_likelihood_grid.npz",
        grid_edges=edges, L=L.astype(np.float32),
        source_id=df["source_id"].to_numpy(dtype=np.int64),
        n_samples=nsamp, vgrf_point=point,
        vgrf_mean=mean, vgrf_std=std, p_below_25=pbel, inside_grid=inside,
    )
    summary = {
        "n_stars": int(len(df)),
        "grid_max_kms": GRID_MAX,
        "grid_step_kms": GRID_STEP,
        "mc_schedule": MC,
        "n_refined": int(refine.sum()),
        "median_vgrf_sigma_kms": float(np.nanmedian(std)),
        "vgrf_sigma_percentiles": {p: float(np.nanpercentile(std, p)) for p in (5, 16, 50, 84, 95)},
        "mean_mass_inside_grid": float(np.mean(inside)),
        "min_mass_inside_grid": float(np.min(inside)),
        "p_below_25_sum": float(np.nansum(pbel)),
        "p_below_25_vs_phase0e_max_abs_diff": float(
            np.nanmax(np.abs(pbel - df["P_vgrf_below_25"].to_numpy(dtype=float)))),
        "elapsed_seconds": round(time.time() - t0, 1),
    }
    (out_dir / "vgrf_likelihood_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
