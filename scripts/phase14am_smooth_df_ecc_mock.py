"""Regenerate the single-ellipsoid smooth-DF eccentricity mock and its
Anderson-Darling comparisons from released bundle products only.

This backs the manuscript sentence that the Tier A+B+C eccentricity
distribution is rejected against both the adjacent 25-50 km/s matched
control band and the single-ellipsoid smooth-continuation mock at
asymptotic p < 0.001.

Recipe (identical to the referee-readiness WP-4 construction):
1. Fit a single trivariate Gaussian to the Cartesian Galactocentric
   velocities of the full 25-260 km/s matched-control library
   (phase14/control_orbits.fits, columns vx_kms/vy_kms/vz_kms).
2. Draw 2,000,000 velocities from that ellipsoid (seed 20260523) and keep
   the |v| < 25 km/s tail.
3. Assign each low-speed draw an eccentricity borrowed from observed
   control orbits by inverse-distance-weighted sampling over the 15
   nearest control stars in standardised velocity space, keeping the mock
   tied to the same orbit library as the manuscript control comparison.
4. Truncate to the Tier A+B+C sample size and write the realisation to
   phase14/smooth_df_eccentricity_mock.csv plus an Anderson-Darling
   summary JSON.

Run from the bundle root:

    python scripts/phase14am_smooth_df_ecc_mock.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.spatial import cKDTree
from scipy.stats import anderson_ksamp

BUNDLE = Path(__file__).resolve().parents[1]
CONTROL = BUNDLE / "phase14" / "control_orbits.fits"
ORBITS = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.fits"
OUT_CSV = BUNDLE / "phase14" / "smooth_df_eccentricity_mock.csv"
OUT_JSON = BUNDLE / "phase14" / "smooth_df_ecc_mock_summary.json"

SEED = 20260523
N_DRAWS = 2_000_000
V_LOW = 25.0


def ad(sample_a: np.ndarray, sample_b: np.ndarray) -> dict:
    res = anderson_ksamp([sample_a, sample_b])
    return {
        "A2": float(res.statistic),
        "asymptotic_p_capped": float(res.significance_level),
        "p_lt_0p001": bool(res.significance_level <= 0.001),
        "n_a": int(len(sample_a)),
        "n_b": int(len(sample_b)),
    }


def main() -> int:
    ctrl = Table.read(CONTROL).to_pandas()
    for col in ("vx_kms", "vy_kms", "vz_kms"):
        if col not in ctrl.columns:
            raise RuntimeError(
                f"{CONTROL} lacks {col}; this bundle predates the "
                "self-contained mock (v1.0.7+)."
            )
    band = ctrl["band"]
    if band.dtype == object:
        band = band.str.decode("utf-8", errors="ignore")

    vel = ctrl[["vx_kms", "vy_kms", "vz_kms"]].to_numpy(float)
    ecc_ctrl = ctrl["ecc"].to_numpy(float)
    mu = vel.mean(axis=0)
    cov = np.cov(vel, rowvar=False) + np.eye(3) * 1e-6

    slow = Table.read(ORBITS)
    ecc_slow = np.asarray(slow["static_ecc"], dtype=float)
    n_slow = len(ecc_slow)

    rng = np.random.default_rng(SEED)
    mock_vel = rng.multivariate_normal(mu, cov, size=N_DRAWS)
    mock_speed = np.linalg.norm(mock_vel, axis=1)
    low_vel = mock_vel[mock_speed < V_LOW]
    if len(low_vel) < n_slow:
        raise RuntimeError("not enough low-speed smooth-DF draws")

    std = vel.std(axis=0, ddof=1)
    vel_std = (vel - mu) / std
    low_std = (low_vel[: 8 * n_slow] - mu) / std
    tree = cKDTree(vel_std)
    dist, idx = tree.query(low_std, k=15)
    w = 1.0 / np.maximum(dist, 1e-6)
    w = w / w.sum(axis=1, keepdims=True)
    choices = np.array(
        [rng.choice(idx.shape[1], p=w[i]) for i in range(idx.shape[0])]
    )
    borrowed = idx[np.arange(idx.shape[0]), choices]
    ecc_mock = ecc_ctrl[borrowed[:n_slow]]

    pd.DataFrame({"eccentricity_smooth_df_mock": ecc_mock}).to_csv(
        OUT_CSV, index=False
    )

    ecc_2550 = ecc_ctrl[(band == "vgrf_25_50").to_numpy()]
    summary = {
        "seed": SEED,
        "n_draws": N_DRAWS,
        "recipe": (
            "single trivariate-Gaussian ellipsoid fit to the 25-260 km/s "
            "matched-control velocities; |v|<25 km/s draws assigned "
            "eccentricities by inverse-distance-weighted sampling over the "
            "15 nearest control orbits in standardised velocity space"
        ),
        "mock_ecc_median": float(np.median(ecc_mock)),
        "slow_ecc_median": float(np.median(ecc_slow)),
        "ad_slow_vs_control_25_50": ad(ecc_slow, ecc_2550),
        "ad_slow_vs_smooth_df_mock": ad(ecc_slow, ecc_mock),
        "note": (
            "scipy anderson_ksamp caps the asymptotic significance level "
            "at 0.001; p_lt_0p001=true means rejection at asymptotic "
            "p<0.001 as quoted in the manuscript"
        ),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["ad_slow_vs_smooth_df_mock"], indent=2))
    print(json.dumps(summary["ad_slow_vs_control_25_50"], indent=2))
    print("wrote", OUT_CSV)
    print("wrote", OUT_JSON)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
