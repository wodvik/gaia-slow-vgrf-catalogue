"""Phase 16L -- population-prior tier counts at alternative velocity thresholds.

tab_threshold reports how the catalogue responds to moving the selection
threshold to 15 and 20 km/s. Those rows were computed under the forward
tiering and must be restated on the adopted definition.

This needs no new deconvolution. The latent Vgrf distribution P(v) is a
property of the candidate pool, not of the threshold; only the integration
limit changes. So the same converged reconstruction yields
P_pop(Vgrf < X) = sum_{v_j < X} w_ij for any X on the grid.

Orbit statistics at each threshold come from the existing orbit products: a
lower threshold gives a strict subset of the 25 km/s Tier A+B+C, so every
star required is already integrated.

Outputs: phase14/latent_deconvolution/threshold_variation_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.stats import beta

BUNDLE = Path(__file__).resolve().parents[1]
IN_DIR = BUNDLE / "phase14" / "latent_deconvolution"
NPZ = IN_DIR / "vgrf_likelihood_grid.npz"
CAT = BUNDLE / "catalogues"

VMAX = 60.0
SMOOTH_SIGMA = 1.5
MAX_ITER = 500
THRESHOLDS = [15.0, 20.0, 25.0]


def pct(k, n):
    f = 100 * k / n
    lo = 100 * beta.ppf(0.16, k + 0.5, n - k + 0.5) if k > 0 else 0.0
    hi = 100 * beta.ppf(0.84, k + 0.5, n - k + 0.5) if k < n else 100.0
    return f, f - lo, hi - f


def main() -> None:
    d = np.load(NPZ)
    edges = d["grid_edges"]
    centres_full = 0.5 * (edges[:-1] + edges[1:])
    keep = centres_full <= VMAX
    centres = centres_full[keep]
    L = d["L"].astype(np.float64)[:, keep]
    row = L.sum(axis=1)
    good = row > 0
    L = L[good] / row[good, None]
    sid = d["source_id"][good]
    point = d["vgrf_point"][good]
    n, m = L.shape

    j = np.arange(m)
    K = np.exp(-0.5 * ((j[:, None] - j[None, :]) / SMOOTH_SIGMA) ** 2)
    K /= K.sum(axis=1, keepdims=True)
    P = np.full(m, 1.0 / m)
    for _ in range(MAX_ITER):
        num = L * P[None, :]
        den = num.sum(axis=1, keepdims=True)
        den[den <= 0] = 1.0
        P = (num / den).sum(axis=0) / n
        P = K @ P
        P /= P.sum()

    num = L * P[None, :]
    den = num.sum(axis=1)
    den[den <= 0] = 1.0

    orb = Table.read(CAT / "catalogue_expanded_orbits_tierABC.fits")
    osid = np.asarray(orb["source_id"]).astype(np.int64)
    oe = np.asarray(orb["static_ecc"], dtype=float)
    orb_map = dict(zip(osid.tolist(), oe.tolist()))

    out = {"smoothing_sigma_bins": SMOOTH_SIGMA, "thresholds": {}}
    for X in THRESHOLDS:
        below = centres < X
        ppop = np.clip(num[:, below].sum(axis=1) / den, 0.0, 1.0)
        tier_abc = ppop > 0.50
        ids = sid[tier_abc]
        ecc = np.array([orb_map[int(s)] for s in ids if int(s) in orb_map])
        f, lo, hi = pct(int((ecc > 0.95).sum()), len(ecc)) if len(ecc) else (np.nan,) * 3
        out["thresholds"][str(X)] = {
            "point_estimate_below": int((point < X).sum()),
            "tier_ABC": int(tier_abc.sum()),
            "n_with_orbits": int(len(ecc)),
            "median_ecc": float(np.median(ecc)) if len(ecc) else float("nan"),
            "frac_e_gt_0p95": f, "err_lo": lo, "err_hi": hi,
        }
        print(f"  <{X:.0f}: point={int((point<X).sum()):5d}  TierABC={int(tier_abc.sum()):5d}  "
              f"orbits={len(ecc):5d}  ecc={np.median(ecc) if len(ecc) else float('nan'):.3f}  "
              f"e>0.95={f:.1f}+{hi:.1f}-{lo:.1f}%")

    (IN_DIR / "threshold_variation_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
