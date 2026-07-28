"""Phase 16E -- corrected purities, rank preservation, and retiering options.

Settles two questions the latent reconstruction raises:

  1. The abstract quotes "score-implied purity" for each tier, computed as the
     mean FORWARD probability of its members. If the forward scores are
     biased high, those purities are biased high too, and that is an error
     rather than a presentational choice. This script recomputes them under
     the population prior.

  2. Whether retiering on the corrected scores would reorder the catalogue or
     merely move the thresholds. If the correction is monotonic, the released
     ranking is intact and retiering is a relabelling; if not, membership
     genuinely changes and the released tiers are misordered.

Outputs: phase14/latent_deconvolution/purity_and_ranking_summary.json
         phase14/latent_deconvolution/latent_vgrf_per_star_regularised.csv
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kendalltau

BUNDLE = Path(__file__).resolve().parents[1]
IN_DIR = BUNDLE / "phase14" / "latent_deconvolution"
NPZ = IN_DIR / "vgrf_likelihood_grid.npz"

VGRF_CUT = 25.0
VMAX = 60.0
SMOOTH_SIGMA = 1.5
MAX_ITER = 500
THRESH = {"A": 0.95, "B": 0.84, "C": 0.50}


def gaussian_kernel_matrix(m: int, sigma: float) -> np.ndarray:
    j = np.arange(m)
    K = np.exp(-0.5 * ((j[:, None] - j[None, :]) / sigma) ** 2)
    return K / K.sum(axis=1, keepdims=True)


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
    p_fwd = d["p_below_25"][good]
    below = centres < VGRF_CUT
    n, m = L.shape

    K = gaussian_kernel_matrix(m, SMOOTH_SIGMA)
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
    p_lat = num[:, below].sum(axis=1) / den
    # The EM normalisation can leave a value an epsilon above unity. Clip so the
    # released products never carry a probability outside [0, 1].
    p_lat = np.clip(p_lat, 0.0, 1.0)

    df = pd.DataFrame({"source_id": sid, "P_forward": p_fwd, "P_latent": p_lat})
    df["delta"] = df.P_latent - df.P_forward
    df.to_csv(IN_DIR / "latent_vgrf_per_star_regularised.csv", index=False)

    out = {"smoothing_sigma_bins": SMOOTH_SIGMA, "n_stars": int(n)}

    # --- 1. Purity, on the RELEASED (forward-defined) tiers ---
    purity = {}
    for name, lo in [("A", THRESH["A"]), ("A+B", THRESH["B"]), ("A+B+C", THRESH["C"])]:
        sel = p_fwd > lo
        purity[name] = {
            "n_members": int(sel.sum()),
            "score_implied_purity_forward": float(p_fwd[sel].mean()),
            "score_implied_purity_population_prior": float(p_lat[sel].mean()),
            "expected_true_slow_forward": float(p_fwd[sel].sum()),
            "expected_true_slow_population_prior": float(p_lat[sel].sum()),
        }
    out["purity_on_released_tiers"] = purity

    # --- 2. Rank preservation ---
    rho = spearmanr(p_fwd, p_lat).statistic
    tau = kendalltau(p_fwd, p_lat).statistic
    # Pair inversions restricted to stars that matter (either score > 0.1)
    m2 = (p_fwd > 0.1) | (p_lat > 0.1)
    rho2 = spearmanr(p_fwd[m2], p_lat[m2]).statistic
    out["rank_preservation"] = {
        "spearman_all": float(rho),
        "kendall_tau_all": float(tau),
        "n_relevant": int(m2.sum()),
        "spearman_relevant": float(rho2),
        "monotone_violations_frac": float(1.0 - abs(tau)),
    }

    # --- 3. Retiering at the same thresholds ---
    retier = {}
    for name, lo in [("A", THRESH["A"]), ("A+B", THRESH["B"]), ("A+B+C", THRESH["C"])]:
        old = p_fwd > lo
        new = p_lat > lo
        retier[name] = {
            "n_old": int(old.sum()), "n_new": int(new.sum()),
            "n_retained": int((old & new).sum()),
            "n_dropped": int((old & ~new).sum()),
            "n_added": int((~old & new).sum()),
            "new_purity": float(p_lat[new].mean()) if new.any() else float("nan"),
        }
    out["retiering_same_thresholds"] = retier

    # --- 4. What forward threshold reproduces the corrected tier sizes? ---
    equiv = {}
    for name, lo in [("A", THRESH["A"]), ("A+B", THRESH["B"]), ("A+B+C", THRESH["C"])]:
        k = int((p_lat > lo).sum())
        if k > 0:
            equiv[name] = {"n": k,
                           "equivalent_forward_threshold": float(np.sort(p_fwd)[::-1][k - 1])}
    out["forward_threshold_matching_corrected_size"] = equiv

    (IN_DIR / "purity_and_ranking_summary.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
