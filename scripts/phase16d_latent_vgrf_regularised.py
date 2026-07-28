"""Phase 16D -- regularised latent Vgrf reconstruction + iteration diagnostics.

Phase 16C ran plain expectation-maximisation. Unregularised EM (equivalently
Richardson-Lucy) is a maximum-likelihood deconvolution, and at high iteration
counts it amplifies noise into spurious bin-to-bin oscillations. The Phase 16C
figure shows exactly that: adjacent 1 km/s bins swinging between 0.02 and 0.115
across 30-50 km/s. Banik et al. (2026) hit the same problem and added a
gradient-smoothness penalty to their MCMC reconstruction (their Eq. 10; their
Fig. 15 shows the unsmoothed variant oscillating in the same way).

This script does two things:

  1. ITERATION DIAGNOSTIC. Tracks the recovered N(Vgrf < 25) and a roughness
     statistic as a function of EM iteration, so we can see whether the
     answer is an artefact of running to the iteration cap.

  2. REGULARISED RECONSTRUCTION (EMS; Silverman et al. 1990). A Gaussian
     smoothing kernel is applied to the latent distribution after each M
     step, which is the direct analogue of Banik's gradient penalty. We scan
     the kernel width so the sensitivity of the answer to the amount of
     smoothing is visible rather than hidden.

Outputs: phase14/latent_deconvolution/latent_vgrf_regularised_summary.json
         figures/latent_vgrf_regularised.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
IN_DIR = BUNDLE / "phase14" / "latent_deconvolution"
NPZ = IN_DIR / "vgrf_likelihood_grid.npz"
FIG_DIRS = [BUNDLE / "figures", REPO / "release" / "figures"]

VGRF_CUT = 25.0
VMAX = 60.0
CHECKPOINTS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
SMOOTH_SIGMAS = [0.0, 0.75, 1.5, 3.0]   # bins (1 km/s each)
MAX_ITER = 2000


def gaussian_kernel_matrix(m: int, sigma: float) -> np.ndarray:
    """Row-stochastic Gaussian smoothing operator on an m-bin grid."""
    if sigma <= 0:
        return np.eye(m)
    j = np.arange(m)
    K = np.exp(-0.5 * ((j[:, None] - j[None, :]) / sigma) ** 2)
    return K / K.sum(axis=1, keepdims=True)


def roughness(P: np.ndarray) -> float:
    """Mean squared second difference, normalised by mean(P)^2."""
    d2 = np.diff(P, n=2)
    return float(np.mean(d2 ** 2) / max(np.mean(P) ** 2, 1e-300))


def run_em(L: np.ndarray, below: np.ndarray, sigma: float, max_iter: int,
           checkpoints: list[int] | None = None):
    n, m = L.shape
    K = gaussian_kernel_matrix(m, sigma)
    P = np.full(m, 1.0 / m)
    track = []
    cps = set(checkpoints or [])
    for it in range(1, max_iter + 1):
        num = L * P[None, :]
        den = num.sum(axis=1, keepdims=True)
        den[den <= 0] = 1.0
        P = (num / den).sum(axis=0) / n
        if sigma > 0:
            P = K @ P
            P /= P.sum()
        if it in cps:
            track.append({"iter": it,
                          "count_below_25": float(n * P[below].sum()),
                          "roughness": roughness(P)})
    return P, track


def tier_counts(L: np.ndarray, P: np.ndarray, below: np.ndarray) -> dict:
    num = L * P[None, :]
    den = num.sum(axis=1)
    den[den <= 0] = 1.0
    post = num[:, below].sum(axis=1) / den
    return ({k: int((post > t).sum()) for k, t in
             [("A", 0.95), ("A+B", 0.84), ("A+B+C", 0.50)]}, post)


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
    below = centres < VGRF_CUT
    n = L.shape[0]
    print(f"{n:,} stars x {L.shape[1]} bins; cut at {VGRF_CUT} km/s")

    results = {}
    curves = {}
    for sigma in SMOOTH_SIGMAS:
        P, track = run_em(L, below, sigma, MAX_ITER, CHECKPOINTS)
        counts, _ = tier_counts(L, P, below)
        results[str(sigma)] = {
            "smoothing_sigma_bins": sigma,
            "count_below_25": float(n * P[below].sum()),
            "tier_counts": counts,
            "roughness": roughness(P),
            "iteration_track": track,
        }
        curves[sigma] = P
        print(f"sigma={sigma:4.2f}  N(<25)={n*P[below].sum():7.1f}  "
              f"roughness={roughness(P):9.4f}  tiers={counts}")

    summary = {
        "n_stars": int(n),
        "vmax_kms": VMAX,
        "vgrf_cut_kms": VGRF_CUT,
        "max_iter": MAX_ITER,
        "note": ("EMS (Silverman et al. 1990) regularised reconstruction; "
                 "sigma=0 reproduces the unregularised Phase 16C run."),
        "by_smoothing": results,
    }
    (IN_DIR / "latent_vgrf_regularised_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for sigma in SMOOTH_SIGMAS:
        lab = "no smoothing" if sigma == 0 else rf"$\sigma={sigma}$ bins"
        ax[0].step(centres, curves[sigma], where="mid", lw=1.4, label=lab)
    ax[0].axvline(VGRF_CUT, color="k", ls="--", lw=1.0)
    ax[0].set_xlim(0, VMAX)
    ax[0].set_xlabel(r"$V_{\rm grf}$ (km s$^{-1}$)")
    ax[0].set_ylabel("Probability per 1 km s$^{-1}$ bin")
    ax[0].legend(fontsize=8, frameon=False)
    ax[0].set_title("Latent reconstruction vs smoothing", fontsize=10)

    for sigma in SMOOTH_SIGMAS:
        tr = results[str(sigma)]["iteration_track"]
        lab = "no smoothing" if sigma == 0 else rf"$\sigma={sigma}$ bins"
        ax[1].plot([t["iter"] for t in tr], [t["count_below_25"] for t in tr],
                   marker="o", ms=3, lw=1.2, label=lab)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("EM iteration")
    ax[1].set_ylabel(r"Recovered $N(V_{\rm grf}<25)$")
    ax[1].legend(fontsize=8, frameon=False)
    ax[1].set_title("Iteration dependence of the recovered count", fontsize=10)
    fig.tight_layout()
    for fd in FIG_DIRS:
        if fd.parent.exists():
            fd.mkdir(parents=True, exist_ok=True)
            fig.savefig(fd / "latent_vgrf_regularised.png", dpi=180)
    plt.close(fig)

    print("\nIteration tracks:")
    for sigma in SMOOTH_SIGMAS:
        print(f" sigma={sigma}:")
        for t in results[str(sigma)]["iteration_track"]:
            print(f"   it={t['iter']:5d}  N(<25)={t['count_below_25']:7.1f}  "
                  f"rough={t['roughness']:9.4f}")


if __name__ == "__main__":
    main()
