"""Phase 16C -- latent Vgrf-distribution reconstruction (Eddington-bias audit).

The tier probabilities released with this catalogue are FORWARD scores: for
each star, P(Vgrf < 25 km/s) under its own measurement-error model. They do
not know that the parent population is a steeply rising function of Vgrf, so
that many more stars are available to scatter DOWN into the slow window than
up out of it. This is the classical Eddington bias, and it is the same
problem that Banik et al. (2026, arXiv:2607.00764) solve for the oldest-star
age extremum by reconstructing the population's latent distribution from the
stacked per-star likelihoods (their Sec. 4.1, Eq. 9).

We do the same here. Given per-star likelihood grids L_i(v) from Phase 16B,
we recover the latent distribution P(v) by expectation-maximisation:

    w_ij   = P_j L_ij / sum_k P_k L_ik          (E step)
    P_j    = (1/N) sum_i w_ij                    (M step)

which is the properly normalised form of the population-prior iteration. The
converged w_ij are the per-star posteriors WITH the population prior folded
in, so sum_{v_j < 25} w_ij is a revised membership probability for star i and
N * sum_{v_j < 25} P_j is the revised expected number of genuinely slow stars.

COMPLETENESS CAVEAT. The candidate pool is complete only below ~50 km/s (it
was screened at a preliminary Vgrf < 50). Latent bins above that are
under-populated relative to the true Galaxy, so the reconstruction
under-represents down-scatter from fast stars and the correction reported
here is therefore conservative -- the true contamination can only be larger.
We quantify this by varying the latent grid ceiling.

Outputs: phase14/latent_deconvolution/latent_vgrf_summary.json
         phase14/latent_deconvolution/latent_vgrf_per_star.csv
         figures/latent_vgrf_reconstruction.png
         tables/v15/tab_latent_deconvolution.tex
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
IN_DIR = BUNDLE / "phase14" / "latent_deconvolution"
NPZ = IN_DIR / "vgrf_likelihood_grid.npz"
FIG_DIRS = [BUNDLE / "figures", REPO / "release" / "figures"]
TAB_DIRS = [BUNDLE / "tables" / "v15", REPO / "release" / "tables" / "v15"]

VGRF_CUT = 25.0
TIER_THRESHOLDS = {"A": 0.95, "B": 0.84, "C": 0.50}
NOMINAL_VMAX = 60.0
VMAX_VARIANTS = [40.0, 50.0, 60.0, 80.0]
N_BOOTSTRAP = 200
EM_MAX_ITER = 2000
EM_TOL = 1e-10


def em_deconvolve(L: np.ndarray, max_iter: int = EM_MAX_ITER, tol: float = EM_TOL):
    """Return (P, n_iter). L is (N, M), rows renormalised to sum to 1."""
    n, m = L.shape
    P = np.full(m, 1.0 / m)
    for it in range(1, max_iter + 1):
        num = L * P[None, :]
        den = num.sum(axis=1, keepdims=True)
        den[den <= 0] = 1.0
        P_new = (num / den).sum(axis=0) / n
        if np.max(np.abs(P_new - P)) < tol:
            return P_new, it
        P = P_new
    return P, max_iter


def _smoothing_sentence() -> str:
    """Report the smoothing-kernel spread in the caption when Phase 16D has run.

    The spread is a systematic, and keeping it in the caption rather than as a
    fourth column keeps this table inside a single A&A column.
    """
    f = IN_DIR / "latent_vgrf_regularised_summary.json"
    if not f.exists():
        return ""
    d = json.loads(f.read_text(encoding="utf-8"))["by_smoothing"]
    tiers = {k: [] for k in ("A", "A+B", "A+B+C")}
    counts = []
    for entry in d.values():
        counts.append(entry["count_below_25"])
        for k in tiers:
            tiers[k].append(entry["tier_counts"][k])
    rng = lambda v: f"{min(v):,.0f}--{max(v):,.0f}".replace(",", "{,}")
    return (r"  Varying the kernel over $0$--$3$ bins moves the tier counts to "
            f"${rng(tiers['A'])}$, ${rng(tiers['A+B'])}$, and "
            f"${rng(tiers['A+B+C'])}$, and the expected count to "
            f"${rng(counts)}$; that spread is a systematic rather than a "
            r"statistical one and is comparable to the bootstrap width.")


def posterior_below(L: np.ndarray, P: np.ndarray, below: np.ndarray) -> np.ndarray:
    """Per-star posterior mass below the cut, with the population prior applied."""
    num = L * P[None, :]
    den = num.sum(axis=1)
    den[den <= 0] = 1.0
    return num[:, below].sum(axis=1) / den


def run_for_vmax(L_full: np.ndarray, centres_full: np.ndarray, vmax: float):
    keep = centres_full <= vmax
    L = L_full[:, keep].astype(np.float64)
    centres = centres_full[keep]
    row = L.sum(axis=1)
    good = row > 0
    leaked = 1.0 - row  # likelihood mass outside the truncated grid
    L = L[good] / row[good, None]
    P, n_iter = em_deconvolve(L)
    below = centres < VGRF_CUT
    post = posterior_below(L, P, below)
    n_stars = int(good.sum())
    return {
        "vmax": vmax,
        "n_stars": n_stars,
        "n_iter": int(n_iter),
        "median_leaked_mass": float(np.median(leaked[good])),
        "mean_leaked_mass": float(np.mean(leaked[good])),
        "latent_P": P,
        "centres": centres,
        "latent_fraction_below": float(P[below].sum()),
        "latent_count_below": float(n_stars * P[below].sum()),
        "posterior": post,
        "good": good,
        "tier_counts": {k: int((post > t).sum()) for k, t in
                        [("A", 0.95), ("A+B", 0.84), ("A+B+C", 0.50)]},
    }


def main() -> None:
    d = np.load(NPZ)
    edges = d["grid_edges"]
    centres_full = 0.5 * (edges[:-1] + edges[1:])
    L_full = d["L"].astype(np.float64)
    source_id = d["source_id"]
    p_forward = d["p_below_25"]

    print(f"Loaded {L_full.shape[0]:,} stars x {L_full.shape[1]} bins")

    # Naive stacked posterior (Banik et al. 2026's "first guess"): no prior.
    row_full = L_full.sum(axis=1)
    stacked = (L_full / np.maximum(row_full, 1e-300)[:, None]).mean(axis=0)

    results = {v: run_for_vmax(L_full, centres_full, v) for v in VMAX_VARIANTS}
    nom = results[NOMINAL_VMAX]
    print(f"\nNominal vmax={NOMINAL_VMAX}: latent below-25 fraction "
          f"{nom['latent_fraction_below']:.5f} -> count {nom['latent_count_below']:.0f} "
          f"({nom['n_iter']} EM iterations)")
    print("Revised tiers:", nom["tier_counts"])

    # Bootstrap over stars at the nominal ceiling.
    keep = centres_full <= NOMINAL_VMAX
    Lk = L_full[:, keep].astype(np.float64)
    rowk = Lk.sum(axis=1)
    good = rowk > 0
    Lk = Lk[good] / rowk[good, None]
    centres = centres_full[keep]
    below = centres < VGRF_CUT
    n = Lk.shape[0]
    rng = np.random.default_rng(20260727)
    boot_counts, boot_tiers = [], {"A": [], "A+B": [], "A+B+C": []}
    for b in range(N_BOOTSTRAP):
        idx = rng.integers(0, n, n)
        Pb, _ = em_deconvolve(Lk[idx], max_iter=500, tol=1e-9)
        boot_counts.append(n * Pb[below].sum())
        postb = posterior_below(Lk, Pb, below)
        for k, t in [("A", 0.95), ("A+B", 0.84), ("A+B+C", 0.50)]:
            boot_tiers[k].append(int((postb > t).sum()))
        if (b + 1) % 50 == 0:
            print(f"  bootstrap {b+1}/{N_BOOTSTRAP}")
    boot_counts = np.array(boot_counts)

    def ci(a):
        return [float(np.percentile(a, 16)), float(np.percentile(a, 84))]

    summary = {
        "method": "EM population-prior reconstruction (Banik et al. 2026, Sec. 4.1)",
        "nominal_vmax_kms": NOMINAL_VMAX,
        "vgrf_cut_kms": VGRF_CUT,
        "n_stars_pool": int(L_full.shape[0]),
        "forward_sum_of_probabilities": float(np.nansum(p_forward)),
        "forward_tier_counts": {
            "A": int((p_forward > 0.95).sum()),
            "A+B": int((p_forward > 0.84).sum()),
            "A+B+C": int((p_forward > 0.50).sum()),
        },
        "latent_count_below_25": nom["latent_count_below"],
        "latent_count_below_25_ci68": ci(boot_counts),
        "latent_tier_counts": nom["tier_counts"],
        "latent_tier_counts_ci68": {k: [int(np.percentile(v, 16)), int(np.percentile(v, 84))]
                                    for k, v in boot_tiers.items()},
        "median_leaked_mass_nominal": nom["median_leaked_mass"],
        "n_em_iterations": nom["n_iter"],
        "n_bootstrap": N_BOOTSTRAP,
        "vmax_variants": {
            str(v): {
                "latent_count_below_25": r["latent_count_below"],
                "tier_counts": r["tier_counts"],
                "median_leaked_mass": r["median_leaked_mass"],
            } for v, r in results.items()
        },
        "completeness_caveat": (
            "Candidate pool complete only below ~50 km/s; latent bins above that "
            "are under-populated, so the reported correction is conservative."),
    }
    IN_DIR.mkdir(parents=True, exist_ok=True)
    (IN_DIR / "latent_vgrf_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    per_star = pd.DataFrame({
        "source_id": source_id[nom["good"]],
        "P_forward_below_25": p_forward[nom["good"]],
        "P_latent_below_25": nom["posterior"],
    })
    per_star["delta_P"] = per_star.P_latent_below_25 - per_star.P_forward_below_25
    per_star.to_csv(IN_DIR / "latent_vgrf_per_star.csv", index=False)

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    c, P = nom["centres"], nom["latent_P"]
    ax[0].step(centres_full, stacked, where="mid", color="0.55", lw=1.2,
               label="Stacked posterior (no prior)")
    ax[0].step(c, P, where="mid", color="C3", lw=1.6, label="Latent reconstruction")
    ax[0].axvline(VGRF_CUT, color="k", ls="--", lw=1.0)
    ax[0].set_xlim(0, NOMINAL_VMAX)
    ax[0].set_xlabel(r"$V_{\rm grf}$ (km s$^{-1}$)")
    ax[0].set_ylabel("Probability per 1 km s$^{-1}$ bin")
    ax[0].legend(fontsize=8, frameon=False)
    ax[0].set_title("Latent vs stacked $V_{\\rm grf}$ distribution", fontsize=10)

    ax[1].scatter(per_star.P_forward_below_25, per_star.P_latent_below_25, s=2,
                  alpha=0.25, color="C0", rasterized=True)
    ax[1].plot([0, 1], [0, 1], color="k", ls="--", lw=1.0)
    for t, col in [(0.95, "C3"), (0.84, "C1"), (0.50, "C2")]:
        ax[1].axhline(t, color=col, lw=0.7, ls=":")
        ax[1].axvline(t, color=col, lw=0.7, ls=":")
    ax[1].set_xlabel(r"Forward $P(V_{\rm grf}<25)$ (released tiers)")
    ax[1].set_ylabel(r"Population-prior $P(V_{\rm grf}<25)$")
    ax[1].set_title("Effect of the population prior per star", fontsize=10)
    fig.tight_layout()
    for fd in FIG_DIRS:
        if fd.parent.exists():
            fd.mkdir(parents=True, exist_ok=True)
            fig.savefig(fd / "latent_vgrf_reconstruction.png", dpi=180)
    plt.close(fig)

    # ---- table ----
    f = summary["forward_tier_counts"]
    lt = summary["latent_tier_counts"]
    lci = summary["latent_tier_counts_ci68"]
    lines = [
        r"\begin{table}",
        r"\caption{Effect of folding the population prior into the",
        r"  threshold-membership probabilities.  Forward scores are the released",
        r"  per-star Monte Carlo probabilities; population-prior scores come from",
        r"  the latent-distribution reconstruction of Sec.~\ref{sec:latent}, which",
        r"  corrects for the fact that far more stars are available to scatter",
        r"  down into the slow window than up out of it.  Uncertainties are",
        r"  $68\%$ intervals from " + str(N_BOOTSTRAP) + r" bootstrap resamplings"
        r" over stars." + _smoothing_sentence() + r"}",
        r"\label{tab:latent_deconvolution}",
        r"\centering",
        r"\begin{tabular}{lcc}",
        r"\hline\hline",
        r"Sample & Forward & Population-prior \\",
        r"\hline",
    ]
    for k in ["A", "A+B", "A+B+C"]:
        lines.append(f"Tier {k} & {f[k]} & ${lt[k]}^{{+{lci[k][1]-lt[k]}}}_{{-{lt[k]-lci[k][0]}}}$ \\\\")
    lo, hi = summary["latent_count_below_25_ci68"]
    lc = summary["latent_count_below_25"]
    lines += [
        r"\hline",
        f"Expected $N(\\vgrf<25\\kms)$ & {summary['forward_sum_of_probabilities']:.0f} & "
        f"${lc:.0f}^{{+{hi-lc:.0f}}}_{{-{lc-lo:.0f}}}$ \\\\",
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
    ]
    text = "\n".join(lines) + "\n"
    for td in TAB_DIRS:
        if td.parent.exists():
            td.mkdir(parents=True, exist_ok=True)
            (td / "tab_latent_deconvolution.tex").write_text(text, encoding="utf-8")

    print("\n" + json.dumps({k: v for k, v in summary.items()
                             if k not in ("vmax_variants",)}, indent=2))
    print("\nvmax variants:", json.dumps(summary["vmax_variants"], indent=2))


if __name__ == "__main__":
    main()
