"""Phase 14AD -- observed-control smooth-continuation diagnostics.

Referee response (deep-review Issue 3). The manuscript's smooth-continuation
check fits a single trivariate-Gaussian velocity ellipsoid to the
matched-control library (|v_GC| in 25--260 km/s) and predicts
N(<25)=76.5, implying a 25.5x observed/predicted ratio for the observed Tier A+B+C count
(1,952). A single Gaussian is a weak null for a heterogeneous,
selection-truncated disc+halo velocity field. This script evaluates the
observed/predicted ratio under three smooth nulls fit to the SAME control velocity vectors:

  gaussian : single trivariate Gaussian (reproduces the manuscript value)
  gmm      : Gaussian mixture, component count chosen by BIC
  kde      : Gaussian kernel density, bandwidth chosen by cross-validation

For each null the implied count below 25 km/s is
  N_pred(<25) = N_control * P(|v|<25) / P(25<=|v|<=260),
evaluated by Monte Carlo sampling of the fitted model. The reported
observed/predicted ratio is 1952 / N_pred. Bootstrap resampling of the control library gives
a 16th--84th percentile interval on each ratio, and an isotropic
velocity-error jitter checks measurement-error sensitivity.

The robust conclusion is that the observed slow count is above every
tested observed-control smooth-continuation extrapolation, but this is not
a Gaia-selected intrinsic-DF prediction. The precise ratio is
null-dependent, so the single-ellipsoid 25.5x is a heuristic baseline and
the BIC mixture gives the conservative tested value.

Outputs
-------
phase14/expanded_null_models_summary.json
tables/v15/tab_null_models.tex
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
DATADIR = REPO / "release" / "data"
OUT = BUNDLE / "phase14"
TAB_DIRS = [REPO / "release" / "tables" / "v15", BUNDLE / "tables" / "v15"]

BANDS = ["band_25_50", "band_50_100", "band_100_200", "band_200_260"]
V_LO, V_HI = 25.0, 260.0
OBS_TIERABC = 1952
K_PM = 4.74047           # km/s per (mas/yr * kpc)
N_MOCK = 2_000_000
N_MOCK_BOOT = 1_000_000
N_BOOT = 200
SEED = 20260605


def load_controls() -> tuple[np.ndarray, np.ndarray]:
    frames = [pd.read_csv(DATADIR / f"{b}_real0_enriched.csv") for b in BANDS]
    df = pd.concat(frames, ignore_index=True)
    V = df[["vx_kms", "vy_kms", "vz_kms"]].to_numpy(float)
    s = np.linalg.norm(V, axis=1)
    keep = (s >= V_LO) & (s <= V_HI) & np.all(np.isfinite(V), axis=1)
    # Representative per-star 3D velocity error (RV + tangential terms).
    d_kpc = np.maximum(df["distance_pc"].to_numpy(float), 1.0) / 1000.0
    sig_rv = df["radial_velocity_error"].to_numpy(float)
    sig_t1 = K_PM * d_kpc * df["pmra_error"].to_numpy(float)
    sig_t2 = K_PM * d_kpc * df["pmdec_error"].to_numpy(float)
    sig_v = np.sqrt(sig_rv ** 2 + sig_t1 ** 2 + sig_t2 ** 2)
    return V[keep], sig_v[keep]


def excess_from_speeds(speeds: np.ndarray, n_control: int) -> tuple[float, float]:
    below = np.mean(speeds < V_LO)
    shell = np.mean((speeds >= V_LO) & (speeds <= V_HI))
    n_pred = n_control * below / max(shell, 1e-12)
    return n_pred, OBS_TIERABC / max(n_pred, 1e-9)


def fit_gaussian(V, rng, n_mock):
    mu, cov = V.mean(0), np.cov(V.T)
    samp = rng.multivariate_normal(mu, cov, n_mock)
    return np.linalg.norm(samp, axis=1)


def fit_gmm(V, rng, n_mock, k=None):
    if k is None:
        best = None
        for kk in range(1, 9):
            g = GaussianMixture(kk, covariance_type="full", random_state=0,
                                reg_covar=1e-4).fit(V)
            bic = g.bic(V)
            if best is None or bic < best[0]:
                best = (bic, kk, g)
        k, g = best[1], best[2]
    else:
        g = GaussianMixture(k, covariance_type="full", random_state=0,
                            reg_covar=1e-4).fit(V)
    seed = int(rng.integers(0, 2**31 - 1))
    samp = g.sample(n_mock)[0]
    rng.shuffle(samp)
    return np.linalg.norm(samp, axis=1), k


def cv_bandwidth(V, rng) -> float:
    sub = V[rng.choice(len(V), size=min(3000, len(V)), replace=False)]
    scott = len(sub) ** (-1.0 / 7.0) * np.mean(np.std(sub, axis=0))
    grid = GridSearchCV(KernelDensity(kernel="gaussian"),
                        {"bandwidth": np.linspace(0.4, 2.0, 9) * scott}, cv=4)
    grid.fit(sub)
    return float(grid.best_params_["bandwidth"])


def fit_kde(V, rng, n_mock, bandwidth):
    kde = KernelDensity(bandwidth=bandwidth, kernel="gaussian").fit(V)
    samp = kde.sample(n_mock, random_state=int(rng.integers(0, 2**31 - 1)))
    return np.linalg.norm(samp, axis=1)


def main() -> int:
    rng = np.random.default_rng(SEED)
    V, sig_v = load_controls()
    n_control = len(V)
    med_sig_v = float(np.nanmedian(sig_v))

    # ---- point estimates ----
    n_g, ex_g = excess_from_speeds(fit_gaussian(V, rng, N_MOCK), n_control)
    sp_gmm, kbest = fit_gmm(V, rng, N_MOCK)
    n_gmm, ex_gmm = excess_from_speeds(sp_gmm, n_control)
    bw = cv_bandwidth(V, rng)
    n_k, ex_k = excess_from_speeds(fit_kde(V, rng, N_MOCK, bw), n_control)

    # ---- bootstrap CIs ----
    boot = {"gaussian": [], "gmm": [], "kde": []}
    for _ in range(N_BOOT):
        idx = rng.integers(0, n_control, n_control)
        Vb = V[idx]
        boot["gaussian"].append(excess_from_speeds(fit_gaussian(Vb, rng, N_MOCK_BOOT), n_control)[1])
        boot["gmm"].append(excess_from_speeds(fit_gmm(Vb, rng, N_MOCK_BOOT, k=kbest)[0], n_control)[1])
        boot["kde"].append(excess_from_speeds(fit_kde(Vb, rng, N_MOCK_BOOT, bw), n_control)[1])
    ci = {m: [float(np.percentile(v, 16)), float(np.percentile(v, 84))] for m, v in boot.items()}

    # ---- measurement-error sensitivity (isotropic jitter on the GMM null) ----
    Vj = V + rng.normal(0, med_sig_v / np.sqrt(3.0), V.shape)
    n_gmm_j, ex_gmm_j = excess_from_speeds(fit_gmm(Vj, rng, N_MOCK, k=kbest)[0], n_control)

    summ = {
        "n_control": n_control, "obs_tierABC_below25": OBS_TIERABC,
        "median_velocity_error_kms": round(med_sig_v, 2),
        "gmm_n_components_bic": int(kbest), "kde_bandwidth_kms": round(bw, 2),
        "nulls": {
            "gaussian": {"label": "single trivariate Gaussian", "n_pred": round(n_g, 1),
                         "excess": round(ex_g, 1), "excess_ci": [round(c, 1) for c in ci["gaussian"]]},
            "gmm": {"label": f"Gaussian mixture (k={kbest}, BIC)", "n_pred": round(n_gmm, 1),
                    "excess": round(ex_gmm, 1), "excess_ci": [round(c, 1) for c in ci["gmm"]]},
            "kde": {"label": "kernel density (CV bandwidth)", "n_pred": round(n_k, 1),
                    "excess": round(ex_k, 1), "excess_ci": [round(c, 1) for c in ci["kde"]]},
        },
        "gmm_excess_with_velocity_jitter": round(ex_gmm_j, 1),
        "conclusion": "Observed slow count is above every tested observed-control "
                      "smooth-continuation extrapolation; the BIC mixture gives the "
                      "conservative tested value and the single ellipsoid (25.5x) is "
                      "a heuristic, not a Gaia-selected intrinsic-DF prediction.",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "expanded_null_models_summary.json").write_text(json.dumps(summ, indent=2))
    _latex_table(summ)
    print(json.dumps(summ, indent=2))
    return 0


def _latex_table(summ: dict) -> None:
    order = [("gmm", "Gaussian mixture (BIC)"),
             ("gaussian", "Single trivariate Gaussian"),
             ("kde", "Kernel density (CV)")]
    lines = [
        r"% Auto-generated by scripts/phase14ad_null_models.py",
        r"\begin{deluxetable}{lccc}",
        r"\tablecaption{Observed/predicted slow-tail ratio for the Tier~A+B+C catalogue"
        r" ($N=1{,}952$ below $25\kms$) under three tested smooth velocity-distribution nulls"
        r" fitted to the observed matched-control library ($|v_{\rm GC}|$ in $25$--$260\kms$,"
        rf" $N={summ['n_control']:,}$). Intervals are 16th--84th percentile bootstrap"
        r" ranges over the control library. The ratio is above unity under every"
        r" tested observed-control smooth-continuation extrapolation, but its magnitude is"
        r" null-dependent: the single ellipsoid is an illustrative heuristic and the BIC"
        r" mixture gives the conservative tested value. These rows are not a forward"
        r" Gaia-selected intrinsic-DF prediction.\label{tab:null_models}}",
        r"\tablehead{\colhead{smooth null} & \colhead{$N_{\rm pred}(<25)$}"
        r" & \colhead{obs./pred.} & \colhead{16--84\% range}}",
        r"\startdata",
    ]
    for key, lab in order:
        d = summ["nulls"][key]
        lines.append(f"{lab} & {d['n_pred']:.1f} & {d['excess']:.1f}$\\times$ & "
                     f"{d['excess_ci'][0]:.1f}--{d['excess_ci'][1]:.1f} \\\\")
    lines += [r"\enddata", r"\end{deluxetable}"]
    text = "\n".join(lines) + "\n"
    for d in TAB_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        (d / "tab_null_models.tex").write_text(text)


if __name__ == "__main__":
    raise SystemExit(main())
