"""Phase 14AC -- matched-control balance diagnostics and weighted inference.

Referee response (deep-review Issue 4). The original matched-control
comparison (phase14t_expanded_control_reweight.py) kernel-weights the four
velocity-stratified control bands to the Tier A+B+C slow catalogue on three
observability marginals (G, log d, |sin b|). A balance audit shows that the
product-of-1D-kernel weights do NOT achieve covariate balance: the standardized
mean differences (SMDs) of the matched marginals are not reduced (and can grow)
once joint covariate correlations are accounted for.

This script therefore evaluates the monotonic median-pericentre trend under an
ENSEMBLE of weighting schemes that span very different assumptions:

  raw      : no weighting
  kernel   : the paper's product-of-1D observability kernel (sensitivity)
  ipw      : propensity-score ATT weights, e(x)/(1-e(x)) from L2 logistic
  ebal     : entropy balancing -- weights that match all covariate means
             exactly (|SMD| ~ 1e-12) at maximum entropy (Hainmueller 2012)

The scientific claim -- median pericentre increases monotonically with Vgrf,
with the slow sample at the compact end -- is robust to ALL four schemes. The
entropy-balanced member is the principled weighted-inference layer: it achieves
exact covariate balance with a reported effective sample size, and the trend
survives it.

No Gaia query and no orbit re-integration: joins released control orbit products
(band, R_peri, kernel weights) to the enriched control-band CSVs that
carry the full covariate set (l, b, BP-RP, RUWE, parallax SNR, G_RVS).

Outputs
-------
phase14/expanded_control_balance.csv            per band x covariate SMDs (all schemes)
phase14/expanded_control_balance_summary.json   ESS + R_peri/ecc per scheme + monotonic flags
figures/fig_control_balance_love.pdf            Love plot (before / kernel / entropy-balanced)
tables/v15/tab_control_balance.tex              manuscript balance table
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.special import logsumexp
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
MASTER = BUNDLE / "catalogues" / "catalogue_expanded_master.fits"
CONTROLS = REPO / "release" / "_iterations" / "v2" / "phase5" / "control_orbits.fits"
WEIGHTS = BUNDLE / "phase14" / "expanded_control_weights.csv"
DATADIR = REPO / "release" / "data"
OUT = BUNDLE / "phase14"
FIG = BUNDLE / "figures"
TAB = REPO / "release" / "tables" / "v15"

SLOW_RPERI_PC = 115.516  # Tier A+B+C point-estimate slow median (fixed anchor)

BANDS = [
    ("vgrf_25_50", "25--50", "band_25_50_real0_enriched.csv"),
    ("vgrf_50_100", "50--100", "band_50_100_real0_enriched.csv"),
    ("vgrf_100_200", "100--200", "band_100_200_real0_enriched.csv"),
    ("vgrf_200_260", "200--260", "band_200_260_real0_enriched.csv"),
]

COVS = [  # (key, manuscript label)
    ("G", r"$G$"),
    ("logd", r"$\log_{10}(d/\mathrm{pc})$"),
    ("abs_sin_b", r"$|\sin b|$"),
    ("sin_l", r"$\sin l$"),
    ("cos_l", r"$\cos l$"),
    ("bp_rp", r"$BP-RP$"),
    ("log_ruwe", r"$\ln\,\mathrm{RUWE}$"),
    ("log_plx_snr", r"$\log_{10}(\varpi/\sigma_\varpi)$"),
    ("grvs", r"$G_\mathrm{RVS}$"),
]
COV_KEYS = [c[0] for c in COVS]


def decode(s: pd.Series) -> pd.Series:
    if s.dtype == object and len(s) and isinstance(s.iloc[0], (bytes, bytearray)):
        return s.str.decode("utf-8").str.strip()
    return s.astype(str).str.strip()


def feat_frame(df: pd.DataFrame) -> pd.DataFrame:
    dist = df["dist_pc"] if "dist_pc" in df.columns else df["distance_pc"]
    return pd.DataFrame({
        "G": df["phot_g_mean_mag"].to_numpy(float),
        "logd": np.log10(np.maximum(dist.to_numpy(float), 1.0)),
        "abs_sin_b": np.abs(np.sin(np.deg2rad(df["b"].to_numpy(float)))),
        "sin_l": np.sin(np.deg2rad(df["l"].to_numpy(float))),
        "cos_l": np.cos(np.deg2rad(df["l"].to_numpy(float))),
        "bp_rp": df["bp_rp"].to_numpy(float),
        "log_ruwe": np.log(np.maximum(df["ruwe"].to_numpy(float), 0.5)),
        "log_plx_snr": np.log10(np.maximum(df["parallax_over_error"].to_numpy(float), 0.1)),
        "grvs": df["grvs_mag"].to_numpy(float),
    })


def ess(w: np.ndarray) -> float:
    ok = np.isfinite(w) & (w > 0)
    w = w[ok]
    return float(w.sum() ** 2 / np.sum(w * w)) if len(w) else 0.0


def weighted_quantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[ok], w[ok]
    if len(x) == 0:
        return float("nan")
    order = np.argsort(x)
    x, w = x[order], w[order]
    cdf = np.cumsum(w) / np.sum(w)
    return float(np.interp(q, cdf, x))


def smd(slow_x: np.ndarray, ctrl_x: np.ndarray, w: np.ndarray) -> float:
    """Standardized mean difference, control(weighted) - slow, with the
    pre-match pooled-SD denominator fixed (Austin 2011)."""
    ms = np.nanmean(slow_x)
    vs = np.nanvar(slow_x, ddof=1)
    vc = np.nanvar(ctrl_x, ddof=1)
    ok = np.isfinite(ctrl_x) & np.isfinite(w) & (w > 0)
    mc = np.sum(w[ok] * ctrl_x[ok]) / np.sum(w[ok])
    sd = np.sqrt((vs + vc) / 2.0)
    return float((mc - ms) / sd) if sd > 0 else float("nan")


def ipw_att_weights(slow_feat: pd.DataFrame, ctrl_feat: pd.DataFrame) -> np.ndarray:
    Xs, Xc = slow_feat.to_numpy(float), ctrl_feat.to_numpy(float)
    X = np.vstack([Xs, Xc])
    y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xc))])
    ok = np.all(np.isfinite(X), axis=1)
    scaler = StandardScaler().fit(X[ok])
    Xz = scaler.transform(np.nan_to_num(X, nan=0.0))
    lr = LogisticRegression(max_iter=2000, C=1.0).fit(Xz[ok], y[ok])
    e = lr.predict_proba(Xz)[:, 1]
    w = np.zeros(len(Xc))
    e_ctrl = e[len(Xs):]
    finite = np.all(np.isfinite(Xc), axis=1)
    w[finite] = e_ctrl[finite] / np.clip(1.0 - e_ctrl[finite], 1e-6, None)
    pos = w > 0
    if pos.any():
        w = np.clip(w, 0, np.nanpercentile(w[pos], 99.5))
        w = w / np.mean(w[w > 0])
    return w


def entropy_balance(ctrl_feat: pd.DataFrame, target: np.ndarray, sd: np.ndarray) -> np.ndarray:
    """Entropy-balancing weights matching all covariate means to `target`
    exactly (Hainmueller 2012), solved by Newton with backtracking line search
    on the convex dual L(lam)=logsumexp(-X lam). Returns weights over ALL rows
    (0 where covariates are non-finite)."""
    C = ctrl_feat.to_numpy(float)
    finite = np.all(np.isfinite(C), axis=1)
    X = (C[finite] - target) / sd
    lam = np.zeros(X.shape[1])

    def loss(l):
        return logsumexp(-X @ l)

    cur = loss(lam)
    for _ in range(500):
        v = -X @ lam
        w = np.exp(v - logsumexp(v))
        g = -(w @ X)               # gradient of L
        if np.max(np.abs(w @ X)) < 1e-10:
            break
        H = (X * w[:, None]).T @ X - np.outer(w @ X, w @ X)
        d = -np.linalg.solve(H + 1e-10 * np.eye(X.shape[1]), g)
        t = 1.0
        while t > 1e-8 and loss(lam + t * d) > cur - 1e-4 * t * (g @ d):
            t *= 0.5
        lam = lam + t * d
        cur = loss(lam)
    v = -X @ lam
    wfin = np.exp(v - logsumexp(v))
    out = np.zeros(len(C))
    out[finite] = wfin / wfin.mean()   # mean-1 normalization
    return out


def main() -> int:
    master = Table.read(MASTER).to_pandas()
    master["tier"] = decode(master["tier"])
    slow = master[master["tier"].isin(["A", "B", "C"])].copy()
    slow_feat = feat_frame(slow)
    target = np.nanmean(slow_feat.to_numpy(float), axis=0)

    hw = pd.read_csv(WEIGHTS)
    hw["band"] = hw["band"].astype(str).str.strip()
    ctrl_orb = Table.read(CONTROLS).to_pandas()
    ctrl_orb["band"] = decode(ctrl_orb["band"])

    bal_rows = []
    summ = {"n_slow_tierABC": int(len(slow)), "slow_median_Rperi_pc": SLOW_RPERI_PC, "bands": {}}
    for key, label, csv in BANDS:
        cov = pd.read_csv(DATADIR / csv).replace([np.inf, -np.inf], np.nan)
        orb = (ctrl_orb[ctrl_orb["band"] == key][["source_id", "R_peri_kpc", "ecc"]]
               .rename(columns={"R_peri_kpc": "R_peri_orb_kpc", "ecc": "ecc_orb"}))
        band = orb.merge(cov, on="source_id", how="left")
        band = band.merge(hw[hw["band"] == key][["source_id", "expanded_kw_weight"]],
                          on="source_id", how="left")
        cf = feat_frame(band)
        sd = np.sqrt((np.nanvar(slow_feat.to_numpy(float), axis=0, ddof=1)
                      + np.nanvar(cf.to_numpy(float), axis=0, ddof=1)) / 2.0)

        w = {
            "raw": np.ones(len(band)),
            "kernel": np.where(np.isfinite(band["expanded_kw_weight"]), band["expanded_kw_weight"], 0.0).astype(float),
            "ipw": ipw_att_weights(slow_feat, cf),
            "ebal": entropy_balance(cf, target, sd),
        }

        for ckey, clabel in COVS:
            row = {"band": key, "covariate": ckey, "label": clabel}
            for scheme in w:
                row[f"smd_{scheme}"] = smd(slow_feat[ckey].to_numpy(), cf[ckey].to_numpy(), w[scheme])
            bal_rows.append(row)

        rperi = band["R_peri_orb_kpc"].to_numpy(float) * 1000.0
        ecc = band["ecc_orb"].to_numpy(float)
        bsum = {"label": label, "n_control": int(len(band))}
        for scheme, wv in w.items():
            bsum[f"ess_{scheme}"] = ess(wv)
            bsum[f"median_Rperi_pc_{scheme}"] = weighted_quantile(rperi, wv, 0.5)
            bsum[f"median_ecc_{scheme}"] = weighted_quantile(ecc, wv, 0.5)
            band_rows = [r for r in bal_rows if r["band"] == key]
            bsum[f"max_abs_smd_{scheme}"] = float(np.nanmax(np.abs([r[f"smd_{scheme}"] for r in band_rows])))
        summ["bands"][key] = bsum

    bal = pd.DataFrame(bal_rows)
    OUT.mkdir(parents=True, exist_ok=True)
    bal.to_csv(OUT / "expanded_control_balance.csv", index=False)

    for scheme in ["raw", "kernel", "ipw", "ebal"]:
        rseq = [SLOW_RPERI_PC] + [summ["bands"][b[0]][f"median_Rperi_pc_{scheme}"] for b in BANDS]
        eseq = [summ["bands"][BANDS[0][0]][f"median_ecc_{scheme}"]]  # placeholder filled below
        summ[f"rperi_sequence_{scheme}_pc"] = [round(v, 1) for v in rseq]
        summ[f"rperi_monotonic_increasing_{scheme}"] = bool(np.all(np.diff(rseq) > 0))
        eseq = [summ["bands"][b[0]][f"median_ecc_{scheme}"] for b in BANDS]
        summ[f"ecc_sequence_{scheme}"] = [round(v, 3) for v in eseq]
        summ[f"ecc_monotonic_decreasing_{scheme}"] = bool(np.all(np.diff(eseq) < 0))
    (OUT / "expanded_control_balance_summary.json").write_text(json.dumps(summ, indent=2))

    _love_plot(bal)
    _latex_table(bal, summ)

    # Mirror the manuscript-facing products into both the working tree
    # (release/, used by release/main.tex) and the bundle snapshot.
    import shutil
    rel_fig = REPO / "release" / "figures"
    bun_tab = BUNDLE / "tables" / "v15"
    rel_fig.mkdir(parents=True, exist_ok=True)
    bun_tab.mkdir(parents=True, exist_ok=True)
    shutil.copy(FIG / "fig_control_balance_love.pdf", rel_fig / "fig_control_balance_love.pdf")
    shutil.copy(TAB / "tab_control_balance.tex", bun_tab / "tab_control_balance.tex")

    print(json.dumps(summ, indent=2))
    return 0


def _love_plot(bal: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    labels = [c[1] for c in COVS]
    y = np.arange(len(COV_KEYS))[::-1]
    # mean over the four bands of |SMD| for each scheme
    def mean_abs(scheme):
        return np.array([np.nanmean(np.abs(bal[bal.covariate == k][f"smd_{scheme}"])) for k in COV_KEYS])
    fig, ax = plt.subplots(figsize=(5.0, 3.7))
    ax.scatter(mean_abs("raw"), y, facecolors="none", edgecolors="#555555", s=34, lw=1.0, label="unweighted")
    ax.scatter(mean_abs("kernel"), y, marker="s", color="#ff7f0e", s=24, label="kernel sensitivity")
    ax.scatter(mean_abs("ebal"), y, marker="o", color="#2ca02c", s=26, label="entropy-balanced")
    ax.axvline(0.1, color="grey", ls="--", lw=0.8)
    ax.axvline(0.25, color="grey", ls=":", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel(r"mean $|\mathrm{standardized\ mean\ difference}|$ over the four control bands")
    ax.set_xlim(left=-0.005)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "fig_control_balance_love.pdf", bbox_inches="tight", dpi=250)
    plt.close(fig)


def _latex_table(bal: pd.DataFrame, summ: dict) -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    bands = [b[0] for b in BANDS]
    lines = [
        r"% Auto-generated by scripts/phase14ac_control_balance.py",
        r"\begin{deluxetable*}{lcccccccc}",
        r"\tablecaption{Covariate balance of the velocity-stratified matched controls"
        r" relative to the Tier~A+B+C slow catalogue, expressed as standardized mean"
        r" differences (SMD). For each band the columns give the SMD before weighting"
        r" (unweighted) and after entropy balancing; $|\mathrm{SMD}|<0.1$ indicates good"
        r" balance. The paper's product-of-1D observability kernel does not reduce these"
        r" SMDs (Fig.~\ref{fig:control_balance}), whereas entropy balancing matches every"
        r" covariate mean to $|\mathrm{SMD}|\lesssim10^{-3}$. The footer gives the"
        r" entropy-balanced effective sample size $N_{\rm eff}$ and weighted median"
        r" pericentre, which preserves its monotonic increase with $\vgrf$.\label{tab:control_balance}}",
        r"\tablehead{\colhead{} & \multicolumn{2}{c}{25--50} & \multicolumn{2}{c}{50--100}"
        r" & \multicolumn{2}{c}{100--200} & \multicolumn{2}{c}{200--260} \\"
        r" \colhead{covariate} & \colhead{before} & \colhead{ebal.} & \colhead{before} & \colhead{ebal.}"
        r" & \colhead{before} & \colhead{ebal.} & \colhead{before} & \colhead{ebal.}}",
        r"\startdata",
    ]
    for ckey, clabel in COVS:
        cells = []
        for band in bands:
            r = bal[(bal["band"] == band) & (bal["covariate"] == ckey)].iloc[0]
            cells.append(f"{r['smd_raw']:+.2f}")
            cells.append(f"{abs(r['smd_ebal']):.3f}")
        lines.append(f"{clabel} & " + " & ".join(cells) + r" \\")
    lines.append(r"\hline")
    ess_cells = " & ".join(
        [f"\\multicolumn{{2}}{{c}}{{{summ['bands'][b]['ess_ebal']:.0f}}}" for b in bands])
    lines.append(r"$N_{\rm eff}$ (ebal.) & " + ess_cells + r" \\")
    rp_cells = " & ".join(
        [f"\\multicolumn{{2}}{{c}}{{{summ['bands'][b]['median_Rperi_pc_ebal']:.0f}}}" for b in bands])
    lines.append(r"med. $R_{\rm peri}$ (pc) & " + rp_cells + r" \\")
    lines += [r"\enddata", r"\end{deluxetable*}"]
    (TAB / "tab_control_balance.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
