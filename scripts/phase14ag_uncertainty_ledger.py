"""Phase 14AG -- uncertainty ledger for primary numerical claims.

Referee response (deep-review Issue 11). Consolidates the primary numbers
and attaches an uncertainty interval to each, computing bootstrap confidence
intervals for the key medians (which previously carried none) and Wilson
intervals for the pipeline fractions. The excess-factor intervals are taken
from phase14ad_null_models.py.

Outputs
-------
phase14/expanded_uncertainty_ledger.json
tables/v15/tab_uncertainty_ledger.tex
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
ORBITS = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.fits"
MASTER = BUNDLE / "catalogues" / "catalogue_expanded_master.fits"
CONTROLS = REPO / "release" / "_iterations" / "v2" / "phase5" / "control_orbits.fits"
WEIGHTS = BUNDLE / "phase14" / "expanded_control_weights.csv"
NULLS = BUNDLE / "phase14" / "expanded_null_models_summary.json"
OUT = BUNDLE / "phase14"
TAB_DIRS = [REPO / "release" / "tables" / "v15", BUNDLE / "tables" / "v15"]

B = 4000
SEED = 20260605
# Injection-recovery counts from the GeDR3mock test (Sec. 4.7 / tab_recovery_summary).
REC_SLOW_NUM, REC_SLOW_DEN = 6086, 6944          # conditional recovery 87.6%
REC_E2E_NUM, REC_E2E_DEN = 6086, 10000           # end-to-end 60.9%
LEAK_NUM, LEAK_DEN = 282, 10000                  # 25-50 leakage 2.82%


def boot_median_ci(x, rng, weights=None):
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    x = x[ok]
    if weights is not None:
        w = np.asarray(weights, float)[ok]
    meds = np.empty(B)
    n = len(x)
    for i in range(B):
        idx = rng.integers(0, n, n)
        if weights is None:
            meds[i] = np.median(x[idx])
        else:
            xs, ws = x[idx], w[idx]
            o = np.argsort(xs)
            cdf = np.cumsum(ws[o]) / np.sum(ws[o])
            meds[i] = np.interp(0.5, cdf, xs[o])
    return float(np.percentile(meds, 16)), float(np.percentile(meds, 84))


def wilson(k, n, z=1.0):
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return float((c - h) / d), float((c + h) / d)


def main() -> int:
    rng = np.random.default_rng(SEED)
    orb = Table.read(ORBITS).to_pandas()
    ecc = orb["static_ecc"].to_numpy(float)
    rperi = orb["static_R_peri_kpc"].to_numpy(float) * 1000.0
    rapo = orb["static_R_apo_kpc"].to_numpy(float)

    ctrl = Table.read(CONTROLS).to_pandas()
    ctrl["band"] = ctrl["band"].str.decode("utf-8").str.strip() if ctrl["band"].dtype == object and isinstance(ctrl["band"].iloc[0], bytes) else ctrl["band"].astype(str).str.strip()
    w = pd.read_csv(WEIGHTS)
    w["band"] = w["band"].astype(str).str.strip()
    c25 = ctrl[ctrl["band"] == "vgrf_25_50"].merge(
        w[w["band"] == "vgrf_25_50"][["source_id", "expanded_kw_weight"]], on="source_id", how="left")
    c25_rp = c25["R_peri_kpc"].to_numpy(float) * 1000.0
    c25_w = c25["expanded_kw_weight"].to_numpy(float)

    nulls = json.loads(NULLS.read_text())["nulls"]

    master = Table.read(MASTER).to_pandas()
    tier = master["tier"].str.decode("utf-8").str.strip() if master["tier"].dtype == object and isinstance(master["tier"].iloc[0], bytes) else master["tier"].astype(str).str.strip()
    prob = master["P_vgrf_below_25"].to_numpy(float)

    rows = []

    def add(name, value, lo, hi, method):
        rows.append({"quantity": name, "value": value, "lo": lo, "hi": hi, "method": method})

    def wmed(x, w):
        ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
        xs, ws = x[ok], w[ok]
        o = np.argsort(xs)
        return float(np.interp(0.5, np.cumsum(ws[o]) / np.sum(ws[o]), xs[o]))

    e_lo, e_hi = boot_median_ci(ecc, rng)
    add(r"Median eccentricity (Tier A+B+C, point)", f"{np.nanmedian(ecc):.3f}", f"{e_lo:.3f}", f"{e_hi:.3f}", "bootstrap")
    rp_lo, rp_hi = boot_median_ci(rperi, rng)
    add(r"Median $R_{\rm peri}$ (pc, point)", f"{np.nanmedian(rperi):.1f}", f"{rp_lo:.0f}", f"{rp_hi:.0f}", "bootstrap")
    ra_lo, ra_hi = boot_median_ci(rapo, rng)
    add(r"Median $R_{\rm apo}$ (kpc, point)", f"{np.nanmedian(rapo):.2f}", f"{ra_lo:.2f}", f"{ra_hi:.2f}", "bootstrap")
    c_lo, c_hi = boot_median_ci(c25_rp, rng, weights=c25_w)
    add(r"Control 25--50 median $R_{\rm peri}$ (pc)", f"{wmed(c25_rp, c25_w):.1f}", f"{c_lo:.0f}", f"{c_hi:.0f}", "weighted bootstrap")

    def add_purity(name, mask):
        p = prob[mask]
        purity = float(np.mean(p))
        se = float(np.sqrt(np.sum(p * (1.0 - p))) / len(p))
        lo = max(0.0, purity - se)
        hi = min(1.0, purity + se)
        add(name, f"{100*purity:.0f}", f"{100*lo:.0f}", f"{100*hi:.0f}", "score Bernoulli SE")

    add_purity(r"Nominal purity (Tier A+B, \%)", tier.isin(["A", "B"]).to_numpy())
    add_purity(r"Nominal purity (Tier A+B+C, \%)", tier.isin(["A", "B", "C"]).to_numpy())

    add(r"Excess factor (BIC mixture)", f"{nulls['gmm']['excess']:.1f}",
        f"{nulls['gmm']['excess_ci'][0]:.1f}", f"{nulls['gmm']['excess_ci'][1]:.1f}", "bootstrap")
    add(r"Excess factor (single ellipsoid)", f"{nulls['gaussian']['excess']:.1f}",
        f"{nulls['gaussian']['excess_ci'][0]:.1f}", f"{nulls['gaussian']['excess_ci'][1]:.1f}", "bootstrap")

    for name, k, n, val in [
        (r"Pipeline recovery, conditional (\%)", REC_SLOW_NUM, REC_SLOW_DEN, "87.6"),
        (r"Pipeline recovery, end-to-end (\%)", REC_E2E_NUM, REC_E2E_DEN, "60.9"),
        (r"25--50 leakage fraction (\%)", LEAK_NUM, LEAK_DEN, "2.82"),
    ]:
        lo, hi = wilson(k, n)
        add(name, val, f"{lo*100:.1f}", f"{hi*100:.1f}", "Wilson")

    ledger = {"n_bootstrap": B, "rows": rows}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "expanded_uncertainty_ledger.json").write_text(json.dumps(ledger, indent=2))
    _latex(rows)
    print(json.dumps(ledger, indent=2))
    return 0


def _latex(rows) -> None:
    lines = [
        r"% Auto-generated by scripts/phase14ag_uncertainty_ledger.py",
        r"\begin{deluxetable*}{lccl}",
        r"\tablecaption{Uncertainty ledger for the primary numerical claims. Intervals are"
        r" central 68\% ranges: bootstrap resampling of the catalogue for the medians and"
        r" excess factors, Wilson score intervals for the binomial pipeline fractions,"
        r" and score-implied Bernoulli standard errors for nominal purity. Rows labelled"
        r" ``point'' use deterministic point-estimate orbit products; the"
        r" uncertainty-propagated eccentricity summary is Table~\ref{tab:orbit_summary}.\label{tab:uncertainty_ledger}}",
        r"\tablehead{\colhead{quantity} & \colhead{value} & \colhead{68\% interval}"
        r" & \colhead{method}}",
        r"\startdata",
    ]
    for r in rows:
        lines.append(f"{r['quantity']} & {r['value']} & {r['lo']}--{r['hi']} & {r['method']} \\\\")
    lines += [r"\enddata", r"\end{deluxetable*}"]
    text = "\n".join(lines) + "\n"
    for d in TAB_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        (d / "tab_uncertainty_ledger.tex").write_text(text)


if __name__ == "__main__":
    raise SystemExit(main())
