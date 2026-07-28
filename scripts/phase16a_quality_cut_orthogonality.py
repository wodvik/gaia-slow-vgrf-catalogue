"""Phase 16A -- quality-cut orthogonality and direction-of-bias audit.

A quality cut applied to a sample that is itself selected at a kinematic
extreme is only defensible if the cut axis is distinct from the selection
axis. Following the argument made by Banik et al. (2026, MNRAS, submitted;
arXiv:2607.00764) for age-extremum samples, we make that explicit here in
two ways:

  1. ORTHOGONALITY. For every observational quality cut we report the pass
     rate as a function of point-estimate Vgrf across the full candidate
     pool. A cut whose pass rate is flat in Vgrf cannot manufacture the
     slow sample, because it is blind to the quantity used to select it.

  2. DIRECTION OF BIAS. For every cut we report how the surviving counts
     and the median Vgrf move. Cuts that only ever remove stars, and that
     leave the surviving Vgrf distribution unchanged, cannot inflate the
     catalogue.

This pass also adds the Gaia image-parameter-determination cut
`ipd_frac_multi_peak <= 1` (Banik et al. 2026 Sec. 2.2, after Banik et al.
2024; Pace, Erkal & Li 2022 adopt <= 2), which targets marginally RESOLVED
companions and is therefore complementary to RUWE rather than redundant
with it.

Outputs: phase14/quality_cut_orthogonality_summary.json
         tables/v15/tab_cut_orthogonality.tex
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
INP = BUNDLE / "private_inputs" / "expanded_candidates_mc_tiered.csv"
OUT_JSON = BUNDLE / "phase14" / "quality_cut_orthogonality_summary.json"
TAB_DIRS = [BUNDLE / "tables" / "v15", REPO / "release" / "tables" / "v15"]

VGRF_COL = "vgrf_default_exact"
P_COL = "P_vgrf_below_25"
TIERS = {"A": 0.95, "A+B": 0.84, "A+B+C": 0.50}

# Vgrf bins spanning the region where the candidate pool is complete.
VGRF_BINS = [0.0, 10.0, 20.0, 25.0, 30.0, 40.0, 50.0]


def rv_quality_ok(d: pd.DataFrame) -> pd.Series:
    """The released rvs_quality_ok flag, joined from the catalogue.

    Deriving this from the raw rv_chisq_pvalue / rv_expected_sig_to_noise
    columns only approximates the released screen (it omits the marginal
    classes) and disagreed with the Gold subset by one star. Use the shipped
    flag so this table and the Gold numbers in the text cannot drift apart.
    """
    from astropy.table import Table as _T
    m = _T.read(BUNDLE / "catalogues" / "catalogue_retier_master.fits")
    flag = dict(zip(np.asarray(m["source_id"]).astype(np.int64).tolist(),
                    np.asarray(m["rvs_quality_ok"]).astype(bool).tolist()))
    return pd.Series([flag.get(int(s), False) for s in d["source_id"].to_numpy()],
                     index=d.index)


def build_cuts(d: pd.DataFrame) -> dict[str, pd.Series]:
    sigma_d = (d["dist_hi_pc_final_screen"] - d["dist_lo_pc_final_screen"]) / 2.0
    frac_d = sigma_d / d["dist_pc_final_screen"]
    return {
        "RUWE $<1.4$": d["ruwe"] < 1.4,
        r"\texttt{ipd\_frac\_multi\_peak} $\leq 1$": d["ipd_frac_multi_peak"] <= 1,
        "RVS quality ok": rv_quality_ok(d),
        r"$\varpi/\sigma_\varpi>10$": d["parallax_over_error"] > 10,
        r"$\sigma_d/d<0.15$": frac_d < 0.15,
    }


def wilson(k: int, n: int, z: float = 1.0) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return ctr - half, ctr + half


def main() -> None:
    d = pd.read_csv(INP, low_memory=False)
    # Survival is reported on the ADOPTED (population-prior) tiers; the
    # orthogonality statistic itself is pool-wide and tier-independent.
    lat = pd.read_csv(BUNDLE / "phase14" / "latent_deconvolution"
                      / "latent_vgrf_per_star_regularised.csv")
    pmap = dict(zip(lat.source_id.to_numpy(), lat.P_latent.to_numpy()))
    d[P_COL] = [pmap.get(int(s), 0.0) for s in d["source_id"].to_numpy()]
    cuts = build_cuts(d)
    v = d[VGRF_COL]

    summary: dict = {
        "n_candidate_pool": int(len(d)),
        "vgrf_bins_kms": VGRF_BINS,
        "cuts": {},
    }

    for name, mask in cuts.items():
        mask = mask.fillna(False)
        entry: dict = {"pool_pass_rate": float(mask.mean())}

        # 1. Orthogonality: pass rate vs Vgrf across the whole pool.
        rates = []
        for lo, hi in zip(VGRF_BINS[:-1], VGRF_BINS[1:]):
            inb = (v >= lo) & (v < hi)
            n = int(inb.sum())
            k = int((inb & mask).sum())
            lo_ci, hi_ci = wilson(k, n)
            rates.append({"lo": lo, "hi": hi, "n": n, "pass": k,
                          "rate": k / n if n else float("nan"),
                          "ci68": [lo_ci, hi_ci]})
        entry["pass_rate_vs_vgrf"] = rates
        finite = [r["rate"] for r in rates if r["n"] > 0]
        entry["pass_rate_spread"] = float(max(finite) - min(finite)) if finite else None

        # 2. Direction of bias: per-tier survival and median-Vgrf shift.
        per_tier = {}
        for tname, thr in TIERS.items():
            sel = d[P_COL] > thr
            n = int(sel.sum())
            k = int((sel & mask).sum())
            lo_ci, hi_ci = wilson(k, n)
            med_before = float(np.nanmedian(v[sel]))
            med_after = float(np.nanmedian(v[sel & mask]))
            per_tier[tname] = {
                "n_before": n, "n_after": k,
                "survival": k / n if n else float("nan"),
                "survival_ci68": [lo_ci, hi_ci],
                "median_vgrf_before": med_before,
                "median_vgrf_after": med_after,
                "median_vgrf_shift": med_after - med_before,
            }
        entry["per_tier"] = per_tier
        summary["cuts"][name] = entry

    # Joint application of all cuts (the Gold-style stack, plus the new IPD cut).
    joint = np.logical_and.reduce([m.fillna(False).to_numpy() for m in cuts.values()])
    joint_tier = {}
    for tname, thr in TIERS.items():
        sel = (d[P_COL] > thr).to_numpy()
        n = int(sel.sum())
        k = int((sel & joint).sum())
        joint_tier[tname] = {"n_before": n, "n_after": k,
                             "survival": k / n if n else float("nan"),
                             "median_vgrf_before": float(np.nanmedian(v[sel])),
                             "median_vgrf_after": float(np.nanmedian(v[sel & joint]))}
    summary["all_cuts_joint"] = joint_tier

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {OUT_JSON}")

    # LaTeX table. table* : the cut labels make this too wide for one A&A column.
    lines = [
        r"\begin{table*}",
        r"\caption{Observational quality cuts are orthogonal to the \vgrf{}",
        r"  selection axis.  ``Pass-rate spread'' is the largest difference in",
        r"  pass rate between any two \vgrf{} bins spanning $0$--$50\kms$ across",
        r"  the full 20{,}829-source candidate pool; a small spread means the cut",
        r"  is effectively blind to \vgrf{} and therefore cannot manufacture the",
        r"  slow sample.  Survival columns give the fraction of each tier retained.}",
        r"\label{tab:cut_orthogonality}",
        r"\centering",
        r"\begin{tabular}{lccc}",
        r"\hline\hline",
        r"Cut & Pass-rate spread & \multicolumn{2}{c}{Survival} \\",
        r" & ($0$--$50\kms$) & Tier A+B & Tier A+B+C \\",
        r"\hline",
    ]
    for name, entry in summary["cuts"].items():
        ab = entry["per_tier"]["A+B"]
        abc = entry["per_tier"]["A+B+C"]
        lines.append(
            f"{name} & {entry['pass_rate_spread']:.3f} & "
            f"{ab['n_after']}/{ab['n_before']} ({100*ab['survival']:.1f}\\%) & "
            f"{abc['n_after']}/{abc['n_before']} ({100*abc['survival']:.1f}\\%) \\\\"
        )
    jab, jabc = joint_tier["A+B"], joint_tier["A+B+C"]
    lines += [
        r"\hline",
        f"All cuts jointly & --- & {jab['n_after']}/{jab['n_before']} "
        f"({100*jab['survival']:.1f}\\%) & {jabc['n_after']}/{jabc['n_before']} "
        f"({100*jabc['survival']:.1f}\\%) \\\\",
        r"\hline",
        r"\end{tabular}",
        r"\end{table*}",
    ]
    text = "\n".join(lines) + "\n"
    for td in TAB_DIRS:
        if td.parent.exists():
            td.mkdir(parents=True, exist_ok=True)
            (td / "tab_cut_orthogonality.tex").write_text(text, encoding="utf-8")
            print(f"wrote {td / 'tab_cut_orthogonality.tex'}")

    # Console report.
    for name, entry in summary["cuts"].items():
        print(f"\n{name}: pool pass {100*entry['pool_pass_rate']:.1f}%, "
              f"spread {entry['pass_rate_spread']:.3f}")
        for r in entry["pass_rate_vs_vgrf"]:
            if r["n"]:
                print(f"   {r['lo']:5.1f}-{r['hi']:5.1f}: {r['rate']:.3f}  (n={r['n']})")
        for tname, t in entry["per_tier"].items():
            print(f"   {tname:6s} {t['n_after']}/{t['n_before']} "
                  f"({100*t['survival']:.1f}%)  med Vgrf {t['median_vgrf_before']:.2f}"
                  f" -> {t['median_vgrf_after']:.2f} ({t['median_vgrf_shift']:+.2f})")
    print("\nJOINT:", json.dumps(joint_tier, indent=2))


if __name__ == "__main__":
    main()
