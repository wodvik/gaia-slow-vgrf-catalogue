"""Phase 5C — Splash / Aurora / GSE population decomposition.

Uses the chemistry catalogue from 5D (APOGEE+GALAH crossmatch, n=141)
and the XP/GSP-Phot metallicity from 5E (n=1963) to classify v2 Tier
A+B+C stars into three populations:

  - **Aurora** (in-situ metal-poor proto-disk):
      [Fe/H] < -1.0  AND  [alpha/Fe] >= +0.25
  - **Splash** (heated thick disk):
      -1.0 <= [Fe/H] < +0.0  AND  [alpha/Fe] >= +0.15
  - **GSE / Sausage** (accreted, low-alpha):
      [Fe/H] < -0.7  AND  [alpha/Fe] < +0.20
  - **Disk** (contamination):
      [Fe/H] >= 0  OR  ([Fe/H] >= -0.3 AND [alpha/Fe] < +0.15)
  - **unclassified**: anything else

Cuts are standard literature definitions (Belokurov+2020, Bonaca+2020,
Conroy+2022 for Aurora; Bonaca+2017 for Splash; Helmi+2018 for GSE).

For stars with chemistry only from GSP-Phot (no [alpha/Fe]), we use
[M/H] alone to flag low-MH (Aurora-or-GSE bucket).

Outputs
-------
release/v2/phase5/populations_v2.fits         (per-star labels)
release/v2/phase5/gate5C_populations.json     (tier-resolved fractions)
release/v2/phase5/gate5C_chem_pop.png         (alpha-FeH plane with
                                               population colouring)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase5"


def log(m): print(f"[5C t={time.time()-T0:5.1f}s] {m}", flush=True)


def classify(feh, alpha):
    """Population label from ([Fe/H], [alpha/Fe]). NaN-safe."""
    if not np.isfinite(feh):
        return "unclassified"
    if feh >= 0.0:
        return "disk"
    if not np.isfinite(alpha):
        # GSP-Phot only: bucket by [Fe/H]
        if feh < -1.0:
            return "low_MH_no_alpha"  # Aurora-or-GSE
        else:
            return "intermediate_MH_no_alpha"  # Splash-or-disk
    # Both [Fe/H] and [alpha/Fe] available
    if feh < -1.0 and alpha >= 0.25:
        return "Aurora"
    if feh < -0.7 and alpha < 0.20:
        return "GSE"
    if -1.0 <= feh < 0.0 and alpha >= 0.15:
        return "Splash"
    if feh >= -0.3 and alpha < 0.15:
        return "disk"
    return "unclassified"


def main():
    global T0; T0 = time.time()

    log("loading 5D chemistry + 5E XP MH...")
    chem = Table.read(OUT / "chemistry_v2.fits").to_pandas()
    xp   = Table.read(OUT / "xp_metallicity_v2.fits").to_pandas()
    if chem["tier"].dtype == object and isinstance(chem["tier"].iloc[0], (bytes, bytearray)):
        chem["tier"] = chem["tier"].str.decode("utf-8")
        chem["survey"] = chem["survey"].str.decode("utf-8")
    if xp["tier"].dtype == object and isinstance(xp["tier"].iloc[0], (bytes, bytearray)):
        xp["tier"] = xp["tier"].str.decode("utf-8")

    # Build base table: one row per star, with best available [Fe/H] and
    # [alpha/Fe]. A few sources have both APOGEE and GALAH rows in the
    # chemistry table; use a deterministic one-row-per-source view for
    # population bookkeeping so tier counts cannot be inflated by joins.
    chem["survey_rank"] = chem["survey"].map({"APOGEE": 0, "GALAH": 1}).fillna(9)
    chem["finite_rank"] = ~(
        np.isfinite(chem["FeH"].astype(float))
        & np.isfinite(chem["alpha_proxy"].astype(float))
    )
    chem_one = (chem.sort_values(["source_id", "finite_rank", "survey_rank"])
                    .drop_duplicates("source_id", keep="first"))

    base = xp[["source_id", "tier", "P_vgrf_below_25",
               "vgrf_default", "MH_xp"]].copy()
    base = base.merge(
        chem_one[["source_id", "FeH", "MgFe", "alpha_proxy", "survey"]],
        on="source_id", how="left", suffixes=("", "_chem")
    )
    # Coalesce: prefer APOGEE/GALAH FeH if present, else GSP-Phot MH
    base["best_FeH"]   = base["FeH"].where(base["FeH"].notna(), base["MH_xp"])
    base["best_alpha"] = base["alpha_proxy"]
    base["chem_source"] = np.where(
        base["FeH"].notna(),
        np.where(base["survey"].notna(), base["survey"], "spec_unknown"),
        "GSP-Phot",
    )

    base["population"] = base.apply(
        lambda r: classify(r["best_FeH"], r["best_alpha"]), axis=1)
    Table.from_pandas(base).write(OUT / "populations_v2.fits", overwrite=True)
    log(f"wrote populations_v2.fits ({len(base)} rows)")

    # Tier-resolved fractions
    summary = {"by_tier": {}, "global": {}}
    pops = ["Aurora", "Splash", "GSE", "low_MH_no_alpha",
            "intermediate_MH_no_alpha", "disk", "unclassified"]
    for tier_set, lab in [(["A", "B"], "tier_AB"),
                          (["A", "B", "C"], "tier_ABC"),
                          (["X"], "tier_X")]:
        m = base["tier"].isin(tier_set)
        n = int(m.sum())
        if n == 0: continue
        counts = {p: int(((base["population"] == p) & m).sum()) for p in pops}
        fracs = {p: round(counts[p] / n, 4) for p in pops}
        summary["by_tier"][lab] = {"n": n, "counts": counts, "fracs": fracs}

    # Headline interpretation: one-row-per-source alpha-chemistry subset
    # within Tier A+B+C. These counts are the only population-region counts
    # used in the manuscript.
    AB_C = summary["by_tier"]["tier_ABC"]
    alpha_mask = (base["tier"].isin(["A", "B", "C"]) &
                  base["FeH"].notna() &
                  base["alpha_proxy"].notna())
    alpha_counts = {
        p: int(((base["population"] == p) & alpha_mask).sum())
        for p in ["Aurora", "Splash", "GSE", "disk", "unclassified"]
    }
    n_alpha = int(alpha_mask.sum())
    summary["headline"] = {
        "n_TierABC": AB_C["n"],
        "n_with_alpha_chem": n_alpha,
        "alpha_subsample_counts": alpha_counts,
        "alpha_subsample_fracs": {
            p: (alpha_counts[p] / max(n_alpha, 1))
            for p in alpha_counts
        },
        "n_GSPPhot_only_no_alpha": (AB_C["counts"]["low_MH_no_alpha"] +
                                     AB_C["counts"]["intermediate_MH_no_alpha"]),
        "Aurora_fraction_of_alphasubsample": alpha_counts["Aurora"] / max(n_alpha, 1),
        "Splash_fraction_of_alphasubsample": alpha_counts["Splash"] / max(n_alpha, 1),
        "GSE_fraction_of_alphasubsample": alpha_counts["GSE"] / max(n_alpha, 1),
    }

    (OUT / "gate5C_populations.json").write_text(json.dumps(summary, indent=2))
    log("wrote gate5C_populations.json")
    print(json.dumps(summary["headline"], indent=2))
    print(json.dumps(summary["by_tier"], indent=2))

    # Plot: alpha-FeH plane coloured by population (alpha-subsample only)
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5))
    pal = {"Aurora": "darkviolet", "Splash": "C1", "GSE": "C2",
           "disk": "C0", "unclassified": "0.6",
           "low_MH_no_alpha": "C3", "intermediate_MH_no_alpha": "C4"}
    tier_ABC = base["tier"].isin(["A", "B", "C"])
    for pop, color in pal.items():
        m = tier_ABC & (base["population"] == pop) & base["alpha_proxy"].notna()
        if m.sum() == 0: continue
        ax.scatter(base.loc[m, "best_FeH"], base.loc[m, "best_alpha"],
                   color=color, s=22, alpha=0.85,
                   label=f"{pop} (n={int(m.sum())})", edgecolors="white")
    ax.axhline(0.15, ls=":", color="black", alpha=0.5)
    ax.axhline(0.25, ls=":", color="black", alpha=0.5)
    ax.axvline(-1.0, ls=":", color="black", alpha=0.5)
    ax.axvline(-0.7, ls=":", color="black", alpha=0.5)
    ax.set_xlabel("[Fe/H]"); ax.set_ylabel(r"[$\alpha$/Fe]")
    ax.set_title("Phase 5C — chemistry plane decomposition (Tier A+B+C, alpha-subsample)")
    ax.legend(); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(OUT / "gate5C_chem_pop.png", dpi=140); plt.close(fig)
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
