"""Phase 5D — APOGEE + GALAH chemistry recompute under v2 inputs.

The cross-match files release/data/crossmatch_apogee.csv (77 stars) and
release/data/crossmatch_galah.csv (66 stars) are pre-cross-matched
against the input candidate list. The chemical abundances themselves
are not affected by ZP correction (those come from the spectroscopy);
what changes in v2 is *which* stars are in the slow-Vgrf sample and
which tier each falls into.

This script:
1. Cross-matches the chemistry tables to catalogue_v2 (all 2,859
   candidates carrying the release tier).
2. Reports tier-resolved chemistry distributions:
     [Fe/H], [a/Fe] / [Mg/Fe], and the (Mg, Al) "Splash plane".
3. Outputs a per-star chemistry table for downstream Phase 5C use.

No new VizieR queries — uses the cached on-disk crossmatches.

Outputs
-------
release/v2/phase5/chemistry_v2.fits
release/v2/phase5/gate5D_chemistry.json
release/v2/phase5/gate5D_alpha_FeH.png
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
OUT.mkdir(parents=True, exist_ok=True)


def log(m): print(f"[5D t={time.time()-T0:5.1f}s] {m}", flush=True)


def main():
    global T0; T0 = time.time()

    log("loading APOGEE + GALAH crossmatches and v2 catalogue...")
    apo = pd.read_csv(REPO / "release/data/crossmatch_apogee.csv")
    gal = pd.read_csv(REPO / "release/data/crossmatch_galah.csv")
    v2 = Table.read(REPO / "release/v2/phase1/catalogue_v2.fits").to_pandas()
    if v2["tier"].dtype == object and isinstance(v2["tier"].iloc[0], (bytes, bytearray)):
        v2["tier"] = v2["tier"].str.decode("utf-8")

    # --- Standardise APOGEE columns ---
    apo_clean = apo.rename(columns={
        "__Fe_H_":   "FeH",  "e__Fe_H_":  "e_FeH",
        "__a_M_":    "alphaM", "e__a_M_": "e_alphaM",
        "__Mg_Fe_":  "MgFe", "e__Mg_Fe_": "e_MgFe",
        "__Al_Fe_":  "AlFe", "e__Al_Fe_": "e_AlFe",
        "__Mn_Fe_":  "MnFe", "e__Mn_Fe_": "e_MnFe",
    })[["source_id", "Teff", "logg", "FeH", "e_FeH", "alphaM", "MgFe",
        "AlFe", "MnFe"]]
    apo_clean["survey"] = "APOGEE"

    gal_clean = gal.rename(columns={
        "__Fe_H_":     "FeH",   "e__Fe_H_":     "e_FeH",
        "__alpha_Fe_": "alphaFe", "e__alpha_Fe_": "e_alphaFe",
        "__Mg_Fe_":    "MgFe",  "e__Mg_Fe_":    "e_MgFe",
        "__Al_Fe_":    "AlFe",  "e__Al_Fe_":    "e_AlFe",
        "__Mn_Fe_":    "MnFe",  "e__Mn_Fe_":    "e_MnFe",
    })[["source_id", "Teff", "logg", "FeH", "e_FeH", "alphaFe", "MgFe",
        "AlFe", "MnFe"]]
    gal_clean["survey"] = "GALAH"

    # --- Merge with v2 tiers ---
    keep = ["source_id", "tier", "P_vgrf_below_25", "vgrf_default", "rv_quality"]
    v2_min = v2[keep].copy()
    apo_v2 = apo_clean.merge(v2_min, on="source_id", how="inner")
    gal_v2 = gal_clean.merge(v2_min, on="source_id", how="inner")

    log(f"APOGEE matched to v2: {len(apo_v2)} / {len(apo_clean)} ; "
        f"GALAH matched: {len(gal_v2)} / {len(gal_clean)}")

    # Standardize an "alpha-like" column: APOGEE has alphaM, GALAH has alphaFe.
    # Approximate alphaFe ~ alphaM (decent for most stars; documented limitation)
    apo_v2 = apo_v2.assign(alpha_proxy=apo_v2["alphaM"])
    gal_v2 = gal_v2.assign(alpha_proxy=gal_v2["alphaFe"])
    keep_cols = ["source_id", "tier", "P_vgrf_below_25", "vgrf_default",
                 "rv_quality", "Teff", "logg", "FeH", "e_FeH",
                 "MgFe", "AlFe", "MnFe", "alpha_proxy", "survey"]
    chem = pd.concat([apo_v2[keep_cols], gal_v2[keep_cols]],
                     ignore_index=True)
    Table.from_pandas(chem).write(OUT / "chemistry_v2.fits", overwrite=True)
    log(f"wrote chemistry_v2.fits ({len(chem)} rows)")

    # --- Tier-resolved distributions ---
    def stat(arr):
        a = np.asarray(arr, dtype=float)
        a = a[np.isfinite(a)]
        if not len(a): return None
        return {"n": int(len(a)),
                "p16": float(np.percentile(a, 16)),
                "p50": float(np.percentile(a, 50)),
                "p84": float(np.percentile(a, 84))}

    summary = {"by_tier": {}, "by_survey": {}}
    for tier_set, label in [(["A", "B"], "tier_AB"),
                             (["A", "B", "C"], "tier_ABC"),
                             (["A", "B", "C", "D"], "tier_ABCD"),
                             (["X"], "tier_X_excluded")]:
        m = chem["tier"].isin(tier_set)
        s = chem.loc[m]
        summary["by_tier"][label] = {
            "n": int(m.sum()),
            "FeH":         stat(s["FeH"]),
            "alpha_proxy": stat(s["alpha_proxy"]),
            "MgFe":        stat(s["MgFe"]),
            "AlFe":        stat(s["AlFe"]),
            "MnFe":        stat(s["MnFe"]),
        }
    for surv in ("APOGEE", "GALAH"):
        m = chem["survey"] == surv
        summary["by_survey"][surv] = {
            "n": int(m.sum()),
            "FeH":         stat(chem.loc[m, "FeH"]),
            "alpha_proxy": stat(chem.loc[m, "alpha_proxy"]),
        }

    (OUT / "gate5D_chemistry.json").write_text(json.dumps(summary, indent=2))
    log("wrote gate5D_chemistry.json")
    print(json.dumps(summary, indent=2))

    # --- Plots: alpha vs [Fe/H], coloured by tier ---
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5))
    color_map = {"A": "C3", "B": "C1", "C": "C2", "D": "C0", "X": "0.6"}
    for tier in ("X", "D", "C", "B", "A"):
        m = chem["tier"] == tier
        if m.sum() == 0: continue
        ax.scatter(chem.loc[m, "FeH"], chem.loc[m, "alpha_proxy"],
                   color=color_map[tier], s=20, alpha=0.8,
                   label=f"Tier {tier} (n={int(m.sum())})", edgecolors="white")
    ax.set_xlabel("[Fe/H]")
    ax.set_ylabel(r"[$\alpha$/Fe] (proxy: APOGEE [a/M] or GALAH [$\alpha$/Fe])")
    ax.set_title("Phase 5D — chemistry vs v2 tier (APOGEE + GALAH crossmatch)")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "gate5D_alpha_FeH.png", dpi=140); plt.close(fig)
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
