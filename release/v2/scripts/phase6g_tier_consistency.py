"""Phase 6G — tier-consistency master table.

Aggregate every Phase 1-6 quantity at Tier A / B / C separately.
Flag any quantity whose Tier-A value differs from Tier-C by more than
the within-tier MC uncertainty (measurement-driven vs population-driven).

Inputs:
  release/v2/phase1/catalogue_v2.fits           (tiers, P, vgrf)
  release/v2/phase4/catalogue_v2_orbits.fits    (point-est orbits)
  release/v2/phase5/populations_v2.fits         (chemistry/pops)
  release/v2/phase5/xp_metallicity_v2.fits      (GSP-Phot MH)
  release/v2/phase6/catalogue_mc_orbits.fits    (MC percentiles)

Output:
  release/v2/phase6/tier_consistency_master.csv
  release/v2/phase6/gate6G_tier_consistency.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
OUT = REPO / "release/v2/phase6"


def decode_bytes(df, *cols):
    for c in cols:
        if c in df.columns and df[c].dtype == object \
                and len(df) and isinstance(df[c].iloc[0], (bytes, bytearray)):
            df[c] = df[c].str.decode("utf-8")


def stat(arr):
    a = np.asarray(arr, dtype=float); a = a[np.isfinite(a)]
    if not len(a): return None
    return {"n": int(len(a)),
            "p16": float(np.percentile(a, 16)),
            "p50": float(np.percentile(a, 50)),
            "p84": float(np.percentile(a, 84))}


def main():
    v2     = Table.read(REPO / "release/v2/phase1/catalogue_v2.fits").to_pandas()
    orbits = Table.read(REPO / "release/v2/phase4/catalogue_v2_orbits.fits").to_pandas()
    pops   = Table.read(REPO / "release/v2/phase5/populations_v2.fits").to_pandas()
    xp     = Table.read(REPO / "release/v2/phase5/xp_metallicity_v2.fits").to_pandas()
    mc_path = REPO / "release/v2/phase6/catalogue_mc_orbits.fits"
    mc = Table.read(mc_path).to_pandas() if mc_path.exists() else None

    decode_bytes(v2, "tier", "dist_source", "rv_quality")
    decode_bytes(orbits, "tier", "rv_quality")
    decode_bytes(pops, "tier", "population", "chem_source")
    decode_bytes(xp, "tier")
    if mc is not None:
        decode_bytes(mc, "tier")

    rows = []
    for tier in ("A", "B", "C"):
        v_t = v2[v2["tier"] == tier]
        o_t = orbits[orbits["tier"] == tier]
        p_t = pops[pops["tier"] == tier]
        xp_t = xp[xp["tier"] == tier]
        n = len(v_t)

        d = {"tier": tier, "n": n}
        # Catalogue counts
        d["frac_of_total_2859"] = round(n / 2859, 4)
        # Vgrf
        s = stat(v_t["vgrf_default"]); d["vgrf_p50"] = s["p50"] if s else None
        d["vgrf_p16"] = s["p16"] if s else None
        d["vgrf_p84"] = s["p84"] if s else None
        # Point-estimate orbits
        if len(o_t):
            d["static_R_peri_pc_p50_pt"] = round(o_t["static_R_peri_kpc"].median() * 1000, 1)
            d["static_R_apo_kpc_p50_pt"] = round(o_t["static_R_apo_kpc"].median(), 2)
            d["static_z_max_kpc_p50_pt"] = round(o_t["static_z_max_kpc"].median(), 3)
            d["static_ecc_p50_pt"]       = round(o_t["static_ecc"].median(), 3)
            d["barred_R_peri_pc_p50_pt"] = round(o_t["barred_R_peri_kpc"].median() * 1000, 1)
            d["barred_R_apo_kpc_p50_pt"] = round(o_t["barred_R_apo_kpc"].median(), 2)
            d["n_R_peri_lt_100pc_static"] = int((o_t["static_R_peri_kpc"] < 0.100).sum())
            d["n_R_peri_lt_100pc_barred"] = int((o_t["barred_R_peri_kpc"] < 0.100).sum())
            d["n_min_rsph_lt_10pc_static"] = int((o_t["static_min_r_sph_kpc"] < 0.010).sum())
            d["n_min_rsph_lt_10pc_barred"] = int((o_t["barred_min_r_sph_kpc"] < 0.010).sum())
            d["n_bridgers_static"] = int(((o_t["static_R_peri_kpc"] < 2.0) & (o_t["static_R_apo_kpc"] > 15.0)).sum())
            d["n_bridgers_barred"] = int(((o_t["barred_R_peri_kpc"] < 2.0) & (o_t["barred_R_apo_kpc"] > 15.0)).sum())
            # Resonance
            r = o_t["res_ratio_OmegaR_over_dPhi"].dropna()
            d["frac_OLR_pm0p3"] = round((np.abs(r - 2.0) < 0.3).mean(), 4) if len(r) else None
            d["frac_4_1_pm0p3"] = round((np.abs(r - 4.0) < 0.3).mean(), 4) if len(r) else None
            d["frac_ILR_pm0p3"] = round((np.abs(r + 2.0) < 0.3).mean(), 4) if len(r) else None
        # Chemistry
        if len(p_t):
            d["n_with_alpha_chem"] = int(p_t["alpha_proxy"].notna().sum())
            d["frac_Splash"]       = int((p_t["population"] == "Splash").sum())
            d["frac_GSE"]          = int((p_t["population"] == "GSE").sum())
            d["frac_Aurora"]       = int((p_t["population"] == "Aurora").sum())
            d["frac_disk"]         = int((p_t["population"] == "disk").sum())
        # XP MH
        if len(xp_t):
            s = stat(xp_t.get("MH_xp", xp_t.get("mh_gspphot")))
            d["xp_MH_p16"] = s["p16"] if s else None
            d["xp_MH_p50"] = s["p50"] if s else None
            d["xp_MH_p84"] = s["p84"] if s else None
        # MC orbits (Phase 6A)
        if mc is not None:
            mc_t = mc[mc["tier"] == tier]
            if len(mc_t):
                d["mc_R_peri_pc_p16_p50"] = round(mc_t["R_peri_kpc_p16"].median() * 1000, 1)
                d["mc_R_peri_pc_p50_p50"] = round(mc_t["R_peri_kpc_p50"].median() * 1000, 1)
                d["mc_R_peri_pc_p84_p50"] = round(mc_t["R_peri_kpc_p84"].median() * 1000, 1)
                d["mc_ecc_p16_p50"] = round(mc_t["ecc_p16"].median(), 3)
                d["mc_ecc_p50_p50"] = round(mc_t["ecc_p50"].median(), 3)
                d["mc_ecc_p84_p50"] = round(mc_t["ecc_p84"].median(), 3)
        rows.append(d)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tier_consistency_master.csv", index=False)
    print(df.to_string(index=False))

    # Tier consistency check: A vs C
    flags = []
    quantities = ["static_R_peri_pc_p50_pt", "static_R_apo_kpc_p50_pt",
                   "static_ecc_p50_pt", "frac_OLR_pm0p3", "xp_MH_p50"]
    for q in quantities:
        if q in df.columns:
            try:
                vals = df.set_index("tier")[q]
                ratio = abs((vals.get("A") - vals.get("C")) / max(abs(vals.get("C", 0)), 1e-9))
                flags.append({"quantity": q,
                               "tier_A": vals.get("A"),
                               "tier_C": vals.get("C"),
                               "rel_diff": float(ratio),
                               "measurement_driven_flag": bool(ratio > 0.3)})
            except Exception:
                pass
    summary = {"per_tier_master": rows,
               "tier_consistency_flags": flags}
    (OUT / "gate6G_tier_consistency.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(flags, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
