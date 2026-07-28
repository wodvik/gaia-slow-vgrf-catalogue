"""Phase 16G -- recompute every manuscript-facing number on the retiered sample.

Collects, side by side, the forward-tier value currently in the manuscript and
the population-prior-tier replacement, so the text sweep can be done against
computed values rather than by hand.

Outputs: phase14/latent_deconvolution/manuscript_numbers.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
CAT = BUNDLE / "catalogues"
MC = BUNDLE / "phase14" / "expanded_orbit_mc"
OUT = BUNDLE / "phase14" / "latent_deconvolution" / "manuscript_numbers.json"


def jeffreys(k: int, n: int) -> tuple[float, float]:
    from scipy.stats import beta
    if n == 0:
        return float("nan"), float("nan")
    lo = beta.ppf(0.16, k + 0.5, n - k + 0.5) if k > 0 else 0.0
    hi = beta.ppf(0.84, k + 0.5, n - k + 0.5) if k < n else 1.0
    return float(lo), float(hi)


def kfrac(k, n):
    lo, hi = jeffreys(int(k), int(n))
    return {"k": int(k), "n": int(n), "pct": 100.0 * k / n if n else float("nan"),
            "pct_lo": 100 * lo, "pct_hi": 100 * hi}


def main() -> None:
    m = Table.read(CAT / "catalogue_retier_master.fits")
    orb_all = Table.read(CAT / "catalogue_expanded_orbits_tierABC.fits")
    mc_all = Table.read(MC / "expanded_catalogue_mc_orbits.fits")

    tier = np.asarray(m["tier"]).astype(str)
    tier_f = np.asarray(m["tier_forward"]).astype(str)
    sid = np.asarray(m["source_id"]).astype(np.int64)

    newA = tier == "A"
    newAB = np.isin(tier, ["A", "B"])
    newABC = np.isin(tier, ["A", "B", "C"])
    oldAB = np.isin(tier_f, ["A", "B"])
    oldABC = np.isin(tier_f, ["A", "B", "C"])

    ids_AB = set(sid[newAB].tolist())
    ids_ABC = set(sid[newABC].tolist())
    o_ab = np.array([int(s) in ids_AB for s in np.asarray(orb_all["source_id"])])
    o_abc = np.array([int(s) in ids_ABC for s in np.asarray(orb_all["source_id"])])
    c_abc = np.array([int(s) in ids_ABC for s in np.asarray(mc_all["source_id"])])
    orb = orb_all[o_abc]
    mc = mc_all[c_abc]

    R = {}
    R["counts"] = {"A": int(newA.sum()), "A+B": int(newAB.sum()), "A+B+C": int(newABC.sum()),
                   "C_only": int((tier == "C").sum()), "D": int((tier == "D").sum()),
                   "forward": {"A+B": int(oldAB.sum()), "A+B+C": int(oldABC.sum())}}

    # Parallax SNR split (manuscript Sec. 2.x)
    snr = np.asarray(m["parallax_over_error"], dtype=float)
    for lbl, sel in [("A+B", newAB), ("A+B+C", newABC)]:
        s = snr[sel]
        R[f"snr_{lbl}"] = {
            "in_5_10": kfrac(int(((s > 5) & (s <= 10)).sum()), int(sel.sum())),
            "gt_10": kfrac(int((s > 10).sum()), int(sel.sum())),
        }

    # RUWE / IPD / NSS
    ruwe = np.asarray(m["ruwe"], dtype=float)
    ipd = np.asarray(m["ipd_frac_multi_peak"], dtype=float)
    for lbl, sel in [("A+B", newAB), ("A+B+C", newABC)]:
        R[f"ruwe_gt14_{lbl}"] = kfrac(int((ruwe[sel] > 1.4).sum()), int(sel.sum()))
        R[f"ipd_le1_{lbl}"] = kfrac(int((ipd[sel] <= 1).sum()), int(sel.sum()))
    if "nss_two_body" in orb_all.colnames:
        nss = np.asarray(orb_all["nss_two_body"]).astype(bool)
        R["nss_A+B"] = kfrac(int(nss[o_ab].sum()), int(o_ab.sum()))
        R["nss_A+B+C"] = kfrac(int(nss[o_abc].sum()), int(o_abc.sum()))

    # Action reliability
    if "action_reliability_flag" in orb_all.colnames:
        fl = np.asarray(orb["action_reliability_flag"]).astype(str)
        R["action_sampled_poor_A+B+C"] = kfrac(int((fl == "sampled_poor").sum()), len(orb))
        R["action_flag_counts"] = {k: int((fl == k).sum()) for k in np.unique(fl)}

    # MC eccentricity statistics
    e = np.asarray(mc["ecc_p50"], dtype=float)
    R["ecc_MC"] = {
        "median": float(np.median(e)),
        "frac_gt_0.8": kfrac(int((e > 0.8).sum()), len(e)),
        "frac_gt_0.95": kfrac(int((e > 0.95).sum()), len(e)),
        "frac_gt_0.99": kfrac(int((e > 0.99).sum()), len(e)),
    }
    rp = np.asarray(mc["R_peri_kpc_p50"], dtype=float)
    R["MC_Rperi_pc_median"] = float(np.median(rp) * 1000)
    R["MC_Rperi_all_inside_1kpc"] = bool(np.all(rp < 1.0))
    R["MC_Rapo_kpc_median"] = float(np.median(np.asarray(mc["R_apo_kpc_p50"], dtype=float)))

    # Static point-estimate reach
    srp = np.asarray(orb["static_R_peri_kpc"], dtype=float)
    R["static_Rperi_lt_100pc"] = kfrac(int((srp < 0.1).sum()), len(srp))
    R["static_Rperi_all_lt_1kpc"] = bool(np.all(srp < 1.0))
    R["static_ecc_median"] = float(np.median(np.asarray(orb["static_ecc"], dtype=float)))
    R["static_Rperi_pc_median"] = float(np.median(srp) * 1000)
    R["static_Rapo_kpc_median"] = float(
        np.median(np.asarray(orb["static_R_apo_kpc"], dtype=float)))
    brp = np.asarray(orb["barred_R_peri_kpc"], dtype=float)
    R["barred_Rperi_lt_100pc"] = kfrac(int((brp < 0.1).sum()), len(brp))

    # Chemistry, both samples
    for lbl, sel in [("popprior_A+B+C", o_abc), ("forward_A+B+C", np.ones(len(orb_all), bool))]:
        cs = np.asarray(orb_all["chem_survey"]).astype(str)[sel]
        fe = np.asarray(orb_all["feh_spec"], dtype=float)[sel]
        R[f"chem_{lbl}"] = {"APOGEE": int((cs == "APOGEE").sum()),
                            "GALAH": int((cs == "GALAH").sum()),
                            "with_feh": int(np.isfinite(fe).sum()),
                            "n": int(sel.sum())}

    # Gold subset
    sd = (np.asarray(m["dist_hi_pc"], dtype=float) - np.asarray(m["dist_lo_pc"], dtype=float)) / 2
    fd = sd / np.asarray(m["dist_pc"], dtype=float)
    gold = (ruwe < 1.4) & np.asarray(m["rvs_quality_ok"]).astype(bool) & (snr > 10) & (fd < 0.15)
    R["gold"] = {"popprior_A+B": int((newAB & gold).sum()),
                 "popprior_A+B_with_ipd": int((newAB & gold & (ipd <= 1)).sum()),
                 "forward_A+B": int((oldAB & gold).sum())}
    R["frac_dist_gt_015_A+B+C"] = kfrac(int((fd[newABC] > 0.15).sum()), int(newABC.sum()))

    OUT.write_text(json.dumps(R, indent=2), encoding="utf-8")
    print(json.dumps(R, indent=2))


if __name__ == "__main__":
    main()
