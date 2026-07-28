"""Phase 16F -- retier the catalogue on population-prior probabilities.

The released tiers were defined on the FORWARD Monte Carlo score
P(Vgrf<25 | data_i), which carries no information about the parent
population. Phase 16C/16D showed that this biases the implied purity high by
a large factor, and Phase 16E showed the population-prior correction is
perfectly NESTED: no star enters any tier, so retiering moves the threshold
without reordering the catalogue.

This script rebuilds the catalogue products on the corrected score:

  P_pop = sum_{v<25} w_ij   from the regularised EM reconstruction

Both scores ship. P_forward is retained under its original column name so
the released products stay backward-compatible with the deposited catalogue,
and P_pop becomes the tier-defining quantity.

Because the new tiers are nested subsets of the old ones, every per-star
orbit product already exists and is simply re-aggregated -- no orbit is
re-integrated.

Outputs (written alongside the existing products, suffixed _retier):
  catalogues/catalogue_retier_master.fits
  catalogues/catalogue_retier_tier{A,AB,ABC}.fits
  catalogues/catalogue_retier_orbits_tierABC.fits
  catalogues/retier_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table, join

BUNDLE = Path(__file__).resolve().parents[1]
CAT = BUNDLE / "catalogues"
LAT = BUNDLE / "phase14" / "latent_deconvolution"
MC = BUNDLE / "phase14" / "expanded_orbit_mc"

THRESH = {"A": 0.95, "B": 0.84, "C": 0.50}


def assign_tier(p: np.ndarray, point_below: np.ndarray) -> np.ndarray:
    tier = np.full(len(p), "X", dtype=object)
    tier[point_below] = "D"
    tier[p > THRESH["C"]] = "C"
    tier[p > THRESH["B"]] = "B"
    tier[p > THRESH["A"]] = "A"
    tier[(~point_below) & (p <= THRESH["C"])] = "X"
    return tier


def pct(a, name, scale=1.0):
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)] * scale
    if a.size == 0:
        return {"n": 0}
    return {"n": int(a.size),
            "p16": float(np.percentile(a, 16)),
            "p50": float(np.percentile(a, 50)),
            "p84": float(np.percentile(a, 84)),
            "min": float(a.min()), "max": float(a.max())}


def main() -> None:
    master = Table.read(CAT / "catalogue_expanded_master.fits")
    lat = pd.read_csv(LAT / "latent_vgrf_per_star_regularised.csv")
    lat_map = dict(zip(lat.source_id.to_numpy(), lat.P_latent.to_numpy()))

    # The master FITS does not carry the image-parameter column; pull it from
    # the candidate table so the Gold/IPD variants can be evaluated here.
    aux = pd.read_csv(BUNDLE / "private_inputs" / "expanded_candidates_mc_tiered.csv",
                      usecols=["source_id", "ipd_frac_multi_peak"], low_memory=False)
    ipd_map = dict(zip(aux.source_id.to_numpy(), aux.ipd_frac_multi_peak.to_numpy()))
    master["ipd_frac_multi_peak"] = np.array(
        [ipd_map.get(int(s), np.nan) for s in np.asarray(master["source_id"])], dtype=float)

    sid = np.asarray(master["source_id"])
    p_pop = np.array([lat_map.get(int(s), np.nan) for s in sid])
    n_missing = int(np.isnan(p_pop).sum())
    # Stars absent from the likelihood run (incomplete 6D) keep a zero score:
    # they were never tier members under either definition.
    p_pop_filled = np.clip(np.where(np.isfinite(p_pop), p_pop, 0.0), 0.0, 1.0)
    assert (p_pop_filled >= 0.0).all() and (p_pop_filled <= 1.0).all()

    master["P_vgrf_below_25_forward"] = np.asarray(master["P_vgrf_below_25"], dtype=float)
    master["P_vgrf_below_25_popprior"] = p_pop_filled
    master["tier_forward"] = np.asarray(master["tier"]).astype(str)

    point_below = np.asarray(master["vgrf_default"], dtype=float) < 25.0
    new_tier = assign_tier(p_pop_filled, point_below)
    master["tier"] = np.array(new_tier, dtype=str)

    is_A = master["tier"] == "A"
    is_AB = np.isin(master["tier"], ["A", "B"])
    is_ABC = np.isin(master["tier"], ["A", "B", "C"])

    master.write(CAT / "catalogue_retier_master.fits", overwrite=True)
    master[is_A].write(CAT / "catalogue_retier_tierA.fits", overwrite=True)
    master[is_AB].write(CAT / "catalogue_retier_tierAB.fits", overwrite=True)
    master[is_ABC].write(CAT / "catalogue_retier_tierABC.fits", overwrite=True)

    # ---- verify nesting against the released tiers ----
    old_A = master["tier_forward"] == "A"
    old_AB = np.isin(master["tier_forward"], ["A", "B"])
    old_ABC = np.isin(master["tier_forward"], ["A", "B", "C"])
    nesting = {
        "A": {"added": int((is_A & ~old_A).sum()), "dropped": int((~is_A & old_A).sum())},
        "A+B": {"added": int((is_AB & ~old_AB).sum()), "dropped": int((~is_AB & old_AB).sum())},
        "A+B+C": {"added": int((is_ABC & ~old_ABC).sum()), "dropped": int((~is_ABC & old_ABC).sum())},
    }

    # ---- re-aggregate orbit products over the nested subsets ----
    orb = Table.read(CAT / "catalogue_expanded_orbits_tierABC.fits")
    mc = Table.read(MC / "expanded_catalogue_mc_orbits.fits")
    new_ids = set(int(s) for s in np.asarray(master["source_id"])[is_ABC])
    new_ab_ids = set(int(s) for s in np.asarray(master["source_id"])[is_AB])

    orb_sel = np.array([int(s) in new_ids for s in np.asarray(orb["source_id"])])
    orb_ab = np.array([int(s) in new_ab_ids for s in np.asarray(orb["source_id"])])
    mc_sel = np.array([int(s) in new_ids for s in np.asarray(mc["source_id"])])

    orb_new = orb[orb_sel]
    # The `tier` column must carry the POPULATION-PRIOR label, not the forward
    # one. Subsetting the rows alone leaves forward labels inside the retiered
    # row set, which silently mis-styles every downstream per-tier plot and
    # legend (the totals look right; the A/B/C split does not).
    tier_map = dict(zip(np.asarray(master["source_id"]).astype(np.int64).tolist(),
                        np.asarray(master["tier"]).astype(str).tolist()))
    new_labels = np.array([tier_map[int(s)] for s in np.asarray(orb_new["source_id"])],
                          dtype=str)
    orb_new["tier_forward"] = np.asarray(orb_new["tier"]).astype(str)
    orb_new["tier"] = new_labels
    orb_new["tier_popprior"] = new_labels
    assert set(np.unique(new_labels)) <= {"A", "B", "C"}, "unexpected tier label"
    orb_new.write(CAT / "catalogue_retier_orbits_tierABC.fits", overwrite=True)

    summary = {
        "method": "population-prior retiering (Phase 16C/16D/16E)",
        "thresholds": THRESH,
        "n_pool": int(len(master)),
        "n_missing_popprior_score": n_missing,
        "tier_counts_forward": {k: int((master["tier_forward"] == k).sum())
                                for k in ["A", "B", "C", "D", "X"]},
        "tier_counts_popprior": {k: int((master["tier"] == k).sum())
                                 for k in ["A", "B", "C", "D", "X"]},
        "headline_counts": {
            "forward": {"A": int(old_A.sum()), "A+B": int(old_AB.sum()),
                        "A+B+C": int(old_ABC.sum())},
            "popprior": {"A": int(is_A.sum()), "A+B": int(is_AB.sum()),
                         "A+B+C": int(is_ABC.sum())},
        },
        "nesting_check": nesting,
        "purity": {
            "forward_tiers_forward_score": {
                "A": float(np.asarray(master["P_vgrf_below_25_forward"])[old_A].mean()),
                "A+B": float(np.asarray(master["P_vgrf_below_25_forward"])[old_AB].mean()),
                "A+B+C": float(np.asarray(master["P_vgrf_below_25_forward"])[old_ABC].mean()),
            },
            "forward_tiers_popprior_score": {
                "A": float(p_pop_filled[old_A].mean()),
                "A+B": float(p_pop_filled[old_AB].mean()),
                "A+B+C": float(p_pop_filled[old_ABC].mean()),
            },
            "popprior_tiers_popprior_score": {
                "A": float(p_pop_filled[is_A].mean()),
                "A+B": float(p_pop_filled[is_AB].mean()),
                "A+B+C": float(p_pop_filled[is_ABC].mean()),
            },
        },
        "expected_true_slow": {
            "forward_tiers": float(np.asarray(master["P_vgrf_below_25_forward"])[old_ABC].sum()),
            "popprior_tiers": float(p_pop_filled[is_ABC].sum()),
        },
    }

    # MC posterior-median headline quantities (the abstract's numbers).
    mc_new = mc[mc_sel]
    summary["mc_headline_popprior_tierABC"] = {
        "R_peri_p50_pc": pct(mc_new["R_peri_kpc_p50"], "rp", 1000.0),
        "R_apo_p50_kpc": pct(mc_new["R_apo_kpc_p50"], "ra"),
        "ecc_p50": pct(mc_new["ecc_p50"], "e"),
        "z_max_p50_kpc": pct(mc_new["z_max_kpc_p50"], "z"),
    }
    summary["mc_headline_forward_tierABC"] = {
        "R_peri_p50_pc": pct(mc["R_peri_kpc_p50"], "rp", 1000.0),
        "R_apo_p50_kpc": pct(mc["R_apo_kpc_p50"], "ra"),
        "ecc_p50": pct(mc["ecc_p50"], "e"),
        "z_max_p50_kpc": pct(mc["z_max_kpc_p50"], "z"),
    }

    # Point-estimate (static/barred) summaries.
    summary["static_popprior"] = {
        "R_peri_kpc_tierAB": pct(orb[orb_ab]["static_R_peri_kpc"], "x"),
        "R_peri_kpc_tierABC": pct(orb_new["static_R_peri_kpc"], "x"),
        "ecc_tierABC": pct(orb_new["static_ecc"], "x"),
        "R_apo_kpc_tierABC": pct(orb_new["static_R_apo_kpc"], "x"),
        "barred_R_peri_kpc_tierABC": pct(orb_new["barred_R_peri_kpc"], "x"),
        "n_static_R_peri_lt_100pc": int((np.asarray(orb_new["static_R_peri_kpc"]) < 0.1).sum()),
        "n_barred_R_peri_lt_100pc": int((np.asarray(orb_new["barred_R_peri_kpc"]) < 0.1).sum()),
    }

    # Quality-flag fractions on the new tiers.
    def frac(tbl, mask, col, cond):
        v = np.asarray(tbl[col], dtype=float)[mask]
        v = v[np.isfinite(v)]
        return {"n": int(v.size), "k": int(cond(v).sum()),
                "frac": float(cond(v).mean()) if v.size else float("nan")}

    summary["quality_popprior"] = {
        "ruwe_gt_1p4_tierAB": frac(master, is_AB, "ruwe", lambda v: v > 1.4),
        "ruwe_gt_1p4_tierABC": frac(master, is_ABC, "ruwe", lambda v: v > 1.4),
        "ipd_multi_peak_gt1_tierAB": frac(master, is_AB, "ipd_frac_multi_peak", lambda v: v > 1),
        "ipd_multi_peak_gt1_tierABC": frac(master, is_ABC, "ipd_frac_multi_peak", lambda v: v > 1),
    }
    if "nss_two_body" in orb.colnames:
        summary["quality_popprior"]["nss_tierABC"] = int(
            np.asarray(orb_new["nss_two_body"]).astype(bool).sum())
        summary["quality_popprior"]["nss_tierAB"] = int(
            np.asarray(orb[orb_ab]["nss_two_body"]).astype(bool).sum())

    # Gold subset under the new tiers.
    sigma_d = (np.asarray(master["dist_hi_pc"], dtype=float)
               - np.asarray(master["dist_lo_pc"], dtype=float)) / 2.0
    frac_d = sigma_d / np.asarray(master["dist_pc"], dtype=float)
    gold_cuts = (
        (np.asarray(master["ruwe"], dtype=float) < 1.4)
        & (np.asarray(master["rvs_quality_ok"]).astype(bool)
           if "rvs_quality_ok" in master.colnames else True)
        & (np.asarray(master["parallax_over_error"], dtype=float) > 10)
        & (frac_d < 0.15)
    )
    summary["gold"] = {
        "forward_tierAB": int((old_AB & gold_cuts).sum()),
        "popprior_tierAB": int((is_AB & gold_cuts).sum()),
        "popprior_tierAB_with_ipd": int(
            (is_AB & gold_cuts
             & (np.asarray(master["ipd_frac_multi_peak"], dtype=float) <= 1)).sum()),
    }

    # Chemistry overlap under the new tiers.
    if "chem_survey" in orb.colnames:
        cs = np.asarray(orb_new["chem_survey"]).astype(str)
        summary["chemistry_popprior_tierABC"] = {
            k: int((cs == k).sum()) for k in np.unique(cs)}

    (CAT / "retier_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
