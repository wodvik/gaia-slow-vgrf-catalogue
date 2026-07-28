"""Phase 16K -- recompute the hand-authored table values on the adopted tiers.

The tables under tables/v15/ are hand-authored rather than generated, and were
written against the forward-score tiers. This script recomputes every value
they contain from the retiered products so the edits can be made against
computed numbers rather than by hand.

Prints a block per table. Nothing is written; the tables are edited from this
output so the provenance of each change stays visible in the diff.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from astropy.table import Table
from pathlib import Path
from scipy.stats import beta

BUNDLE = Path(__file__).resolve().parents[1]
CAT = BUNDLE / "catalogues"
MC = BUNDLE / "phase14" / "expanded_orbit_mc"


def jeff(k, n):
    lo = beta.ppf(0.16, k + 0.5, n - k + 0.5) if k > 0 else 0.0
    hi = beta.ppf(0.84, k + 0.5, n - k + 0.5) if k < n else 1.0
    return 100 * lo, 100 * hi


def pct(k, n):
    f = 100 * k / n
    lo, hi = jeff(k, n)
    return f"${f:.1f}^{{+{hi-f:.1f}}}_{{-{f-lo:.1f}}}$"


def main():
    m = Table.read(CAT / "catalogue_retier_master.fits")
    orb = Table.read(CAT / "catalogue_retier_orbits_tierABC.fits")
    mcall = Table.read(MC / "expanded_catalogue_mc_orbits.fits")

    tier = np.asarray(m["tier"]).astype(str)
    sid = np.asarray(m["source_id"]).astype(np.int64)
    ruwe = np.asarray(m["ruwe"], dtype=float)
    snr = np.asarray(m["parallax_over_error"], dtype=float)
    rvok = np.asarray(m["rvs_quality_ok"]).astype(bool)
    sd = (np.asarray(m["dist_hi_pc"], dtype=float) - np.asarray(m["dist_lo_pc"], dtype=float)) / 2
    fd = sd / np.asarray(m["dist_pc"], dtype=float)
    vgrf = np.asarray(m["vgrf_default"], dtype=float)
    gmag = np.asarray(m["phot_g_mean_mag"], dtype=float)
    dist = np.asarray(m["dist_pc"], dtype=float)

    AB = np.isin(tier, ["A", "B"])
    print("=" * 70, "\ntab_gold_empirical")
    steps = [("Tier A+B base", AB),
             ("& RUWE<1.4", AB & (ruwe < 1.4)),
             ("& RVS ok", AB & (ruwe < 1.4) & rvok),
             ("& parallax SNR>10", AB & (ruwe < 1.4) & rvok & (snr > 10)),
             ("& sigma_d/d<0.15", AB & (ruwe < 1.4) & rvok & (snr > 10) & (fd < 0.15))]
    for lbl, msk in steps:
        print(f"  {lbl:<22} {int(msk.sum())}")
    gold = steps[-1][1]
    print(f"  median Vgrf  {np.nanmedian(vgrf[gold]):.1f}")
    if gmag is not None:
        print(f"  median G      {np.nanmedian(gmag[gold]):.1f}")
    print(f"  median dist  {np.nanmedian(dist[gold]):.0f}")

    # --- gold dynamical: MC posterior medians on Tier A+B and Gold ---
    gold_ids = set(sid[gold].tolist())
    ab_ids = set(sid[AB].tolist())
    msid = np.asarray(mcall["source_id"]).astype(np.int64)
    sel_ab = np.array([int(s) in ab_ids for s in msid])
    sel_gd = np.array([int(s) in gold_ids for s in msid])
    # tab_gold_dynamical states "quantities here are point-estimate outputs of
    # the adopted static potential" -- so it must be built from the orbit table,
    # NOT from the Monte Carlo posterior medians (which run systematically less
    # eccentric because the error convolution regresses e away from unity).
    print("=" * 70, "\ntab_gold_dynamical  (POINT ESTIMATE, static potential)")
    osid = np.asarray(orb["source_id"]).astype(np.int64)
    o_ab = np.array([int(s) in ab_ids for s in osid])
    o_gd = np.array([int(s) in gold_ids for s in osid])
    for lbl, s in [("Tier A+B", o_ab), ("Gold", o_gd)]:
        e = np.asarray(orb["static_ecc"], dtype=float)[s]
        rp = np.asarray(orb["static_R_peri_kpc"], dtype=float)[s] * 1000
        ra = np.asarray(orb["static_R_apo_kpc"], dtype=float)[s]
        zm = np.asarray(orb["static_z_max_kpc"], dtype=float)[s]
        n = len(e)
        print(f"  {lbl:<10} N={n}  ecc={np.median(e):.3f}  "
              f"e>0.8={100*np.mean(e>0.8):.1f}  e>0.95={pct(int((e>0.95).sum()),n)}  "
              f"e>0.99={pct(int((e>0.99).sum()),n)}")
        print(f"             Rperi={np.median(rp):.1f} pc  Rapo={np.median(ra):.2f} kpc  "
              f"zmax={np.median(zm):.2f} kpc")

    # --- RV quality sensitivity (point-estimate medians on orbit table) ---
    print("=" * 70, "\ntab_rv_quality_sensitivity")
    otier = np.asarray(orb["tier"]).astype(str)
    orp = np.asarray(orb["static_R_peri_kpc"], dtype=float) * 1000
    oe = np.asarray(orb["static_ecc"], dtype=float)
    ook = np.asarray(orb["rvs_quality_ok"]).astype(bool)
    for lbl, msk in [("A+B+C all", np.isin(otier, ["A", "B", "C"])),
                     ("A+B+C excl poor", np.isin(otier, ["A", "B", "C"]) & ook),
                     ("A+B all", np.isin(otier, ["A", "B"])),
                     ("A+B excl poor", np.isin(otier, ["A", "B"]) & ook)]:
        print(f"  {lbl:<18} N={int(msk.sum()):5d}  Rperi={np.median(orp[msk]):.0f} pc  "
              f"e={np.median(oe[msk]):.3f}")

    # --- band comparison slow row ---
    print("=" * 70, "\ntab_band_comparison (slow row)")
    abc = np.isin(otier, ["A", "B", "C"])
    print(f"  slow A+B+C  N={int(abc.sum())}  Rperi={np.median(orp[abc]):.0f} pc  "
          f"e={np.median(oe[abc]):.3f}")

    # --- reach table ---
    print("=" * 70, "\ntab_reach")
    brp = np.asarray(orb["barred_R_peri_kpc"], dtype=float) * 1000
    srs = np.asarray(orb["static_min_r_sph_kpc"], dtype=float) * 1000
    brs = np.asarray(orb["barred_min_r_sph_kpc"], dtype=float) * 1000
    for lbl, arr, thr in [("Rperi<100pc static", orp, 100), ("Rperi<50pc static", orp, 50),
                          ("Rperi<10pc static", orp, 10),
                          ("Rperi<100pc barred", brp, 100), ("Rperi<50pc barred", brp, 50),
                          ("Rperi<10pc barred", brp, 10),
                          ("min r_sph<100pc static", srs, 100), ("min r_sph<10pc static", srs, 10),
                          ("min r_sph<100pc barred", brs, 100), ("min r_sph<10pc barred", brs, 10)]:
        print(f"  {lbl:<24} {int((arr[abc] < thr).sum())}")

    # --- threshold table (point-estimate counts at alternative cuts) ---
    print("=" * 70, "\ntab_threshold (population-prior tiers at each Vgrf cut)")
    print("  NOTE: alternative-threshold tiers require re-running the MC at that")
    print("  cut; only the <25 row is recomputable from shipped products.")
    print(f"  <25: point-estimate pool={int((vgrf<25).sum())}, Tier A+B+C={int(abc.sum())}, "
          f"e={np.median(oe[abc]):.3f}, e>0.95={pct(int((oe[abc]>0.95).sum()), int(abc.sum()))}")


if __name__ == "__main__":
    main()
