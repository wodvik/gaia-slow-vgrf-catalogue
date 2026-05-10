"""Phase 6A — adaptive MC orbits for Tier A+B+C (the deferred 4D).

For each of the 632 Tier A+B+C stars:
  1. Sample 500 ICs from full uncertainty (3x3 (parallax,pmra,pmdec)
     covariance + RV + asymmetric Bailer-Jones distance).
  2. Integrate each in pot_static_full for 4 Gyr at trajsize=1001
     (4 Myr per step; sufficient for R_peri/R_apo).
  3. Compute Stäckel-fudge actions on each sample IC.
  4. Reduce to 16/50/84 percentiles per star.
  5. Adaptive refinement to 5,000 realisations if R_peri 16th < 200 pc
     OR R_peri 84th > 5 kpc OR star is in the 10 barred-bridger list
     from Phase 5F.
  6. Convergence test: random 100-star sub-sample re-run at 10,000.

Output: release/v2/phase6/catalogue_mc_orbits.fits
        release/v2/phase6/convergence_test.json
        release/v2/phase6/gate6A_mc.json

Throughput: with agama OpenMP across 24 cores, ~5,000 orbits/s. Total
orbit count for 6A is ~316K base + ~200×5000 refine + 100×10000 test
~= 1.7M orbits = ~6 minutes wall clock.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import agama
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase6"
OUT.mkdir(parents=True, exist_ok=True)
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"

R0 = float(CONFIG["solar_variants"]["default"]["R0_kpc"])
V0 = float(CONFIG["solar_variants"]["default"]["Vc_kms"])
SEED = int(CONFIG["mc"]["random_seed"])
GYR = 1.0 / 0.9778

agama.setUnits(length=1, mass=1, velocity=1)


def log(m):
    print(f"[6A t={time.time()-T0:6.1f}s] {m}", flush=True)


def galcen_frame():
    sv = CONFIG["solar_variants"]["default"]
    return coord.Galactocentric(
        galcen_distance=sv["R0_kpc"] * u.kpc,
        z_sun=sv["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            sv["U_kms"] * u.km/u.s,
            (sv["Vc_kms"] + sv["V_kms"]) * u.km/u.s,
            sv["W_kms"] * u.km/u.s,
        ),
    )


def cov3(plx_e, pmra_e, pmdec_e, plx_pmra, plx_pmdec, pmra_pmdec):
    n = len(plx_e)
    sig = np.stack([plx_e, pmra_e, pmdec_e], axis=-1)
    rho = np.zeros((n, 3, 3))
    for i in (0, 1, 2):
        rho[:, i, i] = 1.0
    rho[:, 0, 1] = rho[:, 1, 0] = plx_pmra
    rho[:, 0, 2] = rho[:, 2, 0] = plx_pmdec
    rho[:, 1, 2] = rho[:, 2, 1] = pmra_pmdec
    return rho * (sig[:, :, None] * sig[:, None, :])


def chol3(C):
    out = np.zeros_like(C)
    for i in range(C.shape[0]):
        try:
            out[i] = np.linalg.cholesky(C[i])
        except np.linalg.LinAlgError:
            out[i] = np.linalg.cholesky(C[i] + 1e-9 * np.eye(3))
    return out


def split_normal(med, lo, hi, n_samp, rng):
    sl = np.maximum(med - lo, 1.0)
    sh = np.maximum(hi - med, 1.0)
    u_ = rng.standard_normal((len(med), n_samp))
    sig = np.where(u_ < 0, sl[:, None], sh[:, None])
    return med[:, None] + u_ * sig


def sample_ics(df, n_samp, rng):
    """Sample (n_stars, n_samp, 6) galactocentric ICs."""
    n = len(df)
    mu = np.stack([df["parallax_zpcorr"].to_numpy(),
                   df["pmra"].to_numpy(),
                   df["pmdec"].to_numpy()], axis=-1)
    C = cov3(df["parallax_error"].to_numpy(),
             df["pmra_error"].to_numpy(),
             df["pmdec_error"].to_numpy(),
             df["parallax_pmra_corr"].to_numpy(),
             df["parallax_pmdec_corr"].to_numpy(),
             df["pmra_pmdec_corr"].to_numpy())
    L = chol3(C)
    z = rng.standard_normal((n, n_samp, 3))
    s = mu[:, None, :] + np.einsum("nij,nsj->nsi", L, z)
    pmra_s = s[:, :, 1]; pmdec_s = s[:, :, 2]
    rv_mu = df["radial_velocity"].to_numpy()
    rv_e = df["radial_velocity_error"].to_numpy()
    rv_s = rv_mu[:, None] + rng.standard_normal((n, n_samp)) * rv_e[:, None]
    dist_s = split_normal(df["dist_pc"].to_numpy(),
                          df["dist_lo_pc"].to_numpy(),
                          df["dist_hi_pc"].to_numpy(), n_samp, rng)
    dist_s = np.maximum(dist_s, 1.0)
    ra = np.broadcast_to(df["ra"].to_numpy()[:, None], (n, n_samp))
    dec = np.broadcast_to(df["dec"].to_numpy()[:, None], (n, n_samp))
    galcen = galcen_frame()
    icrs = coord.SkyCoord(ra=ra.flatten() * u.deg, dec=dec.flatten() * u.deg,
                          distance=dist_s.flatten() * u.pc,
                          pm_ra_cosdec=pmra_s.flatten() * u.mas/u.yr,
                          pm_dec=pmdec_s.flatten() * u.mas/u.yr,
                          radial_velocity=rv_s.flatten() * u.km/u.s,
                          frame="icrs")
    g = icrs.transform_to(galcen)
    ic = np.column_stack([
        g.x.to_value(u.kpc), g.y.to_value(u.kpc), g.z.to_value(u.kpc),
        g.v_x.to_value(u.km/u.s), g.v_y.to_value(u.km/u.s), g.v_z.to_value(u.km/u.s),
    ])  # shape (n*n_samp, 6)
    return ic.reshape(n, n_samp, 6)


def integrate_orbits(ic_flat, pot, n_steps=1001):
    """Integrate (M, 6) ICs and return per-orbit summary stats."""
    res = agama.orbit(potential=pot, ic=ic_flat, time=4.0 * GYR,
                      trajsize=n_steps)
    n = len(ic_flat)
    R_peri = np.zeros(n); R_apo = np.zeros(n)
    z_max = np.zeros(n);  ecc = np.zeros(n); rsph_min = np.zeros(n)
    for i in range(n):
        traj = np.asarray(res[i, 1])
        Rcyl = np.sqrt(traj[:, 0]**2 + traj[:, 1]**2)
        z = traj[:, 2]
        rsph = np.sqrt(Rcyl**2 + z**2)
        R_peri[i] = Rcyl.min(); R_apo[i] = Rcyl.max()
        z_max[i] = np.abs(z).max()
        ecc[i] = (R_apo[i] - R_peri[i]) / (R_apo[i] + R_peri[i])
        rsph_min[i] = rsph.min()
    return R_peri, R_apo, z_max, ecc, rsph_min


def percentiles(arr_n_x_s):
    """For a (n_stars, n_samples) array, return (n_stars, 3) of p16/50/84."""
    return np.percentile(arr_n_x_s, [16, 50, 84], axis=1).T  # shape (n,3)


def main():
    global T0; T0 = time.time()
    log(f"agama {agama.__version__}")

    # Load inputs
    log("loading catalogue + Phase 1 outputs ...")
    src = pd.read_csv(REPO / CONFIG["input"]["source_csv"])
    v2 = Table.read(REPO / "release/v2/phase1/catalogue_v2.fits").to_pandas()
    if v2["tier"].dtype == object and isinstance(v2["tier"].iloc[0], (bytes, bytearray)):
        v2["tier"] = v2["tier"].str.decode("utf-8")
    dist = Table.read(REPO / "release/v2/phase1/catalogue_dist.fits").to_pandas()
    zp = Table.read(REPO / "release/v2/phase1/catalogue_zpcorr.fits").to_pandas()

    # Tier A+B+C only
    tABC = v2["tier"].isin(["A", "B", "C"])
    sids = v2.loc[tABC, "source_id"].astype(int).to_numpy()
    log(f"Tier A+B+C: {len(sids)} stars")

    df = (src.merge(zp[["source_id", "parallax_zpcorr",
                         "parallax_zpcorr_error", "zpcorr_value_uas"]],
                     on="source_id")
              .merge(dist[["source_id", "dist_pc", "dist_lo_pc",
                           "dist_hi_pc"]],
                     on="source_id"))
    df = df[df["source_id"].isin(sids)].reset_index(drop=True)
    log(f"merged: {len(df)} rows")

    # Bridger source_ids from Phase 5F
    bridgers = Table.read(REPO / "release/v2/phase5/bridger_candidates.fits").to_pandas()
    bridger_ids = set(bridgers["source_id"].astype(int).tolist())
    log(f"bridger candidates from 5F: {len(bridger_ids)}")

    # Load static_full
    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    af = agama.ActionFinder(pot_axi, interp=True)

    # ----- Pass 1: 500 realisations for all 632 stars -----
    rng = np.random.default_rng(SEED + 6000)
    n0 = 500
    log(f"pass 1: sampling {n0} ICs/star for {len(df)} stars...")
    ic_3d = sample_ics(df, n0, rng)  # (n,500,6)
    log(f"  built ICs in {time.time()-T0:.1f}s; "
        f"{ic_3d.shape[0]*ic_3d.shape[1]:,} total orbits to integrate")
    ic_flat = ic_3d.reshape(-1, 6)
    log("  integrating...")
    Rp, Ra, zm, ec, rs = integrate_orbits(ic_flat, pot_axi, 1001)
    log(f"  integration done in {time.time()-T0:.1f}s")

    # Reshape per-star
    n_stars = len(df)
    Rp = Rp.reshape(n_stars, n0); Ra = Ra.reshape(n_stars, n0)
    zm = zm.reshape(n_stars, n0); ec = ec.reshape(n_stars, n0)
    rs = rs.reshape(n_stars, n0)

    # Actions on the same samples
    log("  computing actions...")
    acts = af(ic_flat)  # (n*n0, 3)
    JR = acts[:, 0].reshape(n_stars, n0)
    Jz = acts[:, 1].reshape(n_stars, n0)
    Jphi = acts[:, 2].reshape(n_stars, n0)
    log(f"  actions done at {time.time()-T0:.1f}s")

    # Reduce per-star to percentiles
    pcts = {
        "R_peri_kpc":   percentiles(Rp),
        "R_apo_kpc":    percentiles(Ra),
        "z_max_kpc":    percentiles(zm),
        "ecc":          percentiles(ec),
        "min_r_sph_kpc": percentiles(rs),
        "J_R":          percentiles(JR),
        "J_z":          percentiles(Jz),
        "J_phi":        percentiles(Jphi),
    }
    n_real_per_star = np.full(n_stars, n0, dtype=int)

    # ----- Identify refinement set -----
    Rp_p16 = pcts["R_peri_kpc"][:, 0]
    Rp_p84 = pcts["R_peri_kpc"][:, 2]
    refine_mask = (
        (Rp_p16 < 0.200) |
        (Rp_p84 > 5.0)   |
        df["source_id"].isin(bridger_ids).to_numpy()
    )
    n_refine = int(refine_mask.sum())
    log(f"refinement set: {n_refine} / {n_stars} stars (criteria union)")

    # ----- Pass 2: 5000 realisations for refinement set, CHUNKED -----
    if n_refine > 0:
        log("pass 2: 5000 realisations for refinement subset (chunked)...")
        df_ref = df.loc[refine_mask].reset_index(drop=True)
        idx_global = np.where(refine_mask)[0]
        rng2 = np.random.default_rng(SEED + 6001)

        # Chunk by 50 stars at a time => 50*5000=250K orbits/chunk;
        # ~12 GB transient trajectory memory.
        CHUNK = 50
        for start in range(0, n_refine, CHUNK):
            stop = min(start + CHUNK, n_refine)
            log(f"  chunk {start}-{stop} of {n_refine}")
            df_c = df_ref.iloc[start:stop].reset_index(drop=True)
            ic_c_3d = sample_ics(df_c, 5000, rng2)  # (m,5000,6)
            ic_c_flat = ic_c_3d.reshape(-1, 6)
            Rp_c, Ra_c, zm_c, ec_c, rs_c = integrate_orbits(ic_c_flat, pot_axi, 1001)
            acts_c = af(ic_c_flat)
            m = stop - start
            Rp_c = Rp_c.reshape(m, 5000); Ra_c = Ra_c.reshape(m, 5000)
            zm_c = zm_c.reshape(m, 5000); ec_c = ec_c.reshape(m, 5000)
            rs_c = rs_c.reshape(m, 5000)
            JR_c = acts_c[:, 0].reshape(m, 5000)
            Jz_c = acts_c[:, 1].reshape(m, 5000)
            Jphi_c = acts_c[:, 2].reshape(m, 5000)
            chunk_idx = idx_global[start:stop]
            pcts["R_peri_kpc"][chunk_idx]    = percentiles(Rp_c)
            pcts["R_apo_kpc"][chunk_idx]     = percentiles(Ra_c)
            pcts["z_max_kpc"][chunk_idx]     = percentiles(zm_c)
            pcts["ecc"][chunk_idx]           = percentiles(ec_c)
            pcts["min_r_sph_kpc"][chunk_idx] = percentiles(rs_c)
            pcts["J_R"][chunk_idx]           = percentiles(JR_c)
            pcts["J_z"][chunk_idx]           = percentiles(Jz_c)
            pcts["J_phi"][chunk_idx]         = percentiles(Jphi_c)
            n_real_per_star[chunk_idx] = 5000
            # Free memory
            del Rp_c, Ra_c, zm_c, ec_c, rs_c, JR_c, Jz_c, Jphi_c, ic_c_3d, ic_c_flat, acts_c
        log(f"  pass 2 done at {time.time()-T0:.1f}s")

    # ----- Convergence test: 100 random stars at 10,000, CHUNKED -----
    rng3 = np.random.default_rng(SEED + 6002)
    n_test = 100
    test_idx = rng3.choice(n_stars, size=n_test, replace=False)
    df_test = df.iloc[test_idx].reset_index(drop=True)
    log(f"convergence test: {n_test} stars × 10000 realisations (chunked)...")
    Rp_t_p = np.zeros((n_test, 3)); ec_t_p = np.zeros((n_test, 3))
    CHUNK_T = 25  # 25 * 10000 = 250K orbits/chunk = ~12 GB transient
    for start in range(0, n_test, CHUNK_T):
        stop = min(start + CHUNK_T, n_test)
        df_tc = df_test.iloc[start:stop].reset_index(drop=True)
        ic_tc_3d = sample_ics(df_tc, 10000, rng3)
        ic_tc_flat = ic_tc_3d.reshape(-1, 6)
        Rp_tc, Ra_tc, zm_tc, ec_tc, rs_tc = integrate_orbits(ic_tc_flat, pot_axi, 1001)
        m = stop - start
        Rp_t_p[start:stop] = percentiles(Rp_tc.reshape(m, 10000))
        ec_t_p[start:stop] = percentiles(ec_tc.reshape(m, 10000))
        del Rp_tc, Ra_tc, zm_tc, ec_tc, rs_tc, ic_tc_3d, ic_tc_flat
    # Compare to current pcts
    cur_Rp = pcts["R_peri_kpc"][test_idx]
    cur_ec = pcts["ecc"][test_idx]
    dRp = np.abs(Rp_t_p - cur_Rp) * 1000  # pc
    dec = np.abs(ec_t_p - cur_ec)
    convergence = {
        "n_test_stars": n_test,
        "Rp_pc_max_diff_p50": float(dRp[:, 1].max()),
        "Rp_pc_p84_diff_p50": float(np.percentile(dRp[:, 1], 84)),
        "ecc_max_diff_p50": float(dec[:, 1].max()),
        "ecc_p84_diff_p50": float(np.percentile(dec[:, 1], 84)),
        "spec_R_peri_50pc":  bool(dRp[:, 1].max() < 50),
        "spec_ecc_0p005":    bool(dec[:, 1].max() < 0.005),
    }
    (OUT / "convergence_test.json").write_text(json.dumps(convergence, indent=2))
    log(f"convergence: max |dR_peri p50|={dRp[:,1].max():.1f} pc, "
        f"max |de p50|={dec[:,1].max():.4f}")

    # ----- Assemble master MC catalogue -----
    out = Table()
    out["source_id"] = df["source_id"].to_numpy()
    out["tier"] = df.merge(v2[["source_id", "tier"]],
                            on="source_id")["tier"].to_numpy().astype(str)
    out["n_realisations"] = n_real_per_star
    for k, p in pcts.items():
        out[f"{k}_p16"] = p[:, 0]
        out[f"{k}_p50"] = p[:, 1]
        out[f"{k}_p84"] = p[:, 2]
    out.write(OUT / "catalogue_mc_orbits.fits", overwrite=True)
    log(f"wrote catalogue_mc_orbits.fits ({len(out)} rows)")

    # Gate6A summary
    def stat(v):
        v = np.asarray(v, dtype=float); v = v[np.isfinite(v)]
        return ({"n": int(len(v)),
                 "p16": float(np.percentile(v, 16)),
                 "p50": float(np.percentile(v, 50)),
                 "p84": float(np.percentile(v, 84))} if len(v) else None)
    summary = {
        "n_TierABC": int(n_stars),
        "n_pass1_realisations": n0,
        "n_refined_to_5000": n_refine,
        "n_convergence_test_at_10000": n_test,
        "convergence": convergence,
        "headline": {
            "median_R_peri_p50_pc": stat(pcts["R_peri_kpc"][:, 1] * 1000),
            "median_R_peri_p16_pc": stat(pcts["R_peri_kpc"][:, 0] * 1000),
            "median_R_peri_p84_pc": stat(pcts["R_peri_kpc"][:, 2] * 1000),
            "median_ecc_p50":       stat(pcts["ecc"][:, 1]),
            "median_R_apo_p50_kpc": stat(pcts["R_apo_kpc"][:, 1]),
            "median_z_max_p50_kpc": stat(pcts["z_max_kpc"][:, 1]),
        },
    }
    # Per-tier breakdown
    tiers_df = df.merge(v2[["source_id", "tier"]], on="source_id")
    summary["per_tier"] = {}
    for t in ("A", "B", "C"):
        m = (tiers_df["tier"] == t).to_numpy()
        if m.sum() == 0: continue
        summary["per_tier"][t] = {
            "n": int(m.sum()),
            "median_R_peri_p50_pc": stat(pcts["R_peri_kpc"][m, 1] * 1000),
            "median_ecc_p50":       stat(pcts["ecc"][m, 1]),
            "median_R_apo_p50_kpc": stat(pcts["R_apo_kpc"][m, 1]),
        }
    (OUT / "gate6A_mc.json").write_text(json.dumps(summary, indent=2))
    log("wrote gate6A_mc.json")
    print(json.dumps(summary, indent=2))
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
