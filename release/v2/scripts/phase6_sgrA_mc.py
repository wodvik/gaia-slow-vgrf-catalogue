"""Phase 6 supplemental — MC propagation of 4 Sgr-A* candidate approachers.

Runs 5,000-realisation MC on the 4 Phase-6B candidates:
  static_full Tier A:  1846633734516771840
  barred Tier B:        431014850227388672
  barred Tier A:       3738499345877271808
  barred Tier B:       4287640292264861056

For each candidate, integrate 5000 sample orbits at trajsize=40001
(0.1 Myr) and compute:
  - P(min r_sph < 10 pc) -- the strict approach test
  - P(min r_sph < 100 pc) -- a softer claim
  - 16/50/84 percentiles of min r_sph

Static MC integrates in pot_static_full; barred MC in pot_rotating
with the same Omega applied at integration time.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd, yaml, agama
import astropy.coordinates as coord, astropy.units as u
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase6"
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"
OMEGA_P = -float(CONFIG["bar_pattern_speeds_kms_kpc"]["default"])
SEED = int(CONFIG["mc"]["random_seed"])
GYR = 1.0 / 0.9778
agama.setUnits(length=1, mass=1, velocity=1)

CANDIDATES = [
    {"source_id": 1846633734516771840, "potential": "static",  "tier": "A"},
    {"source_id":  431014850227388672, "potential": "barred",  "tier": "B"},
    {"source_id": 3738499345877271808, "potential": "barred",  "tier": "A"},
    {"source_id": 4287640292264861056, "potential": "barred",  "tier": "B"},
]
N_SAMP = 5000


def log(m): print(f"[sgrA t={time.time()-T0:5.1f}s] {m}", flush=True)


def galcen():
    sv = CONFIG["solar_variants"]["default"]
    return coord.Galactocentric(
        galcen_distance=sv["R0_kpc"] * u.kpc, z_sun=sv["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            sv["U_kms"] * u.km/u.s,
            (sv["Vc_kms"] + sv["V_kms"]) * u.km/u.s,
            sv["W_kms"] * u.km/u.s))


def cov3(plx_e, pmra_e, pmdec_e, p_pmra, p_pmdec, pmra_pmdec):
    n = len(plx_e); sig = np.stack([plx_e, pmra_e, pmdec_e], axis=-1)
    rho = np.zeros((n, 3, 3))
    for i in (0,1,2): rho[:, i, i] = 1.0
    rho[:, 0, 1] = rho[:, 1, 0] = p_pmra
    rho[:, 0, 2] = rho[:, 2, 0] = p_pmdec
    rho[:, 1, 2] = rho[:, 2, 1] = pmra_pmdec
    return rho * (sig[:, :, None] * sig[:, None, :])


def chol(C):
    out = np.zeros_like(C)
    for i in range(C.shape[0]):
        try: out[i] = np.linalg.cholesky(C[i])
        except np.linalg.LinAlgError:
            out[i] = np.linalg.cholesky(C[i] + 1e-9 * np.eye(3))
    return out


def split_normal(med, lo, hi, n_samp, rng):
    sl = np.maximum(med - lo, 1.0); sh = np.maximum(hi - med, 1.0)
    u_ = rng.standard_normal((len(med), n_samp))
    sig = np.where(u_ < 0, sl[:, None], sh[:, None])
    return med[:, None] + u_ * sig


def sample_ics(df_one_row, n_samp, rng):
    df = df_one_row
    n = 1
    mu = np.array([[df["parallax_zpcorr"], df["pmra"], df["pmdec"]]])
    C = cov3(np.array([df["parallax_error"]]), np.array([df["pmra_error"]]),
             np.array([df["pmdec_error"]]),
             np.array([df["parallax_pmra_corr"]]),
             np.array([df["parallax_pmdec_corr"]]),
             np.array([df["pmra_pmdec_corr"]]))
    L = chol(C)
    z = rng.standard_normal((1, n_samp, 3))
    s = mu[:, None, :] + np.einsum("nij,nsj->nsi", L, z)
    pmra_s = s[0, :, 1]; pmdec_s = s[0, :, 2]
    rv_s = df["radial_velocity"] + rng.standard_normal(n_samp) * df["radial_velocity_error"]
    dist_s = split_normal(np.array([df["dist_pc"]]), np.array([df["dist_lo_pc"]]),
                           np.array([df["dist_hi_pc"]]), n_samp, rng)[0]
    dist_s = np.maximum(dist_s, 1.0)
    icrs = coord.SkyCoord(
        ra=np.full(n_samp, df["ra"]) * u.deg,
        dec=np.full(n_samp, df["dec"]) * u.deg,
        distance=dist_s * u.pc,
        pm_ra_cosdec=pmra_s * u.mas/u.yr,
        pm_dec=pmdec_s * u.mas/u.yr,
        radial_velocity=rv_s * u.km/u.s, frame="icrs")
    g = icrs.transform_to(galcen())
    return np.column_stack([
        g.x.to_value(u.kpc), g.y.to_value(u.kpc), g.z.to_value(u.kpc),
        g.v_x.to_value(u.km/u.s), g.v_y.to_value(u.km/u.s), g.v_z.to_value(u.km/u.s)])


def main():
    global T0; T0 = time.time()
    log(f"agama {agama.__version__}")

    src = pd.read_csv(REPO / CONFIG["input"]["source_csv"])
    zp = Table.read(REPO / "release/v2/phase1/catalogue_zpcorr.fits").to_pandas()
    dist = Table.read(REPO / "release/v2/phase1/catalogue_dist.fits").to_pandas()
    df_all = (src.merge(zp[["source_id", "parallax_zpcorr"]], on="source_id")
                  .merge(dist[["source_id", "dist_pc", "dist_lo_pc", "dist_hi_pc"]],
                         on="source_id"))

    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    pot_rot = agama.Potential(file=str(WORK / "MWPotentialHunter24_rot_default.ini"))

    rng = np.random.default_rng(SEED + 6700)
    results = []
    for cand in CANDIDATES:
        sid = cand["source_id"]
        log(f"--- candidate {sid} ({cand['potential']} Tier {cand['tier']}) ---")
        row = df_all[df_all["source_id"] == sid].iloc[0]
        ic = sample_ics(row, N_SAMP, rng)  # (5000, 6)
        log(f"  built {len(ic)} ICs; integrating dt=0.1 Myr trajsize=40001...")
        if cand["potential"] == "static":
            res = agama.orbit(potential=pot_axi, ic=ic, time=4.0 * GYR,
                               trajsize=40001)
            pot_used = "static"
        else:
            res = agama.orbit(potential=pot_rot, ic=ic, time=4.0 * GYR,
                               trajsize=40001, Omega=OMEGA_P)
            pot_used = "barred"
        log(f"  integration done in {time.time()-T0:.1f}s")
        min_rsph = np.zeros(N_SAMP)
        for k in range(N_SAMP):
            traj = np.asarray(res[k, 1])
            r_s = np.sqrt(np.sum(traj[:, :3]**2, axis=1))
            min_rsph[k] = r_s.min()
        P_lt_10pc = float(np.mean(min_rsph < 0.010))
        P_lt_100pc = float(np.mean(min_rsph < 0.100))
        results.append({
            "source_id": int(sid),
            "tier": cand["tier"],
            "potential": pot_used,
            "n_realisations": N_SAMP,
            "min_rsph_pc_p16": float(np.percentile(min_rsph, 16) * 1000),
            "min_rsph_pc_p50": float(np.percentile(min_rsph, 50) * 1000),
            "min_rsph_pc_p84": float(np.percentile(min_rsph, 84) * 1000),
            "min_rsph_pc_min": float(min_rsph.min() * 1000),
            "min_rsph_pc_max": float(min_rsph.max() * 1000),
            "P_lt_10pc": P_lt_10pc,
            "P_lt_100pc": P_lt_100pc,
            "P_lt_1kpc":  float(np.mean(min_rsph < 1.0)),
        })
        log(f"  P(<10pc)={P_lt_10pc:.3f}  P(<100pc)={P_lt_100pc:.3f}  "
            f"min_rsph p50={results[-1]['min_rsph_pc_p50']:.0f} pc")

    summary = {
        "n_candidates": len(CANDIDATES),
        "candidates": results,
        "n_confirmed_lt_10pc_P_gt_0p5":   sum(r["P_lt_10pc"] > 0.5 for r in results),
        "n_confirmed_lt_100pc_P_gt_0p5":  sum(r["P_lt_100pc"] > 0.5 for r in results),
    }
    (OUT / "sgrA_candidate_mc.json").write_text(json.dumps(summary, indent=2))
    log("wrote sgrA_candidate_mc.json")
    print(json.dumps(summary, indent=2))
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
