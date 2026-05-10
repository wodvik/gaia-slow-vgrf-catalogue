"""Phase 6C — solar parameter sensitivity (full-sample re-tier and re-orbit).

For each of the four solar variants from config.yml, re-derive per-star
Vgrf using ZP-corrected parallax + Bailer-Jones distance, re-tier, and
integrate Tier A+B+C in pot_static_full at the *appropriate* solar
parameters for the Galactocentric frame.

Note: the potential pot_static_full was built with R0=8.178, V0=232 (the
default). When evaluating in a different solar frame, the orbit
geometry shifts because the ICs (galactocentric XYZ/VxVyVz) are
recomputed under the new frame. We do NOT rebuild the potential per
variant — that would be a different test (halo-mass effect).
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
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"

VGRF_CUT = float(CONFIG["vgrf_cutoff_kms"])
TIERS = CONFIG["tiers"]
SEED = int(CONFIG["mc"]["random_seed"])
GYR = 1.0 / 0.9778
agama.setUnits(length=1, mass=1, velocity=1)


def log(m): print(f"[6C t={time.time()-T0:5.1f}s] {m}", flush=True)


def galcen_frame(name, params):
    if name == "rb20":
        Vy = params["Vc_kms"]
    else:
        Vy = params["Vc_kms"] + params["V_kms"]
    return coord.Galactocentric(
        galcen_distance=params["R0_kpc"] * u.kpc,
        z_sun=params["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            params["U_kms"] * u.km/u.s,
            Vy * u.km/u.s,
            params["W_kms"] * u.km/u.s,
        ),
    )


def assign_tier(P, point_below):
    tier = np.where(point_below, "D", "X").astype("U2")
    tier[P > TIERS["C"]] = "C"
    tier[P > TIERS["B"]] = "B"
    tier[P > TIERS["A"]] = "A"
    tier[(~point_below) & (P <= TIERS["C"])] = "X"
    return tier


def cov3(plx_e, pmra_e, pmdec_e, p_pmra, p_pmdec, pmra_pmdec):
    n = len(plx_e); sig = np.stack([plx_e, pmra_e, pmdec_e], axis=-1)
    rho = np.zeros((n, 3, 3))
    for i in (0, 1, 2): rho[:, i, i] = 1.0
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


def variant_run(df, variant_name, variant_params, n_samp, rng):
    """Compute P(Vgrf<25) per star + point Vgrf + galactocentric ICs under
    the given solar variant."""
    n = len(df)
    galcen = galcen_frame(variant_name, variant_params)
    plx = df["parallax_zpcorr"].to_numpy()
    plx_e = df["parallax_error"].to_numpy()
    mu = np.stack([plx, df["pmra"].to_numpy(), df["pmdec"].to_numpy()], axis=-1)
    C = cov3(plx_e, df["pmra_error"].to_numpy(), df["pmdec_error"].to_numpy(),
             df["parallax_pmra_corr"].to_numpy(),
             df["parallax_pmdec_corr"].to_numpy(),
             df["pmra_pmdec_corr"].to_numpy())
    L = chol(C)
    z = rng.standard_normal((n, n_samp, 3))
    s = mu[:, None, :] + np.einsum("nij,nsj->nsi", L, z)
    plx_s = s[:, :, 0]; pmra_s = s[:, :, 1]; pmdec_s = s[:, :, 2]
    rv_s = (df["radial_velocity"].to_numpy()[:, None] +
            rng.standard_normal((n, n_samp)) *
            df["radial_velocity_error"].to_numpy()[:, None])
    dist_s = split_normal(df["dist_pc"].to_numpy(),
                          df["dist_lo_pc"].to_numpy(),
                          df["dist_hi_pc"].to_numpy(), n_samp, rng)
    dist_s = np.maximum(dist_s, 1.0)
    ra = np.broadcast_to(df["ra"].to_numpy()[:, None], (n, n_samp)).flatten()
    dec = np.broadcast_to(df["dec"].to_numpy()[:, None], (n, n_samp)).flatten()
    icrs = coord.SkyCoord(ra=ra * u.deg, dec=dec * u.deg,
                          distance=dist_s.flatten() * u.pc,
                          pm_ra_cosdec=pmra_s.flatten() * u.mas/u.yr,
                          pm_dec=pmdec_s.flatten() * u.mas/u.yr,
                          radial_velocity=rv_s.flatten() * u.km/u.s, frame="icrs")
    g = icrs.transform_to(galcen)
    vgrf = np.sqrt(g.v_x.to_value(u.km/u.s)**2 +
                    g.v_y.to_value(u.km/u.s)**2 +
                    g.v_z.to_value(u.km/u.s)**2).reshape(n, n_samp)
    P = np.mean(vgrf < VGRF_CUT, axis=1)
    # Point estimate ICs (no MC)
    icrs_pt = coord.SkyCoord(
        ra=df["ra"].to_numpy() * u.deg, dec=df["dec"].to_numpy() * u.deg,
        distance=df["dist_pc"].to_numpy() * u.pc,
        pm_ra_cosdec=df["pmra"].to_numpy() * u.mas/u.yr,
        pm_dec=df["pmdec"].to_numpy() * u.mas/u.yr,
        radial_velocity=df["radial_velocity"].to_numpy() * u.km/u.s, frame="icrs"
    )
    g_pt = icrs_pt.transform_to(galcen)
    ic = np.column_stack([
        g_pt.x.to_value(u.kpc), g_pt.y.to_value(u.kpc), g_pt.z.to_value(u.kpc),
        g_pt.v_x.to_value(u.km/u.s), g_pt.v_y.to_value(u.km/u.s),
        g_pt.v_z.to_value(u.km/u.s),
    ])
    vgrf_pt = np.sqrt(np.sum(ic[:, 3:6]**2, axis=1))
    return P, vgrf_pt, ic


def main():
    global T0; T0 = time.time()
    log(f"agama {agama.__version__}")

    src = pd.read_csv(REPO / CONFIG["input"]["source_csv"])
    zp = Table.read(REPO / "release/v2/phase1/catalogue_zpcorr.fits").to_pandas()
    dist = Table.read(REPO / "release/v2/phase1/catalogue_dist.fits").to_pandas()
    df = (src.merge(zp[["source_id", "parallax_zpcorr"]], on="source_id")
              .merge(dist[["source_id", "dist_pc", "dist_lo_pc", "dist_hi_pc"]],
                     on="source_id"))
    log(f"merged {len(df)} stars")

    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))

    rows = []
    for vname, vparams in CONFIG["solar_variants"].items():
        log(f"--- variant {vname} ---")
        rng = np.random.default_rng(SEED + abs(hash(vname)) % 1_000_000)
        P, vgrf_pt, ic = variant_run(df, vname, vparams, 5000, rng)
        point_below = vgrf_pt < VGRF_CUT
        tier = assign_tier(P, point_below)

        # Tier A+B+C orbits
        tABC_mask = np.isin(tier, ["A", "B", "C"])
        log(f"  Tier counts: A={np.sum(tier=='A')}, B={np.sum(tier=='B')}, "
            f"C={np.sum(tier=='C')}, D={np.sum(tier=='D')}, X={np.sum(tier=='X')}")
        ic_ABC = ic[tABC_mask]
        log(f"  integrating {len(ic_ABC)} Tier A+B+C orbits in static_full ...")
        t0 = time.time()
        res = agama.orbit(potential=pot_axi, ic=ic_ABC, time=4.0 * GYR,
                           trajsize=1001)
        log(f"  agama.orbit done in {time.time()-t0:.1f}s")
        Rp = np.zeros(len(ic_ABC)); Ra = np.zeros(len(ic_ABC))
        zm = np.zeros(len(ic_ABC)); rsph = np.zeros(len(ic_ABC))
        for k in range(len(ic_ABC)):
            traj = np.asarray(res[k, 1])
            Rcyl = np.sqrt(traj[:, 0]**2 + traj[:, 1]**2)
            r_s = np.sqrt(Rcyl**2 + traj[:, 2]**2)
            Rp[k] = Rcyl.min(); Ra[k] = Rcyl.max()
            zm[k] = np.abs(traj[:, 2]).max(); rsph[k] = r_s.min()
        ecc = (Ra - Rp) / (Ra + Rp)
        rows.append({
            "variant": vname,
            **{k: vparams[k] for k in ("R0_kpc","z_sun_pc","Vc_kms","U_kms","V_kms","W_kms")},
            "n_TierA": int(np.sum(tier == "A")),
            "n_TierB": int(np.sum(tier == "B")),
            "n_TierC": int(np.sum(tier == "C")),
            "n_TierD": int(np.sum(tier == "D")),
            "median_vgrf_kms": float(np.nanmedian(vgrf_pt)),
            "median_R_peri_pc_TierABC":  float(np.median(Rp) * 1000),
            "median_R_apo_kpc_TierABC":  float(np.median(Ra)),
            "median_z_max_kpc_TierABC":  float(np.median(zm)),
            "median_ecc_TierABC":        float(np.median(ecc)),
            "n_R_peri_lt_100pc_TierABC": int(np.sum(Rp < 0.100)),
            "n_min_rsph_lt_10pc_TierABC": int(np.sum(rsph < 0.010)),
        })
    pd.DataFrame(rows).to_csv(OUT / "solar_variant_table.csv", index=False)
    print(json.dumps(rows, indent=2))
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
