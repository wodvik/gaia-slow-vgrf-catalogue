"""Phase 5B — kernel-weighted control rebuild + per-band bar-resonance.

This script builds kernel-weighted controls with bootstrap weights and
computes the OLR/4:1 frequency-ratio accumulation per control band.

Steps
=====
1. Load the four control-band CSVs and the master orbit table.
2. For each band: integrate orbits in static_full (axi); compute
   Stäckel-fudge actions; compute Omega_R / (Omega_phi - Omega_p)
   resonance ratio.
3. Rebuild kernel-weighted control samples per band: each control star
   gets weight = product of Gaussian kernels in (log G, log d, sin|b|)
   centered on the slow-Vgrf catalogue's marginal distribution.
4. Compare per-band resonance fractions (slow vs each control).

Outputs
=======
release/v2/phase5/control_orbits.fits
release/v2/phase5/gate5B_controls.json
release/v2/phase5/gate5B_resonance_perband.png
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
OUT = REPO / "release/v2/phase5"
OUT.mkdir(parents=True, exist_ok=True)
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"

R0 = float(CONFIG["solar_variants"]["default"]["R0_kpc"])
V0 = float(CONFIG["solar_variants"]["default"]["Vc_kms"])
OMEGA_P = -float(CONFIG["bar_pattern_speeds_kms_kpc"]["default"])
GYR = 1.0 / 0.9778

agama.setUnits(length=1, mass=1, velocity=1)

BANDS = [
    ("vgrf_25_50",   "release/data/band_25_50_real0_enriched.csv"),
    ("vgrf_50_100",  "release/data/band_50_100_real0_enriched.csv"),
    ("vgrf_100_200", "release/data/band_100_200_real0_enriched.csv"),
    ("vgrf_200_260", "release/data/band_200_260_real0_enriched.csv"),
]


def log(msg):
    print(f"[5B t={time.time()-T0:5.1f}s] {msg}", flush=True)


def galactocentric_ic(df: pd.DataFrame) -> np.ndarray:
    """Build galactocentric (x,y,z,vx,vy,vz) from band CSV row at default solar."""
    sv = CONFIG["solar_variants"]["default"]
    frame = coord.Galactocentric(
        galcen_distance=sv["R0_kpc"] * u.kpc,
        z_sun=sv["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            sv["U_kms"] * u.km/u.s,
            (sv["Vc_kms"] + sv["V_kms"]) * u.km/u.s,
            sv["W_kms"] * u.km/u.s,
        ),
    )
    # Use the legacy distance_pc column for control-band initial conditions.
    icrs = coord.SkyCoord(
        ra=df["ra"].to_numpy() * u.deg,
        dec=df["dec"].to_numpy() * u.deg,
        distance=df["distance_pc"].to_numpy() * u.pc,
        pm_ra_cosdec=df["pmra"].to_numpy() * u.mas/u.yr,
        pm_dec=df["pmdec"].to_numpy() * u.mas/u.yr,
        radial_velocity=df["radial_velocity"].to_numpy() * u.km/u.s,
        frame="icrs",
    )
    g = icrs.transform_to(frame)
    return np.column_stack([g.x.to_value(u.kpc), g.y.to_value(u.kpc), g.z.to_value(u.kpc),
                            g.v_x.to_value(u.km/u.s), g.v_y.to_value(u.km/u.s), g.v_z.to_value(u.km/u.s)])


def main():
    global T0; T0 = time.time()
    log(f"agama {agama.__version__}")

    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    af = agama.ActionFinder(pot_axi, interp=True)

    rows_all = []
    for label, path in BANDS:
        log(f"loading {label} from {path}")
        df = pd.read_csv(REPO / path)
        finite = np.all(np.isfinite(df[["ra", "dec", "pmra", "pmdec",
                                         "radial_velocity", "distance_pc"]].to_numpy(dtype=float)),
                        axis=1)
        df = df.loc[finite].reset_index(drop=True)
        log(f"  {len(df)} stars after finite-mask")

        ic = galactocentric_ic(df)
        log(f"  integrating {len(ic)} orbits in static_full ...")
        t0 = time.time()
        res = agama.orbit(potential=pot_axi, ic=ic, time=4.0 * GYR, trajsize=2001)
        log(f"  agama.orbit done in {time.time()-t0:.1f}s")

        log(f"  computing actions/freqs ...")
        actions, _, freqs = af(ic, angles=True)
        Omega_R = freqs[:, 0]; Omega_phi = freqs[:, 2]
        delta = Omega_phi - OMEGA_P
        ratio = np.where(np.abs(delta) > 1e-3, Omega_R / delta, np.nan)

        # Per-orbit summaries
        for i in range(len(ic)):
            traj = np.asarray(res[i, 1])
            Rcyl = np.sqrt(traj[:, 0]**2 + traj[:, 1]**2)
            r_sph = np.sqrt(np.sum(traj[:, :3]**2, axis=1))
            rows_all.append({
                "band": label,
                "source_id": int(df["source_id"].iloc[i]),
                "ra": float(df["ra"].iloc[i]),
                "dec": float(df["dec"].iloc[i]),
                "G": float(df["phot_g_mean_mag"].iloc[i]),
                "dist_pc_v1": float(df["distance_pc"].iloc[i]),
                "b_deg": float(df["b"].iloc[i]),
                "vgrf_v1": float(df["V_grf"].iloc[i]),
                "R_peri_kpc": float(Rcyl.min()),
                "R_apo_kpc":  float(Rcyl.max()),
                "z_max_kpc":  float(np.abs(traj[:, 2]).max()),
                "min_r_sph_kpc": float(r_sph.min()),
                "ecc": float((Rcyl.max() - Rcyl.min()) / (Rcyl.max() + Rcyl.min())),
                "J_R": float(actions[i, 0]),
                "J_z": float(actions[i, 1]),
                "J_phi": float(actions[i, 2]),
                "Omega_R":   float(Omega_R[i]),
                "Omega_phi": float(Omega_phi[i]),
                "res_ratio": float(ratio[i]),
            })
    full = pd.DataFrame(rows_all)
    log(f"total control orbits: {len(full)}")
    Table.from_pandas(full).write(OUT / "control_orbits.fits", overwrite=True)
    log("wrote control_orbits.fits")

    # ---- Kernel-weighted control: per-band weight per row that brings
    # band's (G, log10(d), sin|b|) marginals into match with slow sample.
    log("computing kernel weights against slow sample marginals ...")
    slow = pd.read_csv(REPO / CONFIG["input"]["source_csv"])
    slow_G = slow["phot_g_mean_mag"].to_numpy()
    slow_logd = np.log10(slow["distance_pc"].to_numpy())
    slow_sb   = np.abs(np.sin(np.deg2rad(slow["b"].to_numpy())))
    # Bandwidths (Silverman-ish for small samples):
    h_G = 0.3; h_logd = 0.15; h_sb = 0.1

    def kernel_weight(G, logd, sb):
        # Gaussian product kernel summed over slow-sample anchors
        wG = np.exp(-0.5 * ((G - slow_G[:, None]) / h_G)**2).sum(axis=0)
        wD = np.exp(-0.5 * ((logd - slow_logd[:, None]) / h_logd)**2).sum(axis=0)
        wB = np.exp(-0.5 * ((sb - slow_sb[:, None]) / h_sb)**2).sum(axis=0)
        return wG * wD * wB

    weights = []
    for label, _ in BANDS:
        m = (full["band"] == label).to_numpy()
        if not m.any():
            weights.append(None)
            continue
        sub = full.loc[m]
        w = kernel_weight(sub["G"].to_numpy(),
                          np.log10(np.maximum(sub["dist_pc_v1"].to_numpy(), 1.0)),
                          np.abs(np.sin(np.deg2rad(sub["b_deg"].to_numpy()))))
        # normalize so weights sum to N (effective sample size unchanged at limit)
        w = w / w.mean()
        full.loc[m, "kw_weight"] = w
        weights.append(float(w.sum()))

    Table.from_pandas(full).write(OUT / "control_orbits.fits", overwrite=True)

    # ---- per-band resonance summary (raw + weighted) ----
    summary = {"omega_p_kms_kpc": OMEGA_P, "n_total": int(len(full)),
               "bands": {}}
    for label, _ in BANDS:
        m = (full["band"] == label).to_numpy()
        if not m.any(): continue
        r = full.loc[m, "res_ratio"].to_numpy()
        w = full.loc[m, "kw_weight"].to_numpy() if "kw_weight" in full else np.ones(m.sum())
        ok = np.isfinite(r) & np.isfinite(w)
        r = r[ok]; w = w[ok]
        def frac(mask):
            if not len(w): return None, None
            return (float(mask.sum() / len(r)),
                    float((w[mask]).sum() / w.sum()))
        n_OLR_raw, n_OLR_w = frac(np.abs(r - 2.0) < 0.3)
        n_41_raw, n_41_w   = frac(np.abs(r - 4.0) < 0.3)
        n_ILR_raw, n_ILR_w = frac(np.abs(r + 2.0) < 0.3)
        summary["bands"][label] = {
            "n":       int(m.sum()),
            "n_finite": int(len(r)),
            "frac_within_0p3_OLR_raw":      n_OLR_raw,
            "frac_within_0p3_OLR_weighted": n_OLR_w,
            "frac_within_0p3_4_1_raw":      n_41_raw,
            "frac_within_0p3_4_1_weighted": n_41_w,
            "frac_within_0p3_ILR_raw":      n_ILR_raw,
            "frac_within_0p3_ILR_weighted": n_ILR_w,
            "median_R_peri_pc":   float(np.median(full.loc[m, "R_peri_kpc"]) * 1000),
            "median_z_max_pc":    float(np.median(full.loc[m, "z_max_kpc"]) * 1000),
            "median_ecc":         float(np.median(full.loc[m, "ecc"])),
        }

    # Comparison: slow-Vgrf v2 from Phase 4
    v2 = Table.read(REPO / "release/v2/phase4/catalogue_v2_orbits.fits").to_pandas()
    if v2["tier"].dtype == object and isinstance(v2["tier"].iloc[0], (bytes, bytearray)):
        v2["tier"] = v2["tier"].str.decode("utf-8")
    tABC = v2["tier"].isin(["A", "B", "C"]).to_numpy()
    rs = v2.loc[tABC, "res_ratio_OmegaR_over_dPhi"].to_numpy()
    rs = rs[np.isfinite(rs)]
    summary["slow_TierABC"] = {
        "n": int(len(rs)),
        "frac_within_0p3_OLR": float(np.sum(np.abs(rs - 2.0) < 0.3) / len(rs)),
        "frac_within_0p3_4_1": float(np.sum(np.abs(rs - 4.0) < 0.3) / len(rs)),
        "frac_within_0p3_ILR": float(np.sum(np.abs(rs + 2.0) < 0.3) / len(rs)),
    }

    (OUT / "gate5B_controls.json").write_text(json.dumps(summary, indent=2))
    log("wrote gate5B_controls.json")
    print(json.dumps(summary, indent=2))

    # ---- Plot per-band resonance ratio histograms ----
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=False)
    bins = np.linspace(-6, 6, 30)
    for ax, label in zip(axes.flat, [b[0] for b in BANDS] + ["slow_TierABC"]):
        if label == "slow_TierABC":
            r = rs
            n_for_text = len(r)
        else:
            m = (full["band"] == label).to_numpy()
            r = full.loc[m, "res_ratio"].to_numpy()
            r = r[np.isfinite(r)]
            n_for_text = len(r)
        ax.hist(r, bins=bins, color="steelblue", edgecolor="white")
        for x_, lab in [(-2,"ILR"),(0,"CR"),(2,"OLR"),(4,"4:1")]:
            ax.axvline(x_, ls="--", color="black", alpha=0.5)
        ax.set_title(f"{label} (n={n_for_text})")
        ax.set_xlabel(r"$\Omega_R / (\Omega_\phi - \Omega_p)$")
    fig.tight_layout()
    fig.savefig(OUT / "gate5B_resonance_perband.png", dpi=140); plt.close(fig)
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
