"""Phase 4 — point-estimate orbit integration + Stackel actions + bar
resonance gating.

Inputs
------
- release/v2/phase1/catalogue_v2.fits     -- master Phase 1 table
- release/v2/phase1/catalogue_vgrf.fits   -- per-star (x,y,z,vx,vy,vz)
                                             under default solar
- release/v2/phase3_agama/_hunter24_workdir/MWPotentialHunter24_axi.ini
- release/v2/phase3_agama/_hunter24_workdir/MWPotentialHunter24_rot_default.ini

Compute (point-estimate, no MC)
-------------------------------
For each of 2,859 stars, in static_full (axi) and barred_default (rotating
@ Omega_p = -37.5 km/s/kpc):
  R_peri, R_apo, |z|_max, eccentricity, min r_sph, n_peri (peri count over 4 Gyr)
  energy_conservation_rel  (static)
  jacobi_conservation_rel  (barred — E_J = E - Omega * L_z)
Stackel-fudge actions on static_full: J_R, J_z, J_phi, Omega_R, Omega_phi, Omega_z
Bar-resonance frequency ratio: Omega_R / (Omega_phi - Omega_p), histogrammed
to look for accumulation at OLR (=2), CR (->inf), and the 4:1 ultra-harmonic.

Outputs
-------
- release/v2/phase4/catalogue_v2_orbits.fits
- release/v2/phase4/gate4.json
- release/v2/phase4/gate4_rperi_hist.png
- release/v2/phase4/gate4_resonance_ratio.png
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
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase4"
OUT.mkdir(parents=True, exist_ok=True)
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"

R0 = float(CONFIG["solar_variants"]["default"]["R0_kpc"])
V0 = float(CONFIG["solar_variants"]["default"]["Vc_kms"])
OMEGA_P_DEFAULT = -float(CONFIG["bar_pattern_speeds_kms_kpc"]["default"])  # -37.5
GYR_TO_AGAMA = 1.0 / 0.9778  # 1 Gyr in agama time units (kpc/(km/s))
T_INTEGRATE = 4.0 * GYR_TO_AGAMA   # 4 Gyr forward
TRAJ_STEPS = 2001                  # ~2 Myr resolution; sufficient for peri/apo

agama.setUnits(length=1, mass=1, velocity=1)


def log(msg: str) -> None:
    print(f"[phase4 t={time.time()-T0:6.1f}s] {msg}", flush=True)


def load_inputs():
    log("loading catalogue + ICs ...")
    src = pd.read_csv(REPO / CONFIG["input"]["source_csv"])
    v2 = Table.read(REPO / "release/v2/phase1/catalogue_v2.fits").to_pandas()
    if v2["tier"].dtype == object and isinstance(v2["tier"].iloc[0],
                                                 (bytes, bytearray)):
        v2["tier"] = v2["tier"].str.decode("utf-8")
    vgrf = Table.read(REPO / "release/v2/phase1/catalogue_vgrf.fits").to_pandas()

    # vgrf already carries vgrf_default; pull only the non-duplicate cols from v2.
    df = (vgrf.merge(v2[["source_id", "tier", "P_vgrf_below_25",
                         "rv_quality"]],
                     on="source_id"))
    log(f"merged {len(df)} stars (tier counts: "
        f"{df['tier'].value_counts().to_dict()})")
    ic = np.column_stack([
        df["x_kpc"].to_numpy(),
        df["y_kpc"].to_numpy(),
        df["z_kpc"].to_numpy(),
        df["vx_kms"].to_numpy(),
        df["vy_kms"].to_numpy(),
        df["vz_kms"].to_numpy(),
    ])
    return df, ic


def orbit_summary(times: np.ndarray, traj: np.ndarray, pot, omega_p: float):
    """Reduce one orbit trajectory to scalar summary stats."""
    x, y, z, vx, vy, vz = (traj[:, i] for i in range(6))
    R_cyl = np.sqrt(x*x + y*y)
    r_sph = np.sqrt(x*x + y*y + z*z)
    Lz = x*vy - y*vx                        # specific z-angular momentum
    KE = 0.5 * (vx*vx + vy*vy + vz*vz)
    PE = pot.potential(traj[:, :3])
    E = KE + PE
    E_J = E - omega_p * Lz                  # Jacobi (= E if omega_p=0)

    # Find peri count (zero-crossings of dR/dt).
    dR = np.diff(R_cyl)
    n_peri = int(np.sum((dR[:-1] < 0) & (dR[1:] > 0)))

    return {
        "R_peri_kpc": float(R_cyl.min()),
        "R_apo_kpc":  float(R_cyl.max()),
        "z_max_kpc":  float(np.abs(z).max()),
        "min_r_sph_kpc": float(r_sph.min()),
        "ecc": float((R_cyl.max() - R_cyl.min()) / (R_cyl.max() + R_cyl.min())),
        "n_peri": n_peri,
        "E_mean":         float(np.mean(E)),
        "E_drift_rel":    float(abs((E[-1] - E[0]) / E[0])),
        "E_range_rel":    float((E.max() - E.min()) / abs(E.mean())),
        "EJ_drift_rel":   float(abs((E_J[-1] - E_J[0]) / E_J[0])),
        "EJ_range_rel":   float((E_J.max() - E_J.min()) / abs(E_J.mean())),
        "Lz_kpc_kms":     float(np.mean(Lz)),
    }


def integrate_one_potential(ic, pot, omega_p, label):
    log(f"integrating {len(ic)} orbits in {label} (Omega={omega_p}) "
        f"for {T_INTEGRATE*0.9778:.1f} Gyr ({TRAJ_STEPS} steps)...")
    t0 = time.time()
    kwargs = dict(potential=pot, ic=ic, time=T_INTEGRATE,
                  trajsize=TRAJ_STEPS)
    if omega_p != 0.0:
        kwargs["Omega"] = omega_p
    result = agama.orbit(**kwargs)
    log(f"  agama.orbit returned shape {np.asarray(result).shape} in "
        f"{time.time()-t0:.1f}s")
    n = len(ic)
    rows = []
    for i in range(n):
        times = np.asarray(result[i, 0])
        traj  = np.asarray(result[i, 1])
        s = orbit_summary(times, traj, pot, omega_p)
        s["star_idx"] = i
        rows.append(s)
    return pd.DataFrame(rows)


def stackel_actions(ic, pot_axi):
    log(f"computing Stackel-fudge actions for {len(ic)} stars in static_full...")
    t0 = time.time()
    af = agama.ActionFinder(pot_axi, interp=True)
    # ActionFinder returns (J_R, J_z, J_phi, Omega_R, Omega_z, Omega_phi)
    actions, angles, freqs = af(ic, angles=True)
    log(f"  done in {time.time()-t0:.1f}s")
    return pd.DataFrame({
        "J_R":     actions[:, 0],
        "J_z":     actions[:, 1],
        "J_phi":   actions[:, 2],
        "Omega_R":   freqs[:, 0],
        "Omega_z":   freqs[:, 1],
        "Omega_phi": freqs[:, 2],
    })


def main():
    global T0
    T0 = time.time()
    log(f"agama {agama.__version__}")

    df, ic = load_inputs()

    # --- Load potentials ---
    log("loading Hunter+2024 axi + rotating potentials ...")
    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    pot_rot = agama.Potential(file=str(WORK / "MWPotentialHunter24_rot_default.ini"))

    # --- 4A: orbits in static_full ---
    orb_static = integrate_one_potential(ic, pot_axi, 0.0, "static_full")
    orb_static = orb_static.add_prefix("static_")

    # --- 4B: orbits in barred_default ---
    orb_barred = integrate_one_potential(ic, pot_rot, OMEGA_P_DEFAULT,
                                         "barred_default")
    orb_barred = orb_barred.add_prefix("barred_")

    # --- 4C: Stackel-fudge actions (static_full) ---
    try:
        actions = stackel_actions(ic, pot_axi)
    except Exception as e:
        log(f"WARN: action computation failed: {e}; skipping 4C")
        actions = pd.DataFrame({k: np.full(len(ic), np.nan)
                                for k in ("J_R", "J_z", "J_phi",
                                          "Omega_R", "Omega_z", "Omega_phi")})

    # --- 4E: bar-resonance frequency ratio ---
    # Resonance condition (in static frequencies):
    #   bar resonance at  m*(Omega_phi - Omega_p) - n*Omega_R = 0
    #   CR    : Omega_phi = Omega_p          (m=1, n=0)
    #   OLR   : Omega_R = 2*(Omega_phi - Omega_p)  (m=2, n=2)
    #   ILR   : Omega_R = -2*(Omega_phi - Omega_p) (m=2, n=2, opposite sign)
    #   4:1   : Omega_R = 4*(Omega_phi - Omega_p)  ultraharmonic
    # We compute the rate ratio  Omega_R / (Omega_phi - Omega_p)
    om_p_aga = OMEGA_P_DEFAULT  # -37.5 km/s/kpc; agama freqs are in same units.
    delta = actions["Omega_phi"] - om_p_aga
    ratio = np.where(np.abs(delta) > 1e-3,
                     actions["Omega_R"] / delta, np.nan)
    actions["res_ratio_OmegaR_over_dPhi"] = ratio

    # --- assemble per-star master orbit table ---
    log("assembling master output ...")
    df = df.reset_index(drop=True)
    out = pd.concat([
        df[["source_id", "tier", "P_vgrf_below_25", "vgrf_default",
            "rv_quality", "x_kpc", "y_kpc", "z_kpc",
            "vx_kms", "vy_kms", "vz_kms"]].reset_index(drop=True),
        orb_static.reset_index(drop=True),
        orb_barred.reset_index(drop=True),
        actions.reset_index(drop=True),
    ], axis=1)
    # de-duplicate any star_idx columns from prefix
    out = out.loc[:, ~out.columns.duplicated()]
    tab = Table.from_pandas(out)
    tab.write(OUT / "catalogue_v2_orbits.fits", overwrite=True)
    log(f"wrote {OUT / 'catalogue_v2_orbits.fits'}")

    # --- Gate 4 summary ---
    tier_AB = out["tier"].isin(["A", "B"])
    tier_ABC = out["tier"].isin(["A", "B", "C"])

    def stat(mask, col):
        v = out.loc[mask, col].to_numpy()
        v = v[np.isfinite(v)]
        if not len(v):
            return None
        return {"n": int(len(v)),
                "p16": float(np.percentile(v, 16)),
                "p50": float(np.percentile(v, 50)),
                "p84": float(np.percentile(v, 84)),
                "min": float(np.min(v)),
                "max": float(np.max(v))}

    summary = {
        "n_processed": int(len(out)),
        "Omega_p_kms_kpc": float(om_p_aga),
        "T_integrate_Gyr": 4.0,
        "headline": {
            "median_R_peri_static_pc_TierAB":  stat(tier_AB,  "static_R_peri_kpc")["p50"] * 1000,
            "median_R_peri_static_pc_TierABC": stat(tier_ABC, "static_R_peri_kpc")["p50"] * 1000,
            "median_R_peri_barred_pc_TierABC": stat(tier_ABC, "barred_R_peri_kpc")["p50"] * 1000,
            "v1_median_R_peri_pc": 55,
            "n_R_peri_lt_100pc_static_TierABC":  int(((out["static_R_peri_kpc"] < 0.100) & tier_ABC).sum()),
            "n_R_peri_lt_100pc_barred_TierABC":  int(((out["barred_R_peri_kpc"] < 0.100) & tier_ABC).sum()),
            "n_min_r_sph_lt_10pc_static_TierABC": int(((out["static_min_r_sph_kpc"] < 0.010) & tier_ABC).sum()),
            "n_min_r_sph_lt_10pc_barred_TierABC": int(((out["barred_min_r_sph_kpc"] < 0.010) & tier_ABC).sum()),
            "n_bridgers_static_TierABC":  int(((out["static_R_peri_kpc"] < 2.0) & (out["static_R_apo_kpc"] > 15.0) & tier_ABC).sum()),
            "n_bridgers_barred_TierABC":  int(((out["barred_R_peri_kpc"] < 2.0) & (out["barred_R_apo_kpc"] > 15.0) & tier_ABC).sum()),
        },
        "energy_conservation": {
            "static_E_drift_p84_TierABC":   stat(tier_ABC, "static_E_drift_rel")["p84"],
            "barred_E_drift_p84_TierABC":   stat(tier_ABC, "barred_E_drift_rel")["p84"],
            "barred_EJ_drift_p84_TierABC":  stat(tier_ABC, "barred_EJ_drift_rel")["p84"],
        },
        "stat_static_R_peri_kpc": {"all": stat(np.ones(len(out), dtype=bool), "static_R_peri_kpc"),
                                   "tierAB": stat(tier_AB, "static_R_peri_kpc"),
                                   "tierABC": stat(tier_ABC, "static_R_peri_kpc")},
        "stat_static_R_apo_kpc": {"tierABC": stat(tier_ABC, "static_R_apo_kpc")},
        "stat_static_z_max_kpc": {"tierABC": stat(tier_ABC, "static_z_max_kpc")},
        "stat_static_ecc":       {"tierABC": stat(tier_ABC, "static_ecc")},
        "stat_barred_R_peri_kpc": {"tierABC": stat(tier_ABC, "barred_R_peri_kpc")},
        "stat_barred_R_apo_kpc":  {"tierABC": stat(tier_ABC, "barred_R_apo_kpc")},
    }

    # --- 4E: resonance accumulation evidence ---
    r = out.loc[tier_ABC, "res_ratio_OmegaR_over_dPhi"].to_numpy()
    r = r[np.isfinite(r)]
    summary["bar_resonance"] = {
        "n_TierABC": int(len(r)),
        "histogram_bins":   list(np.linspace(-6, 6, 25)),
        "histogram_counts": np.histogram(r, bins=np.linspace(-6, 6, 25))[0].tolist(),
        "n_within_0p3_of_OLR_(ratio_2)":  int(np.sum(np.abs(r - 2.0) < 0.3)),
        "n_within_0p3_of_4_1_(ratio_4)":  int(np.sum(np.abs(r - 4.0) < 0.3)),
        "n_within_0p3_of_minus2_ILR":     int(np.sum(np.abs(r + 2.0) < 0.3)),
        "n_high_CR_like_ratio_lt_0p1":    int(np.sum(np.abs(r) < 0.1)),
    }

    (OUT / "gate4.json").write_text(json.dumps(summary, indent=2))
    log("wrote gate4.json")
    print(json.dumps(summary["headline"], indent=2))
    print(json.dumps(summary["energy_conservation"], indent=2))
    print(json.dumps(summary["bar_resonance"], indent=2))

    # --- Plots ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.logspace(-3, 1.5, 60)  # 1 pc to ~30 kpc
    for label, col, color in [
        ("static_full Tier A+B+C",  "static_R_peri_kpc", "C0"),
        ("barred_default Tier A+B+C", "barred_R_peri_kpc", "C3"),
    ]:
        v = out.loc[tier_ABC, col].to_numpy()
        v = v[np.isfinite(v) & (v > 0)]
        ax.hist(v, bins=bins, alpha=0.6, label=label, color=color)
    ax.axvline(0.055, ls="--", color="black", label="reference median 55 pc")
    ax.set_xscale("log"); ax.set_xlabel("R_peri (kpc)"); ax.set_ylabel("count")
    ax.set_title("Phase 4 — R_peri distribution (Tier A+B+C)")
    ax.legend(); fig.tight_layout()
    fig.savefig(OUT / "gate4_rperi_hist.png", dpi=140); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(r, bins=np.linspace(-6, 6, 60), color="steelblue", edgecolor="white")
    for x_, lab in [(2, "OLR (2:1)"), (4, "4:1 UH"), (-2, "ILR"), (0, "CR")]:
        ax.axvline(x_, ls="--", color="black", alpha=0.6)
        ax.text(x_, ax.get_ylim()[1]*0.92, lab, rotation=90,
                ha="right", va="top", fontsize=8)
    ax.set_xlabel(r"$\Omega_R / (\Omega_\phi - \Omega_p)$")
    ax.set_ylabel("count (Tier A+B+C)")
    ax.set_title("Phase 4 — bar-resonance ratio distribution")
    fig.tight_layout()
    fig.savefig(OUT / "gate4_resonance_ratio.png", dpi=140); plt.close(fig)
    log(f"DONE in {time.time()-T0:.1f}s")


if __name__ == "__main__":
    sys.exit(main() or 0)
