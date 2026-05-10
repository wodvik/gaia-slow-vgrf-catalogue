"""Phase 6B — re-integrate barred_default at dt=0.1 Myr.

Phase 4B used trajsize=2001 (~2 Myr/step). Now re-run at trajsize=40001
(0.1 Myr/step over 4 Gyr).

Targets:
  - Drop Jacobi-drift p84 to <=1e-5 (or document why not)
  - Resolve whether the 10 surviving bridgers are real or step-artifact
  - Resolve whether N(R_peri<100pc) = 423 holds at finer dt
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
OUT = REPO / "release/v2/phase6"
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"
OMEGA_P = -float(CONFIG["bar_pattern_speeds_kms_kpc"]["default"])
GYR = 1.0 / 0.9778
agama.setUnits(length=1, mass=1, velocity=1)


def log(m):
    print(f"[6B t={time.time()-T0:6.1f}s] {m}", flush=True)


def main():
    global T0; T0 = time.time()
    log(f"agama {agama.__version__}")

    log("loading Phase 4 orbit catalogue...")
    o4 = Table.read(REPO / "release/v2/phase4/catalogue_v2_orbits.fits").to_pandas()
    if o4["tier"].dtype == object and isinstance(o4["tier"].iloc[0], (bytes, bytearray)):
        o4["tier"] = o4["tier"].str.decode("utf-8")
    n = len(o4)
    log(f"{n} stars from Phase 4")

    pot_rot = agama.Potential(file=str(WORK / "MWPotentialHunter24_rot_default.ini"))
    ic = o4[["x_kpc", "y_kpc", "z_kpc",
             "vx_kms", "vy_kms", "vz_kms"]].to_numpy()

    # 0.1 Myr step over 4 Gyr -> 40001 steps. Per-chunk memory at 150 stars:
    # 150 × 40001 × 6 × 8 = 2.9 GB transient.
    chunk = 150
    rows = []
    for i in range(0, n, chunk):
        ic_c = ic[i:i+chunk]
        log(f"chunk {i}-{i+len(ic_c)} ({len(ic_c)} stars) at dt=0.1 Myr...")
        t0 = time.time()
        res = agama.orbit(potential=pot_rot, ic=ic_c, time=4.0 * GYR,
                           trajsize=40001, Omega=OMEGA_P)
        log(f"  agama.orbit done in {time.time()-t0:.1f}s")
        for k in range(len(ic_c)):
            traj = np.asarray(res[k, 1])
            x, y, z, vx, vy, vz = (traj[:, j] for j in range(6))
            Rcyl = np.sqrt(x*x + y*y); rsph = np.sqrt(Rcyl**2 + z*z)
            Lz = x*vy - y*vx
            KE = 0.5*(vx*vx + vy*vy + vz*vz)
            PE = pot_rot.potential(traj[:, :3])
            E = KE + PE
            E_J = E - OMEGA_P * Lz
            rows.append({
                "star_idx": i + k,
                "source_id": int(o4["source_id"].iloc[i + k]),
                "tier": str(o4["tier"].iloc[i + k]),
                "fine_R_peri_kpc": float(Rcyl.min()),
                "fine_R_apo_kpc":  float(Rcyl.max()),
                "fine_z_max_kpc":  float(np.abs(z).max()),
                "fine_min_r_sph_kpc": float(rsph.min()),
                "fine_ecc": float((Rcyl.max() - Rcyl.min()) / (Rcyl.max() + Rcyl.min())),
                "fine_EJ_drift_rel": float(abs((E_J[-1] - E_J[0]) / E_J[0])),
                "fine_EJ_range_rel": float((E_J.max() - E_J.min()) / abs(E_J.mean())),
                # Original Phase 4 values (2 Myr step) for diff
                "old_R_peri_kpc": float(o4["barred_R_peri_kpc"].iloc[i + k]),
                "old_R_apo_kpc":  float(o4["barred_R_apo_kpc"].iloc[i + k]),
                "old_z_max_kpc":  float(o4["barred_z_max_kpc"].iloc[i + k]),
                "old_ecc":        float(o4["barred_ecc"].iloc[i + k]),
                "old_EJ_drift":   float(o4["barred_EJ_drift_rel"].iloc[i + k]),
            })
    df = pd.DataFrame(rows)
    Table.from_pandas(df).write(OUT / "orbits_barred_dt0p1.fits", overwrite=True)
    log(f"wrote orbits_barred_dt0p1.fits ({len(df)} rows)")

    # Gate 6B summary
    tABC = df["tier"].isin(["A", "B", "C"])
    bridgers_old = (df["old_R_peri_kpc"] < 2.0) & (df["old_R_apo_kpc"] > 15.0) & tABC
    bridgers_new = (df["fine_R_peri_kpc"] < 2.0) & (df["fine_R_apo_kpc"] > 15.0) & tABC
    summary = {
        "n_total": int(len(df)),
        "n_TierABC": int(tABC.sum()),
        "EJ_drift_p50_old":  float(df.loc[tABC, "old_EJ_drift"].median()),
        "EJ_drift_p84_old":  float(df.loc[tABC, "old_EJ_drift"].quantile(0.84)),
        "EJ_drift_p50_fine": float(df.loc[tABC, "fine_EJ_drift_rel"].median()),
        "EJ_drift_p84_fine": float(df.loc[tABC, "fine_EJ_drift_rel"].quantile(0.84)),
        "EJ_drift_p99_fine": float(df.loc[tABC, "fine_EJ_drift_rel"].quantile(0.99)),
        "median_R_peri_pc_old":  float(df.loc[tABC, "old_R_peri_kpc"].median() * 1000),
        "median_R_peri_pc_fine": float(df.loc[tABC, "fine_R_peri_kpc"].median() * 1000),
        "n_R_peri_lt_100pc_old":  int(((df["old_R_peri_kpc"] < 0.100) & tABC).sum()),
        "n_R_peri_lt_100pc_fine": int(((df["fine_R_peri_kpc"] < 0.100) & tABC).sum()),
        "n_min_rsph_lt_10pc_fine": int(((df["fine_min_r_sph_kpc"] < 0.010) & tABC).sum()),
        "n_bridgers_old":  int(bridgers_old.sum()),
        "n_bridgers_fine": int(bridgers_new.sum()),
        "bridgers_old_source_ids":  df.loc[bridgers_old, "source_id"].tolist(),
        "bridgers_fine_source_ids": df.loc[bridgers_new, "source_id"].tolist(),
    }
    (OUT / "gate6B_barred_dt.json").write_text(json.dumps(summary, indent=2))
    log("wrote gate6B_barred_dt.json")
    print(json.dumps(summary, indent=2))
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
