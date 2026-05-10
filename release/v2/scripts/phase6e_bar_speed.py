"""Phase 6E — bar pattern speed sensitivity.

500-star stratified sub-sample of Tier A+B+C, drawn to span (R_apo, tier)
bins. Re-integrate sub-sample in pot_barred_slow (Omega=33) and
pot_barred_fast (Omega=41) at dt=0.1 Myr (40001 steps over 4 Gyr).
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
SEED = int(CONFIG["mc"]["random_seed"])
BARS = CONFIG["bar_pattern_speeds_kms_kpc"]
GYR = 1.0 / 0.9778
agama.setUnits(length=1, mass=1, velocity=1)


def log(m): print(f"[6E t={time.time()-T0:5.1f}s] {m}", flush=True)


def main():
    global T0; T0 = time.time()

    o4 = Table.read(REPO / "release/v2/phase4/catalogue_v2_orbits.fits").to_pandas()
    if o4["tier"].dtype == object and isinstance(o4["tier"].iloc[0], (bytes, bytearray)):
        o4["tier"] = o4["tier"].str.decode("utf-8")
    tABC = o4["tier"].isin(["A", "B", "C"]).to_numpy()
    sub = o4.loc[tABC].reset_index(drop=True)

    # Stratify by (tier, R_apo bin)
    rng = np.random.default_rng(SEED + 6500)
    bins_apo = np.linspace(0, 30, 7)  # 5 kpc bins
    sub["apo_bin"] = np.digitize(sub["barred_R_apo_kpc"], bins_apo)
    n_target = 500
    # Draw proportionally; floor at 1 per non-empty cell
    cells = sub.groupby(["tier", "apo_bin"])
    counts = cells.size()
    n_cells = len(counts)
    base = max(1, n_target // n_cells)
    picks = []
    for (t, b), grp in cells:
        n_take = min(base, len(grp))
        picks.append(grp.sample(n_take, random_state=int(rng.integers(1<<31))))
    chosen = pd.concat(picks).drop_duplicates("source_id").reset_index(drop=True)
    log(f"stratified subsample: {len(chosen)} / {len(sub)} stars")

    ic = chosen[["x_kpc", "y_kpc", "z_kpc",
                  "vx_kms", "vy_kms", "vz_kms"]].to_numpy()

    # For each bar speed, build a custom rotating .ini (Phase 3 style).
    rows = []
    for label, omega in [("slow", BARS["slow"]),
                          ("default", BARS["default"]),
                          ("fast", BARS["fast"])]:
        omega_signed = -abs(omega)
        ini = WORK / f"MWPotentialHunter24_rot_phase6_{label}.ini"
        full_ini = WORK / "MWPotentialHunter24_full.ini"
        ini.write_text(
            f"# Phase 6E custom Omega = {omega_signed} km/s/kpc\n"
            f"[Potential]\nfile={full_ini}\n"
            f"rotation=[[0,-0.44],[1,{-0.44 + omega_signed}]]\n"
        )
        pot = agama.Potential(file=str(ini))
        log(f"--- {label} (Omega={omega_signed}) ---  trajsize=40001 (dt~0.1 Myr) ---")
        t0 = time.time()
        res = agama.orbit(potential=pot, ic=ic, time=4.0 * GYR,
                           trajsize=40001, Omega=omega_signed)
        log(f"  agama.orbit done in {time.time()-t0:.1f}s")
        for k in range(len(ic)):
            traj = np.asarray(res[k, 1])
            x, y, z, vx, vy, vz = (traj[:, j] for j in range(6))
            Rcyl = np.sqrt(x*x + y*y); rsph = np.sqrt(Rcyl**2 + z*z)
            Lz = x*vy - y*vx
            E = 0.5*(vx*vx+vy*vy+vz*vz) + pot.potential(traj[:, :3])
            E_J = E - omega_signed * Lz
            rows.append({
                "source_id": int(chosen["source_id"].iloc[k]),
                "tier": str(chosen["tier"].iloc[k]),
                "Omega_p": omega_signed,
                "speed_label": label,
                "R_peri_kpc": float(Rcyl.min()),
                "R_apo_kpc":  float(Rcyl.max()),
                "z_max_kpc":  float(np.abs(z).max()),
                "min_r_sph_kpc": float(rsph.min()),
                "ecc": float((Rcyl.max() - Rcyl.min()) / (Rcyl.max() + Rcyl.min())),
                "EJ_drift_rel": float(abs((E_J[-1] - E_J[0]) / E_J[0])),
            })
    df = pd.DataFrame(rows)
    Table.from_pandas(df).write(OUT / "bar_speed_subsample.fits", overwrite=True)

    # Per-speed median summary
    summary = {"n_stars_subsample": int(len(chosen)),
               "by_speed": {}}
    for label in ("slow", "default", "fast"):
        m = df["speed_label"] == label
        summary["by_speed"][label] = {
            "Omega_p": float(df.loc[m, "Omega_p"].iloc[0]),
            "median_R_peri_pc": float(df.loc[m, "R_peri_kpc"].median() * 1000),
            "median_R_apo_kpc": float(df.loc[m, "R_apo_kpc"].median()),
            "median_ecc":       float(df.loc[m, "ecc"].median()),
            "EJ_drift_p84":     float(df.loc[m, "EJ_drift_rel"].quantile(0.84)),
            "n_R_peri_lt_100pc": int((df.loc[m, "R_peri_kpc"] < 0.100).sum()),
            "n_bridgers": int(((df.loc[m, "R_peri_kpc"] < 2.0) &
                                (df.loc[m, "R_apo_kpc"] > 15.0)).sum()),
        }
    (OUT / "gate6E_bar_speed.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
