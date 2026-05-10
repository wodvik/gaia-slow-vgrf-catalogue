"""Phase 6D — halo mass sensitivity (full-sample point-estimate orbits).

Re-integrate Tier A+B+C in pot_static_full_halo_0.85 and
pot_static_full_halo_1.15. Phase 3 already built these as
delta-NFW additions on top of the Hunter+2024 axi base.

We rebuild the variant potentials inline here (Phase 3's pickle wasn't
loaded back; cleanest is to recreate them at runtime).
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
GYR = 1.0 / 0.9778
agama.setUnits(length=1, mass=1, velocity=1)


def log(m): print(f"[6D t={time.time()-T0:5.1f}s] {m}", flush=True)


def make_halo_variant(axi_pot, factor):
    delta_halo = agama.Potential(type="NFW", scaleRadius=19.6,
                                 mass=1.27e12 * (factor - 1.0))
    return agama.Potential(axi_pot, delta_halo)


def reduce_orbit(traj):
    Rcyl = np.sqrt(traj[:, 0]**2 + traj[:, 1]**2)
    rsph = np.sqrt(Rcyl**2 + traj[:, 2]**2)
    return {
        "R_peri_kpc": float(Rcyl.min()),
        "R_apo_kpc":  float(Rcyl.max()),
        "z_max_kpc":  float(np.abs(traj[:, 2]).max()),
        "min_r_sph_kpc": float(rsph.min()),
        "ecc": float((Rcyl.max() - Rcyl.min()) / (Rcyl.max() + Rcyl.min())),
    }


def main():
    global T0; T0 = time.time()
    log(f"agama {agama.__version__}")

    log("loading Phase 4 orbits + base potential ...")
    o4 = Table.read(REPO / "release/v2/phase4/catalogue_v2_orbits.fits").to_pandas()
    if o4["tier"].dtype == object and isinstance(o4["tier"].iloc[0], (bytes, bytearray)):
        o4["tier"] = o4["tier"].str.decode("utf-8")
    tABC = o4["tier"].isin(["A", "B", "C"]).to_numpy()
    sub = o4.loc[tABC].reset_index(drop=True)
    log(f"Tier A+B+C: {len(sub)}")
    ic = sub[["x_kpc", "y_kpc", "z_kpc",
              "vx_kms", "vy_kms", "vz_kms"]].to_numpy()

    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    pot_h085 = make_halo_variant(pot_axi, 0.85)
    pot_h115 = make_halo_variant(pot_axi, 1.15)

    rows = []
    for label, pot in [("halo_0p85", pot_h085), ("default", pot_axi),
                        ("halo_1p15", pot_h115)]:
        log(f"--- variant {label} ---")
        t0 = time.time()
        res = agama.orbit(potential=pot, ic=ic, time=4.0 * GYR, trajsize=1001)
        log(f"  agama.orbit done in {time.time()-t0:.1f}s")
        Rp = np.zeros(len(ic)); Ra = np.zeros(len(ic))
        zm = np.zeros(len(ic)); rsph = np.zeros(len(ic))
        for k in range(len(ic)):
            s = reduce_orbit(np.asarray(res[k, 1]))
            Rp[k] = s["R_peri_kpc"]; Ra[k] = s["R_apo_kpc"]
            zm[k] = s["z_max_kpc"]; rsph[k] = s["min_r_sph_kpc"]
        ecc = (Ra - Rp) / (Ra + Rp)
        rows.append({
            "variant": label, "n_TierABC": int(len(ic)),
            "median_R_peri_pc":  float(np.median(Rp) * 1000),
            "median_R_apo_kpc":  float(np.median(Ra)),
            "median_z_max_kpc":  float(np.median(zm)),
            "median_ecc":        float(np.median(ecc)),
            "n_R_peri_lt_100pc": int(np.sum(Rp < 0.100)),
            "n_min_rsph_lt_10pc": int(np.sum(rsph < 0.010)),
        })
    pd.DataFrame(rows).to_csv(OUT / "halo_mass_table.csv", index=False)

    # Per-star delta vs default
    default_idx = next(i for i, r in enumerate(rows) if r["variant"] == "default")
    print(json.dumps(rows, indent=2))
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
