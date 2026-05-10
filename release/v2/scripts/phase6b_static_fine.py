"""Phase 6B supplemental — static_full at trajsize=40001 (sampling check).

Phase 6B revealed that the median R_peri shifted dramatically (66 pc -> 16 pc
in barred) when the OUTPUT trajectory sampling went from 2 Myr to 0.1 Myr —
the integrator itself was accurate but coarse sampling missed deep peri
passages. This script repeats the test for static_full to see if the
manuscript's headline static R_peri also moves under finer sampling.
"""
from __future__ import annotations

import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd, yaml, agama
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase6"
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"
GYR = 1.0 / 0.9778
agama.setUnits(length=1, mass=1, velocity=1)


def log(m): print(f"[6B-static t={time.time()-T0:5.1f}s] {m}", flush=True)


def main():
    global T0; T0 = time.time()
    o4 = Table.read(REPO / "release/v2/phase4/catalogue_v2_orbits.fits").to_pandas()
    if o4["tier"].dtype == object and isinstance(o4["tier"].iloc[0], (bytes, bytearray)):
        o4["tier"] = o4["tier"].str.decode("utf-8")
    tABC = o4["tier"].isin(["A", "B", "C"]).to_numpy()
    sub = o4.loc[tABC].reset_index(drop=True)
    ic = sub[["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]].to_numpy()
    log(f"static_full, Tier A+B+C: {len(ic)} stars")

    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))

    # Chunk 150 to stay under 3 GB transient.
    rows = []
    for i in range(0, len(ic), 150):
        ic_c = ic[i:i+150]
        log(f"chunk {i}-{i+len(ic_c)} at trajsize=40001 (0.1 Myr)...")
        t0 = time.time()
        res = agama.orbit(potential=pot_axi, ic=ic_c, time=4.0 * GYR,
                           trajsize=40001)
        log(f"  done in {time.time()-t0:.1f}s")
        for k in range(len(ic_c)):
            traj = np.asarray(res[k, 1])
            x, y, z = traj[:, 0], traj[:, 1], traj[:, 2]
            Rcyl = np.sqrt(x*x + y*y); rsph = np.sqrt(Rcyl**2 + z*z)
            rows.append({
                "source_id": int(sub["source_id"].iloc[i+k]),
                "tier": str(sub["tier"].iloc[i+k]),
                "fine_R_peri_kpc":   float(Rcyl.min()),
                "fine_R_apo_kpc":    float(Rcyl.max()),
                "fine_z_max_kpc":    float(np.abs(z).max()),
                "fine_min_r_sph_kpc": float(rsph.min()),
                "old_R_peri_kpc":    float(sub["static_R_peri_kpc"].iloc[i+k]),
                "old_min_r_sph_kpc": float(sub["static_min_r_sph_kpc"].iloc[i+k]),
            })
    df = pd.DataFrame(rows)
    Table.from_pandas(df).write(OUT / "orbits_static_dt0p1.fits", overwrite=True)
    summary = {
        "n_TierABC": int(len(ic)),
        "median_R_peri_pc_old":  float(df["old_R_peri_kpc"].median() * 1000),
        "median_R_peri_pc_fine": float(df["fine_R_peri_kpc"].median() * 1000),
        "n_R_peri_lt_100pc_old":  int((df["old_R_peri_kpc"] < 0.100).sum()),
        "n_R_peri_lt_100pc_fine": int((df["fine_R_peri_kpc"] < 0.100).sum()),
        "n_min_rsph_lt_10pc_old":  int((df["old_min_r_sph_kpc"] < 0.010).sum()),
        "n_min_rsph_lt_10pc_fine": int((df["fine_min_r_sph_kpc"] < 0.010).sum()),
    }
    (OUT / "gate6B_static_fine.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
