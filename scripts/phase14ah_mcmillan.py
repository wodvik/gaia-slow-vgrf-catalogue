"""Phase 14AH -- McMillan (2017) comparison-potential sensitivity.

Referee response (deep-review Issue 7). Adds a widely used, independently
calibrated axisymmetric Milky Way potential (McMillan 2017, MNRAS 465, 76;
agama's bundled McMillan17.ini) as a standard benchmark, to show the
radial/compact orbit result is not an artefact of the less-standard,
non-Gaia-DR3-calibrated Hunter+2024 model family.

Same conventions as phase14v_expanded_potential_sensitivity.py: the fixed
Tier A+B+C orbit initial conditions, kpc/Msun/(km/s) units, 4 Gyr,
40001-point readout, parabola-interpolated cylindrical pericentre. Runs
under WSL/agama.

Outputs
-------
phase14/expanded_mcmillan_per_star.csv
phase14/expanded_mcmillan_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import agama
import numpy as np
import pandas as pd

BUNDLE = Path(__file__).resolve().parents[1]
OUT = BUNDLE / "phase14"
IN = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.csv"
MCM_INI = Path(agama.__file__).parent / "data" / "McMillan17.ini"

GYR = 1.0 / 0.9778
T_INT = 4.0 * GYR
TRAJSIZE = 40001
R0 = 8.178

agama.setUnits(length=1, mass=1, velocity=1)


def pmin(R):
    dR = np.diff(R)
    idx = np.where((dR[:-1] < 0) & (dR[1:] > 0))[0] + 1
    best = float(R.min())
    for k in idx:
        y0, y1, y2 = R[k - 1], R[k], R[k + 1]
        den = y0 - 2 * y1 + y2
        if den > 0:
            ym = y1 - (y0 - y2) ** 2 / (8.0 * den)
            if ym < best:
                best = ym
    return best


def reduce_orbits(result):
    rows = []
    for _, traj in result:
        arr = np.asarray(traj)
        x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
        rcyl = np.sqrt(x * x + y * y)
        rsph = np.sqrt(rcyl * rcyl + z * z)
        rperi = float(pmin(rcyl))
        rapo = float(rcyl.max())
        rows.append({
            "R_peri_kpc": rperi, "R_apo_kpc": rapo,
            "z_max_kpc": float(np.abs(z).max()),
            "min_r_sph_kpc": float(rsph.min()),
            "ecc": float((rapo - rperi) / (rapo + rperi)),
        })
    return pd.DataFrame(rows)


def main() -> int:
    pot = agama.Potential(file=str(MCM_INI))
    vc = float((-R0 * pot.force(R0, 0.0, 0.0)[0]) ** 0.5)
    print(f"[14AH] McMillan17 loaded; circular velocity at R0={R0} kpc: {vc:.1f} km/s", flush=True)

    stars = pd.read_csv(IN)
    ic = stars[["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]].to_numpy(float)
    print(f"[14AH] integrating McMillan17 for {len(ic)} stars", flush=True)
    res = agama.orbit(potential=pot, ic=ic, time=T_INT, trajsize=TRAJSIZE)
    df = reduce_orbits(res)
    df.insert(0, "variant", "mcmillan17")
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "expanded_mcmillan_per_star.csv", index=False)

    bridger = (df["R_peri_kpc"] < 2.0) & (df["R_apo_kpc"] > 15.0)
    summary = {
        "variant": "McMillan+2017",
        "potential_file": str(MCM_INI),
        "circular_velocity_R0_kms": round(vc, 1),
        "n": int(len(df)),
        "median_R_peri_pc": round(float(df["R_peri_kpc"].median() * 1000.0), 1),
        "median_R_apo_kpc": round(float(df["R_apo_kpc"].median()), 2),
        "median_ecc": round(float(df["ecc"].median()), 3),
        "frac_e_gt_0p8_pct": round(float((df["ecc"] > 0.8).mean() * 100.0), 1),
        "n_R_peri_lt_100pc": int((df["R_peri_kpc"] < 0.1).sum()),
        "n_R_peri_lt_50pc": int((df["R_peri_kpc"] < 0.05).sum()),
        "n_R_peri_lt_10pc": int((df["R_peri_kpc"] < 0.01).sum()),
        "n_min_rsph_lt_100pc": int((df["min_r_sph_kpc"] < 0.1).sum()),
        "n_min_rsph_lt_10pc": int((df["min_r_sph_kpc"] < 0.01).sum()),
        "n_bridgers": int(bridger.sum()),
        "trajsize": TRAJSIZE,
    }
    (OUT / "expanded_mcmillan_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
