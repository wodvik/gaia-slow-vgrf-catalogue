"""Phase 14V -- expanded potential sensitivity table.

Runs compact point-estimate sensitivity checks from the expanded
Tier A+B+C orbit initial conditions:

* static Hunter+2024 default
* static halo-mass variants at 0.85x and 1.15x
* barred Hunter/Sormani variants at Omega_p = 33, 37.5, 41 km/s/kpc

This is intentionally a table-level sensitivity pass, not a replacement
for the full Phase 6 posterior orbit battery.
"""
from __future__ import annotations

import json
from pathlib import Path

import agama
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2].parent
BUNDLE = Path(__file__).resolve().parents[1]
OUT = BUNDLE / "phase14"
WORK = REPO / "release" / "_iterations" / "v2" / "phase3_agama" / "_hunter24_workdir"
IN = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.csv"
OUT.mkdir(parents=True, exist_ok=True)

GYR = 1.0 / 0.9778
T_INT = 4.0 * GYR
TRAJSIZE = 40001
ANGLE_BAR = -0.44

agama.setUnits(length=1, mass=1, velocity=1)


def pmin(R):
    """Exact pericentre via parabolic interpolation at dR sign-changes."""
    dR = np.diff(R)
    idx = np.where((dR[:-1] < 0) & (dR[1:] > 0))[0] + 1
    best = float(R.min())
    for k in idx:
        y0, y1, y2 = R[k-1], R[k], R[k+1]
        den = y0 - 2*y1 + y2
        if den > 0:
            ym = y1 - (y0 - y2)**2 / (8.0*den)
            if ym < best:
                best = ym
    return best


def reduce_orbits(result, pot=None, omega=None):
    rows = []
    for _, traj in result:
        arr = np.asarray(traj)
        x, y, z, vx, vy, vz = (arr[:, j] for j in range(6))
        rcyl = np.sqrt(x * x + y * y)
        rsph = np.sqrt(rcyl * rcyl + z * z)
        rperi = float(pmin(rcyl))
        rapo = float(rcyl.max())
        row = {
            "R_peri_kpc": rperi,
            "R_apo_kpc": rapo,
            "z_max_kpc": float(np.abs(z).max()),
            "min_r_sph_kpc": float(rsph.min()),
            "ecc": float((rapo - rperi) / (rapo + rperi)),
        }
        if pot is not None and omega is not None:
            lz = x * vy - y * vx
            energy = 0.5 * (vx * vx + vy * vy + vz * vz) + pot.potential(arr[:, :3])
            ej = energy - omega * lz
            row["EJ_drift_rel"] = float(abs((ej[-1] - ej[0]) / ej[0]))
        rows.append(row)
    return pd.DataFrame(rows)


def summarise(label, kind, df):
    bridger = (df["R_peri_kpc"] < 2.0) & (df["R_apo_kpc"] > 15.0)
    row = {
        "variant": label,
        "kind": kind,
        "n": int(len(df)),
        "median_R_peri_pc": float(df["R_peri_kpc"].median() * 1000.0),
        "median_R_apo_kpc": float(df["R_apo_kpc"].median()),
        "median_ecc": float(df["ecc"].median()),
        "n_R_peri_lt_100pc": int((df["R_peri_kpc"] < 0.1).sum()),
        "n_bridgers": int(bridger.sum()),
    }
    if "EJ_drift_rel" in df:
        row["EJ_drift_p84"] = float(df["EJ_drift_rel"].quantile(0.84))
    return row


def main() -> int:
    stars = pd.read_csv(IN)
    ic = stars[["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]].to_numpy(float)
    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    pot_full = agama.Potential(file=str(WORK / "MWPotentialHunter24_full.ini"))
    pot_bar = agama.Potential(potential=pot_full, rotation=ANGLE_BAR)
    delta_low = agama.Potential(type="NFW", scaleRadius=19.6, mass=1.27e12 * (0.85 - 1.0))
    delta_high = agama.Potential(type="NFW", scaleRadius=19.6, mass=1.27e12 * (1.15 - 1.0))
    variants = []

    for label, pot in [
        ("static_halo_0p85", agama.Potential(pot_axi, delta_low)),
        ("static_default", pot_axi),
        ("static_halo_1p15", agama.Potential(pot_axi, delta_high)),
    ]:
        print(f"[14V] integrating {label} for {len(ic)} stars", flush=True)
        res = agama.orbit(potential=pot, ic=ic, time=T_INT, trajsize=TRAJSIZE)
        tab = reduce_orbits(res)
        tab.insert(0, "variant", label)
        variants.append(tab)

    for label, omega_abs in [
        ("bar_omega33", 33.0),
        ("bar_omega37p5", 37.5),
        ("bar_omega41", 41.0),
    ]:
        omega = -omega_abs
        print(f"[14V] integrating {label} for {len(ic)} stars", flush=True)
        res = agama.orbit(potential=pot_bar, ic=ic, time=T_INT, trajsize=TRAJSIZE, Omega=omega)
        tab = reduce_orbits(res, pot=pot_bar, omega=omega)
        tab.insert(0, "variant", label)
        variants.append(tab)

    all_rows = pd.concat(variants, ignore_index=True)
    all_rows.to_csv(OUT / "expanded_potential_sensitivity_per_star.csv", index=False)

    summary_rows = []
    for variant, group in all_rows.groupby("variant", sort=False):
        kind = "barred" if variant.startswith("bar_") else "static"
        row = summarise(variant, kind, group)
        row["trajsize"] = TRAJSIZE
        row["pericentre_method"] = "parabola-interpolated cylindrical R minimum at dR=0"
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "expanded_potential_sensitivity_summary.csv", index=False)
    (OUT / "expanded_potential_sensitivity_summary.json").write_text(
        json.dumps(summary_rows, indent=2)
    )
    print(json.dumps(summary_rows, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
