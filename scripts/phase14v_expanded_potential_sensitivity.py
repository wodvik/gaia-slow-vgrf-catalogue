"""Phase 14V -- expanded potential sensitivity table.

Runs compact point-estimate sensitivity checks from the expanded Tier A+B+C
orbit initial conditions. The shipped v1.0.6-review script regenerates the
18-variant table used by the manuscript:

* static Hunter+2024 default plus halo-mass variants
* barred Hunter/Sormani pattern-speed variants
* barred angle variants at Omega_p = 37.5 km/s/kpc
* static halo-flattening q_z variants

The q_z sweep rebuilds Hunter+2024 static potentials from the local Hunter
workdir under release/_iterations. This is a table-level sensitivity pass, not
a replacement for the full Phase 6 posterior orbit battery.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import agama
import numpy as np
import pandas as pd

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parent.parent
OUT = BUNDLE / "phase14"
POTENTIALS = BUNDLE / "potentials"
HUNTER_WORK = REPO / "release" / "_iterations" / "v2" / "phase3_agama" / "_hunter24_workdir"
IN = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.csv"

GYR = 1.0 / 0.9778
T_INT = 4.0 * GYR
TRAJSIZE = 40001
ANGLE_BAR = -0.44

STATIC_MASS_VARIANTS = [
    ("static_halo_0p85", 0.85),
    ("static_default", 1.00),
    ("static_halo_1p15", 1.15),
]
OMEGA_VARIANTS = [
    ("bar_omega33", 33.0),
    ("bar_omega37p5", 37.5),
    ("bar_omega41", 41.0),
    ("bar_omega24", 24.0),
    ("bar_omega28", 28.0),
]
BAR_ANGLE_VARIANTS = [
    ("bar_omega37p5_ang20", 20.0),
    ("bar_omega37p5_ang25", 25.0),
    ("bar_omega37p5_ang30", 30.0),
    ("bar_omega37p5_ang35", 35.0),
]
QZ_VARIANTS = [
    ("static_halo_qz0p80", 0.80),
    ("static_halo_qz0p85", 0.85),
    ("static_halo_qz0p90", 0.90),
    ("static_halo_qz0p95", 0.95),
    ("static_halo_qz1p00", 1.00),
    ("static_halo_qz1p05", 1.05),
]

PER_STAR_COLUMNS = [
    "variant",
    "R_peri_kpc",
    "R_apo_kpc",
    "z_max_kpc",
    "min_r_sph_kpc",
    "ecc",
    "EJ_drift_rel",
    "source_id",
]

agama.setUnits(length=1, mass=1, velocity=1)


def log(message: str) -> None:
    print(f"[14V] {message}", flush=True)


def pmin(r: np.ndarray) -> float:
    """Pericentre via parabolic interpolation at dR sign changes."""
    d_r = np.diff(r)
    idx = np.where((d_r[:-1] < 0) & (d_r[1:] > 0))[0] + 1
    best = float(r.min())
    for k in idx:
        y0, y1, y2 = r[k - 1], r[k], r[k + 1]
        den = y0 - 2 * y1 + y2
        if den > 0:
            ym = y1 - (y0 - y2) ** 2 / (8.0 * den)
            if ym < best:
                best = float(ym)
    return best


def reduce_orbits(result, source_ids: np.ndarray, pot=None, omega=None) -> pd.DataFrame:
    rows = []
    for source_id, (_, traj) in zip(source_ids, result):
        arr = np.asarray(traj)
        x, y, z, vx, vy, vz = (arr[:, j] for j in range(6))
        rcyl = np.sqrt(x * x + y * y)
        rsph = np.sqrt(rcyl * rcyl + z * z)
        rperi = pmin(rcyl)
        rapo = float(rcyl.max())
        row = {
            "R_peri_kpc": rperi,
            "R_apo_kpc": rapo,
            "z_max_kpc": float(np.abs(z).max()),
            "min_r_sph_kpc": float(rsph.min()),
            "ecc": float((rapo - rperi) / (rapo + rperi)),
            "source_id": source_id,
        }
        if pot is not None and omega is not None:
            lz = x * vy - y * vx
            energy = 0.5 * (vx * vx + vy * vy + vz * vz) + pot.potential(arr[:, :3])
            ej = energy - omega * lz
            row["EJ_drift_rel"] = float(abs((ej[-1] - ej[0]) / ej[0]))
        rows.append(row)
    return pd.DataFrame(rows)


def summarise(label: str, kind: str, df: pd.DataFrame) -> dict:
    bridger = (df["R_peri_kpc"] < 2.0) & (df["R_apo_kpc"] > 15.0)
    row = {
        "variant": label,
        "kind": kind,
        "n": int(len(df)),
        "median_R_peri_pc": float(df["R_peri_kpc"].median() * 1000.0),
        "median_R_apo_kpc": float(df["R_apo_kpc"].median()),
        "median_ecc": float(df["ecc"].median()),
        "frac_ecc_gt_0p95": float((df["ecc"] > 0.95).mean() * 100.0),
        "n_R_peri_lt_100pc": int((df["R_peri_kpc"] < 0.1).sum()),
        "n_bridgers": int(bridger.sum()),
        "trajsize": TRAJSIZE,
        "pericentre_method": "parabola-interpolated cylindrical R minimum at dR=0",
    }
    if "EJ_drift_rel" in df.columns and df["EJ_drift_rel"].notna().any():
        row["EJ_drift_p84"] = float(df["EJ_drift_rel"].quantile(0.84))
    return row


def hunter_axi_with_halo_mass(scale: float):
    pot_axi = agama.Potential(file=str(POTENTIALS / "MWPotentialHunter24_axi.ini"))
    if scale == 1.0:
        return pot_axi
    delta = agama.Potential(type="NFW", scaleRadius=19.6, mass=1.27e12 * (scale - 1.0))
    return agama.Potential(pot_axi, delta)


def hunter_static_with_halo_qz(q_z: float):
    """Rebuild Hunter+2024 static potential while varying only dark-halo q_z."""
    if not HUNTER_WORK.exists():
        raise FileNotFoundError(
            f"Hunter+2024 workdir required for q_z sweep is missing: {HUNTER_WORK}"
        )
    sys.path.insert(0, str(HUNTER_WORK))
    from example_mw_bar_potential import makeBarDensity

    params_bh = dict(type="Plummer", mass=4.1e6, scaleRadius=1e-3)
    params_nsc = dict(
        type="Spheroid",
        mass=6.1e7,
        gamma=0.71,
        beta=4,
        alpha=1,
        axisRatioZ=0.73,
        scaleRadius=0.0059,
        outerCutoffRadius=0.1,
    )
    params_nsd = [
        dict(
            type="Spheroid",
            densityNorm=2.00583e12,
            gamma=0,
            beta=0,
            alpha=1,
            axisRatioZ=0.37,
            outerCutoffRadius=0.00506,
            cutoffStrength=0.72,
        ),
        dict(
            type="Spheroid",
            densityNorm=1.53e12,
            gamma=0,
            beta=0,
            alpha=1,
            axisRatioZ=0.37,
            outerCutoffRadius=0.0246,
            cutoffStrength=0.79,
        ),
    ]
    params_dark = dict(
        type="Spheroid",
        densityNorm=2.774e11,
        gamma=0,
        beta=0,
        alpha=1,
        axisRatioZ=q_z,
        outerCutoffRadius=8.682e-6,
        cutoffStrength=0.1704,
    )
    params_disk = [
        dict(
            type="Disk",
            surfaceDensity=1.332e9,
            scaleRadius=2.0,
            scaleHeight=0.3,
            innerCutoffRadius=2.7,
            sersicIndex=1,
        ),
        dict(
            type="Disk",
            surfaceDensity=8.97e8,
            scaleRadius=2.8,
            scaleHeight=0.9,
            innerCutoffRadius=2.7,
            sersicIndex=1,
        ),
    ]
    params_gas = [
        dict(
            type="Disk",
            surfaceDensity=5.81e7,
            scaleRadius=7,
            scaleHeight=-0.085,
            innerCutoffRadius=4,
            sersicIndex=1,
        ),
        dict(
            type="Disk",
            surfaceDensity=2.68e9,
            scaleRadius=1.5,
            scaleHeight=-0.045,
            innerCutoffRadius=12,
            sersicIndex=1,
        ),
    ]
    pot_mul = agama.Potential(
        type="Multipole",
        density=agama.Density(params_dark, params_bh, params_nsc, *params_nsd),
        lmax=12,
        gridSizeR=36,
        rmin=1e-4,
        rmax=1000,
    )
    pot_cyl_axi = agama.Potential(
        type="CylSpline",
        density=agama.Density(makeBarDensity(), *(params_disk + params_gas)),
        gridSizeR=30,
        gridSizez=32,
        Rmin=0.1,
        Rmax=200,
        zmin=0.05,
        zmax=200,
        mmax=0,
    )
    return agama.Potential(pot_mul, pot_cyl_axi)


def write_range_summary(summary: pd.DataFrame) -> None:
    ranges = {}
    for prefix, label in [
        ("bar_omega37p5_ang", "bar_angle_sweep"),
        ("static_halo_qz", "halo_flattening_sweep"),
    ]:
        sub = summary[summary["variant"].str.startswith(prefix)]
        if sub.empty:
            continue
        ranges[label] = {
            "variants": sub["variant"].tolist(),
            "median_R_peri_pc": [
                float(sub["median_R_peri_pc"].min()),
                float(sub["median_R_peri_pc"].max()),
            ],
            "median_R_apo_kpc": [
                float(sub["median_R_apo_kpc"].min()),
                float(sub["median_R_apo_kpc"].max()),
            ],
            "median_ecc": [
                float(sub["median_ecc"].min()),
                float(sub["median_ecc"].max()),
            ],
            "frac_ecc_gt_0p95": [
                float(sub["frac_ecc_gt_0p95"].min()),
                float(sub["frac_ecc_gt_0p95"].max()),
            ],
            "n_bridgers": [
                int(sub["n_bridgers"].min()),
                int(sub["n_bridgers"].max()),
            ],
        }
    (OUT / "expanded_potential_sensitivity_range_summary.json").write_text(
        json.dumps(ranges, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    stars = pd.read_csv(
        IN,
        usecols=["source_id", "x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"],
    )
    source_ids = stars["source_id"].to_numpy()
    ic = stars[["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]].to_numpy(float)
    log(f"loaded {len(ic)} orbit initial conditions")

    variants = []
    for label, halo_scale in STATIC_MASS_VARIANTS:
        log(f"integrating {label}")
        pot = hunter_axi_with_halo_mass(halo_scale)
        res = agama.orbit(potential=pot, ic=ic, time=T_INT, trajsize=TRAJSIZE)
        tab = reduce_orbits(res, source_ids)
        tab.insert(0, "variant", label)
        variants.append(tab)

    pot_full = agama.Potential(file=str(POTENTIALS / "MWPotentialHunter24_full.ini"))
    pot_bar_default = agama.Potential(potential=pot_full, rotation=ANGLE_BAR)
    for label, omega_abs in OMEGA_VARIANTS:
        omega = -omega_abs
        log(f"integrating {label}")
        res = agama.orbit(
            potential=pot_bar_default,
            ic=ic,
            time=T_INT,
            trajsize=TRAJSIZE,
            Omega=omega,
        )
        tab = reduce_orbits(res, source_ids, pot=pot_bar_default, omega=omega)
        tab.insert(0, "variant", label)
        variants.append(tab)

    for label, ang_deg in BAR_ANGLE_VARIANTS:
        omega = -37.5
        angle = -math.radians(ang_deg)
        log(f"integrating {label} angle={ang_deg:g} deg")
        pot_bar_angle = agama.Potential(potential=pot_full, rotation=angle)
        res = agama.orbit(
            potential=pot_bar_angle,
            ic=ic,
            time=T_INT,
            trajsize=TRAJSIZE,
            Omega=omega,
        )
        tab = reduce_orbits(res, source_ids, pot=pot_bar_angle, omega=omega)
        tab.insert(0, "variant", label)
        variants.append(tab)

    for label, q_z in QZ_VARIANTS:
        log(f"building static flattened-halo potential {label} q_z={q_z:g}")
        pot_qz = hunter_static_with_halo_qz(q_z)
        log(f"integrating {label}")
        res = agama.orbit(potential=pot_qz, ic=ic, time=T_INT, trajsize=TRAJSIZE)
        tab = reduce_orbits(res, source_ids)
        tab.insert(0, "variant", label)
        variants.append(tab)

    all_rows = pd.concat(variants, ignore_index=True)
    all_rows = all_rows.reindex(columns=PER_STAR_COLUMNS)
    all_rows.to_csv(OUT / "expanded_potential_sensitivity_per_star.csv", index=False)

    summary_rows = []
    for variant, group in all_rows.groupby("variant", sort=False):
        kind = "barred" if variant.startswith("bar_") else "static"
        summary_rows.append(summarise(variant, kind, group))
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "expanded_potential_sensitivity_summary.csv", index=False)
    (OUT / "expanded_potential_sensitivity_summary.json").write_text(
        json.dumps(summary_rows, indent=2) + "\n",
        encoding="utf-8",
    )
    write_range_summary(summary)
    print(json.dumps(summary_rows, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
