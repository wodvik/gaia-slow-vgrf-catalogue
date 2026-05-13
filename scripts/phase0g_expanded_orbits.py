"""
Expanded-catalogue point-estimate orbit integration.

Run under WSL, where AGAMA is installed:

    wsl -- python3 /mnt/c/Users/humbl/GAIA2026/release/v2/scripts/phase0g_expanded_orbits.py

This integrates the expanded Tier A+B+C sample from Phase 0E in the same
Hunter+2024 static and barred potentials used by the manuscript's v2 Phase 4.
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table
import agama


if platform.system().lower() == "linux":
    REPO = Path("/mnt/c/Users/humbl/GAIA2026")
    DEFAULT_INPUT = Path("/mnt/d/GAIA/parent_scan/expanded_candidates_mc_tiered.csv")
else:
    REPO = Path("C:/Users/humbl/GAIA2026")
    DEFAULT_INPUT = Path("D:/GAIA/parent_scan/expanded_candidates_mc_tiered.csv")

CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase0_expanded"
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"

DEFAULT_SOLAR = CONFIG["solar_variants"]["default"]
OMEGA_P_DEFAULT = -float(CONFIG["bar_pattern_speeds_kms_kpc"]["default"])
ANGLE_BAR = -0.44
GYR_TO_AGAMA = 1.0 / 0.9778
T_INTEGRATE = 4.0 * GYR_TO_AGAMA
TRAJ_STEPS = 2001

agama.setUnits(length=1, mass=1, velocity=1)


def log(msg: str) -> None:
    print(f"[phase0g t={time.time() - T0:7.1f}s] {msg}", flush=True)


def galcen_frame() -> coord.Galactocentric:
    return coord.Galactocentric(
        galcen_distance=DEFAULT_SOLAR["R0_kpc"] * u.kpc,
        z_sun=DEFAULT_SOLAR["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            DEFAULT_SOLAR["U_kms"] * u.km / u.s,
            (DEFAULT_SOLAR["Vc_kms"] + DEFAULT_SOLAR["V_kms"]) * u.km / u.s,
            DEFAULT_SOLAR["W_kms"] * u.km / u.s,
        ),
    )


def initial_conditions(df: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    frame = galcen_frame()
    icrs = coord.SkyCoord(
        ra=df["ra"].to_numpy(dtype=float) * u.deg,
        dec=df["dec"].to_numpy(dtype=float) * u.deg,
        distance=df["dist_pc_final_screen"].to_numpy(dtype=float) * u.pc,
        pm_ra_cosdec=df["pmra"].to_numpy(dtype=float) * u.mas / u.yr,
        pm_dec=df["pmdec"].to_numpy(dtype=float) * u.mas / u.yr,
        radial_velocity=df["radial_velocity"].to_numpy(dtype=float) * u.km / u.s,
        frame="icrs",
    )
    g = icrs.transform_to(frame)
    df = df.copy()
    df["x_kpc"] = g.x.to_value(u.kpc)
    df["y_kpc"] = g.y.to_value(u.kpc)
    df["z_kpc"] = g.z.to_value(u.kpc)
    df["vx_kms"] = g.v_x.to_value(u.km / u.s)
    df["vy_kms"] = g.v_y.to_value(u.km / u.s)
    df["vz_kms"] = g.v_z.to_value(u.km / u.s)
    df["vgrf_default_orbit_ic"] = np.sqrt(
        df["vx_kms"] ** 2 + df["vy_kms"] ** 2 + df["vz_kms"] ** 2
    )
    ic = np.column_stack(
        [
            df["x_kpc"].to_numpy(),
            df["y_kpc"].to_numpy(),
            df["z_kpc"].to_numpy(),
            df["vx_kms"].to_numpy(),
            df["vy_kms"].to_numpy(),
            df["vz_kms"].to_numpy(),
        ]
    )
    return df, ic


def classify_rvs_quality(df: pd.DataFrame) -> np.ndarray:
    n = len(df)
    quality = np.full(n, "ok", dtype=object)
    pval = pd.to_numeric(df["rv_chisq_pvalue"], errors="coerce").to_numpy(dtype=float)
    snr_rv = pd.to_numeric(df["rv_expected_sig_to_noise"], errors="coerce").to_numpy(dtype=float)
    amp = pd.to_numeric(df["rv_amplitude_robust"], errors="coerce").to_numpy(dtype=float)
    n_trans = pd.to_numeric(df["rv_nb_transits"], errors="coerce").to_numpy(dtype=float)
    teff = pd.to_numeric(df["rv_template_teff"], errors="coerce").to_numpy(dtype=float)
    grvs = pd.to_numeric(df["grvs_mag"], errors="coerce").to_numpy(dtype=float)
    poor = (
        (np.isfinite(pval) & (pval < 0.01))
        | (np.isfinite(snr_rv) & (snr_rv < 2.0))
    )
    marginal = (
        (np.isfinite(amp) & np.isfinite(n_trans) & (amp > 10.0) & (n_trans < 10))
        | (
            np.isfinite(teff)
            & np.isfinite(grvs)
            & ((teff < 3500.0) | (teff > 7000.0))
            & (grvs > 12.0)
        )
    )
    quality[marginal] = "marginal"
    quality[poor] = "poor"
    return quality


def orbit_summary(traj: np.ndarray, pot, omega_p: float) -> dict:
    x, y, z, vx, vy, vz = (traj[:, i] for i in range(6))
    R_cyl = np.sqrt(x * x + y * y)
    r_sph = np.sqrt(x * x + y * y + z * z)
    Lz = x * vy - y * vx
    KE = 0.5 * (vx * vx + vy * vy + vz * vz)
    PE = pot.potential(traj[:, :3])
    E = KE + PE
    E_J = E - omega_p * Lz
    dR = np.diff(R_cyl)
    n_peri = int(np.sum((dR[:-1] < 0) & (dR[1:] > 0)))
    return {
        "R_peri_kpc": float(R_cyl.min()),
        "R_apo_kpc": float(R_cyl.max()),
        "z_max_kpc": float(np.abs(z).max()),
        "min_r_sph_kpc": float(r_sph.min()),
        "ecc": float((R_cyl.max() - R_cyl.min()) / (R_cyl.max() + R_cyl.min())),
        "n_peri": n_peri,
        "E_drift_rel": float(abs((E[-1] - E[0]) / E[0])),
        "E_range_rel": float((E.max() - E.min()) / abs(E.mean())),
        "EJ_drift_rel": float(abs((E_J[-1] - E_J[0]) / E_J[0])),
        "EJ_range_rel": float((E_J.max() - E_J.min()) / abs(E_J.mean())),
        "Lz_kpc_kms": float(np.mean(Lz)),
    }


def integrate(ic: np.ndarray, pot, omega_p: float, label: str) -> pd.DataFrame:
    log(f"integrating {len(ic)} {label} orbits for 4 Gyr")
    kwargs = {"potential": pot, "ic": ic, "time": T_INTEGRATE, "trajsize": TRAJ_STEPS}
    if omega_p != 0.0:
        kwargs["Omega"] = omega_p
    result = agama.orbit(**kwargs)
    rows = []
    for i in range(len(ic)):
        traj = np.asarray(result[i, 1])
        row = orbit_summary(traj, pot, omega_p)
        row["star_idx"] = i
        rows.append(row)
    out = pd.DataFrame(rows).add_prefix(f"{label}_")
    log(f"finished {label}")
    return out


def stackel_actions(ic: np.ndarray, pot_axi) -> pd.DataFrame:
    log(f"computing Stackel actions for {len(ic)} stars")
    af = agama.ActionFinder(pot_axi, interp=True)
    actions, angles, freqs = af(ic, angles=True)
    out = pd.DataFrame(
        {
            "J_R": actions[:, 0],
            "J_z": actions[:, 1],
            "J_phi": actions[:, 2],
            "Omega_R": freqs[:, 0],
            "Omega_z": freqs[:, 1],
            "Omega_phi": freqs[:, 2],
        }
    )
    delta = out["Omega_phi"] - OMEGA_P_DEFAULT
    out["res_ratio_OmegaR_over_dPhi"] = np.where(
        np.abs(delta) > 1e-3, out["Omega_R"] / delta, np.nan
    )
    return out


def pct(df: pd.DataFrame, mask: np.ndarray, col: str) -> dict:
    v = df.loc[mask, col].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    return {
        "n": int(len(v)),
        "p16": float(np.percentile(v, 16)),
        "p50": float(np.percentile(v, 50)),
        "p84": float(np.percentile(v, 84)),
        "min": float(v.min()),
        "max": float(v.max()),
    }


def run(args: argparse.Namespace) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    src = pd.read_csv(args.input_csv)
    src["tier"] = src["tier"].astype(str)
    src["source_id"] = src["source_id"].astype("int64")
    src["source_in_old_preselection"] = src["source_in_old_preselection"].astype(bool)
    src["rv_quality"] = classify_rvs_quality(src)
    src["rvs_quality_ok"] = src["rv_quality"] == "ok"
    df = src[src["tier"].isin(["A", "B", "C"])].copy().reset_index(drop=True)
    log(f"loaded {len(src)} expanded rows; integrating {len(df)} Tier A+B+C rows")
    df, ic = initial_conditions(df)

    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    pot_full = agama.Potential(file=str(WORK / "MWPotentialHunter24_full.ini"))
    pot_bar = agama.Potential(potential=pot_full, rotation=ANGLE_BAR)

    static = integrate(ic, pot_axi, 0.0, "static")
    barred = integrate(ic, pot_bar, OMEGA_P_DEFAULT, "barred")
    try:
        actions = stackel_actions(ic, pot_axi)
    except Exception as exc:
        log(f"WARN action computation failed: {exc}")
        actions = pd.DataFrame(
            {
                c: np.full(len(df), np.nan)
                for c in ["J_R", "J_z", "J_phi", "Omega_R", "Omega_z", "Omega_phi",
                          "res_ratio_OmegaR_over_dPhi"]
            }
        )

    keep = [
        "source_id",
        "tier",
        "P_vgrf_below_25",
        "rv_quality",
        "rvs_quality_ok",
        "source_in_old_preselection",
        "vgrf_default_exact",
        "vgrf_default_orbit_ic",
        "dist_pc_final_screen",
        "dist_source_final_screen",
        "x_kpc",
        "y_kpc",
        "z_kpc",
        "vx_kms",
        "vy_kms",
        "vz_kms",
    ]
    out = pd.concat([df[keep].reset_index(drop=True), static, barred, actions], axis=1)
    out = out.loc[:, ~out.columns.duplicated()]
    out_path = OUT / "catalogue_expanded_orbits_tierABC.fits"
    Table.from_pandas(out).write(out_path, overwrite=True)
    csv_path = OUT / "catalogue_expanded_orbits_tierABC.csv"
    out.to_csv(csv_path, index=False)

    tier_ab = out["tier"].isin(["A", "B"]).to_numpy()
    tier_abc = out["tier"].isin(["A", "B", "C"]).to_numpy()
    old = out["source_in_old_preselection"].to_numpy(dtype=bool)
    summary = {
        "input_csv": str(args.input_csv),
        "output_fits": str(out_path),
        "output_csv": str(csv_path),
        "n_integrated_tierABC": int(len(out)),
        "n_tierAB": int(tier_ab.sum()),
        "n_tierABC": int(tier_abc.sum()),
        "n_old_preselection_tierABC": int(old.sum()),
        "n_outside_old_tierABC": int((~old).sum()),
        "static_R_peri_kpc_tierAB": pct(out, tier_ab, "static_R_peri_kpc"),
        "static_R_peri_kpc_tierABC": pct(out, tier_abc, "static_R_peri_kpc"),
        "barred_R_peri_kpc_tierABC": pct(out, tier_abc, "barred_R_peri_kpc"),
        "static_ecc_tierABC": pct(out, tier_abc, "static_ecc"),
        "static_R_apo_kpc_tierABC": pct(out, tier_abc, "static_R_apo_kpc"),
        "n_static_R_peri_lt_100pc_tierABC": int(((out["static_R_peri_kpc"] < 0.1) & tier_abc).sum()),
        "n_barred_R_peri_lt_100pc_tierABC": int(((out["barred_R_peri_kpc"] < 0.1) & tier_abc).sum()),
        "n_static_bridgers_tierABC": int(
            ((out["static_R_peri_kpc"] < 2.0) & (out["static_R_apo_kpc"] > 15.0) & tier_abc).sum()
        ),
        "n_barred_bridgers_tierABC": int(
            ((out["barred_R_peri_kpc"] < 2.0) & (out["barred_R_apo_kpc"] > 15.0) & tier_abc).sum()
        ),
    }
    summary_path = OUT / "expanded_orbit_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", default=str(DEFAULT_INPUT))
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    T0 = time.time()
    main()
