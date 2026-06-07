"""WP-7 GeDR3mock injection-recovery diagnostic.

This public copy lives in the review bundle and reproduces the
GeDR3mock injection-recovery diagnostic products:

  tables/v15/tab_recovery_summary.tex
  figures/fig_recovery_heatmap.pdf

The synthetic stars use GeDR3mock sightlines, photometry, stellar
parameters, and Gaia-like astrometric/RVS uncertainties.  Controlled
Galactocentric velocities are injected into the four requested true-Vgrf
bins so the rare 0-25 km/s regime has enough statistics for a recovery
measurement.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyvo
from astropy.coordinates import SkyCoord, Galactocentric, CartesianDifferential
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table


BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
OUT_REVIEW = BUNDLE / "phase14" / "injection_recovery"
PHASE15 = OUT_REVIEW
RELEASE_TABLE = BUNDLE / "tables/v15/tab_recovery_summary.tex"
RELEASE_FIG = BUNDLE / "figures/fig_recovery_heatmap.pdf"
BUNDLE_TABLE = RELEASE_TABLE
BUNDLE_FIG = RELEASE_FIG
PHASE0E = BUNDLE / "scripts/phase0e_expanded_mc_tiering.py"

VGRF_DECADES = [(0.0, 25.0), (25.0, 50.0), (50.0, 100.0), (100.0, 200.0)]
RNG_SEED = 20260523
SF_FLOOR = 0.02

SOLAR = {
    "R0_kpc": 8.178,
    "z_sun_pc": 25.0,
    "Vc_kms": 229.0,  # WP1/O1 (Phase-2): adopt 229 (Eilers 2019 / Hunter axisym); was 232.0
    "U_kms": 11.1,
    "V_kms": 12.24,
    "W_kms": 7.25,
}


def log(message: str) -> None:
    print(f"[wp7 t={time.time() - T0:8.1f}s] {message}", flush=True)


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO)).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def load_phase0e_module():
    spec = importlib.util.spec_from_file_location("phase0e_expanded_mc_tiering", PHASE0E)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {PHASE0E}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def galcen_frame() -> Galactocentric:
    return Galactocentric(
        galcen_distance=SOLAR["R0_kpc"] * u.kpc,
        z_sun=SOLAR["z_sun_pc"] * u.pc,
        galcen_v_sun=CartesianDifferential(
            SOLAR["U_kms"] * u.km / u.s,
            (SOLAR["Vc_kms"] + SOLAR["V_kms"]) * u.km / u.s,
            SOLAR["W_kms"] * u.km / u.s,
        ),
    )


def ensure_dirs() -> None:
    OUT_REVIEW.mkdir(parents=True, exist_ok=True)
    PHASE15.mkdir(parents=True, exist_ok=True)
    RELEASE_TABLE.parent.mkdir(parents=True, exist_ok=True)
    RELEASE_FIG.parent.mkdir(parents=True, exist_ok=True)
    BUNDLE_TABLE.parent.mkdir(parents=True, exist_ok=True)
    BUNDLE_FIG.parent.mkdir(parents=True, exist_ok=True)


def fetch_gedr3mock_seed(target_rows: int, force: bool) -> pd.DataFrame:
    cache = PHASE15 / f"gedr3mock_seed_top{target_rows}.fits"
    if cache.exists() and not force:
        log(f"using cached GeDR3mock seed {rel(cache)}")
        return Table.read(cache).to_pandas()

    log(f"querying GAVO TAP for {target_rows:,} GeDR3mock observable seeds")
    query = f"""
    SELECT TOP {target_rows}
      source_id, ra, dec, l, b,
      parallax, parallax_error,
      pmra_error, pmdec_error, radial_velocity_error,
      phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
      phot_rvs_mean_mag, bp_rp, logg, teff_val,
      astrometric_params_solved, random_index
    FROM gedr3mock.main
    WHERE random_index < 200000000
      AND phot_g_mean_mag < 17
      AND phot_rvs_mean_mag < 14
      AND teff_val BETWEEN 3500 AND 7000
      AND parallax > 0
      AND parallax_error > 0
      AND pmra_error IS NOT NULL
      AND pmdec_error IS NOT NULL
      AND radial_velocity_error IS NOT NULL
      AND bp_rp IS NOT NULL
      AND logg IS NOT NULL
    """
    service = pyvo.dal.TAPService("https://dc.g-vo.org/tap")
    job = service.submit_job(query, maxrec=target_rows)
    job.run()
    job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=1800)
    if job.phase != "COMPLETED":
        raise RuntimeError(f"GeDR3mock TAP query failed with phase={job.phase}")
    table = job.fetch_result().to_table()
    job.delete()
    table.write(cache, overwrite=True)
    shutil.copy2(cache, OUT_REVIEW / cache.name)
    log(f"cached {len(table):,} GeDR3mock rows to {rel(cache)}")
    return table.to_pandas()


def clean_seed(seed: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "ra",
        "dec",
        "l",
        "b",
        "parallax",
        "parallax_error",
        "pmra_error",
        "pmdec_error",
        "radial_velocity_error",
        "phot_g_mean_mag",
        "phot_rvs_mean_mag",
        "bp_rp",
        "logg",
        "teff_val",
    ]
    out = seed.copy()
    for col in cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    finite = np.all(np.isfinite(out[cols].to_numpy(float)), axis=1)
    positive_err = (
        (out["parallax_error"] > 0)
        & (out["pmra_error"] > 0)
        & (out["pmdec_error"] > 0)
        & (out["radial_velocity_error"] > 0)
    )
    observable = (
        (out["phot_g_mean_mag"] < 17)
        & (out["phot_rvs_mean_mag"] < 14)
        & out["teff_val"].between(3500, 7000)
        & (out["parallax"] > 0)
    )
    out = out.loc[finite & positive_err & observable].reset_index(drop=True)
    if out.empty:
        raise RuntimeError("no usable GeDR3mock rows after observability cuts")
    return out


def sample_base_rows(seed: pd.DataFrame, n_per_bin: int, rng: np.random.Generator) -> pd.DataFrame:
    total = n_per_bin * len(VGRF_DECADES)
    replace = len(seed) < total
    idx = rng.choice(len(seed), size=total, replace=replace)
    out = seed.iloc[idx].reset_index(drop=True).copy()
    out["wp7_mock_id"] = np.arange(total, dtype=np.int64)
    decade_labels: list[str] = []
    decade_lo: list[float] = []
    decade_hi: list[float] = []
    for lo, hi in VGRF_DECADES:
        decade_labels.extend([f"{int(lo)}-{int(hi)}"] * n_per_bin)
        decade_lo.extend([lo] * n_per_bin)
        decade_hi.extend([hi] * n_per_bin)
    out["vgrf_bin"] = decade_labels
    out["vgrf_true_lo"] = np.array(decade_lo)
    out["vgrf_true_hi"] = np.array(decade_hi)
    return out


def inject_true_velocities(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    frame = galcen_frame()
    n = len(df)
    d_pc = 1000.0 / df["parallax"].to_numpy(float)
    icrs_pos = SkyCoord(
        ra=df["ra"].to_numpy(float) * u.deg,
        dec=df["dec"].to_numpy(float) * u.deg,
        distance=d_pc * u.pc,
        frame="icrs",
    )
    gal_pos = icrs_pos.transform_to(frame)

    speeds = np.empty(n, dtype=float)
    for label, (lo, hi) in zip([f"{int(a)}-{int(b)}" for a, b in VGRF_DECADES], VGRF_DECADES):
        mask = df["vgrf_bin"].eq(label).to_numpy()
        speeds[mask] = rng.uniform(lo, hi, size=int(mask.sum()))

    direction = rng.normal(size=(n, 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    vel = direction * speeds[:, None]

    gal_full = SkyCoord(
        x=gal_pos.x,
        y=gal_pos.y,
        z=gal_pos.z,
        v_x=vel[:, 0] * u.km / u.s,
        v_y=vel[:, 1] * u.km / u.s,
        v_z=vel[:, 2] * u.km / u.s,
        frame=frame,
        representation_type="cartesian",
        differential_type="cartesian",
    )
    icrs_full = gal_full.transform_to(coord.ICRS())

    out = df.copy()
    out["dist_true_pc"] = d_pc
    out["vgrf_true_kms"] = speeds
    out["pmra_true"] = icrs_full.pm_ra_cosdec.to_value(u.mas / u.yr)
    out["pmdec_true"] = icrs_full.pm_dec.to_value(u.mas / u.yr)
    out["radial_velocity_true"] = icrs_full.radial_velocity.to_value(u.km / u.s)
    out["vx_gal_true_kms"] = vel[:, 0]
    out["vy_gal_true_kms"] = vel[:, 1]
    out["vz_gal_true_kms"] = vel[:, 2]
    return out


def approx_nu_eff(bp_rp: np.ndarray) -> np.ndarray:
    # Smooth monotonic approximation sufficient for L21 validity gating.
    nu = 1.76 - 0.23 * bp_rp + 0.015 * bp_rp * bp_rp
    return np.clip(nu, 1.11, 1.89)


def compute_zero_point_columns(df: pd.DataFrame) -> pd.DataFrame:
    from zero_point import zpt

    zpt.load_tables()
    out = df.copy()
    coords = SkyCoord(ra=out["ra"].to_numpy(float) * u.deg, dec=out["dec"].to_numpy(float) * u.deg)
    ecl = coords.barycentrictrueecliptic.lat.to_value(u.deg)
    g = out["phot_g_mean_mag"].to_numpy(float)
    bp_rp = out["bp_rp"].to_numpy(float)
    nu = approx_nu_eff(bp_rp)
    pseudocolour = np.full(len(out), 1.5)
    if "astrometric_params_solved" in out:
        asp_raw = pd.to_numeric(out["astrometric_params_solved"], errors="coerce").to_numpy(float)
    else:
        asp_raw = np.full(len(out), 31.0)
    asp = np.where(np.isfinite(asp_raw), asp_raw, 31).astype(int)
    asp = np.where(np.isin(asp, [31, 95]), asp, 31)

    zpcorr = zpt.get_zpt(g, nu, pseudocolour, ecl, asp, _warnings=False)
    valid = (
        np.isin(asp, [31, 95])
        & (g > 6.0)
        & (g < 21.0)
        & np.where(asp == 31, (nu > 1.1) & (nu < 1.9), (pseudocolour > 1.24) & (pseudocolour < 1.72))
    )
    out["nu_eff_used_in_astrometry"] = nu
    out["pseudocolour"] = pseudocolour
    out["ecl_lat"] = ecl
    out["astrometric_params_solved"] = asp
    out["zpcorr_value_mas"] = zpcorr
    out["zpcorr_value_uas"] = zpcorr * 1000.0
    out["zpcorr_valid"] = valid
    return out


def add_measurement_noise(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    out = compute_zero_point_columns(df)
    n = len(out)
    parallax_true = out["parallax"].to_numpy(float)
    parallax_error = out["parallax_error"].to_numpy(float)
    pmra_error = out["pmra_error"].to_numpy(float)
    pmdec_error = out["pmdec_error"].to_numpy(float)
    rv_error = out["radial_velocity_error"].to_numpy(float)

    out["parallax_true_mas"] = parallax_true
    out["parallax"] = parallax_true + out["zpcorr_value_mas"].to_numpy(float) + rng.normal(0.0, parallax_error)
    out["parallax_zpcorr"] = np.where(
        out["zpcorr_valid"].to_numpy(bool),
        out["parallax"].to_numpy(float) - out["zpcorr_value_mas"].to_numpy(float),
        out["parallax"].to_numpy(float),
    )
    out["parallax_zpcorr_error"] = parallax_error
    out["pmra"] = out["pmra_true"].to_numpy(float) + rng.normal(0.0, pmra_error, size=n)
    out["pmdec"] = out["pmdec_true"].to_numpy(float) + rng.normal(0.0, pmdec_error, size=n)
    out["radial_velocity"] = (
        out["radial_velocity_true"].to_numpy(float) + rng.normal(0.0, rv_error, size=n)
    )
    out["grvs_mag"] = out["phot_rvs_mean_mag"].to_numpy(float)
    out["rv_template_teff"] = out["teff_val"].to_numpy(float)
    out["rv_template_logg"] = out["logg"].to_numpy(float)
    out["parallax_over_error"] = out["parallax"].to_numpy(float) / parallax_error
    out["parallax_zpcorr_over_error"] = out["parallax_zpcorr"].to_numpy(float) / parallax_error
    for col in ["parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr"]:
        out[col] = 0.0
    return out


def add_distance_posterior(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    plx = out["parallax_zpcorr"].to_numpy(float)
    plx_err = out["parallax_zpcorr_error"].to_numpy(float)
    med = np.where(plx > 0, 1000.0 / plx, np.nan)
    sigma = np.abs(med * (plx_err / plx))
    lo = np.maximum(1.0, med - sigma)
    hi = med + sigma
    out["distance_inv_zpcorr_pc"] = med
    out["rpgeo"] = med
    out["rpgeo_lo"] = lo
    out["rpgeo_hi"] = hi
    out["dist_pc_final_screen"] = med
    out["dist_lo_pc_final_screen"] = lo
    out["dist_hi_pc_final_screen"] = hi
    out["dist_source_final_screen"] = "bailer_jones_2021_high_snr_surrogate"
    return out


def add_rvs_availability(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    from gaiaunlimited.selectionfunctions import DR3RVSSelectionFunction

    out = df.copy()
    sf = DR3RVSSelectionFunction()
    coords = SkyCoord(
        ra=out["ra"].to_numpy(float) * u.deg,
        dec=out["dec"].to_numpy(float) * u.deg,
        frame="icrs",
    )
    grvs = out["grvs_mag"].to_numpy(float)
    bp_rp = out["bp_rp"].to_numpy(float)
    bp_rp_med = float(np.nanmedian(bp_rp))
    bp_rp_filled = np.where(np.isfinite(bp_rp), bp_rp, bp_rp_med)
    sf_value = np.asarray(sf.query(coords, g=grvs, c=bp_rp_filled, fill_nan=True), dtype=float)
    sf_invalid = ~np.isfinite(sf_value) | (sf_value <= 0) | ~np.isfinite(grvs) | ~np.isfinite(bp_rp)
    sf_for_draw = np.clip(np.where(sf_invalid, 0.0, sf_value), 0.0, 1.0)
    out["sf_value_v1convention"] = sf_value
    out["sf_weight_v1convention"] = 1.0 / np.clip(sf_for_draw, SF_FLOOR, 1.0)
    out["sf_invalid"] = sf_invalid
    out["rvs_available_draw"] = rng.uniform(0.0, 1.0, size=len(out)) < sf_for_draw
    return out


def prepare_injection(seed_rows: int, n_per_bin: int, force_tap: bool) -> pd.DataFrame:
    injected_path = PHASE15 / f"wp7_injected_mock_n{n_per_bin}.fits"
    if injected_path.exists() and not force_tap:
        log(f"using cached injected mock {rel(injected_path)}")
        return Table.read(injected_path).to_pandas()

    rng = np.random.default_rng(RNG_SEED)
    seed = clean_seed(fetch_gedr3mock_seed(seed_rows, force_tap))
    log(f"using {len(seed):,} seed rows after observability cuts")
    injected = sample_base_rows(seed, n_per_bin, rng)
    injected = inject_true_velocities(injected, rng)
    injected = add_measurement_noise(injected, rng)
    injected = add_distance_posterior(injected)
    injected = add_rvs_availability(injected, rng)

    table = Table.from_pandas(injected)
    table.write(injected_path, overwrite=True)
    shutil.copy2(injected_path, OUT_REVIEW / injected_path.name)
    log(f"wrote injected mock {rel(injected_path)}")
    return injected


def run_pipeline_mc(injected: pd.DataFrame, n_samples: int, star_batch: int, force_mc: bool) -> pd.DataFrame:
    per_star_path = PHASE15 / f"wp7_recovery_per_star_s{n_samples}.csv"
    if per_star_path.exists() and not force_mc:
        log(f"using cached MC recovery table {rel(per_star_path)}")
        return pd.read_csv(per_star_path)

    phase0e = load_phase0e_module()
    df = injected.copy()
    df["parallax_quality_pass"] = df["parallax_zpcorr_over_error"].to_numpy(float) > 5.0
    df["rvs_quality_pass"] = df["rvs_available_draw"].astype(bool)
    pass_mask = df["parallax_quality_pass"] & df["rvs_quality_pass"]
    work = df.loc[pass_mask].reset_index(drop=True).copy()
    log(
        f"pipeline filters: RVS draw {int(df['rvs_quality_pass'].sum()):,}/{len(df):,}, "
        f"parallax SNR {int(df['parallax_quality_pass'].sum()):,}/{len(df):,}, "
        f"joint {len(work):,}/{len(df):,}"
    )

    point = np.full(len(work), np.nan)
    P = np.full(len(work), np.nan)
    tier = np.full(len(work), "X", dtype=object)
    if len(work):
        dist_med = work["dist_pc_final_screen"].to_numpy(float)
        dist_lo = work["dist_lo_pc_final_screen"].to_numpy(float)
        dist_hi = work["dist_hi_pc_final_screen"].to_numpy(float)
        log(f"running phase0e point_vgrf on {len(work):,} filtered mock stars")
        point = phase0e.point_vgrf(work, dist_med, batch_size=5000)
        work["vgrf_default_exact"] = point
        log(f"running phase0e mc_pass with {n_samples} realisations")
        P = phase0e.mc_pass(
            work,
            dist_med,
            dist_lo,
            dist_hi,
            n_samples,
            np.random.default_rng(RNG_SEED + 1),
            star_batch,
        )
        tier = phase0e.assign_tiers(P, point)
    work["vgrf_default_exact"] = point
    work["P_vgrf_below_25"] = P
    work["mc_realisations"] = n_samples
    work["tier"] = tier
    work["recovered_slow"] = np.isin(tier, ["A", "B", "C"])

    recovered_cols = [
        "wp7_mock_id",
        "vgrf_default_exact",
        "P_vgrf_below_25",
        "mc_realisations",
        "tier",
        "recovered_slow",
    ]
    merged = df.merge(work[recovered_cols], on="wp7_mock_id", how="left")
    merged["recovered_slow"] = merged["recovered_slow"].fillna(False).astype(bool)
    merged["P_vgrf_below_25"] = merged["P_vgrf_below_25"].fillna(0.0)
    merged["tier"] = merged["tier"].fillna("X")
    merged["mc_realisations"] = merged["mc_realisations"].fillna(0).astype(int)
    merged.to_csv(per_star_path, index=False)
    merged.to_csv(OUT_REVIEW / per_star_path.name, index=False)
    log(f"wrote per-star recovery table {rel(per_star_path)}")
    return merged


def summarize_by_decade(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lo, hi in VGRF_DECADES:
        label = f"{int(lo)}-{int(hi)}"
        m = df["vgrf_bin"].eq(label).to_numpy()
        n_true = int(m.sum())
        n_rvs = int((m & df["rvs_quality_pass"].to_numpy(bool)).sum())
        n_pi = int((m & df["parallax_quality_pass"].to_numpy(bool)).sum())
        n_joint = int((m & df["rvs_quality_pass"].to_numpy(bool) & df["parallax_quality_pass"].to_numpy(bool)).sum())
        n_rec = int((m & df["recovered_slow"].to_numpy(bool)).sum())
        rows.append(
            {
                "vgrf_bin": label,
                "true_N": n_true,
                "rvs_available_N": n_rvs,
                "parallax_snr_N": n_pi,
                "pipeline_input_N": n_joint,
                "recovered_N": n_rec,
                "recovery_fraction": n_rec / n_true if n_true else np.nan,
                "pipeline_input_fraction": n_joint / n_true if n_true else np.nan,
                "median_sf_value": float(np.nanmedian(df.loc[m, "sf_value_v1convention"])),
                "median_G": float(np.nanmedian(df.loc[m, "phot_g_mean_mag"])),
                "median_logg": float(np.nanmedian(df.loc[m, "logg"])),
            }
        )
    return pd.DataFrame(rows)


def label_l(l_value: float) -> str:
    if l_value >= 330 or l_value < 30:
        return "inner"
    if l_value < 90:
        return "q1"
    if l_value < 180:
        return "q2"
    if l_value < 270:
        return "q3"
    return "q4"


def binned_slow_recovery(df: pd.DataFrame) -> pd.DataFrame:
    slow = df[df["vgrf_bin"].eq("0-25")].copy()
    slow["l_bin"] = [label_l(float(x)) for x in slow["l"]]
    slow["b_bin"] = pd.cut(
        slow["b"].astype(float),
        bins=[-90, -30, -10, 10, 30, 90],
        labels=["<-30", "-30:-10", "-10:10", "10:30", ">30"],
        include_lowest=True,
    )
    slow["G_bin"] = pd.cut(
        slow["phot_g_mean_mag"].astype(float),
        bins=[0, 10, 12, 14, 15, 16, 17],
        labels=["<10", "10:12", "12:14", "14:15", "15:16", "16:17"],
        include_lowest=True,
    )
    slow["logg_bin"] = pd.cut(
        slow["logg"].astype(float),
        bins=[-1, 2, 3.5, 4.5, 6],
        labels=["<2", "2:3.5", "3.5:4.5", ">4.5"],
        include_lowest=True,
    )
    group_cols = ["l_bin", "b_bin", "G_bin", "logg_bin"]
    grouped = slow.groupby(group_cols, observed=True)
    out = grouped.agg(
        true_N=("wp7_mock_id", "size"),
        recovered_N=("recovered_slow", "sum"),
        rvs_available_N=("rvs_quality_pass", "sum"),
        parallax_snr_N=("parallax_quality_pass", "sum"),
        median_sf_value=("sf_value_v1convention", "median"),
    ).reset_index()
    out["recovery_fraction"] = out["recovered_N"] / out["true_N"]
    return out.sort_values(group_cols).reset_index(drop=True)


def fmt_int(value: int | float) -> str:
    return f"{int(value):,}"


def fmt_frac(value: float) -> str:
    if value == 0:
        return "0.000"
    if value < 0.001:
        return f"{value:.4f}"
    return f"{value:.3f}"


def write_latex_table(summary: pd.DataFrame) -> None:
    lines = [
        r"\begin{center}",
        r"\refstepcounter{table}\label{tab:recovery_summary}",
        r"\begin{minipage}{\columnwidth}",
        r"\small\textbf{Table \thetable.} GeDR3mock injection--recovery summary.  Each true-$V_{\rm grf}$ bin contains controlled kinematic injections into GeDR3mock sightlines, photometry, stellar parameters, and Gaia-like uncertainties after the initial $G<17$, $G_{\rm RVS}<14$, and $3500<T_{\rm eff}<7000$ K observability cuts.  $N_{\rm RVS}$ is the number passing a stochastic GaiaUnlimited DR3-RVS draw evaluated with the same v1 convention used for the catalogue, $g=G_{\rm RVS}$ and $c=BP{-}RP$; $N_{\varpi}$ is the number passing $\varpi_{\rm corr}/\sigma_\varpi>5$; $N_{\rm input}$ is the intersection entering the 500-realisation velocity-threshold Monte Carlo; and $N_{\rm rec}$ is the number recovered as Tier A, B, or C ($P(V_{\rm grf}<25\,{\rm km\,s^{-1}})>0.5$).",
        r"\end{minipage}",
        r"\vspace{0.4ex}",
        "",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{4pt}",
        r"\resizebox{\columnwidth}{!}{%",
        r"\begin{tabular}{lrrrrr}",
        r"\hline\hline",
        r"True $V_{\rm grf}$ bin & $N_{\rm true}$ & $N_{\rm RVS}$ & $N_{\varpi}$ & $N_{\rm input}$ & $N_{\rm rec}/N_{\rm true}$ \\",
        r"\hline",
    ]
    for _, row in summary.iterrows():
        frac = row["recovery_fraction"]
        lines.append(
            f"{row['vgrf_bin']}\\kms & {fmt_int(row['true_N'])} & "
            f"{fmt_int(row['rvs_available_N'])} & {fmt_int(row['parallax_snr_N'])} & "
            f"{fmt_int(row['pipeline_input_N'])} & {fmt_int(row['recovered_N'])}/"
            f"{fmt_int(row['true_N'])} ({fmt_frac(float(frac))}) \\\\"
        )
    lines.extend(
        [
            r"\hline",
            r"\end{tabular}}",
            r"\end{center}",
            "",
        ]
    )
    RELEASE_TABLE.write_text("\n".join(lines), encoding="utf-8")
    shutil.copy2(RELEASE_TABLE, BUNDLE_TABLE)
    log(f"wrote {rel(RELEASE_TABLE)} and bundle mirror")


def make_heatmap(df: pd.DataFrame) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    slow = df[df["vgrf_bin"].eq("0-25")].copy()
    slow["logd_pc"] = np.log10(slow["dist_true_pc"].astype(float))
    g_edges = np.array([6, 10, 12, 13, 14, 15, 16, 17], dtype=float)
    d_edges = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5], dtype=float)
    true_counts, _, _ = np.histogram2d(
        slow["phot_g_mean_mag"].astype(float), slow["logd_pc"], bins=[g_edges, d_edges]
    )
    rec_counts, _, _ = np.histogram2d(
        slow.loc[slow["recovered_slow"], "phot_g_mean_mag"].astype(float),
        slow.loc[slow["recovered_slow"], "logd_pc"],
        bins=[g_edges, d_edges],
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = rec_counts / true_counts
    eff[true_counts == 0] = np.nan

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    mesh = ax.pcolormesh(g_edges, d_edges, eff.T, vmin=0.0, vmax=1.0, cmap="viridis", shading="flat")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(r"recovery efficiency")
    for i in range(len(g_edges) - 1):
        for j in range(len(d_edges) - 1):
            n = int(true_counts[i, j])
            if n > 0:
                ax.text(
                    0.5 * (g_edges[i] + g_edges[i + 1]),
                    0.5 * (d_edges[j] + d_edges[j + 1]),
                    f"{eff[i, j]:.2f}\n({n})",
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="white" if eff[i, j] < 0.55 else "black",
                )
    ax.set_xlabel(r"$G$ (mag)")
    ax.set_ylabel(r"$\log_{10}(d_{\rm true}/{\rm pc})$")
    ax.set_title(r"GeDR3mock injections: true $0<V_{\rm grf}<25$ km s$^{-1}$")
    ax.set_xlim(g_edges[0], g_edges[-1])
    ax.set_ylim(d_edges[0], d_edges[-1])
    fig.tight_layout()
    fig.savefig(RELEASE_FIG)
    plt.close(fig)
    shutil.copy2(RELEASE_FIG, BUNDLE_FIG)
    log(f"wrote {rel(RELEASE_FIG)} and bundle mirror")

    heat_rows = []
    for i in range(len(g_edges) - 1):
        for j in range(len(d_edges) - 1):
            heat_rows.append(
                {
                    "G_lo": float(g_edges[i]),
                    "G_hi": float(g_edges[i + 1]),
                    "logd_lo": float(d_edges[j]),
                    "logd_hi": float(d_edges[j + 1]),
                    "true_N": int(true_counts[i, j]),
                    "recovered_N": int(rec_counts[i, j]),
                    "recovery_fraction": None if true_counts[i, j] == 0 else float(eff[i, j]),
                }
            )
    pd.DataFrame(heat_rows).to_csv(OUT_REVIEW / "wp7_recovery_heatmap_bins.csv", index=False)
    return {
        "G_edges": g_edges.tolist(),
        "logd_edges": d_edges.tolist(),
        "true_counts": true_counts.astype(int).tolist(),
        "recovered_counts": rec_counts.astype(int).tolist(),
        "efficiency": np.where(np.isfinite(eff), eff, np.nan).tolist(),
    }


def write_outputs(df: pd.DataFrame) -> dict[str, Any]:
    summary = summarize_by_decade(df)
    binned = binned_slow_recovery(df)
    summary.to_csv(OUT_REVIEW / "wp7_recovery_summary.csv", index=False)
    binned.to_csv(OUT_REVIEW / "wp7_recovery_l_b_G_logg_bins.csv", index=False)
    write_latex_table(summary)
    heatmap = make_heatmap(df)

    slow = summary[summary["vgrf_bin"].eq("0-25")].iloc[0]
    payload: dict[str, Any] = {
        "seed": RNG_SEED,
        "n_per_true_vgrf_bin": int(summary["true_N"].iloc[0]),
        "mc_realisations": int(df["mc_realisations"].max()),
        "summary_by_true_vgrf_bin": summary.to_dict(orient="records"),
        "slow_0_25_recovery_fraction": float(slow["recovery_fraction"]),
        "slow_0_25_pipeline_input_fraction": float(slow["pipeline_input_fraction"]),
        "rvs_selection_function_convention": "DR3RVSSelectionFunction.query(g=phot_rvs_mean_mag, c=bp_rp), preserving v1 catalogue convention",
        "distance_stage_note": "Mock source_ids have no Bailer-Jones catalogue rows; the script supplies the high-SNR split-normal parallax posterior consumed by the existing v2 MC sampler, which is the BJ21 high-SNR limit after the parallax_over_error>5 filter.",
        "geDR3mock_reference": "Rybizki et al. 2020, PASP, 132, 074501",
        "heatmap": heatmap,
    }
    summary_json = OUT_REVIEW / "wp7_recovery_summary.json"
    summary_json.write_text(json.dumps(payload, indent=2, allow_nan=True), encoding="utf-8")
    shutil.copy2(summary_json, PHASE15 / "wp7_recovery_summary.json")
    log(f"wrote summary {rel(summary_json)}")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-per-bin", type=int, default=10_000)
    parser.add_argument("--seed-rows", type=int, default=60_000)
    parser.add_argument("--mc-realisations", type=int, default=500)
    parser.add_argument("--star-batch", type=int, default=150)
    parser.add_argument("--force-tap", action="store_true")
    parser.add_argument("--force-mc", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_dirs()
    injected = prepare_injection(args.seed_rows, args.n_per_bin, args.force_tap)
    recovered = run_pipeline_mc(injected, args.mc_realisations, args.star_batch, args.force_mc)
    payload = write_outputs(recovered)
    log("WP-7 complete")
    print(json.dumps(payload["summary_by_true_vgrf_bin"], indent=2))
    return 0


if __name__ == "__main__":
    T0 = time.time()
    raise SystemExit(main())
