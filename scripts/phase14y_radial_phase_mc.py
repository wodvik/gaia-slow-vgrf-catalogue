"""Phase 14Y -- present-day Galactocentric radial-phase MC diagnostic.

This diagnostic intentionally uses only the present-day Gaia phase-space
sampling.  It does not integrate orbits and does not depend on the adopted
Galactic potential.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table


BUNDLE = Path(__file__).resolve().parents[1]
INPUT = BUNDLE / "private_inputs" / "expanded_candidates_mc_tiered.csv"
ORBIT_CATALOGUE = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.fits"
TIER_AB_CATALOGUE = BUNDLE / "catalogues" / "catalogue_expanded_tierAB.fits"

# --- population-prior retier switch (Phase 16F) ---
import os as _os
if _os.environ.get("GAIA_RETIER", "").lower() in ("1", "true", "yes"):
    ORBIT_CATALOGUE = BUNDLE / "catalogues" / "catalogue_retier_orbits_tierABC.fits"
    TIER_AB_CATALOGUE = BUNDLE / "catalogues" / "catalogue_retier_tierAB.fits"

CONFIG = yaml.safe_load((BUNDLE / "config.yml").read_text(encoding="utf-8-sig"))
OUT = BUNDLE / "phase14" / "radial_phase"
SUMMARY_JSON = OUT / "radial_phase_mc_summary.json"
TABLE_TEX = BUNDLE / "tables" / "v15" / "tab_radial_phase.tex"

TIERS = ("A", "B", "C")
NEAR_ZERO_KMS = 2.0
DEFAULT_N_SAMP = 5000
DEFAULT_CHUNK_SIZE = 64
SEED = int(CONFIG["mc"]["random_seed"]) + 140_025
T0 = time.time()


def _rel(p) -> str:
    """Bundle-relative path, so released sidecars carry no absolute local path."""
    from pathlib import Path as _P
    try:
        return _P(p).resolve().relative_to(BUNDLE.resolve()).as_posix()
    except ValueError:
        return _P(p).name


def log(message: str) -> None:
    print(f"[14Y-radial-phase t={time.time() - T0:7.1f}s] {message}", flush=True)


def galcen_frame() -> coord.Galactocentric:
    sv = CONFIG["solar_variants"]["default"]
    return coord.Galactocentric(
        galcen_distance=sv["R0_kpc"] * u.kpc,
        z_sun=sv["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            sv["U_kms"] * u.km / u.s,
            (sv["Vc_kms"] + sv["V_kms"]) * u.km / u.s,
            sv["W_kms"] * u.km / u.s,
        ),
    )


def cov3(df: pd.DataFrame) -> np.ndarray:
    sig = np.stack(
        [
            df["parallax_error"].to_numpy(dtype=float),
            df["pmra_error"].to_numpy(dtype=float),
            df["pmdec_error"].to_numpy(dtype=float),
        ],
        axis=-1,
    )
    rho = np.zeros((len(df), 3, 3), dtype=float)
    for i in (0, 1, 2):
        rho[:, i, i] = 1.0
    rho[:, 0, 1] = rho[:, 1, 0] = df["parallax_pmra_corr"].to_numpy(dtype=float)
    rho[:, 0, 2] = rho[:, 2, 0] = df["parallax_pmdec_corr"].to_numpy(dtype=float)
    rho[:, 1, 2] = rho[:, 2, 1] = df["pmra_pmdec_corr"].to_numpy(dtype=float)
    return rho * (sig[:, :, None] * sig[:, None, :])


def chol3(cov: np.ndarray) -> np.ndarray:
    out = np.zeros_like(cov)
    eye = np.eye(3)
    for i in range(cov.shape[0]):
        try:
            out[i] = np.linalg.cholesky(cov[i])
        except np.linalg.LinAlgError:
            out[i] = np.linalg.cholesky(cov[i] + 1e-9 * eye)
    return out


def split_normal(med: np.ndarray, lo: np.ndarray, hi: np.ndarray, n_samp: int, rng: np.random.Generator) -> np.ndarray:
    sig_lo = np.maximum(med - lo, 1.0)
    sig_hi = np.maximum(hi - med, 1.0)
    z = rng.standard_normal((len(med), n_samp))
    sigma = np.where(z < 0.0, sig_lo[:, None], sig_hi[:, None])
    return np.maximum(med[:, None] + z * sigma, 1.0)


def sample_current_phase_space(df: pd.DataFrame, n_samp: int, rng: np.random.Generator) -> np.ndarray:
    mu = np.stack(
        [
            df["parallax_zpcorr"].to_numpy(dtype=float),
            df["pmra"].to_numpy(dtype=float),
            df["pmdec"].to_numpy(dtype=float),
        ],
        axis=-1,
    )
    z = rng.standard_normal((len(df), n_samp, 3))
    ast = mu[:, None, :] + np.einsum("nij,nsj->nsi", chol3(cov3(df)), z)
    rv = (
        df["radial_velocity"].to_numpy(dtype=float)[:, None]
        + rng.standard_normal((len(df), n_samp))
        * df["radial_velocity_error"].to_numpy(dtype=float)[:, None]
    )
    dist = split_normal(
        df["dist_pc_final_screen"].to_numpy(dtype=float),
        df["dist_lo_pc_final_screen"].to_numpy(dtype=float),
        df["dist_hi_pc_final_screen"].to_numpy(dtype=float),
        n_samp,
        rng,
    )
    ra = np.broadcast_to(df["ra"].to_numpy(dtype=float)[:, None], (len(df), n_samp))
    dec = np.broadcast_to(df["dec"].to_numpy(dtype=float)[:, None], (len(df), n_samp))

    icrs = coord.SkyCoord(
        ra=ra.ravel() * u.deg,
        dec=dec.ravel() * u.deg,
        distance=dist.ravel() * u.pc,
        pm_ra_cosdec=ast[:, :, 1].ravel() * u.mas / u.yr,
        pm_dec=ast[:, :, 2].ravel() * u.mas / u.yr,
        radial_velocity=rv.ravel() * u.km / u.s,
        frame="icrs",
    )
    g = icrs.transform_to(galcen_frame())
    ic = np.column_stack(
        [
            g.x.to_value(u.kpc),
            g.y.to_value(u.kpc),
            g.z.to_value(u.kpc),
            g.v_x.to_value(u.km / u.s),
            g.v_y.to_value(u.km / u.s),
            g.v_z.to_value(u.km / u.s),
        ]
    )
    return ic.reshape(len(df), n_samp, 6)


def radial_velocity_sph(ic: np.ndarray) -> np.ndarray:
    x = ic[..., 0]
    y = ic[..., 1]
    z = ic[..., 2]
    vx = ic[..., 3]
    vy = ic[..., 4]
    vz = ic[..., 5]
    r = np.sqrt(x * x + y * y + z * z)
    return (x * vx + y * vy + z * vz) / r


def load_expanded() -> pd.DataFrame:
    cols = [
        "source_id", "tier", "P_vgrf_below_25", "mc_realisations",
        "vgrf_default_exact", "ra", "dec", "parallax_zpcorr",
        "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error",
        "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr",
        "radial_velocity", "radial_velocity_error",
        "dist_pc_final_screen", "dist_lo_pc_final_screen",
        "dist_hi_pc_final_screen",
    ]
    df = pd.read_csv(INPUT, usecols=cols)
    df["tier"] = df["tier"].astype(str).str.strip()
    df["source_id"] = df["source_id"].astype("int64")
    # The candidate CSV carries the FORWARD tier. Under GAIA_RETIER the Monte
    # Carlo half of this table must use the same population-prior sample as the
    # point-estimate half, or the two columns describe different star sets.
    if _os.environ.get("GAIA_RETIER", "").lower() in ("1", "true", "yes"):
        _rt = Table.read(BUNDLE / "catalogues" / "catalogue_retier_master.fits")
        _map = dict(zip(np.asarray(_rt["source_id"]).astype("int64").tolist(),
                        np.asarray(_rt["tier"]).astype(str).tolist()))
        df["tier"] = [_map.get(int(s), "X") for s in df["source_id"].to_numpy()]
    df = df[df["tier"].isin(TIERS)].copy()
    return df.sort_values("source_id").reset_index(drop=True)


def gold_source_ids() -> set[int]:
    df = Table.read(TIER_AB_CATALOGUE).to_pandas()
    frac_dist = (df["dist_hi_pc"].to_numpy(float) - df["dist_lo_pc"].to_numpy(float)) / (
        2.0 * df["dist_pc"].to_numpy(float)
    )
    mask = (
        (df["ruwe"].to_numpy(float) < 1.4)
        & df["rvs_quality_ok"].astype(bool).to_numpy()
        & (df["parallax_over_error"].to_numpy(float) > 10.0)
        & (frac_dist < 0.15)
    )
    return set(df.loc[mask, "source_id"].astype("int64"))


def point_estimate_counts(gold_ids: set[int]) -> dict[str, dict[str, float]]:
    df = Table.read(ORBIT_CATALOGUE).to_pandas()
    df["tier"] = df["tier"].astype(str).str.strip()
    x = df["x_kpc"].to_numpy(float)
    y = df["y_kpc"].to_numpy(float)
    z = df["z_kpc"].to_numpy(float)
    vx = df["vx_kms"].to_numpy(float)
    vy = df["vy_kms"].to_numpy(float)
    vz = df["vz_kms"].to_numpy(float)
    df["v_radial_sph_kms"] = (x * vx + y * vy + z * vz) / np.sqrt(x * x + y * y + z * z)

    specs = {
        "tier_ab": df["tier"].isin(["A", "B"]).to_numpy(),
        "gold_tier_ab": df["source_id"].astype("int64").isin(gold_ids).to_numpy(),
        "tier_abc": df["tier"].isin(["A", "B", "C"]).to_numpy(),
    }
    out: dict[str, dict[str, float]] = {}
    for key, mask in specs.items():
        v = df.loc[mask, "v_radial_sph_kms"].to_numpy(float)
        inward = int(np.sum(v < -NEAR_ZERO_KMS))
        near = int(np.sum(np.abs(v) <= NEAR_ZERO_KMS))
        outward = int(np.sum(v > NEAR_ZERO_KMS))
        nonzero = inward + outward
        out[key] = {
            "n": int(mask.sum()),
            "inward": inward,
            "near_zero": near,
            "outward": outward,
            "nonzero": nonzero,
            "inward_fraction_nonzero": float(inward / nonzero),
        }
    return out


def empty_agg(n_samp: int) -> dict[str, np.ndarray]:
    return {
        "inward": np.zeros(n_samp, dtype=np.int32),
        "near_zero": np.zeros(n_samp, dtype=np.int32),
        "outward": np.zeros(n_samp, dtype=np.int32),
    }


def mc_counts(df: pd.DataFrame, gold_ids: set[int], n_samp: int, chunk_size: int) -> dict[str, dict[str, np.ndarray]]:
    rng = np.random.default_rng(SEED)
    masks = {
        "tier_ab": df["tier"].isin(["A", "B"]).to_numpy(),
        "gold_tier_ab": df["source_id"].isin(gold_ids).to_numpy(),
        "tier_abc": df["tier"].isin(["A", "B", "C"]).to_numpy(),
    }
    agg = {key: empty_agg(n_samp) for key in masks}
    for start in range(0, len(df), chunk_size):
        stop = min(start + chunk_size, len(df))
        chunk = df.iloc[start:stop]
        log(f"sampling current phase space rows {start}-{stop} of {len(df)}")
        vr = radial_velocity_sph(sample_current_phase_space(chunk, n_samp, rng))
        inward = vr < -NEAR_ZERO_KMS
        near = np.abs(vr) <= NEAR_ZERO_KMS
        outward = vr > NEAR_ZERO_KMS
        for key, mask in masks.items():
            local = mask[start:stop]
            if not np.any(local):
                continue
            agg[key]["inward"] += inward[local].sum(axis=0, dtype=np.int32)
            agg[key]["near_zero"] += near[local].sum(axis=0, dtype=np.int32)
            agg[key]["outward"] += outward[local].sum(axis=0, dtype=np.int32)
    return agg


def summarise_mc(agg: dict[str, np.ndarray]) -> dict[str, float | list[float]]:
    inward = agg["inward"]
    near = agg["near_zero"]
    outward = agg["outward"]
    nonzero = inward + outward
    frac = inward / nonzero
    diff = inward - outward
    q = np.percentile(frac, [16, 50, 84])
    return {
        "inward_count_p16_p50_p84": [float(x) for x in np.percentile(inward, [16, 50, 84])],
        "near_zero_count_p16_p50_p84": [float(x) for x in np.percentile(near, [16, 50, 84])],
        "outward_count_p16_p50_p84": [float(x) for x in np.percentile(outward, [16, 50, 84])],
        "inward_fraction_nonzero_p16_p50_p84": [float(x) for x in q],
        "p_inward_count_gt_outward_count": float(np.mean(diff > 0)),
        "p_inward_count_eq_outward_count": float(np.mean(diff == 0)),
    }


def fmt_count(value: float) -> str:
    return f"{int(round(value)):,}".replace(",", r"{,}")


def fmt_int(value: int) -> str:
    return f"{int(value):,}".replace(",", r"{,}")


def fmt_pct_triplet(values: list[float]) -> str:
    lo, med, hi = [100.0 * float(x) for x in values]
    return rf"${med:.1f}^{{+{hi - med:.1f}}}_{{-{med - lo:.1f}}}\%$"


def fmt_prob(value: float) -> str:
    if value >= 0.9995:
        return "$>0.999$"
    if value <= 0.0005:
        return "$<0.001$"
    return f"{value:.3f}"


def write_table(summary: dict) -> None:
    rows = [
        ("tier_ab", r"Tier~A+B"),
        ("gold_tier_ab", r"Gold A+B"),
        ("tier_abc", r"Tier~A+B+C"),
    ]
    lines = [
        r"\begin{center}",
        r"\refstepcounter{table}\label{tab:radial_phase}",
        r"\begin{minipage}{\columnwidth}",
        r"\small\textbf{Table \thetable.} Frame-dependent present-day Galactocentric radial-phase check for the low-\vgrf{} catalogue.",
        r"\end{minipage}",
        r"\vspace{0.4ex}",
        "",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{2.4pt}",
        r"\begin{tabular}{lrrrrcc}",
        r"\hline\hline",
        r"Sample & $N$ & $N_{\rm in}$ & $N_{0}$ & $N_{\rm out}$ & MC $f_{\rm in}$ & $P_{\rm in>out}$ \\",
        r"\hline",
    ]
    for key, label in rows:
        point = summary["point_estimate"][key]
        mc = summary["mc"][key]
        lines.append(
            f"{label} & {fmt_int(point['n'])} & {fmt_int(point['inward'])} & {fmt_int(point['near_zero'])} & "
            f"{fmt_int(point['outward'])} & {fmt_pct_triplet(mc['inward_fraction_nonzero_p16_p50_p84'])} & "
            f"{fmt_prob(mc['p_inward_count_gt_outward_count'])} \\\\"
        )
    lines.extend(
        [
            r"\hline",
            r"\end{tabular}",
            r"\vspace{0.5ex}",
            "",
            r"\begin{minipage}{0.95\columnwidth}",
            r"\footnotesize",
            rf"\textit{{Note.}} $N_{{\rm in}}$, $N_0$, and $N_{{\rm out}}$ are point-estimate counts for spherical Galactocentric $v_r < -{NEAR_ZERO_KMS:g}$\kms, $|v_r|\leq {NEAR_ZERO_KMS:g}$\kms, and $v_r>{NEAR_ZERO_KMS:g}$\kms, respectively.  MC $f_{{\rm in}}$ is the median and 16th--84th percentile interval of $N_{{\rm in}}/(N_{{\rm in}}+N_{{\rm out}})$ from {fmt_int(summary['n_samp'])} present-day Gaia phase-space realisations per star.  This frame-dependent diagnostic uses no orbit integration or Galactic potential.  $P_{{\rm in>out}}$ is conditional on the measured catalogue and Gaia uncertainties; it is not a selection-function-corrected significance test.",
            r"\end{minipage}",
            r"\end{center}",
        ]
    )
    TABLE_TEX.write_text("\n".join(lines) + "\n")
    log(f"wrote {TABLE_TEX}")


def run(n_samp: int, chunk_size: int) -> dict:
    df = load_expanded()
    gold_ids = gold_source_ids()
    log(f"loaded Tier A+B+C N={len(df):,}; Gold Tier A+B N={len(gold_ids):,}")
    point = point_estimate_counts(gold_ids)
    mc_raw = mc_counts(df, gold_ids, n_samp, chunk_size)
    mc = {key: summarise_mc(value) for key, value in mc_raw.items()}
    summary = {
        "phase": "14Y",
        "input_csv": _rel(INPUT),
        "orbit_catalogue": _rel(ORBIT_CATALOGUE),
        "near_zero_kms": NEAR_ZERO_KMS,
        "n_samp": int(n_samp),
        "chunk_size": int(chunk_size),
        "seed": int(SEED),
        "point_estimate": point,
        "mc": mc,
        "outputs": {
            "summary_json": _rel(SUMMARY_JSON),
            "table_tex": _rel(TABLE_TEX),
        },
        "note": "Present-day spherical Galactocentric radial velocity; no orbit integration.",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    log(f"wrote {SUMMARY_JSON}")
    write_table(summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samp", type=int, default=DEFAULT_N_SAMP)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    args = parser.parse_args()
    summary = run(args.n_samp, args.chunk_size)
    for key, label in [("tier_ab", "Tier A+B"), ("gold_tier_ab", "Gold Tier A+B"), ("tier_abc", "Tier A+B+C")]:
        point = summary["point_estimate"][key]
        mc = summary["mc"][key]
        frac = mc["inward_fraction_nonzero_p16_p50_p84"]
        log(
            f"{label}: point in/0/out={point['inward']}/{point['near_zero']}/{point['outward']}; "
            f"MC f_in p16/p50/p84={frac[0]:.4f}/{frac[1]:.4f}/{frac[2]:.4f}; "
            f"P(N_in>N_out)={mc['p_inward_count_gt_outward_count']:.4f}"
        )


if __name__ == "__main__":
    main()
