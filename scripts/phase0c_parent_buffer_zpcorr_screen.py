"""
Screen the full Gaia parent buffer under final parallax-zero-point conventions.

This is a conservative bridge between the broad legacy scan and the expensive
final Bailer-Jones/Monte-Carlo pipeline. It reads the de-duplicated
legacy Vgrf < 200 km/s parent buffer, applies the Lindegren et al. (2021)
parallax zero-point correction, computes a corrected inverse-parallax Vgrf,
and writes a shortlist for final-distance reprocessing.

Default output:
  D:/GAIA/parent_scan/gaia_parent_buffer_zpcorr_inv_screen_vgrf50.csv
  D:/GAIA/parent_scan/gate0_zpcorr_inv_screen_summary.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
DEFAULT_BUFFER = Path("D:/GAIA/parent_scan/gaia_parent_buffer_vgrf200_full_dedup.csv")
DEFAULT_OLD = REPO / "release" / "data" / "slow_stars_enriched_orbits.csv"
DEFAULT_OUT_DIR = Path("D:/GAIA/parent_scan")

K = 4.74047
U_SUN = 11.1
V_SUN = 12.24
W_SUN = 7.25
V_CIRC = 232.0

AG = np.array(
    [
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [+0.4941094279, -0.4448296300, +0.7469822445],
        [-0.8676661490, -0.1980763734, +0.4559837762],
    ]
)

ZPT_VALIDITY = {
    "G_min": 6.0,
    "G_max": 21.0,
    "nu_eff_min": 1.1,
    "nu_eff_max": 1.9,
    "pseudocolour_min": 1.24,
    "pseudocolour_max": 1.72,
}

COUNT_THRESHOLDS = [25, 50, 75, 100, 125, 150, 175, 200]


def load_old_ids(path: Path) -> set[int]:
    old = pd.read_csv(path, usecols=["source_id"])
    return set(old["source_id"].astype("int64").tolist())


def compute_zpcorr(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from zero_point import zpt

    g = df["phot_g_mean_mag"].to_numpy(dtype=float)
    nu = df["nu_eff_used_in_astrometry"].to_numpy(dtype=float)
    psc = df["pseudocolour"].to_numpy(dtype=float)
    ecl = df["ecl_lat"].to_numpy(dtype=float)
    asp = df["astrometric_params_solved"].to_numpy(dtype=int)
    plx = df["parallax"].to_numpy(dtype=float)

    is_5p = asp == 31
    is_6p = asp == 95
    asp_ok = is_5p | is_6p
    g_ok = (g > ZPT_VALIDITY["G_min"]) & (g < ZPT_VALIDITY["G_max"])
    nu_ok = (nu > ZPT_VALIDITY["nu_eff_min"]) & (nu < ZPT_VALIDITY["nu_eff_max"])
    ps_ok = (psc > ZPT_VALIDITY["pseudocolour_min"]) & (
        psc < ZPT_VALIDITY["pseudocolour_max"]
    )
    colour_ok = np.where(is_5p, nu_ok, np.where(is_6p, ps_ok, False))
    valid = asp_ok & g_ok & colour_ok & np.isfinite(plx)

    zpcorr_mas = np.full(len(df), np.nan, dtype=float)
    if asp_ok.any():
        g_in = np.where(np.isfinite(g), g, 12.0)
        nu_in = np.where(np.isfinite(nu), nu, 1.5)
        ps_in = np.where(np.isfinite(psc), psc, 1.5)
        ecl_in = np.where(np.isfinite(ecl), ecl, 0.0)
        asp_in = np.where(asp_ok, asp, 31)
        zpcorr_mas = zpt.get_zpt(
            g_in, nu_in, ps_in, ecl_in, asp_in, _warnings=False
        )
        zpcorr_mas = np.where(asp_ok, zpcorr_mas, np.nan)

    parallax_zpcorr = np.where(valid, plx - zpcorr_mas, plx)
    return zpcorr_mas, parallax_zpcorr, valid


def compute_vgrf_with_parallax(df: pd.DataFrame, parallax_mas: np.ndarray) -> np.ndarray:
    d_kpc = 1.0 / parallax_mas
    ra_rad = np.radians(df["ra"].to_numpy(dtype=float))
    dec_rad = np.radians(df["dec"].to_numpy(dtype=float))
    pmra = df["pmra"].to_numpy(dtype=float)
    pmdec = df["pmdec"].to_numpy(dtype=float)
    rv = df["radial_velocity"].to_numpy(dtype=float)

    v_ra = K * pmra * d_kpc
    v_dec = K * pmdec * d_kpc

    cos_ra = np.cos(ra_rad)
    sin_ra = np.sin(ra_rad)
    cos_dec = np.cos(dec_rad)
    sin_dec = np.sin(dec_rad)

    vx_eq = rv * cos_ra * cos_dec - v_ra * sin_ra - v_dec * cos_ra * sin_dec
    vy_eq = rv * sin_ra * cos_dec + v_ra * cos_ra - v_dec * sin_ra * sin_dec
    vz_eq = rv * sin_dec + v_dec * cos_dec

    v_gal = AG @ np.vstack([vx_eq, vy_eq, vz_eq])
    u_lsr = v_gal[0] + U_SUN
    v_lsr = v_gal[1] + V_SUN
    w_lsr = v_gal[2] + W_SUN
    v_grf = v_lsr + V_CIRC
    return np.sqrt(u_lsr**2 + v_grf**2 + w_lsr**2)


def run(args: argparse.Namespace) -> dict[str, Any]:
    from zero_point import zpt

    zpt.load_tables()
    buffer_csv = Path(args.buffer_csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"gaia_parent_buffer_zpcorr_inv_screen_vgrf{int(args.keep_kms)}.csv"
    if out_csv.exists():
        out_csv.unlink()

    old_ids = load_old_ids(Path(args.old_csv))
    counts = {f"zpcorr_inv_lt_{t}": 0 for t in COUNT_THRESHOLDS}
    counts_outside_old = {f"zpcorr_inv_lt_{t}_outside_old": 0 for t in COUNT_THRESHOLDS}
    n_rows = 0
    n_kept = 0
    n_old_kept = 0
    n_outside_old_kept = 0
    n_zpcorr_valid = 0
    min_v = np.inf
    max_v = -np.inf

    for i, chunk in enumerate(pd.read_csv(buffer_csv, chunksize=args.chunksize), start=1):
        chunk["source_id"] = chunk["source_id"].astype("int64")
        zpcorr_mas, parallax_zpcorr, zpcorr_valid = compute_zpcorr(chunk)
        ok = np.isfinite(parallax_zpcorr) & (parallax_zpcorr > 0)
        vgrf = np.full(len(chunk), np.nan, dtype=float)
        if ok.any():
            vgrf[ok] = compute_vgrf_with_parallax(chunk.loc[ok], parallax_zpcorr[ok])

        old_mask = chunk["source_id"].isin(old_ids).to_numpy()
        finite_v = np.isfinite(vgrf)
        for threshold in COUNT_THRESHOLDS:
            below = finite_v & (vgrf < threshold)
            counts[f"zpcorr_inv_lt_{threshold}"] += int(below.sum())
            counts_outside_old[f"zpcorr_inv_lt_{threshold}_outside_old"] += int(
                (below & ~old_mask).sum()
            )

        keep = old_mask | (finite_v & (vgrf < args.keep_kms))
        out = chunk.loc[keep].copy()
        if not out.empty:
            out["zpcorr_value_uas"] = zpcorr_mas[keep] * 1000.0
            out["zpcorr_valid"] = zpcorr_valid[keep]
            out["parallax_zpcorr"] = parallax_zpcorr[keep]
            out["distance_inv_zpcorr_pc"] = np.where(
                out["parallax_zpcorr"].to_numpy(dtype=float) > 0,
                1000.0 / out["parallax_zpcorr"].to_numpy(dtype=float),
                np.nan,
            )
            out["vgrf_zpcorr_inv"] = vgrf[keep]
            out["source_in_old_preselection"] = old_mask[keep]
            out.to_csv(out_csv, mode="a", index=False, header=not out_csv.exists())

        n_rows += len(chunk)
        n_kept += int(keep.sum())
        n_old_kept += int((keep & old_mask).sum())
        n_outside_old_kept += int((keep & ~old_mask).sum())
        n_zpcorr_valid += int(zpcorr_valid.sum())
        if finite_v.any():
            min_v = min(min_v, float(np.nanmin(vgrf[finite_v])))
            max_v = max(max_v, float(np.nanmax(vgrf[finite_v])))

        if i % 10 == 0:
            print(
                f"chunk {i}: scanned={n_rows:,} kept={n_kept:,} "
                f"outside_old_kept={n_outside_old_kept:,}",
                flush=True,
            )

    summary = {
        "buffer_csv": str(buffer_csv),
        "output_csv": str(out_csv),
        "keep_kms": float(args.keep_kms),
        "n_rows": int(n_rows),
        "n_zpcorr_valid": int(n_zpcorr_valid),
        "n_kept": int(n_kept),
        "n_old_preselection_kept": int(n_old_kept),
        "n_outside_old_kept": int(n_outside_old_kept),
        "vgrf_zpcorr_inv_min": None if not np.isfinite(min_v) else float(min_v),
        "vgrf_zpcorr_inv_max": None if not np.isfinite(max_v) else float(max_v),
        "counts": {**counts, **counts_outside_old},
    }
    summary_path = out_dir / "gate0_zpcorr_inv_screen_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--buffer-csv", default=str(DEFAULT_BUFFER))
    parser.add_argument("--old-csv", default=str(DEFAULT_OLD))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--keep-kms", type=float, default=50.0)
    parser.add_argument("--chunksize", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    summary = run(parse_args())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
