"""
Apply Bailer-Jones distances to the Phase 0C parent-buffer shortlist.

The Phase 0C screen is intentionally conservative and uses L21-corrected
inverse-parallax distances. This script queries Bailer-Jones+2021
photogeometric distances only for that manageable shortlist, recomputes Vgrf
with BJ distances where available, and reports whether stars outside the old
legacy <25 km/s preselection enter the final point-estimate threshold.

Default inputs/outputs live under D:/GAIA/parent_scan.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.table import Table, vstack


REPO = Path(__file__).resolve().parents[3]
DEFAULT_SCREEN = Path("D:/GAIA/parent_scan/gaia_parent_buffer_zpcorr_inv_screen_vgrf50.csv")
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
COUNT_THRESHOLDS = [25, 50, 75, 100]


def query_bailer_jones(source_ids: np.ndarray, cache_path: Path, chunk_size: int) -> Table:
    if cache_path.exists():
        print(f"Using cached BJ table: {cache_path}")
        return Table.read(cache_path)

    import pyvo

    tap = pyvo.dal.TAPService("https://tapvizier.cds.unistra.fr/TAPVizieR/tap")
    chunks = []
    for i in range(0, len(source_ids), chunk_size):
        ids = source_ids[i : i + chunk_size]
        in_clause = ",".join(str(int(x)) for x in ids)
        adql = (
            "SELECT t.Source, t.rgeo, t.\"b_rgeo\" AS rgeo_lo, "
            "t.\"B_rgeo\" AS rgeo_hi, t.rpgeo, "
            "t.\"b_rpgeo\" AS rpgeo_lo, t.\"B_rpgeo\" AS rpgeo_hi, t.Flag "
            f"FROM \"I/352/gedr3dis\" AS t WHERE t.Source IN ({in_clause})"
        )
        print(f"BJ chunk {i}-{i + len(ids)} ({len(ids)} ids)", flush=True)
        res = tap.search(adql).to_table()
        chunks.append(res)
    out = vstack(chunks) if chunks else Table()
    out.write(cache_path, overwrite=True)
    return out


def compute_vgrf_with_distance(df: pd.DataFrame, distance_pc: np.ndarray) -> np.ndarray:
    d_kpc = distance_pc / 1000.0
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
    screen = pd.read_csv(args.screen_csv)
    screen["source_id"] = screen["source_id"].astype("int64")
    source_ids = screen["source_id"].drop_duplicates().to_numpy(dtype=np.int64)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_path = out_dir / "bailer_jones_2021_parent_screen_vgrf50.fits"

    bj = query_bailer_jones(source_ids, cache_path, args.query_chunk_size).to_pandas()
    bj = bj.rename(columns={"Source": "source_id"})
    bj["source_id"] = bj["source_id"].astype("int64")
    df = screen.merge(bj, on="source_id", how="left")

    bj_ok = df["rpgeo"].notna().to_numpy()
    dist_pc = df["distance_inv_zpcorr_pc"].to_numpy(dtype=float).copy()
    dist_source = np.full(len(df), "inv_parallax_zpcorr", dtype=object)
    dist_pc[bj_ok] = df.loc[bj_ok, "rpgeo"].to_numpy(dtype=float)
    dist_source[bj_ok] = "bailer_jones_2021_photogeo"

    finite = np.isfinite(dist_pc) & (dist_pc > 0)
    vgrf = np.full(len(df), np.nan, dtype=float)
    if finite.any():
        vgrf[finite] = compute_vgrf_with_distance(df.loc[finite], dist_pc[finite])

    df["dist_pc_phase0d"] = dist_pc
    df["dist_source_phase0d"] = dist_source
    df["vgrf_bj_or_inv"] = vgrf

    old = df["source_in_old_preselection"].astype(bool).to_numpy()
    finite_v = np.isfinite(vgrf)
    counts = {}
    for threshold in COUNT_THRESHOLDS:
        below = finite_v & (vgrf < threshold)
        counts[f"bj_or_inv_lt_{threshold}"] = int(below.sum())
        counts[f"bj_or_inv_lt_{threshold}_outside_old"] = int((below & ~old).sum())
        counts[f"bj_or_inv_lt_{threshold}_inside_old"] = int((below & old).sum())

    out_csv = out_dir / "gaia_parent_buffer_bj_screen_vgrf50.csv"
    df.to_csv(out_csv, index=False)
    entrants = df[(df["vgrf_bj_or_inv"] < 25.0) & (~df["source_in_old_preselection"].astype(bool))]
    entrants_csv = out_dir / "gate0_bj_vgrf_lt25_outside_old.csv"
    entrants.to_csv(entrants_csv, index=False)

    summary = {
        "screen_csv": str(args.screen_csv),
        "output_csv": str(out_csv),
        "entrants_csv": str(entrants_csv),
        "n_screen_rows": int(len(df)),
        "n_unique_source_id": int(df["source_id"].nunique()),
        "n_bailer_jones_matches": int(bj_ok.sum()),
        "n_inv_parallax_fallback": int((~bj_ok).sum()),
        "counts": counts,
    }
    summary_path = out_dir / "gate0_bj_screen_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen-csv", default=str(DEFAULT_SCREEN))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--query-chunk-size", type=int, default=500)
    return parser.parse_args()


def main() -> None:
    print(json.dumps(run(parse_args()), indent=2))


if __name__ == "__main__":
    main()
