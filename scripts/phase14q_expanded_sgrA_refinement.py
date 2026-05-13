"""Phase 14Q -- expanded-catalogue Sgr A* approacher refinement.

This reruns the Sgr A* close-approach Monte Carlo for the expanded
parent-complete Tier A+B+C catalogue.  Candidate seeds are selected from
the regenerated point-estimate static and barred orbit products by
minimum spherical Galactocentric distance.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import agama
import numpy as np
import pandas as pd
from astropy.table import Table

import phase14d_sgrA_event_refinement as base


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "release/v2/phase14"
MASTER = REPO / "release/v2/phase0_expanded/catalogue_expanded_master.fits"
ORBITS = REPO / "release/v2/phase0_expanded/catalogue_expanded_orbits_tierABC.fits"
SRC = REPO / "release/data/slow_stars_expanded_candidates_vgrf50.csv"


def log(message: str) -> None:
    print(f"[14Q-expanded-SgrA t={time.time() - T0:7.1f}s] {message}", flush=True)


def clean_tier(series: pd.Series) -> pd.Series:
    if series.dtype == object and len(series) and isinstance(series.iloc[0], (bytes, bytearray)):
        return series.str.decode("utf-8").str.strip()
    return series.astype(str).str.strip()


def load_expanded_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    src_cols = [
        "source_id", "ra", "dec", "pmra", "pmdec", "pmra_error", "pmdec_error",
        "parallax_error", "parallax_pmra_corr", "parallax_pmdec_corr",
        "pmra_pmdec_corr", "radial_velocity", "radial_velocity_error",
    ]
    src = pd.read_csv(SRC, usecols=src_cols)
    master = Table.read(MASTER).to_pandas()
    master["tier"] = clean_tier(master["tier"])
    keep = [
        "source_id", "tier", "P_vgrf_below_25", "parallax_zpcorr",
        "dist_pc", "dist_lo_pc", "dist_hi_pc",
    ]
    df_all = src.merge(master[keep], on="source_id", how="inner", validate="one_to_one")

    orbits = Table.read(ORBITS).to_pandas()
    orbits["tier"] = clean_tier(orbits["tier"])
    return df_all, orbits


def select_candidates(orbits: pd.DataFrame, n_per_potential: int) -> list[dict]:
    rows: list[dict] = []
    seen: set[tuple[int, str]] = set()
    for potential in ("static", "barred"):
        col = f"{potential}_min_r_sph_kpc"
        top = orbits.sort_values(col).head(n_per_potential)
        for row in top.itertuples(index=False):
            sid = int(row.source_id)
            key = (sid, potential)
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "source_id": sid,
                "potential": potential,
                "tier": str(row.tier),
                "point_min_rsph_pc": float(getattr(row, col) * 1000.0),
                "point_R_peri_pc": float(getattr(row, f"{potential}_R_peri_kpc") * 1000.0),
            })
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-samp", type=int, default=5000)
    parser.add_argument("--chunk-size", type=int, default=250)
    parser.add_argument("--n-per-potential", type=int, default=5)
    args = parser.parse_args(argv)

    global T0
    T0 = time.time()
    base.T0 = T0
    OUT.mkdir(parents=True, exist_ok=True)
    log(f"agama {agama.__version__}")
    log(f"n_samp={args.n_samp}, chunk_size={args.chunk_size}")

    df_all, orbits = load_expanded_inputs()
    candidates = select_candidates(orbits, args.n_per_potential)
    log(f"selected {len(candidates)} candidate/potential rows")

    pot_axi = agama.Potential(file=str(base.WORK / "MWPotentialHunter24_axi.ini"))
    pot_full = agama.Potential(file=str(base.WORK / "MWPotentialHunter24_full.ini"))
    pot_bar = agama.Potential(potential=pot_full, rotation=base.ANGLE_BAR)
    rng = np.random.default_rng(base.SEED + 14000)

    draw_rows: list[dict] = []
    det_rows: list[dict] = []

    for cand in candidates:
        sid = int(cand["source_id"])
        pot_name = str(cand["potential"])
        tier = str(cand["tier"])
        pot = pot_axi if pot_name == "static" else pot_bar
        log(
            f"candidate {sid} ({pot_name}, Tier {tier}), "
            f"point min={cand['point_min_rsph_pc']:.2f} pc"
        )

        det_ic = base.deterministic_ic(sid, orbits)
        rows = base.integrate_minima(det_ic, pot, pot_name, sid, tier, chunk_size=1)
        for row in rows:
            row["draw"] = -1
            row["point_min_rsph_pc"] = cand["point_min_rsph_pc"]
            row["point_R_peri_pc"] = cand["point_R_peri_pc"]
        det_rows.extend(rows)

        obs = df_all.loc[df_all["source_id"].astype("int64") == sid].iloc[0]
        ic = base.sample_ics(obs, args.n_samp, rng)
        rows = base.integrate_minima(ic, pot, pot_name, sid, tier, args.chunk_size)
        for row in rows:
            row["point_min_rsph_pc"] = cand["point_min_rsph_pc"]
            row["point_R_peri_pc"] = cand["point_R_peri_pc"]
        draw_rows.extend(rows)

    draws_df = pd.DataFrame(draw_rows)
    det_df = pd.DataFrame(det_rows)
    draws_path = OUT / "expanded_sgrA_refinement_draws.csv"
    det_path = OUT / "expanded_sgrA_refinement_deterministic.csv"
    summary_path = OUT / "expanded_sgrA_refinement_summary.csv"
    json_path = OUT / "expanded_sgrA_refinement_summary.json"
    md_path = OUT / "expanded_sgrA_refinement_summary.md"
    draws_df.to_csv(draws_path, index=False)
    det_df.to_csv(det_path, index=False)

    summaries = []
    for cand in candidates:
        sid = int(cand["source_id"])
        pot_name = str(cand["potential"])
        sub = draws_df.loc[
            (draws_df["source_id"].astype("int64") == sid)
            & (draws_df["potential"] == pot_name)
        ].reset_index(drop=True)
        det = det_df.loc[
            (det_df["source_id"].astype("int64") == sid)
            & (det_df["potential"] == pot_name)
        ].reset_index(drop=True)
        row = base.summarize_candidate(sub, det)
        row["point_min_rsph_pc"] = float(cand["point_min_rsph_pc"])
        row["point_R_peri_pc"] = float(cand["point_R_peri_pc"])
        summaries.append(row)

    summary_df = pd.DataFrame(summaries).sort_values(
        ["interp_P_lt_100pc", "interp_min_rsph_pc_p50"], ascending=[False, True]
    )
    summary_df.to_csv(summary_path, index=False)

    counts = {
        "n_candidates": int(len(summary_df)),
        "n_P_lt_10pc_gt_0p5": int(np.sum(summary_df["interp_P_lt_10pc"] > 0.5)),
        "n_P_lt_100pc_gt_0p5": int(np.sum(summary_df["interp_P_lt_100pc"] > 0.5)),
        "n_P_lt_100pc_gt_0p1": int(np.sum(summary_df["interp_P_lt_100pc"] > 0.1)),
    }
    payload = {
        "phase": "14Q",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "expanded_catalogue": str(ORBITS),
        "n_realisations_per_candidate": int(args.n_samp),
        "n_per_potential_seeded": int(args.n_per_potential),
        "trajsize": base.TRAJSIZE,
        "grid_cadence_myr": base.DT_MYR,
        "counts": counts,
        "candidates": summary_df.to_dict(orient="records"),
        "outputs": {
            "draws_csv": str(draws_path),
            "deterministic_csv": str(det_path),
            "summary_csv": str(summary_path),
            "summary_json": str(json_path),
            "summary_md": str(md_path),
        },
    }
    json_path.write_text(json.dumps(payload, indent=2))

    lines = [
        "# Expanded Sgr A* approacher refinement",
        "",
        f"- MC realisations per candidate: {args.n_samp}",
        f"- Candidate seeds: top {args.n_per_potential} by point-estimate minimum spherical distance in each of static and barred integrations.",
        f"- Confirmed P(<10 pc)>0.5: {counts['n_P_lt_10pc_gt_0p5']}",
        f"- Soft P(<100 pc)>0.5: {counts['n_P_lt_100pc_gt_0p5']}",
        "",
        "| source_id | pot | tier | point min pc | interp p16 | interp p50 | interp p84 | P10 | P100 |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in summary_df.itertuples(index=False):
        lines.append(
            f"| {int(r.source_id)} | {r.potential} | {r.tier} | "
            f"{r.point_min_rsph_pc:.2f} | {r.interp_min_rsph_pc_p16:.1f} | "
            f"{r.interp_min_rsph_pc_p50:.1f} | {r.interp_min_rsph_pc_p84:.1f} | "
            f"{r.interp_P_lt_10pc:.3f} | {r.interp_P_lt_100pc:.3f} |"
        )
    md_path.write_text("\n".join(lines))

    log("wrote expanded Sgr A* refinement outputs")
    print(json.dumps({"counts": counts, "summary_csv": str(summary_path)}, indent=2))
    log(f"DONE in {time.time() - T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
