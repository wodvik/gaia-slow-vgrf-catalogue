"""Phase 15A -- rectilinear no-force look-back/look-forward baseline.

This deliberately ignores Galactic gravity.  It propagates the present-day
Galactocentric phase-space vector as

    r(t) = r0 + v0 t

and reports the closest straight-line approach to the Galactic Centre.  The
diagnostic is a sanity baseline for later reverse modelling, not a physical
orbit model.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

import phase14y_radial_phase_mc as phase14y


BUNDLE = Path(__file__).resolve().parents[1]
ORBIT_CATALOGUE = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.fits"
OUT = BUNDLE / "phase14" / "rectilinear_baseline"
SUMMARY_JSON = OUT / "rectilinear_baseline_summary.json"
PER_STAR_CSV = OUT / "rectilinear_baseline_tierAB_per_star.csv"
TABLE_TEX = BUNDLE / "tables" / "v15" / "tab_rectilinear_baseline.tex"

KPC_PER_GYR_PER_KMS = 1.0227121650537077
DEFAULT_N_SAMP = 5000
DEFAULT_CHUNK_SIZE = 64
NEAR_TIME_MYR = 10.0
DIST_THRESHOLDS_PC = (100, 500, 1000)
SEED = phase14y.SEED + 15_001
T0 = time.time()


def log(message: str) -> None:
    print(f"[15A-rectilinear t={time.time() - T0:7.1f}s] {message}", flush=True)


def clean_tier(series: pd.Series) -> pd.Series:
    def one(value) -> str:
        if isinstance(value, (bytes, bytearray)):
            return value.decode("utf-8").strip()
        return str(value).strip()

    return series.map(one)


def rectilinear_metrics(ic: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = ic[..., :3]
    v_kms = ic[..., 3:6]
    v = v_kms * KPC_PER_GYR_PER_KMS
    dot = np.sum(r * v, axis=-1)
    v2 = np.sum(v * v, axis=-1)
    t_gyr = -dot / v2
    closest = r + v * t_gyr[..., None]
    d_pc = np.sqrt(np.sum(closest * closest, axis=-1)) * 1000.0
    rnorm = np.sqrt(np.sum(r * r, axis=-1))
    vr_sph = np.sum(r * v_kms, axis=-1) / rnorm
    return t_gyr * 1000.0, d_pc, vr_sph


def load_point_estimate() -> pd.DataFrame:
    df = Table.read(ORBIT_CATALOGUE).to_pandas()
    df["tier"] = clean_tier(df["tier"])
    ic = df[["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]].to_numpy(float)
    t_myr, d_pc, vr_sph = rectilinear_metrics(ic)
    df["rect_t_min_myr"] = t_myr
    df["rect_d_min_pc"] = d_pc
    df["v_radial_sph_kms"] = vr_sph
    return df


def sample_masks(df: pd.DataFrame, gold_ids: set[int]) -> dict[str, np.ndarray]:
    tier = clean_tier(df["tier"])
    source_id = df["source_id"].astype("int64")
    return {
        "tier_ab": tier.isin(["A", "B"]).to_numpy(),
        "gold_tier_ab": source_id.isin(gold_ids).to_numpy(),
        "tier_abc": tier.isin(["A", "B", "C"]).to_numpy(),
    }


def summarise_point(df: pd.DataFrame, mask: np.ndarray) -> dict[str, float | int]:
    sub = df.loc[mask].copy()
    t = sub["rect_t_min_myr"].to_numpy(float)
    d = sub["rect_d_min_pc"].to_numpy(float)
    out: dict[str, float | int] = {
        "n": int(mask.sum()),
        "d_min_pc_p16": float(np.percentile(d, 16)),
        "d_min_pc_p50": float(np.percentile(d, 50)),
        "d_min_pc_p84": float(np.percentile(d, 84)),
        "abs_t_min_myr_p16": float(np.percentile(np.abs(t), 16)),
        "abs_t_min_myr_p50": float(np.percentile(np.abs(t), 50)),
        "abs_t_min_myr_p84": float(np.percentile(np.abs(t), 84)),
        "past_count": int(np.sum(t < -NEAR_TIME_MYR)),
        "near_now_count": int(np.sum(np.abs(t) <= NEAR_TIME_MYR)),
        "future_count": int(np.sum(t > NEAR_TIME_MYR)),
    }
    for threshold in DIST_THRESHOLDS_PC:
        near = d < threshold
        out[f"d_lt_{threshold}_pc"] = int(np.sum(near))
        out[f"past_d_lt_{threshold}_pc"] = int(np.sum(near & (t < -NEAR_TIME_MYR)))
        out[f"future_d_lt_{threshold}_pc"] = int(np.sum(near & (t > NEAR_TIME_MYR)))
    return out


def empty_count_arrays(n_samp: int) -> dict[str, np.ndarray]:
    out = {
        "past": np.zeros(n_samp, dtype=np.int32),
        "near_now": np.zeros(n_samp, dtype=np.int32),
        "future": np.zeros(n_samp, dtype=np.int32),
    }
    for threshold in DIST_THRESHOLDS_PC:
        out[f"d_lt_{threshold}_pc"] = np.zeros(n_samp, dtype=np.int32)
        out[f"past_d_lt_{threshold}_pc"] = np.zeros(n_samp, dtype=np.int32)
        out[f"future_d_lt_{threshold}_pc"] = np.zeros(n_samp, dtype=np.int32)
    return out


def count_summary(values: np.ndarray) -> list[float]:
    return [float(x) for x in np.percentile(values, [16, 50, 84])]


def summarise_mc_counts(agg: dict[str, np.ndarray]) -> dict[str, list[float]]:
    return {key + "_p16_p50_p84": count_summary(value) for key, value in agg.items()}


def per_star_rows(chunk: pd.DataFrame, t_myr: np.ndarray, d_pc: np.ndarray, vr_sph: np.ndarray) -> list[dict]:
    rows: list[dict] = []
    for i, row in chunk.reset_index(drop=True).iterrows():
        item = {
            "source_id": int(row["source_id"]),
            "tier": str(row["tier"]),
            "vgrf_default_exact": float(row["vgrf_default_exact"]),
            "rect_t_min_myr_p16": float(np.percentile(t_myr[i], 16)),
            "rect_t_min_myr_p50": float(np.percentile(t_myr[i], 50)),
            "rect_t_min_myr_p84": float(np.percentile(t_myr[i], 84)),
            "rect_abs_t_min_myr_p50": float(np.percentile(np.abs(t_myr[i]), 50)),
            "rect_d_min_pc_p16": float(np.percentile(d_pc[i], 16)),
            "rect_d_min_pc_p50": float(np.percentile(d_pc[i], 50)),
            "rect_d_min_pc_p84": float(np.percentile(d_pc[i], 84)),
            "v_radial_sph_kms_p50": float(np.percentile(vr_sph[i], 50)),
            "p_rect_past": float(np.mean(t_myr[i] < -NEAR_TIME_MYR)),
            "p_rect_near_now": float(np.mean(np.abs(t_myr[i]) <= NEAR_TIME_MYR)),
            "p_rect_future": float(np.mean(t_myr[i] > NEAR_TIME_MYR)),
        }
        for threshold in DIST_THRESHOLDS_PC:
            item[f"p_rect_d_lt_{threshold}_pc"] = float(np.mean(d_pc[i] < threshold))
        rows.append(item)
    return rows


def run_mc(df: pd.DataFrame, gold_ids: set[int], n_samp: int, chunk_size: int) -> tuple[dict, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    masks = sample_masks(df, gold_ids)
    agg = {key: empty_count_arrays(n_samp) for key in masks}
    tier_ab_rows: list[dict] = []

    for start in range(0, len(df), chunk_size):
        stop = min(start + chunk_size, len(df))
        chunk = df.iloc[start:stop].copy()
        log(f"sampling rectilinear baseline rows {start}-{stop} of {len(df)}")
        ic = phase14y.sample_current_phase_space(chunk, n_samp, rng)
        t_myr, d_pc, vr_sph = rectilinear_metrics(ic)

        past = t_myr < -NEAR_TIME_MYR
        near_now = np.abs(t_myr) <= NEAR_TIME_MYR
        future = t_myr > NEAR_TIME_MYR
        for key, mask in masks.items():
            local = mask[start:stop]
            if not np.any(local):
                continue
            agg[key]["past"] += past[local].sum(axis=0, dtype=np.int32)
            agg[key]["near_now"] += near_now[local].sum(axis=0, dtype=np.int32)
            agg[key]["future"] += future[local].sum(axis=0, dtype=np.int32)
            for threshold in DIST_THRESHOLDS_PC:
                close = d_pc < threshold
                agg[key][f"d_lt_{threshold}_pc"] += close[local].sum(axis=0, dtype=np.int32)
                agg[key][f"past_d_lt_{threshold}_pc"] += (close[local] & past[local]).sum(axis=0, dtype=np.int32)
                agg[key][f"future_d_lt_{threshold}_pc"] += (close[local] & future[local]).sum(axis=0, dtype=np.int32)

        local_ab = masks["tier_ab"][start:stop]
        if np.any(local_ab):
            tier_ab_rows.extend(per_star_rows(chunk.loc[local_ab], t_myr[local_ab], d_pc[local_ab], vr_sph[local_ab]))

    mc = {key: summarise_mc_counts(value) for key, value in agg.items()}
    per_star = pd.DataFrame(tier_ab_rows).sort_values("rect_d_min_pc_p50").reset_index(drop=True)
    return mc, per_star


def fmt_int(value: int | float) -> str:
    return f"{int(round(float(value))):,}".replace(",", r"{,}")


def fmt_triplet(values: list[float], decimals: int = 0) -> str:
    lo, med, hi = [float(x) for x in values]
    if decimals == 0:
        return rf"{fmt_int(med)} [{fmt_int(lo)}, {fmt_int(hi)}]"
    return rf"{med:.{decimals}f} [{lo:.{decimals}f}, {hi:.{decimals}f}]"


def write_table(summary: dict) -> None:
    rows = [
        ("tier_ab", r"Tier~A+B"),
        ("gold_tier_ab", r"Gold A+B"),
        ("tier_abc", r"Tier~A+B+C"),
    ]
    lines = [
        r"\begin{center}",
        r"\refstepcounter{table}\label{tab:rectilinear_baseline}",
        r"\begin{minipage}{\columnwidth}",
        r"\small\textbf{Table \thetable.} Rectilinear no-force baseline for present-day Galactocentric phase space.",
        r"\end{minipage}",
        r"\vspace{0.4ex}",
        "",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{2.5pt}",
        r"\begin{tabular}{lrrrrr}",
        r"\hline\hline",
        r"Sample & $N$ & $\tilde d_{\rm lin}$ & $\widetilde{|t_{\rm lin}|}$ & MC $N(d<0.5{\rm\,kpc})$ & MC $N(d<1{\rm\,kpc})$ \\",
        r"\hline",
    ]
    for key, label in rows:
        point = summary["point_estimate"][key]
        mc = summary["mc_counts"][key]
        lines.append(
            f"{label} & {fmt_int(point['n'])} & "
            f"{point['d_min_pc_p50'] / 1000.0:.2f} & "
            f"{point['abs_t_min_myr_p50']:.0f} & "
            f"{fmt_triplet(mc['d_lt_500_pc_p16_p50_p84'])} & "
            f"{fmt_triplet(mc['d_lt_1000_pc_p16_p50_p84'])} \\\\"
        )
    lines.extend(
        [
            r"\hline",
            r"\end{tabular}",
            r"\vspace{0.5ex}",
            "",
            r"\begin{minipage}{0.95\columnwidth}",
            r"\footnotesize",
            rf"\textit{{Note.}} The baseline uses $\mathbf{{r}}(t)=\mathbf{{r}}_0+\mathbf{{v}}_0t$ in the adopted Galactocentric frame, with no Galactic potential.  Positive $t_{{\rm lin}}$ is future and negative $t_{{\rm lin}}$ is past; $d_{{\rm lin}}$ is the minimum straight-line Galactocentric distance.  Point-estimate medians of $d_{{\rm lin}}$ are in kpc and $|t_{{\rm lin}}|$ in Myr.  MC count intervals are 16th--84th percentiles from {fmt_int(summary['n_samp'])} present-day Gaia phase-space realisations per star.",
            r"\end{minipage}",
            r"\end{center}",
        ]
    )
    TABLE_TEX.write_text("\n".join(lines) + "\n")
    log(f"wrote {TABLE_TEX}")


def run(n_samp: int, chunk_size: int) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    gold_ids = phase14y.gold_source_ids()

    point_df = load_point_estimate()
    point_masks = sample_masks(point_df, gold_ids)
    point = {key: summarise_point(point_df, mask) for key, mask in point_masks.items()}

    mc_input = phase14y.load_expanded()
    mc_counts, per_star = run_mc(mc_input, gold_ids, n_samp, chunk_size)
    per_star.to_csv(PER_STAR_CSV, index=False)
    log(f"wrote {PER_STAR_CSV}")

    top = point_df.loc[point_masks["tier_ab"]].sort_values("rect_d_min_pc").head(12)
    top_rows = [
        {
            "source_id": int(row["source_id"]),
            "tier": str(row["tier"]),
            "rect_d_min_pc": float(row["rect_d_min_pc"]),
            "rect_t_min_myr": float(row["rect_t_min_myr"]),
            "vgrf_default_exact": float(row["vgrf_default_exact"]),
        }
        for _, row in top.iterrows()
    ]

    summary = {
        "phase": "15A",
        "orbit_catalogue": str(ORBIT_CATALOGUE),
        "input_csv": str(phase14y.INPUT),
        "n_samp": int(n_samp),
        "chunk_size": int(chunk_size),
        "seed": int(SEED),
        "near_time_myr": NEAR_TIME_MYR,
        "kpc_per_gyr_per_kms": KPC_PER_GYR_PER_KMS,
        "point_estimate": point,
        "mc_counts": mc_counts,
        "top_tier_ab_point_estimate": top_rows,
        "outputs": {
            "summary_json": str(SUMMARY_JSON),
            "per_star_csv": str(PER_STAR_CSV),
            "table_tex": str(TABLE_TEX),
        },
        "note": "Rectilinear no-force baseline; not an orbit model.",
    }
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
    for key, label in [("tier_ab", "Tier A+B"), ("gold_tier_ab", "Gold A+B"), ("tier_abc", "Tier A+B+C")]:
        point = summary["point_estimate"][key]
        mc = summary["mc_counts"][key]
        log(
            f"{label}: point median d={point['d_min_pc_p50'] / 1000.0:.2f} kpc, "
            f"median |t|={point['abs_t_min_myr_p50']:.0f} Myr, "
            f"point d<0.5/1.0 kpc={point['d_lt_500_pc']}/{point['d_lt_1000_pc']}, "
            f"MC d<0.5={mc['d_lt_500_pc_p16_p50_p84']}, "
            f"MC d<1.0={mc['d_lt_1000_pc_p16_p50_p84']}"
        )


if __name__ == "__main__":
    main()
