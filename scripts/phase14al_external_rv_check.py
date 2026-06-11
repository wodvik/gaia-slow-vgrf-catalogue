"""Phase 14AL -- external APOGEE/GALAH radial-velocity check.

Queries public VizieR TAP tables for line-of-sight velocities of the
Tier A+B+C catalogue and compares them to Gaia DR3 RVS velocities. The
coverage is selective, but it provides a direct referee-facing check for
the subset with external spectroscopy.
"""
from __future__ import annotations

import json
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from astroquery.utils.tap.core import TapPlus

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
MASTER = BUNDLE / "catalogues" / "catalogue_expanded_master.fits"
OUT = BUNDLE / "phase14" / "external_rv"
TAB_DIRS = [REPO / "release" / "tables" / "v15", BUNDLE / "tables" / "v15"]

VIZIER_TAP_URL = "http://tapvizier.cds.unistra.fr/TAPVizieR/tap"
BATCH_SIZE = 500
PAUSE_SECONDS = 0.5

SURVEYS = {
    "APOGEE": {
        "table": '"III/286/catalog"',
        "id_col": "GaiaEDR3",
        "select_cols": '"GaiaEDR3", "APOGEE", "HRV", "e_HRV", "s_HRV", "Nvis", "FlRV", "SNR"',
        "rv_col": "HRV",
        "rv_err_col": "e_HRV",
    },
    "GALAH": {
        "table": '"J/MNRAS/506/150/stars"',
        "id_col": "GaiaEDR3",
        "select_cols": (
            '"GaiaEDR3", "GALAH", "RVgalah", "e_RVgalah", "redflag", '
            '"snrc1iraf", "snrc2iraf", "snrc3iraf", "snrc4iraf"'
        ),
        "rv_col": "RVgalah",
        "rv_err_col": "e_RVgalah",
    },
}


def clean_tier(s: pd.Series) -> pd.Series:
    if s.dtype == object and len(s) and isinstance(s.iloc[0], bytes):
        return s.str.decode("utf-8").str.strip()
    return s.astype(str).str.strip()


def load_catalogue() -> pd.DataFrame:
    m = Table.read(MASTER).to_pandas()
    m["tier"] = clean_tier(m["tier"])
    cols = ["source_id", "tier", "P_vgrf_below_25", "vgrf_default",
            "radial_velocity", "radial_velocity_error"]
    return m[cols].copy()


def batch_path(survey: str, idx: int) -> Path:
    return OUT / survey.lower() / f"batch_{idx:04d}.csv"


def query_survey(tap: TapPlus, survey: str, ids: np.ndarray) -> pd.DataFrame:
    cfg = SURVEYS[survey]
    sdir = OUT / survey.lower()
    sdir.mkdir(parents=True, exist_ok=True)
    n_batches = int(np.ceil(len(ids) / BATCH_SIZE))
    failures: list[dict[str, str | int]] = []
    for batch_idx in range(n_batches):
        path = batch_path(survey, batch_idx)
        if path.exists():
            continue
        chunk = ids[batch_idx * BATCH_SIZE:min((batch_idx + 1) * BATCH_SIZE, len(ids))]
        id_list = ",".join(str(int(x)) for x in chunk)
        query = f"""
        SELECT {cfg["select_cols"]}
        FROM {cfg["table"]}
        WHERE "{cfg["id_col"]}" IN ({id_list})
        """
        print(f"[14AL] {survey} batch {batch_idx + 1}/{n_batches} ({len(chunk)} IDs)", flush=True)
        try:
            result = tap.launch_job(query).get_results().to_pandas()
            if len(result) and cfg["id_col"] in result.columns:
                result = result.rename(columns={cfg["id_col"]: "source_id"})
            result.to_csv(path, index=False)
            print(f"[14AL]   wrote {len(result)} rows", flush=True)
        except Exception as exc:  # noqa: BLE001
            failures.append({"batch": batch_idx, "error": str(exc)})
            print(f"[14AL]   FAILED: {exc}", flush=True)
        time.sleep(PAUSE_SECONDS)
    (OUT / f"{survey.lower()}_failures.json").write_text(json.dumps(failures, indent=2))

    frames = []
    for path in sorted(sdir.glob("batch_*.csv")):
        try:
            df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            continue
        if len(df):
            frames.append(df)
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if "source_id" in df.columns:
        df["source_id"] = pd.to_numeric(df["source_id"], errors="coerce").astype("Int64")
        df = df.dropna(subset=["source_id"]).drop_duplicates("source_id", keep="first")
        df["source_id"] = df["source_id"].astype(np.int64)
    df.to_csv(OUT / f"{survey.lower()}_rv_crossmatch.csv", index=False)
    return df


def summarize(df: pd.DataFrame, survey: str) -> dict:
    d = df[np.isfinite(df["delta_rv_kms"])]
    ab = d[d["tier"].isin(["A", "B"])]
    if len(d) == 0:
        return {"survey": survey, "n_tierABC": 0, "n_tierAB": 0}
    med = float(np.median(d["delta_rv_kms"]))
    scatter = float(1.4826 * np.median(np.abs(d["delta_rv_kms"] - med)))
    return {
        "survey": survey,
        "n_tierABC": int(len(d)),
        "n_tierAB": int(len(ab)),
        "median_delta_rv_kms": round(med, 3),
        "robust_scatter_kms": round(scatter, 3),
        "median_abs_delta_rv_kms": round(float(np.median(np.abs(d["delta_rv_kms"]))), 3),
        "frac_abs_delta_gt5_pct": round(float(np.mean(np.abs(d["delta_rv_kms"]) > 5.0) * 100), 1),
        "frac_abs_delta_gt10_pct": round(float(np.mean(np.abs(d["delta_rv_kms"]) > 10.0) * 100), 1),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cat = load_catalogue()
    abc = cat[cat["tier"].isin(["A", "B", "C"])].copy()
    ids = np.asarray(sorted(abc["source_id"].astype(np.int64).unique()))
    print(f"[14AL] querying external RVs for {len(ids)} Tier A+B+C source_ids")

    tap = TapPlus(url=VIZIER_TAP_URL)
    rows = []
    summaries = []
    for survey in SURVEYS:
        cfg = SURVEYS[survey]
        x = query_survey(tap, survey, ids)
        if len(x) == 0 or cfg["rv_col"] not in x.columns:
            summaries.append({"survey": survey, "n_tierABC": 0, "n_tierAB": 0,
                              "note": "no external RV rows returned"})
            continue
        x = x.merge(abc, on="source_id", how="inner")
        x["external_survey"] = survey
        x["external_rv_kms"] = pd.to_numeric(x[cfg["rv_col"]], errors="coerce")
        x["external_rv_error_kms"] = pd.to_numeric(x[cfg["rv_err_col"]], errors="coerce")
        x["delta_rv_kms"] = x["radial_velocity"].astype(float) - x["external_rv_kms"]
        rows.append(x)
        summaries.append(summarize(x, survey))

    joined = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(joined):
        keep_cols = [
            "external_survey", "source_id", "tier", "P_vgrf_below_25", "vgrf_default",
            "radial_velocity", "radial_velocity_error", "external_rv_kms",
            "external_rv_error_kms", "delta_rv_kms",
        ]
        extra_cols = [c for c in ["APOGEE", "GALAH", "Nvis", "FlRV", "SNR", "redflag"] if c in joined.columns]
        joined[keep_cols + extra_cols].to_csv(OUT / "external_rv_comparison.csv", index=False)
    else:
        (OUT / "external_rv_comparison.csv").write_text("")

    summary = {
        "timestamp": datetime.now().isoformat(),
        "n_tierABC": int(len(abc)),
        "tap_url": VIZIER_TAP_URL,
        "coverage_note": "External RV coverage is sparse and footprint-selected; it is used as a consistency audit, not as a replacement tiering input.",
        "surveys": summaries,
    }
    (OUT / "external_rv_summary.json").write_text(json.dumps(summary, indent=2))
    write_table(summaries)
    print(json.dumps(summary, indent=2))
    return 0


def write_table(summaries: list[dict]) -> None:
    lines = [
        r"% Auto-generated by scripts/phase14al_external_rv_check.py",
        r"\begin{table}[!htbp]",
        r"\centering",
        r"\begin{minipage}{0.90\columnwidth}",
        r"\caption{External radial-velocity consistency audit for the subset of"
        r" Tier~A+B+C stars with APOGEE DR17 or GALAH DR3 line-of-sight velocities."
        r" $\Delta RV$ is Gaia DR3 RVS minus the external survey velocity; the"
        r" $>5$ column gives the fraction with $|\Delta RV|>5$\kms. The coverage"
        r" is footprint-selected and is used only as a robustness check."
        r"\label{tab:external_rv}}",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{2pt}",
        r"\begin{tabular}{@{}lrrrrr@{}}",
        r"\hline\hline",
        r"survey & $N_{\rm ABC}$ & $N_{\rm AB}$ & med. $\Delta RV$ & "
        r"scatter & $>5$ \\",
        r"\hline",
    ]
    for s in summaries:
        if s.get("n_tierABC", 0) == 0:
            lines.append(f"{s['survey']} & 0 & 0 & -- & -- & -- " + r"\\")
        else:
            lines.append(
                f"{s['survey']} & {s['n_tierABC']} & {s['n_tierAB']} & "
                f"{s['median_delta_rv_kms']:.2f} & {s['robust_scatter_kms']:.2f} & "
                f"{s['frac_abs_delta_gt5_pct']:.1f}\\% \\\\"
            )
    lines += [
        r"\hline",
        r"\end{tabular}",
        r"\end{minipage}",
        r"\end{table}",
    ]
    text = "\n".join(lines) + "\n"
    for d in TAB_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        (d / "tab_external_rv.tex").write_text(text)


if __name__ == "__main__":
    raise SystemExit(main())
