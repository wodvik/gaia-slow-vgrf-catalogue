"""full expanded APOGEE/GALAH spectroscopic cross-match.

Queries VizieR TAP for the full 20,829-source expanded candidate pool
using exact Gaia EDR3/DR3 source identifiers. Results are checkpointed
batch-by-batch so a failed network query can be resumed without losing
completed batches.
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

REPO = Path(__file__).resolve().parents[1]
MASTER = REPO / "catalogues/catalogue_expanded_master.fits"
OUT = REPO / "analysis_products/expanded_spectroscopic_crossmatch"
OUT.mkdir(parents=True, exist_ok=True)

# Use HTTP because the local Windows Python certificate store does not
# validate the CDS chain reliably on this machine. The query payload is
# public catalogue IDs only, and the returned VizieR tables are public.
VIZIER_TAP_URL = "http://tapvizier.cds.unistra.fr/TAPVizieR/tap"
BATCH_SIZE = 500
PAUSE_SECONDS = 1.0


SURVEYS = {
    "apogee_dr17": {
        "table": '"III/286/catalog"',
        "id_col": "GaiaEDR3",
        "select_cols": """
            "GaiaEDR3", "Teff", "e_Teff", "logg", "e_logg",
            "[Fe/H]", "e_[Fe/H]", "[M/H]", "e_[M/H]",
            "[a/M]", "e_[a/M]",
            "[C/Fe]", "e_[C/Fe]", "[N/Fe]", "e_[N/Fe]",
            "[O/Fe]", "e_[O/Fe]", "[Mg/Fe]", "e_[Mg/Fe]",
            "[Al/Fe]", "e_[Al/Fe]", "[Si/Fe]", "e_[Si/Fe]",
            "[Ca/Fe]", "e_[Ca/Fe]", "[Ti/Fe]", "e_[Ti/Fe]",
            "[Mn/Fe]", "e_[Mn/Fe]", "[Ni/Fe]", "e_[Ni/Fe]",
            "[Ce/Fe]", "e_[Ce/Fe]"
        """,
    },
    "galah_dr3": {
        "table": '"J/MNRAS/506/150/stars"',
        "id_col": "GaiaEDR3",
        "select_cols": """
            "GaiaEDR3", "Teff", "e_Teff", "logg", "e_logg",
            "[Fe/H]", "e_[Fe/H]", "[alpha/Fe]", "e_[alpha/Fe]",
            "[C/Fe]", "e_[C/Fe]", "[O/Fe]", "e_[O/Fe]",
            "[Mg/Fe]", "e_[Mg/Fe]", "[Al/Fe]", "e_[Al/Fe]",
            "[Si/Fe]", "e_[Si/Fe]", "[Ca/Fe]", "e_[Ca/Fe]",
            "[Ti/Fe]", "e_[Ti/Fe]", "[Mn/Fe]", "e_[Mn/Fe]",
            "[Ni/Fe]", "e_[Ni/Fe]", "[Ba/Fe]", "e_[Ba/Fe]",
            "[Eu/Fe]", "e_[Eu/Fe]"
        """,
    },
}


def load_ids() -> np.ndarray:
    tab = Table.read(MASTER)
    ids = np.asarray(tab["source_id"], dtype=np.int64)
    return np.unique(ids)


def batch_path(survey: str, idx: int) -> Path:
    return OUT / survey / f"batch_{idx:04d}.csv"


def query_survey(tap: TapPlus, survey: str, source_ids: np.ndarray) -> pd.DataFrame:
    cfg = SURVEYS[survey]
    survey_dir = OUT / survey
    survey_dir.mkdir(exist_ok=True)
    n_batches = int(np.ceil(len(source_ids) / BATCH_SIZE))
    failures: list[dict[str, str | int]] = []

    for batch_idx in range(n_batches):
        path = batch_path(survey, batch_idx)
        if path.exists():
            continue
        start = batch_idx * BATCH_SIZE
        end = min(start + BATCH_SIZE, len(source_ids))
        ids = source_ids[start:end]
        id_list = ",".join(str(int(x)) for x in ids)
        query = f"""
        SELECT {cfg["select_cols"]}
        FROM {cfg["table"]}
        WHERE "{cfg["id_col"]}" IN ({id_list})
        """
        print(f"[14W] {survey} batch {batch_idx + 1}/{n_batches} ({len(ids)} IDs)", flush=True)
        try:
            result = tap.launch_job(query).get_results()
            df = result.to_pandas()
            if len(df) and cfg["id_col"] in df.columns:
                df = df.rename(columns={cfg["id_col"]: "source_id"})
            df.to_csv(path, index=False)
            print(f"[14W]   wrote {len(df)} rows", flush=True)
        except Exception as exc:
            failures.append({"batch": batch_idx, "error": str(exc)})
            print(f"[14W]   FAILED: {exc}", flush=True)
        time.sleep(PAUSE_SECONDS)

    frames = []
    for path in sorted(survey_dir.glob("batch_*.csv")):
        try:
            df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            continue
        if len(df):
            frames.append(df)
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if "source_id" in combined.columns:
        combined["source_id"] = pd.to_numeric(combined["source_id"], errors="coerce").astype("Int64")
        combined = (combined.dropna(subset=["source_id"])
                            .drop_duplicates("source_id", keep="first"))
        combined["source_id"] = combined["source_id"].astype(np.int64)
    combined.to_csv(OUT / f"{survey}_expanded_crossmatch.csv", index=False)
    (OUT / f"{survey}_failures.json").write_text(json.dumps(failures, indent=2))
    return combined


def main() -> int:
    source_ids = load_ids()
    print(f"[14W] loaded {len(source_ids)} expanded source_ids")
    tap = TapPlus(url=VIZIER_TAP_URL)
    summary = {
        "timestamp": datetime.now().isoformat(),
        "n_input": int(len(source_ids)),
        "batch_size": BATCH_SIZE,
        "surveys": {},
    }
    all_ids: set[int] = set()
    for survey in SURVEYS:
        df = query_survey(tap, survey, source_ids)
        ids = set(df["source_id"].astype(int)) if len(df) and "source_id" in df else set()
        all_ids |= ids
        summary["surveys"][survey] = {
            "n_matches": int(len(ids)),
            "match_rate": float(len(ids) / len(source_ids)),
            "table": SURVEYS[survey]["table"].strip('"'),
            "id_column": SURVEYS[survey]["id_col"],
        }
    summary["n_unique_any_apogee_galah"] = int(len(all_ids))
    summary["unique_match_rate"] = float(len(all_ids) / len(source_ids))
    (OUT / "expanded_spectroscopic_crossmatch_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
