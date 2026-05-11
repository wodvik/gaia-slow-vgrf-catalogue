"""
Audit the Phase 0 Gaia parent-buffer scan.

This summarizes the broad legacy-Vgrf buffer and compares it with the
historical slow_stars_enriched_orbits.csv preselection. It is intentionally
lightweight and chunked so it can run on a large buffer CSV.

Outputs:
  release/v2/parent_scan/gate0_parent_buffer_summary.json
  release/v2/parent_scan/gate0_new_legacy_lt25_candidates.csv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
DEFAULT_SCAN_DIR = REPO / "release" / "v2" / "parent_scan"
DEFAULT_BUFFER = DEFAULT_SCAN_DIR / "gaia_parent_buffer_vgrf200_full.csv"
DEFAULT_OLD = REPO / "release" / "data" / "slow_stars_enriched_orbits.csv"
THRESHOLDS = [25, 50, 75, 100, 125, 150, 175, 200]


def load_old_ids(path: Path) -> set[int]:
    old = pd.read_csv(path, usecols=["source_id"])
    return set(old["source_id"].astype("int64").tolist())


def audit_buffer(buffer_csv: Path, old_csv: Path, out_dir: Path, chunksize: int) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    old_ids = load_old_ids(old_csv)

    seen_ids: set[int] = set()
    duplicate_ids: set[int] = set()
    old_recovered: set[int] = set()
    new_legacy_lt25_chunks: list[pd.DataFrame] = []

    counts = {f"legacy_lt_{t}": 0 for t in THRESHOLDS}
    n_rows = 0
    min_v = np.inf
    max_v = -np.inf

    usecols = [
        "source_id",
        "ra",
        "dec",
        "parallax",
        "parallax_over_error",
        "radial_velocity",
        "legacy_v_total_grf",
        "parent_scan_file",
    ]

    for chunk in pd.read_csv(buffer_csv, chunksize=chunksize, usecols=lambda c: c in usecols):
        chunk["source_id"] = chunk["source_id"].astype("int64")
        ids = set(chunk["source_id"].tolist())
        duplicate_ids.update(seen_ids.intersection(ids))
        seen_ids.update(ids)
        old_recovered.update(ids.intersection(old_ids))

        v = chunk["legacy_v_total_grf"].to_numpy(dtype=float)
        n_rows += len(chunk)
        if len(v):
            min_v = min(min_v, float(np.nanmin(v)))
            max_v = max(max_v, float(np.nanmax(v)))
        for threshold in THRESHOLDS:
            counts[f"legacy_lt_{threshold}"] += int(np.sum(v < threshold))

        new_lt25 = chunk[(chunk["legacy_v_total_grf"] < 25.0) & (~chunk["source_id"].isin(old_ids))]
        if not new_lt25.empty:
            new_legacy_lt25_chunks.append(new_lt25.copy())

    missing_old = sorted(old_ids - old_recovered)
    new_lt25_path = out_dir / "gate0_new_legacy_lt25_candidates.csv"
    if new_legacy_lt25_chunks:
        pd.concat(new_legacy_lt25_chunks, ignore_index=True).to_csv(new_lt25_path, index=False)
        n_new_lt25 = int(sum(len(c) for c in new_legacy_lt25_chunks))
    else:
        pd.DataFrame(columns=usecols).to_csv(new_lt25_path, index=False)
        n_new_lt25 = 0

    summary = {
        "buffer_csv": str(buffer_csv),
        "old_preselection_csv": str(old_csv),
        "n_buffer_rows": int(n_rows),
        "n_buffer_unique_source_id": int(len(seen_ids)),
        "n_buffer_duplicate_source_id": int(len(duplicate_ids)),
        "n_old_preselection": int(len(old_ids)),
        "n_old_recovered_in_buffer": int(len(old_recovered)),
        "n_old_missing_from_buffer": int(len(missing_old)),
        "old_missing_source_ids": [int(x) for x in missing_old[:100]],
        "n_new_legacy_lt25_not_in_old": n_new_lt25,
        "new_legacy_lt25_csv": str(new_lt25_path),
        "legacy_v_total_grf_min": None if not np.isfinite(min_v) else float(min_v),
        "legacy_v_total_grf_max": None if not np.isfinite(max_v) else float(max_v),
        "counts": {k: int(v) for k, v in counts.items()},
    }
    out_path = out_dir / "gate0_parent_buffer_summary.json"
    out_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--buffer-csv", default=str(DEFAULT_BUFFER))
    parser.add_argument("--old-csv", default=str(DEFAULT_OLD))
    parser.add_argument("--out-dir", default=str(DEFAULT_SCAN_DIR))
    parser.add_argument("--chunksize", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = audit_buffer(
        Path(args.buffer_csv),
        Path(args.old_csv),
        Path(args.out_dir),
        args.chunksize,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
