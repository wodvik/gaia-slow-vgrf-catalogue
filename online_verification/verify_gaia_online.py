"""Verify local Gaia DR3 catalogue rows against the public Gaia Archive.

The script reads one of the expanded FITS catalogues, queries
``gaiadr3.gaia_source`` by ``source_id`` in batches, and writes a compact JSON
summary plus CSV tables for missing records and column discrepancies.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.table import Table


GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
DEFAULT_OUT_DIR = Path(__file__).resolve().parents[1] / "online_verification"
GAIA_COLUMNS = [
    "source_id",
    "ra",
    "dec",
    "parallax",
    "pmra",
    "pmdec",
    "radial_velocity",
    "phot_g_mean_mag",
]

TOLERANCES = {
    "ra": 1.0e-9,
    "dec": 1.0e-9,
    "parallax": 1.0e-9,
    "pmra": 1.0e-9,
    "pmdec": 1.0e-9,
    "radial_velocity": 1.0e-6,
    "phot_g_mean_mag": 1.0e-6,
}


@dataclass(frozen=True)
class Discrepancy:
    source_id: int
    column: str
    local_value: float | None
    online_value: float | None
    abs_delta: float | None
    tolerance: float


def finite_or_none(value: object) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(number):
        return None
    return number


def chunks(values: list[int], chunk_size: int) -> Iterable[list[int]]:
    for start in range(0, len(values), chunk_size):
        yield values[start : start + chunk_size]


def query_gaia_source_ids(source_ids: list[int], timeout: int, retries: int) -> dict[int, dict[str, str]]:
    ids = ",".join(str(value) for value in source_ids)
    query = (
        "SELECT source_id, ra, dec, parallax, pmra, pmdec, radial_velocity, phot_g_mean_mag "
        "FROM gaiadr3.gaia_source "
        f"WHERE source_id IN ({ids})"
    )
    data = {
        "REQUEST": "doQuery",
        "LANG": "ADQL",
        "FORMAT": "csv",
        "QUERY": query,
    }

    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = requests.post(GAIA_TAP_SYNC, data=data, timeout=timeout)
            response.raise_for_status()
            reader = csv.DictReader(io.StringIO(response.text))
            return {int(row["source_id"]): row for row in reader}
        except Exception as exc:  # noqa: BLE001 - report final network/server error
            last_error = exc
            if attempt < retries:
                time.sleep(min(2**attempt, 30))

    raise RuntimeError(f"Gaia Archive query failed after {retries} attempts: {last_error}")


def compare_rows(local: Table, online: dict[int, dict[str, str]]) -> tuple[list[int], list[Discrepancy]]:
    missing: list[int] = []
    discrepancies: list[Discrepancy] = []

    for row in local:
        source_id = int(row["source_id"])
        online_row = online.get(source_id)
        if online_row is None:
            missing.append(source_id)
            continue

        for column, tolerance in TOLERANCES.items():
            local_value = finite_or_none(row[column])
            online_value = finite_or_none(online_row.get(column))
            if local_value is None and online_value is None:
                continue
            if local_value is None or online_value is None:
                discrepancies.append(
                    Discrepancy(source_id, column, local_value, online_value, None, tolerance)
                )
                continue
            delta = abs(local_value - online_value)
            if delta > tolerance:
                discrepancies.append(
                    Discrepancy(source_id, column, local_value, online_value, delta, tolerance)
                )

    return missing, discrepancies


def write_missing(path: Path, missing: list[int]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source_id"])
        for source_id in missing:
            writer.writerow([source_id])


def write_discrepancies(path: Path, discrepancies: list[Discrepancy]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source_id", "column", "local_value", "online_value", "abs_delta", "tolerance"])
        for item in discrepancies:
            writer.writerow(
                [
                    item.source_id,
                    item.column,
                    item.local_value,
                    item.online_value,
                    item.abs_delta,
                    item.tolerance,
                ]
            )


def bundle_relative(path: Path, bundle_root: Path) -> str:
    try:
        return path.resolve().relative_to(bundle_root.resolve()).as_posix()
    except ValueError:
        return str(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("catalogue", type=Path)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--chunk-size", type=int, default=250)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--sleep", type=float, default=0.15)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    table = Table.read(args.catalogue)
    source_ids = [int(value) for value in table["source_id"]]
    if len(source_ids) != len(set(source_ids)):
        duplicates = len(source_ids) - len(set(source_ids))
        print(f"ERROR: local catalogue contains {duplicates} duplicate source_id values", file=sys.stderr)
        return 2

    online_rows: dict[int, dict[str, str]] = {}
    total_chunks = math.ceil(len(source_ids) / args.chunk_size)
    for index, batch in enumerate(chunks(source_ids, args.chunk_size), start=1):
        online_rows.update(query_gaia_source_ids(batch, timeout=args.timeout, retries=args.retries))
        print(f"queried chunk {index}/{total_chunks}: {len(online_rows)}/{len(source_ids)} rows recovered")
        if args.sleep:
            time.sleep(args.sleep)

    missing, discrepancies = compare_rows(table, online_rows)

    stem = args.catalogue.stem
    missing_path = args.out_dir / f"{stem}_missing_source_ids.csv"
    discrepancy_path = args.out_dir / f"{stem}_discrepancies.csv"
    summary_path = args.out_dir / f"{stem}_online_verification_summary.json"
    bundle_root = args.out_dir.parent

    write_missing(missing_path, missing)
    write_discrepancies(discrepancy_path, discrepancies)

    max_delta_by_column: dict[str, float] = {}
    compared_by_column: dict[str, int] = {column: 0 for column in TOLERANCES}
    for row in table:
        source_id = int(row["source_id"])
        online_row = online_rows.get(source_id)
        if online_row is None:
            continue
        for column in TOLERANCES:
            local_value = finite_or_none(row[column])
            online_value = finite_or_none(online_row.get(column))
            if local_value is None and online_value is None:
                compared_by_column[column] += 1
            elif local_value is not None and online_value is not None:
                compared_by_column[column] += 1
                max_delta_by_column[column] = max(
                    max_delta_by_column.get(column, 0.0),
                    abs(local_value - online_value),
                )

    summary = {
        "catalogue": bundle_relative(args.catalogue, bundle_root),
        "gaia_archive_table": "gaiadr3.gaia_source",
        "gaia_tap_sync_endpoint": GAIA_TAP_SYNC,
        "rows_local": len(table),
        "unique_source_ids_local": len(set(source_ids)),
        "rows_recovered_online": len(online_rows),
        "missing_source_ids": len(missing),
        "discrepancies": len(discrepancies),
        "tolerances": TOLERANCES,
        "compared_by_column": compared_by_column,
        "max_abs_delta_by_column": max_delta_by_column,
        "missing_csv": bundle_relative(missing_path, bundle_root),
        "discrepancies_csv": bundle_relative(discrepancy_path, bundle_root),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 1 if missing or discrepancies else 0


if __name__ == "__main__":
    raise SystemExit(main())
