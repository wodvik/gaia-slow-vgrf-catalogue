"""Smoke/regression checks for the public review bundle.

Run from the release bundle with:

    python tests/smoke_regression.py --bundle-root .
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.table import Table


EXPECTED_COUNTS = {
    "master": 20829,
    "tier_A": 289,
    "tier_AB": 541,
    "tier_ABC": 1952,
    "point_vgrf_lt25": 2755,
}


def fail(message: str) -> None:
    raise AssertionError(message)


def read_table(path: Path) -> Table:
    if not path.exists():
        fail(f"missing required product: {path}")
    return Table.read(path)


def assert_unique_source_ids(table: Table, label: str) -> None:
    values = np.asarray(table["source_id"], dtype=np.int64)
    if len(values) != len(np.unique(values)):
        fail(f"{label} has duplicate source_id values")


def maybe_assert_count(label: str, observed: int) -> None:
    expected = EXPECTED_COUNTS[label]
    if observed != expected:
        fail(f"{label} count mismatch: observed={observed}, expected={expected}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.bundle_root

    master = read_table(root / "catalogues" / "catalogue_expanded_master.fits")
    tier_a = read_table(root / "catalogues" / "catalogue_expanded_tierA.fits")
    tier_ab = read_table(root / "catalogues" / "catalogue_expanded_tierAB.fits")
    tier_abc = read_table(root / "catalogues" / "catalogue_expanded_tierABC.fits")
    orbits = read_table(root / "catalogues" / "catalogue_expanded_orbits_tierABC.fits")

    for label, table in {
        "master": master,
        "tier_A": tier_a,
        "tier_AB": tier_ab,
        "tier_ABC": tier_abc,
        "orbits": orbits,
    }.items():
        assert_unique_source_ids(table, label)

    maybe_assert_count("master", len(master))
    maybe_assert_count("tier_A", len(tier_a))
    maybe_assert_count("tier_AB", len(tier_ab))
    maybe_assert_count("tier_ABC", len(tier_abc))

    required_orbit_columns = {
        "source_id",
        "J_R",
        "J_z",
        "J_phi",
        "action_accuracy_sampled",
        "action_max_fracdiff",
        "action_reliability_flag",
    }
    missing = required_orbit_columns - set(orbits.colnames)
    if missing:
        fail(f"orbit catalogue missing columns: {sorted(missing)}")

    summary_path = root / "catalogues" / "expanded_catalogue_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        maybe_assert_count("point_vgrf_lt25", int(summary["point_estimate_vgrf_lt25"]))

    print("smoke_regression: passed count and structural checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
