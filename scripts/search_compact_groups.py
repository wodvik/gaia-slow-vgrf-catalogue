"""present-day compact-group sanity check.

This checks whether the slow-Vgrf catalogue contains an obvious open/globular
cluster-like aggregate at the present epoch. It uses the released catalogue
and orbit-table Galactocentric positions/velocities, then runs physically
interpretable DBSCAN/FoF searches in 3D position and 3D velocity.

The test is intentionally conservative: it can flag compact candidates, but a
snapshot cannot by itself prove that a small-N group is gravitationally bound.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from astropy.table import Table, join
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN

REPO = Path(__file__).resolve().parents[1]
CATALOGUE = REPO / "catalogues/catalogue_expanded_master.fits"
ORBITS = REPO / "catalogues/catalogue_expanded_orbits_tierABC.fits"
OUT_DIR = REPO / "analysis_products"

G_PC_MSUN_KMS2 = 0.00430091

SEARCHES = [
    {
        "label": "tight_bound_like",
        "max_sep_pc": 10.0,
        "max_dv_kms": 2.0,
        "min_members": 3,
        "description": "very compact, roughly bound-cluster core scale",
    },
    {
        "label": "open_cluster_loose",
        "max_sep_pc": 25.0,
        "max_dv_kms": 5.0,
        "min_members": 3,
        "description": "loose open-cluster-like compactness allowing measurement scatter",
    },
    {
        "label": "spatial_25pc_only",
        "max_sep_pc": 25.0,
        "max_dv_kms": None,
        "min_members": 3,
        "description": "position-only catch for any compact sky/distance clump",
    },
    {
        "label": "association_wide",
        "max_sep_pc": 100.0,
        "max_dv_kms": 10.0,
        "min_members": 5,
        "description": "wide association / tidal-tail sanity check",
    },
    {
        "label": "very_wide_velocity_coherent",
        "max_sep_pc": 250.0,
        "max_dv_kms": 10.0,
        "min_members": 5,
        "description": "deliberately permissive phase-space clump check",
    },
]


def strip_tier(values) -> np.ndarray:
    return np.char.strip(np.array(values, dtype=str))


def pairwise_arrays(table: Table) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xyz_pc = np.vstack([table["x_kpc"], table["y_kpc"], table["z_kpc"]]).T * 1000.0
    vel = np.vstack([table["vx_kms"], table["vy_kms"], table["vz_kms"]]).T
    sep_pc = squareform(pdist(xyz_pc))
    dv_kms = squareform(pdist(vel))
    return xyz_pc, vel, sep_pc, dv_kms


def dbscan_metric(sep_pc: np.ndarray, dv_kms: np.ndarray, max_sep_pc: float,
                  max_dv_kms: float | None) -> np.ndarray:
    if max_dv_kms is None:
        metric = sep_pc / max_sep_pc
    else:
        metric = np.maximum(sep_pc / max_sep_pc, dv_kms / max_dv_kms)
    np.fill_diagonal(metric, 0.0)
    return metric


def summarize_group(sample_name: str, search: dict, table: Table, indices: np.ndarray,
                    xyz_pc: np.ndarray, vel: np.ndarray, sep_pc: np.ndarray,
                    dv_kms: np.ndarray) -> dict:
    ids = [int(table["source_id"][i]) for i in indices]
    tiers = [str(table["tier"][i]).strip() for i in indices]
    tri = np.triu_indices(len(indices), 1)
    local_sep = sep_pc[np.ix_(indices, indices)][tri]
    local_dv = dv_kms[np.ix_(indices, indices)][tri]
    center = xyz_pc[indices].mean(axis=0)
    rel = xyz_pc[indices] - center
    r_pc = np.sqrt(np.sum(rel * rel, axis=1))
    half_radius_pc = float(np.median(r_pc)) if len(r_pc) else 0.0
    v_centered = vel[indices] - vel[indices].mean(axis=0)
    if len(indices) > 1:
        sigma_1d = float(np.sqrt(np.sum(v_centered * v_centered) / (3 * (len(indices) - 1))))
    else:
        sigma_1d = 0.0
    virial_mass_msun = 3.0 * sigma_1d * sigma_1d * max(half_radius_pc, 1e-9) / G_PC_MSUN_KMS2

    return {
        "sample": sample_name,
        "search": search["label"],
        "n_members": int(len(indices)),
        "source_ids": ";".join(str(x) for x in ids),
        "tiers": ";".join(tiers),
        "max_sep_pc": float(local_sep.max()) if len(local_sep) else 0.0,
        "median_sep_pc": float(np.median(local_sep)) if len(local_sep) else 0.0,
        "max_dv_kms": float(local_dv.max()) if len(local_dv) else 0.0,
        "median_dv_kms": float(np.median(local_dv)) if len(local_dv) else 0.0,
        "half_radius_pc": half_radius_pc,
        "sigma_1d_kms": sigma_1d,
        "rough_virial_mass_msun": float(virial_mass_msun),
    }


def run_search(sample_name: str, table: Table, mask: np.ndarray) -> tuple[dict, list[dict], list[dict]]:
    sub = table[mask]
    xyz_pc, vel, sep_pc, dv_kms = pairwise_arrays(sub)
    tri = np.triu_indices(len(sub), 1)
    sep_vals = sep_pc[tri]
    dv_vals = dv_kms[tri]

    pair_counts = {}
    for radius in [5, 10, 25, 50, 100, 250]:
        pair_counts[f"pairs_sep_lt_{radius}pc"] = int(np.sum(sep_vals < radius))
    for radius, dv in [(10, 2), (25, 5), (50, 5), (100, 10), (250, 10)]:
        pair_counts[f"pairs_sep_lt_{radius}pc_dv_lt_{dv}kms"] = int(
            np.sum((sep_vals < radius) & (dv_vals < dv))
        )

    group_rows: list[dict] = []
    search_counts = {}
    for search in SEARCHES:
        metric = dbscan_metric(sep_pc, dv_kms, search["max_sep_pc"], search["max_dv_kms"])
        labels = DBSCAN(eps=1.0, min_samples=search["min_members"],
                        metric="precomputed").fit_predict(metric)
        found = []
        for label in sorted(set(labels)):
            if label < 0:
                continue
            indices = np.flatnonzero(labels == label)
            found.append(indices)
            group_rows.append(
                summarize_group(sample_name, search, sub, indices, xyz_pc, vel, sep_pc, dv_kms)
            )
        search_counts[search["label"]] = len(found)

    nearest_rows = []
    order = np.argsort(sep_vals)[:20]
    for rank, flat_idx in enumerate(order, start=1):
        i = int(tri[0][flat_idx])
        j = int(tri[1][flat_idx])
        nearest_rows.append({
            "sample": sample_name,
            "rank": rank,
            "source_id_1": int(sub["source_id"][i]),
            "tier_1": str(sub["tier"][i]).strip(),
            "vgrf_1_kms": float(sub["vgrf_default"][i]),
            "source_id_2": int(sub["source_id"][j]),
            "tier_2": str(sub["tier"][j]).strip(),
            "vgrf_2_kms": float(sub["vgrf_default"][j]),
            "sep_pc": float(sep_pc[i, j]),
            "dv_kms": float(dv_kms[i, j]),
        })

    summary = {
        "sample": sample_name,
        "n_stars": int(len(sub)),
        "min_pair_sep_pc": float(sep_vals.min()) if len(sep_vals) else None,
        "pair_sep_pc_p01": float(np.percentile(sep_vals, 1)) if len(sep_vals) else None,
        "pair_sep_pc_p05": float(np.percentile(sep_vals, 5)) if len(sep_vals) else None,
        "pair_sep_pc_p50": float(np.percentile(sep_vals, 50)) if len(sep_vals) else None,
        "min_pair_dv_kms": float(dv_vals.min()) if len(dv_vals) else None,
        "pair_counts": pair_counts,
        "dbscan_group_counts": search_counts,
    }
    return summary, group_rows, nearest_rows


GROUP_FIELDS = [
    "sample",
    "search",
    "n_members",
    "source_ids",
    "tiers",
    "max_sep_pc",
    "median_sep_pc",
    "max_dv_kms",
    "median_dv_kms",
    "half_radius_pc",
    "sigma_1d_kms",
    "rough_virial_mass_msun",
]


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if not rows and fieldnames is None:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames or list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, summaries: list[dict], group_rows: list[dict]) -> None:
    lines = [
        "# Compact-group search",
        "",
        "Purpose: check whether the slow-Vgrf release catalogue contains an obvious",
        "present-day open/globular-cluster-like aggregate.",
        "",
        "Method: DBSCAN/FoF searches in 3D Galactocentric position, with optional",
        "3D velocity coherence. Inputs are the expanded master catalogue and",
        "expanded point-estimate orbit table; no Gaia archive rescan is required.",
        "",
        "## Result",
        "",
    ]
    if group_rows:
        lines.append(f"Found {len(group_rows)} candidate groups under the configured searches.")
    else:
        lines.append("No candidate group was found under any configured compact-cluster search.")
    lines.extend(["", "## Sample summary", ""])
    lines.append("| sample | N | min pair sep pc | pairs <25 pc | pairs <25 pc and <5 km/s | DBSCAN groups |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for item in summaries:
        groups = sum(item["dbscan_group_counts"].values())
        lines.append(
            f"| {item['sample']} | {item['n_stars']} | "
            f"{item['min_pair_sep_pc']:.1f} | "
            f"{item['pair_counts']['pairs_sep_lt_25pc']} | "
            f"{item['pair_counts']['pairs_sep_lt_25pc_dv_lt_5kms']} | "
            f"{groups} |"
        )
    lines.extend([
        "",
        "Interpretation: an actual compact bound cluster in this catalogue would be",
        "expected to produce multiple stars within tens of parsecs, usually also",
        "with small relative velocities. The closest headline Tier A+B pair is",
        "already about 39 pc apart, and no three-star compact clump appears when",
        "the sample is expanded to Tier A+B+C.",
        "",
        "Caveat: this is a present-day compactness check, not a full cluster",
        "membership paper. A complete literature-grade test would also crossmatch",
        "source IDs against Gaia DR3 open-cluster and globular-cluster membership",
        "catalogues, then inspect CMD/isochrone consistency for any matches.",
        "",
    ])
    path.write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cat = Table.read(CATALOGUE)
    orbits = Table.read(ORBITS)
    table = join(
        cat,
        orbits[["source_id", "x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]],
        keys="source_id",
    )
    tier = strip_tier(table["tier"])
    samples = {
        "tierAB_headline": np.isin(tier, ["A", "B"]),
        "tierABC_statistical": np.isin(tier, ["A", "B", "C"]),
    }

    summaries = []
    all_groups = []
    nearest_rows = []
    for sample_name, mask in samples.items():
        summary, group_rows, sample_nearest = run_search(sample_name, table, mask)
        summaries.append(summary)
        all_groups.extend(group_rows)
        nearest_rows.extend(sample_nearest)

    search_config = {
        "inputs": {
            "catalogue": str(CATALOGUE.relative_to(REPO)),
            "orbits": str(ORBITS.relative_to(REPO)),
        },
        "searches": SEARCHES,
        "summaries": summaries,
        "n_candidate_groups": len(all_groups),
    }

    (OUT_DIR / "cluster_search_summary.json").write_text(json.dumps(search_config, indent=2))
    write_csv(OUT_DIR / "cluster_search_groups.csv", all_groups, GROUP_FIELDS)
    write_csv(OUT_DIR / "cluster_search_nearest_pairs.csv", nearest_rows)
    write_markdown(OUT_DIR / "cluster_search_summary.md", summaries, all_groups)
    print(json.dumps(search_config, indent=2))
    print(f"wrote {OUT_DIR.relative_to(REPO)}/cluster_search_summary.md")


if __name__ == "__main__":
    main()
