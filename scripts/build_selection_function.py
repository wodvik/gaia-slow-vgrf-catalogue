"""Expanded Gaia DR3 RVS selection-function table.

Evaluates the Castro-Ginard et al. (2023) Gaia DR3 RVS selection
function, via gaiaunlimited, for the parent-complete expanded catalogue.

Outputs are written under analysis_products and are intended to replace
the temporary "selection-function status" note used before the WSL
gaiaunlimited/healpy stack was available.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u


REPO = Path(__file__).resolve().parents[1]
IN_FITS = REPO / "catalogues/catalogue_expanded_master.fits"
OUT = REPO / "analysis_products"
OUT_FITS = OUT / "expanded_selection_function.fits"
OUT_CSV = OUT / "expanded_selection_function_summary.csv"
OUT_JSON = OUT / "expanded_selection_function_summary.json"
STATUS = OUT / "expanded_selection_function_status.md"

SF_FLOOR = 0.02


def log(message: str) -> None:
    print(f"[14Z-sf t={time.time() - T0:7.1f}s] {message}", flush=True)


def clean_tier(series: pd.Series) -> pd.Series:
    if series.dtype == object and len(series):
        first = series.iloc[0]
        if isinstance(first, (bytes, bytearray)):
            return series.str.decode("utf-8").str.strip()
    return series.astype(str).str.strip()


def finite_or_nan(values: pd.Series) -> np.ndarray:
    return pd.to_numeric(values, errors="coerce").to_numpy(float)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cat = Table.read(IN_FITS).to_pandas()
    cat["tier"] = clean_tier(cat["tier"])

    from gaiaunlimited.selectionfunctions import DR3RVSSelectionFunction

    log(f"loaded {IN_FITS} ({len(cat):,} rows)")
    sf = DR3RVSSelectionFunction()
    log("loaded gaiaunlimited DR3RVSSelectionFunction")

    coords = SkyCoord(
        ra=finite_or_nan(cat["ra"]) * u.deg,
        dec=finite_or_nan(cat["dec"]) * u.deg,
        frame="icrs",
    )
    grvs = finite_or_nan(cat["grvs_mag"])
    bp_rp = finite_or_nan(cat["bp_rp"])
    bp_rp_med = float(np.nanmedian(bp_rp))
    bp_rp_filled = np.where(np.isfinite(bp_rp), bp_rp, bp_rp_med)
    sf_value = np.asarray(sf.query(coords, g=grvs, c=bp_rp_filled, fill_nan=True), dtype=float)
    sf_invalid = (
        ~np.isfinite(sf_value)
        | (sf_value <= 0)
        | ~np.isfinite(grvs)
        | ~np.isfinite(bp_rp)
    )
    sf_weight = 1.0 / np.clip(sf_value, SF_FLOOR, 1.0)

    out = Table()
    out["source_id"] = cat["source_id"].to_numpy()
    out["tier"] = cat["tier"].to_numpy(dtype=str)
    out["l"] = finite_or_nan(cat["l"])
    out["b"] = finite_or_nan(cat["b"])
    out["grvs_mag"] = grvs
    out["bp_rp"] = bp_rp
    out["sf_value"] = sf_value
    out["sf_weight"] = sf_weight
    out["sf_invalid"] = sf_invalid
    out["P_vgrf_below_25"] = finite_or_nan(cat["P_vgrf_below_25"])
    out.write(OUT_FITS, overwrite=True)
    log(f"wrote {OUT_FITS}")

    rows: list[dict[str, object]] = []
    masks = {
        "expanded_candidate_pool": np.ones(len(cat), dtype=bool),
        "Tier A": cat["tier"].eq("A").to_numpy(),
        "Tier B": cat["tier"].eq("B").to_numpy(),
        "Tier C": cat["tier"].eq("C").to_numpy(),
        "Tier A+B": cat["tier"].isin(["A", "B"]).to_numpy(),
        "Tier A+B+C": cat["tier"].isin(["A", "B", "C"]).to_numpy(),
    }
    valid = ~sf_invalid
    for label, mask in masks.items():
        m = mask & valid
        vals = sf_value[m]
        weights = sf_weight[m]
        rows.append({
            "sample": label,
            "n_unweighted": int(mask.sum()),
            "n_valid_sf": int(m.sum()),
            "n_invalid_sf": int((mask & sf_invalid).sum()),
            "sf_value_p16": float(np.nanpercentile(vals, 16)) if len(vals) else np.nan,
            "sf_value_p50": float(np.nanpercentile(vals, 50)) if len(vals) else np.nan,
            "sf_value_p84": float(np.nanpercentile(vals, 84)) if len(vals) else np.nan,
            "sum_inverse_sf_weight": float(np.nansum(weights)),
            "median_inverse_sf_weight": float(np.nanmedian(weights)) if len(weights) else np.nan,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_CSV, index=False)
    payload = {
        "input": str(IN_FITS),
        "output_fits": str(OUT_FITS),
        "sf_floor": SF_FLOOR,
        "n_rows": int(len(cat)),
        "n_invalid_sf": int(sf_invalid.sum()),
        "samples": rows,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    STATUS.write_text(
        "# Expanded Selection-Function Status\n\n"
        "Date: 2026-05-10\n\n"
        "The expanded parent-complete catalogue has been evaluated with "
        "the Castro-Ginard et al. (2023) Gaia DR3 RVS selection function "
        "through `gaiaunlimited.selectionfunctions.DR3RVSSelectionFunction` "
        "in a WSL Python environment.  The per-star table is "
        "`expanded_selection_function.fits`; summary statistics are in "
        "`expanded_selection_function_summary.csv` and "
        "`expanded_selection_function_summary.json`.\n\n"
        "These inverse-selection weights are retained as contextual "
        "observability diagnostics, not as a volume-complete Milky Way "
        "deprojection or headline population count.\n"
    )
    log(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    T0 = time.time()
    raise SystemExit(main())
