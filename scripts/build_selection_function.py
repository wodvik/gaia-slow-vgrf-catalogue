"""Expanded Gaia DR3 RVS selection-function table.

Evaluates the Castro-Ginard et al. (2023) Gaia DR3 RVS selection
function, via gaiaunlimited, for the parent-complete expanded catalogue.

Outputs are written under analysis_products/. The
per-source table includes both the GaiaUnlimited selection value and the
underlying parent count used for the low-n reproducibility audit.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import xarray as xr
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u


BUNDLE = Path(__file__).resolve().parents[1]
IN_FITS = BUNDLE / "catalogues" / "catalogue_expanded_master.fits"
OUT = BUNDLE / "analysis_products"
OUT_FITS = OUT / "expanded_selection_function.fits"
OUT_CSV = OUT / "expanded_selection_function_summary.csv"
OUT_JSON = OUT / "expanded_selection_function_summary.json"
LOWN_AUDIT = OUT / "expanded_selection_function_lown_audit.csv"
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


def query_parent_count(sf, coords: SkyCoord, grvs: np.ndarray, bp_rp: np.ndarray) -> np.ndarray:
    """Return the GaiaUnlimited DR3-RVS parent count for each queried cell.

    DR3RVSSelectionFunction stores p=(k+1)/(n+2) on a Galactic nested
    HEALPix x G_RVS x BP-RP grid. The public query method interpolates
    logit(p) with nearest-neighbour coordinates and fills missing cells
    to the prior mean p=0.5. For the low-n audit, missing grid cells are
    therefore parent-count zero cells.
    """
    from gaiaunlimited import utils

    ipix = utils.coord2healpix(coords, "galactic", sf.nside, nest=True)
    parent_raw = sf.ds["n"].interp(
        ipix=xr.DataArray(np.atleast_1d(ipix)),
        g=xr.DataArray(np.atleast_1d(grvs)),
        c=xr.DataArray(np.atleast_1d(bp_rp)),
        method="nearest",
        kwargs=dict(fill_value=None),
    ).to_numpy()
    return np.nan_to_num(parent_raw, nan=0.0).astype(np.int64)


def summarize_sample(
    label: str,
    mask: np.ndarray,
    valid: np.ndarray,
    sf_value: np.ndarray,
    sf_weight: np.ndarray,
    low_parent_count: np.ndarray,
) -> dict[str, object]:
    m = mask & valid
    vals = sf_value[m]
    weights = sf_weight[m]
    low = m & low_parent_count
    low_weights = sf_weight[low]
    sum_weights = float(np.nansum(weights))
    sum_low_weights = float(np.nansum(low_weights))
    return {
        "sample": label,
        "n_unweighted": int(mask.sum()),
        "n_valid_sf": int(m.sum()),
        "n_invalid_sf": int((mask & ~valid).sum()),
        "sf_value_p16": float(np.nanpercentile(vals, 16)) if len(vals) else np.nan,
        "sf_value_p50": float(np.nanpercentile(vals, 50)) if len(vals) else np.nan,
        "sf_value_p84": float(np.nanpercentile(vals, 84)) if len(vals) else np.nan,
        "sum_inverse_sf_weight": sum_weights,
        "median_inverse_sf_weight": float(np.nanmedian(weights)) if len(weights) else np.nan,
        "n_low_parent_count": int(low.sum()),
        "frac_low_parent_count_pct": float(100.0 * low.sum() / m.sum()) if m.sum() else np.nan,
        "sum_inverse_sf_weight_low_parent_count": sum_low_weights,
        "frac_inverse_sf_weight_low_parent_count_pct": (
            float(100.0 * sum_low_weights / sum_weights) if sum_weights else np.nan
        ),
    }


def audit_rows(
    masks: dict[str, np.ndarray],
    valid: np.ndarray,
    sf_weight: np.ndarray,
    low_parent_count: np.ndarray,
    inner: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    specs: Iterable[tuple[str, str, str]] = (
        ("candidate_pool", "expanded_candidate_pool", "full"),
        ("TierA", "Tier A", "full"),
        ("TierB", "Tier B", "full"),
        ("TierC", "Tier C", "full"),
        ("TierAB", "Tier A+B", "full"),
        ("TierABC", "Tier A+B+C", "full"),
        ("TierA", "Tier A", "inner_l330_30"),
        ("TierAB", "Tier A+B", "inner_l330_30"),
        ("TierABC", "Tier A+B+C", "inner_l330_30"),
    )
    for sample, mask_key, region in specs:
        region_mask = inner if region == "inner_l330_30" else np.ones_like(valid, dtype=bool)
        raw_mask = masks[mask_key] & region_mask
        m = raw_mask & valid
        low = m & low_parent_count
        sum_weights = float(np.nansum(sf_weight[m]))
        sum_low_weights = float(np.nansum(sf_weight[low]))
        rows.append({
            "sample": sample,
            "region": region,
            "N": int(raw_mask.sum()),
            "N_valid": int(m.sum()),
            "sum_inv_S": sum_weights,
            "f_n_lt10_pct": float(100.0 * low.sum() / m.sum()) if m.sum() else np.nan,
            "n_lt10_count": int(low.sum()),
            "f_sum_n_lt10_pct": (
                float(100.0 * sum_low_weights / sum_weights) if sum_weights else np.nan
            ),
            "sum_inv_S_n_lt10": sum_low_weights,
        })
    return rows


def write_lown_audit(rows: list[dict[str, object]]) -> None:
    lines = [
        "# Selection-function low-parent-count (n<10) audit summary.",
        "# The per-source parent counts are stored in expanded_selection_function.fits",
        "# as sf_parent_count; sf_prior_dominated_n_lt10 is true when sf_parent_count < 10.",
        "# GaiaUnlimited DR3-RVS cells absent from dr3-rvs-nk.h5 are encoded as",
        "# sf_parent_count = 0, matching the p=0.5 prior-mean fill used by",
        "# DR3RVSSelectionFunction.query(..., fill_nan=True).",
    ]
    csv = pd.DataFrame(rows).to_csv(index=False, float_format="%.6g")
    LOWN_AUDIT.write_text("\n".join(lines) + "\n" + csv)


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
    sf_parent_count = query_parent_count(sf, coords, grvs, bp_rp_filled)
    sf_prior_dominated_n_lt10 = sf_parent_count < 10
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
    out["sf_parent_count"] = sf_parent_count
    out["sf_prior_dominated_n_lt10"] = sf_prior_dominated_n_lt10
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
        rows.append(summarize_sample(
            label,
            mask,
            valid,
            sf_value,
            sf_weight,
            sf_prior_dominated_n_lt10,
        ))
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_CSV, index=False)
    inner = (out["l"].data >= 330.0) | (out["l"].data <= 30.0)
    lown_rows = audit_rows(masks, valid, sf_weight, sf_prior_dominated_n_lt10, inner)
    write_lown_audit(lown_rows)
    payload = {
        "input": str(IN_FITS.relative_to(BUNDLE)),
        "output_fits": str(OUT_FITS.relative_to(BUNDLE)),
        "gaiaunlimited_model": "DR3RVSSelectionFunction",
        "gaiaunlimited_version": getattr(__import__("gaiaunlimited"), "__version__", "unknown"),
        "healpix_order": int(sf.order),
        "healpix_nside": int(sf.nside),
        "parent_count_nan_encoded_as_zero": True,
        "low_parent_count_threshold": 10,
        "sf_floor": SF_FLOOR,
        "n_rows": int(len(cat)),
        "n_invalid_sf": int(sf_invalid.sum()),
        "samples": rows,
        "low_parent_count_audit": lown_rows,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    STATUS.write_text(
        "# Expanded Selection-Function Status\n\n"
        "Release: v1.0.3-review\n\n"
        "The expanded parent-complete catalogue has been evaluated with "
        "the Castro-Ginard et al. (2023) Gaia DR3 RVS selection function "
        "through `gaiaunlimited.selectionfunctions.DR3RVSSelectionFunction` "
        "in a WSL Python environment. The per-star table is "
        "`expanded_selection_function.fits`; it includes `sf_parent_count` "
        "and `sf_prior_dominated_n_lt10` columns for reproducing the "
        "low-parent-count audit. Summary statistics are in "
        "`expanded_selection_function_summary.csv` and "
        "`expanded_selection_function_summary.json`; the table-level "
        "`n<10` reproduction audit is in "
        "`expanded_selection_function_lown_audit.csv`.\n\n"
        "These inverse-selection weights are retained as contextual "
        "observability diagnostics, not as a volume-complete Milky Way "
        "deprojection or headline population count.\n"
    )
    log(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    T0 = time.time()
    raise SystemExit(main())
