"""Phase 2 — gaiaunlimited DR3 selection function.

Builds a HEALPix completeness map for the slow-Vgrf footprint and assigns
per-star numerical-SF weights to the v2 catalogue.

Strategy
========
1. Use `gaiaunlimited.selectionfunctions.DR3RVSSelectionFunction`:
   the official DR3 RVS selection function (Castro-Ginard+2023).
   It returns the probability that a Gaia DR3 source with RV survives,
   given (Galactic l, b) and G_RVS magnitude.

2. For each star in `release/v2/phase1/catalogue_v2.fits`, compute
   sf_value = SF(l, b, G_RVS).
   sf_weight = 1.0 / sf_value (clipped at sf_value > 0.02 to avoid blowups).

3. Aggregate sf_value and sf_weight per HEALPix pixel (Nside=64 -> ~0.84 deg).
   Persist the map and a footprint diagnostic.

Outputs
-------
release/v2/phase2/catalogue_v2_sf.fits          -- per-star SF + weight
release/v2/phase2/sf_healpix_nside64.fits       -- map: pix, n_slow, sum_w
release/v2/phase2/gate2_sf.json                 -- summary stats
release/v2/phase2/gate2_sf_skyplot.png          -- footprint
release/v2/phase2/gate2_sf_weight_histogram.png

Gate 2 expectations
-------------------
- median sf_value in (0.4, 0.95) (most slow stars are bright RVS targets)
- max sf_weight < 50 (the clip)
- pixel-summed weights map should look like a sky distribution, not a
  single hot pixel.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import healpy as hp
from astropy.table import Table
import astropy.units as u
import astropy.coordinates as coord

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT_DIR = REPO / "release/v2/phase2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

V2 = REPO / "release/v2/phase1/catalogue_v2.fits"
SRC = REPO / CONFIG["input"]["source_csv"]

# DR3RVSSelectionFunction factors: g = G_RVS apparent mag, c = G_BP - G_RP.
# Source: gaiaunlimited subsample SF (Castro-Ginard+2023, Zenodo 8300616).

NSIDE = 64
SF_FLOOR = 0.02   # weight clip floor


def main() -> int:
    cat = Table.read(V2).to_pandas()
    # Decode bytes-typed string columns coming from FITS round-trip.
    for c in ("tier", "dist_source", "rv_quality"):
        if c in cat.columns and cat[c].dtype == object:
            sample = cat[c].iloc[0] if len(cat) else b""
            if isinstance(sample, (bytes, bytearray)):
                cat[c] = cat[c].str.decode("utf-8")
    src = pd.read_csv(SRC)[
        ["source_id", "l", "b", "phot_g_mean_mag", "grvs_mag", "bp_rp"]
    ]
    df = cat.merge(src, on="source_id")
    print(f"[2] {len(df)} stars")

    from gaiaunlimited.selectionfunctions import DR3RVSSelectionFunction
    sf = DR3RVSSelectionFunction()

    coords_icrs = coord.SkyCoord(
        ra=df["ra"].to_numpy() * u.deg,
        dec=df["dec"].to_numpy() * u.deg,
        frame="icrs",
    )
    grvs = df["grvs_mag"].to_numpy(dtype=float)
    bp_rp = df["bp_rp"].to_numpy(dtype=float)
    # Replace NaN bp_rp with sample median so the SF can be evaluated;
    # mark the row as sf_invalid for downstream filtering.
    bp_rp_med = float(np.nanmedian(bp_rp))
    bp_rp_filled = np.where(np.isfinite(bp_rp), bp_rp, bp_rp_med)
    bp_rp_is_nan = ~np.isfinite(bp_rp)
    sf_value = sf.query(coords_icrs, g=grvs, c=bp_rp_filled, fill_nan=True)
    sf_value = np.asarray(sf_value, dtype=float)
    print(f"[2] sf_value: min={np.nanmin(sf_value):.3f} "
          f"med={np.nanmedian(sf_value):.3f} max={np.nanmax(sf_value):.3f}")

    sf_value_clipped = np.clip(sf_value, SF_FLOOR, 1.0)
    sf_weight = 1.0 / sf_value_clipped
    sf_invalid = ~np.isfinite(sf_value) | (sf_value <= 0) | bp_rp_is_nan

    # Per-star output
    out_star = Table()
    out_star["source_id"] = df["source_id"].to_numpy()
    out_star["l"] = df["l"].to_numpy()
    out_star["b"] = df["b"].to_numpy()
    out_star["grvs_mag"] = grvs
    out_star["sf_value"] = sf_value
    out_star["sf_weight"] = sf_weight
    out_star["sf_invalid"] = sf_invalid
    out_star["tier"] = df["tier"].to_numpy().astype(str)
    out_star_path = OUT_DIR / "catalogue_v2_sf.fits"
    out_star.write(out_star_path, format="fits", overwrite=True)
    print(f"[2] wrote {out_star_path.relative_to(REPO)}")

    # HEALPix aggregation
    galactic = coord.SkyCoord(l=df["l"].to_numpy() * u.deg,
                              b=df["b"].to_numpy() * u.deg,
                              frame="galactic")
    theta = (90.0 - galactic.b.value) * np.pi / 180.0
    phi = galactic.l.value * np.pi / 180.0
    pix = hp.ang2pix(NSIDE, theta, phi)
    n_pix = hp.nside2npix(NSIDE)

    n_slow = np.bincount(pix, minlength=n_pix)
    sum_w  = np.bincount(pix, weights=sf_weight, minlength=n_pix)
    sum_v  = np.bincount(pix, weights=sf_value, minlength=n_pix)
    out_map = Table()
    out_map["pix"] = np.arange(n_pix, dtype=np.int64)
    out_map["n_slow"] = n_slow
    out_map["sum_sf_weight"] = sum_w
    out_map["sum_sf_value"] = sum_v
    map_path = OUT_DIR / f"sf_healpix_nside{NSIDE}.fits"
    out_map.write(map_path, format="fits", overwrite=True)
    print(f"[2] wrote {map_path.relative_to(REPO)}")

    # ---- Tier-B selection summary ----
    tierAB = (df["tier"] == "A") | (df["tier"] == "B")
    tierABC = tierAB | (df["tier"] == "C")
    summary = {
        "n_total": int(len(df)),
        "sf_value_p16": float(np.nanpercentile(sf_value, 16)),
        "sf_value_p50": float(np.nanpercentile(sf_value, 50)),
        "sf_value_p84": float(np.nanpercentile(sf_value, 84)),
        "sf_value_min": float(np.nanmin(sf_value)),
        "sf_value_max": float(np.nanmax(sf_value)),
        "n_sf_invalid": int(sf_invalid.sum()),
        "sf_floor_clip": SF_FLOOR,
        "sf_weight_p50": float(np.nanmedian(sf_weight)),
        "sf_weight_max": float(np.nanmax(sf_weight)),
        "tier_AB_sum_weight": float(np.sum(sf_weight[tierAB])),
        "tier_AB_n_unweighted": int(tierAB.sum()),
        "tier_ABC_sum_weight": float(np.sum(sf_weight[tierABC])),
        "tier_ABC_n_unweighted": int(tierABC.sum()),
        "n_pixels_occupied": int((n_slow > 0).sum()),
        "nside": NSIDE,
    }
    (OUT_DIR / "gate2_sf.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    # ---- Plots ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(sf_weight[np.isfinite(sf_weight)], bins=80, color="steelblue",
            edgecolor="white")
    ax.set_xlabel("sf_weight (1/SF, clipped)")
    ax.set_ylabel("count")
    ax.set_title("Phase 2 — per-star SF weight distribution")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "gate2_sf_weight_histogram.png", dpi=140)
    plt.close(fig)

    # Mollweide of n_slow
    map_log = np.log10(np.maximum(n_slow, 1))
    map_log[n_slow == 0] = hp.UNSEEN
    fig = plt.figure(figsize=(9, 5))
    hp.mollview(map_log, fig=fig.number, coord="G",
                title=f"Phase 2 — slow-star count map (log10), Nside={NSIDE}",
                cmap="viridis")
    plt.savefig(OUT_DIR / "gate2_sf_skyplot.png", dpi=140,
                bbox_inches="tight")
    plt.close(fig)
    print("[2] wrote gate plots")

    if summary["sf_value_p50"] < 0.1 or summary["sf_value_p50"] > 0.999:
        print("[2] ANOMALY: median SF outside (0.10, 0.999); STOP.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
