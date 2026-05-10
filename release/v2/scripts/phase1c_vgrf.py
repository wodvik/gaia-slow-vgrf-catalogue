"""Phase 1C — Galactocentric coords / V_grf under four solar variants.

Inputs
------
- release/v2/phase1/catalogue_dist.fits       : BJ distance per source
- release/data/slow_stars_enriched_orbits.csv : ra, dec, pm*, RV

Outputs
-------
- release/v2/phase1/catalogue_vgrf.fits :
    source_id,
    vgrf_default, vgrf_grav22, vgrf_lsr6, vgrf_rb20,
    XYZ + VxVyVz under default,
    R_cyl, phi_cyl, z_cyl, vR, vphi, vz under default
- release/v2/phase1/gate1C_vgrf.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT_DIR = REPO / "release/v2/phase1"
SRC_CSV = REPO / CONFIG["input"]["source_csv"]
DIST_PATH = OUT_DIR / "catalogue_dist.fits"
VGRF_CUT = float(CONFIG["vgrf_cutoff_kms"])

VARIANTS = CONFIG["solar_variants"]


def galcen_v_sun(name: str, params: dict) -> coord.CartesianDifferential:
    """Build Galactocentric solar-motion vector for a given variant.

    Convention: galcen_v_sun = (U_sun, V_LSR + V_sun_pec, W_sun)
    where V_LSR = circular speed at solar radius.
    For 'rb20', the published 248.5 km/s already represents
    V_c + V_sun_pec at R0 (the implicit sum), so do not add V again.
    """
    U = params["U_kms"]
    W = params["W_kms"]
    if name == "rb20":
        Vy = params["Vc_kms"]            # already V_c + V_pec
    else:
        Vy = params["Vc_kms"] + params["V_kms"]
    return coord.CartesianDifferential(U * u.km/u.s,
                                       Vy * u.km/u.s,
                                       W * u.km/u.s)


def vgrf_for_variant(icrs: coord.SkyCoord, name: str,
                     params: dict) -> tuple[np.ndarray, coord.SkyCoord]:
    frame = coord.Galactocentric(
        galcen_distance=params["R0_kpc"] * u.kpc,
        z_sun=params["z_sun_pc"] * u.pc,
        galcen_v_sun=galcen_v_sun(name, params),
    )
    g = icrs.transform_to(frame)
    vx = g.v_x.to_value(u.km/u.s)
    vy = g.v_y.to_value(u.km/u.s)
    vz = g.v_z.to_value(u.km/u.s)
    vgrf = np.sqrt(vx**2 + vy**2 + vz**2)
    return vgrf, g


def main() -> int:
    src = pd.read_csv(SRC_CSV)
    dist = Table.read(DIST_PATH).to_pandas()
    df = src.merge(dist[["source_id", "dist_pc", "dist_source",
                         "parallax_zpcorr"]],
                   on="source_id", how="inner")
    print(f"[1C] {len(df)} merged rows")

    # Drop stars without finite distance/RV (essential for Vgrf).
    needed = ["ra", "dec", "pmra", "pmdec", "radial_velocity", "dist_pc"]
    finite = np.all(np.isfinite(df[needed].to_numpy(dtype=float)), axis=1)
    n_drop = int((~finite).sum())
    if n_drop:
        print(f"[1C] dropping {n_drop} stars with non-finite RV/distance/PM")
    df = df.loc[finite].reset_index(drop=True)

    icrs = coord.SkyCoord(
        ra=df["ra"].to_numpy() * u.deg,
        dec=df["dec"].to_numpy() * u.deg,
        distance=df["dist_pc"].to_numpy() * u.pc,
        pm_ra_cosdec=df["pmra"].to_numpy() * u.mas/u.yr,
        pm_dec=df["pmdec"].to_numpy() * u.mas/u.yr,
        radial_velocity=df["radial_velocity"].to_numpy() * u.km/u.s,
        frame="icrs",
    )

    out = Table()
    out["source_id"] = df["source_id"].to_numpy()
    out["dist_pc"] = df["dist_pc"].to_numpy()
    out["v_grf_reference"] = df["V_grf"].to_numpy()  # reference column

    counts: dict[str, int] = {}
    for name, params in VARIANTS.items():
        vgrf, g = vgrf_for_variant(icrs, name, params)
        col = f"vgrf_{name}"
        out[col] = vgrf
        counts[name] = int(np.sum(vgrf < VGRF_CUT))
        print(f"[1C] variant {name}: N(Vgrf<{VGRF_CUT}) = {counts[name]}")
        if name == "default":
            out["x_kpc"] = g.x.to_value(u.kpc)
            out["y_kpc"] = g.y.to_value(u.kpc)
            out["z_kpc"] = g.z.to_value(u.kpc)
            out["vx_kms"] = g.v_x.to_value(u.km/u.s)
            out["vy_kms"] = g.v_y.to_value(u.km/u.s)
            out["vz_kms"] = g.v_z.to_value(u.km/u.s)
            x = g.x.to_value(u.kpc)
            y = g.y.to_value(u.kpc)
            out["R_cyl_kpc"] = np.sqrt(x*x + y*y)
            out["phi_cyl_rad"] = np.arctan2(y, x)
            out["z_cyl_kpc"] = g.z.to_value(u.kpc)

    out_path = OUT_DIR / "catalogue_vgrf.fits"
    out.write(out_path, format="fits", overwrite=True)
    print(f"[1C] wrote {out_path.relative_to(REPO)} ({len(out)} rows)")

    # Shift vs reference column.
    shift = out["vgrf_default"].data - out["v_grf_reference"].data
    pct = np.percentile(shift[np.isfinite(shift)], [16, 50, 84])

    n_reference_below = int(np.sum(out["v_grf_reference"].data < VGRF_CUT))
    summary = {
        "n_processed": int(len(out)),
        "vgrf_cutoff_kms": VGRF_CUT,
        "counts_below_cut_by_variant": counts,
        "n_reference_below_cut": n_reference_below,
        "vgrf_shift_vs_v1_p16": float(pct[0]),
        "vgrf_shift_vs_v1_p50": float(pct[1]),
        "vgrf_shift_vs_v1_p84": float(pct[2]),
        "n_dropped_nonfinite": n_drop,
    }
    (OUT_DIR / "gate1C_vgrf.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    # Plain-language statement (printed only).
    if abs(pct[1]) < 1.0 and max(counts.values()) - min(counts.values()) > 50:
        msg = ("Solar-parameter choice dominates the count change "
               "(median per-star shift small, but variant spread > 50 stars).")
    elif abs(pct[1]) > 5.0:
        msg = ("ZP+BJ-distance correction dominates the count change "
               "(median per-star Vgrf shift > 5 km/s).")
    else:
        msg = ("ZP and solar params contribute roughly comparably.")
    print(f"[1C] interpretation: {msg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
