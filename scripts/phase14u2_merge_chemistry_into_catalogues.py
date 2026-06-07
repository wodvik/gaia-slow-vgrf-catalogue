"""Phase 14U-2 -- merge per-star chemistry into the public catalogue products.

Referee Issue 2 (data bundle): the primary catalogue FITS shipped without
chemistry columns even though the per-star APOGEE/GALAH + Gaia GSP-Phot
products were present in ``phase14/``.  This step joins a compact, clearly
separated set of chemistry columns onto every catalogue FITS by ``source_id``
so the released tables are self-describing for chemodynamic reuse.

Source of truth: ``phase14/expanded_population_context.fits`` produced by
``phase14u_expanded_chemistry.py`` (which itself reproduces the manuscript
chemistry numbers: 1,279 finite GSP-Phot [M/H] and the 117-star spectroscopic
alpha subset in Tier A+B+C).

Design choices:
- Photometric (GSP-Phot) and spectroscopic ([Fe/H]) metallicities are kept in
  SEPARATE columns; they are never blended.  GSP-Phot is a biased low-resolution
  diagnostic for this cool, metal-poor giant population (Andrae et al. 2023).
- ``chem_population`` (Splash/GSE/Aurora/disk/unclassified) is populated ONLY for
  the spectroscopic alpha subset (finite ``feh_spec`` AND ``alpha_spec``); blank
  otherwise, so the chemodynamic class never rests on photometric metallicity.

Idempotent: existing chemistry columns are replaced.  Run from the bundle root:
    python scripts/phase14u2_merge_chemistry_into_catalogues.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
CTX = BUNDLE / "phase14" / "expanded_population_context.fits"
CAT = BUNDLE / "catalogues"

CATALOGUE_FITS = [
    "catalogue_expanded_master.fits",
    "catalogue_expanded_tierA.fits",
    "catalogue_expanded_tierAB.fits",
    "catalogue_expanded_tierABC.fits",
    "catalogue_expanded_orbits_tierABC.fits",
]

FLOAT_COLS = ["mh_gspphot", "mh_gspphot_lo", "mh_gspphot_hi",
              "feh_spec", "feh_spec_err", "alpha_spec"]
STR_COLS = ["chem_survey", "chem_population"]
CHEM_COLS = FLOAT_COLS + STR_COLS

DESCR = {
    "mh_gspphot": "Gaia DR3 GSP-Phot photometric metallicity [M/H] (dex); biased "
                  "toward solar for cool metal-poor giants, contextual only.",
    "mh_gspphot_lo": "Lower (16th percentile) GSP-Phot [M/H] bound (dex).",
    "mh_gspphot_hi": "Upper (84th percentile) GSP-Phot [M/H] bound (dex).",
    "feh_spec": "Spectroscopic [Fe/H] (dex) from the best APOGEE DR17 / GALAH DR3 "
                "exact-source_id match; NaN if no spectroscopic match.",
    "feh_spec_err": "Reported error on feh_spec (dex); NaN if unavailable.",
    "alpha_spec": "Spectroscopic alpha-abundance proxy (dex): [alpha/M] for APOGEE, "
                  "[alpha/Fe] for GALAH (see chem_survey); NaN if unavailable.",
    "chem_survey": "Spectroscopic source of feh_spec/alpha_spec: 'APOGEE', 'GALAH', "
                   "or '' (none). APOGEE preferred when a star is in both.",
    "chem_population": "Chemodynamic class from spectroscopic [Fe/H]+alpha "
                       "(Splash/GSE/Aurora/disk/unclassified) per "
                       "phase14u_expanded_chemistry.classify; '' unless both "
                       "feh_spec and alpha_spec are finite.",
}


def load_chem() -> pd.DataFrame:
    ctx = Table.read(CTX).to_pandas()
    for c in ["survey", "population_context", "tier"]:
        if c in ctx and ctx[c].dtype == object and len(ctx) and isinstance(ctx[c].iloc[0], (bytes, bytearray)):
            ctx[c] = ctx[c].str.decode("utf-8")
    feh = pd.to_numeric(ctx["FeH"], errors="coerce")
    alpha = pd.to_numeric(ctx["alpha_proxy"], errors="coerce")
    spec_ok = feh.notna()
    alpha_ok = feh.notna() & alpha.notna()
    chem = pd.DataFrame({
        "source_id": pd.to_numeric(ctx["source_id"], errors="coerce").astype("int64"),
        "mh_gspphot": pd.to_numeric(ctx["MH_gspphot"], errors="coerce"),
        "mh_gspphot_lo": pd.to_numeric(ctx["MH_gspphot_lo"], errors="coerce"),
        "mh_gspphot_hi": pd.to_numeric(ctx["MH_gspphot_hi"], errors="coerce"),
        "feh_spec": feh,
        "feh_spec_err": pd.to_numeric(ctx["e_FeH"], errors="coerce"),
        "alpha_spec": alpha,
        "chem_survey": np.where(spec_ok, ctx["survey"].astype(str), ""),
        "chem_population": np.where(alpha_ok, ctx["population_context"].astype(str), ""),
    })
    for c in STR_COLS:
        chem[c] = chem[c].astype(str).replace({"nan": "", "None": ""})
    if chem["source_id"].duplicated().any():
        raise RuntimeError("duplicate source_id in chemistry context")
    return chem.set_index("source_id")


def merge_one(path: Path, chem: pd.DataFrame) -> dict:
    t = Table.read(path)
    sid = np.asarray(t["source_id"]).astype("int64")
    sub = chem.reindex(sid)
    for col in CHEM_COLS:
        if col in t.colnames:
            t.remove_column(col)
        vals = sub[col].to_numpy()
        if col in FLOAT_COLS:
            arr = np.array([np.nan if (v is None or (isinstance(v, float) and pd.isna(v))) else float(v)
                            for v in vals], dtype="float64")
        else:
            arr = np.array(["" if (v is None or (isinstance(v, float) and pd.isna(v))) else str(v)
                            for v in vals], dtype="U16")
        t[col] = arr
        t[col].description = DESCR[col]
    t.write(path, overwrite=True)
    n = len(t)
    n_mh = int(np.isfinite(np.asarray(t["mh_gspphot"], dtype=float)).sum())
    n_feh = int(np.isfinite(np.asarray(t["feh_spec"], dtype=float)).sum())
    n_alpha = int((np.asarray(t["chem_population"]) != "").sum())
    return {"file": path.name, "n": n, "finite_mh": n_mh, "finite_feh": n_feh, "alpha_pop": n_alpha}


def main() -> int:
    chem = load_chem()
    print(f"chemistry context: {len(chem)} sources from {CTX.relative_to(BUNDLE)}")
    rows = []
    for name in CATALOGUE_FITS:
        p = CAT / name
        if not p.exists():
            print(f"  SKIP (missing): {name}")
            continue
        rows.append(merge_one(p, chem))
    print(f"\n{'file':42s} {'n':>6} {'finite_mh':>10} {'finite_feh':>11} {'alpha_pop':>10}")
    for r in rows:
        print(f"{r['file']:42s} {r['n']:6d} {r['finite_mh']:10d} {r['finite_feh']:11d} {r['alpha_pop']:10d}")
    print("\nExpected Tier A+B+C: finite_mh=1279, alpha_pop=117 (see expanded_chemistry_summary.json)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
