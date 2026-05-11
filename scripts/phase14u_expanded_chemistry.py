"""Phase 14U -- expanded catalogue chemistry synchronization.

This pass keeps the expanded paper internally synchronized with the full
Phase 14W APOGEE/GALAH cross-match when available, falling back to the
legacy cached matches only if the full expanded files are absent.

Outputs are written under release/v2/phase14 so the legacy Phase 5
products remain intact and auditable.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
OUT = REPO / "release/v2/phase14"
OUT.mkdir(parents=True, exist_ok=True)

MASTER = REPO / "release/v2/phase0_expanded/catalogue_expanded_master.fits"
SRC = REPO / "release/data/slow_stars_expanded_candidates_vgrf50.csv"
XDIR = REPO / "release/v2/phase14/expanded_spectroscopic_crossmatch"
APOGEE = XDIR / "apogee_dr17_expanded_crossmatch.csv"
GALAH = XDIR / "galah_dr3_expanded_crossmatch.csv"
APOGEE_FALLBACK = REPO / "release/data/crossmatch_apogee.csv"
GALAH_FALLBACK = REPO / "release/data/crossmatch_galah.csv"


def decode_strings(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        if df[col].dtype == object and len(df):
            val = df[col].iloc[0]
            if isinstance(val, (bytes, bytearray)):
                df[col] = df[col].str.decode("utf-8")
    return df


def finite_stat(series: pd.Series) -> dict[str, float | int | None]:
    arr = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return {"n": 0, "p16": None, "p50": None, "p84": None}
    return {
        "n": int(len(arr)),
        "p16": float(np.percentile(arr, 16)),
        "p50": float(np.percentile(arr, 50)),
        "p84": float(np.percentile(arr, 84)),
    }


def classify(feh: float, alpha: float) -> str:
    if not np.isfinite(feh):
        return "unclassified"
    if feh >= 0.0:
        return "disk"
    if not np.isfinite(alpha):
        if feh < -1.0:
            return "low_MH_no_alpha"
        return "intermediate_MH_no_alpha"
    if feh < -1.0 and alpha >= 0.25:
        return "Aurora"
    if feh < -0.7 and alpha < 0.20:
        return "GSE"
    if -1.0 <= feh < 0.0 and alpha >= 0.15:
        return "Splash"
    if feh >= -0.3 and alpha < 0.15:
        return "disk"
    return "unclassified"


def load_spectroscopy(master_min: pd.DataFrame) -> pd.DataFrame:
    using_full = APOGEE.exists() and GALAH.exists()
    apo = pd.read_csv(APOGEE if using_full else APOGEE_FALLBACK)
    gal = pd.read_csv(GALAH if using_full else GALAH_FALLBACK)

    apo_clean = apo.rename(columns={
        "__Fe_H_": "FeH",
        "e__Fe_H_": "e_FeH",
        "__a_M_": "alpha_proxy",
        "__Mg_Fe_": "MgFe",
        "__Al_Fe_": "AlFe",
        "__Mn_Fe_": "MnFe",
    })[["source_id", "Teff", "logg", "FeH", "e_FeH", "alpha_proxy",
        "MgFe", "AlFe", "MnFe"]]
    apo_clean["survey"] = "APOGEE"
    apo_clean["survey_rank"] = 0

    gal_clean = gal.rename(columns={
        "__Fe_H_": "FeH",
        "e__Fe_H_": "e_FeH",
        "__alpha_Fe_": "alpha_proxy",
        "__Mg_Fe_": "MgFe",
        "__Al_Fe_": "AlFe",
        "__Mn_Fe_": "MnFe",
    })[["source_id", "Teff", "logg", "FeH", "e_FeH", "alpha_proxy",
        "MgFe", "AlFe", "MnFe"]]
    gal_clean["survey"] = "GALAH"
    gal_clean["survey_rank"] = 1

    spec = pd.concat([apo_clean, gal_clean], ignore_index=True)
    spec["source_id"] = pd.to_numeric(spec["source_id"], errors="coerce").astype("Int64")
    spec = spec.dropna(subset=["source_id"]).copy()
    spec["source_id"] = spec["source_id"].astype(np.int64)
    spec = spec.merge(master_min, on="source_id", how="inner")

    spec["finite_alpha_rank"] = ~(
        np.isfinite(pd.to_numeric(spec["FeH"], errors="coerce"))
        & np.isfinite(pd.to_numeric(spec["alpha_proxy"], errors="coerce"))
    )
    spec_one = (spec.sort_values(["source_id", "finite_alpha_rank", "survey_rank"])
                    .drop_duplicates("source_id", keep="first"))
    return spec, spec_one


def main() -> int:
    master = decode_strings(Table.read(MASTER).to_pandas())
    src = pd.read_csv(SRC)

    keep_master = [
        "source_id", "tier", "P_vgrf_below_25", "vgrf_default",
        "source_in_old_preselection",
    ]
    master_min = master[keep_master].copy()
    master_min["source_id"] = pd.to_numeric(master_min["source_id"], errors="coerce").astype(np.int64)

    gsp_cols = [
        "source_id", "mh_gspphot", "mh_gspphot_lower", "mh_gspphot_upper",
        "teff_gspphot", "logg_gspphot", "azero_gspphot",
    ]
    gsp = src[[c for c in gsp_cols if c in src.columns]].copy()
    gsp["source_id"] = pd.to_numeric(gsp["source_id"], errors="coerce").astype(np.int64)
    gsp = master_min.merge(gsp, on="source_id", how="left")
    gsp = gsp.rename(columns={
        "mh_gspphot": "MH_gspphot",
        "mh_gspphot_lower": "MH_gspphot_lo",
        "mh_gspphot_upper": "MH_gspphot_hi",
    })
    Table.from_pandas(gsp).write(OUT / "expanded_gspphot_metallicity.fits", overwrite=True)

    spec, spec_one = load_spectroscopy(master_min)
    Table.from_pandas(spec).write(OUT / "expanded_spectroscopy_rows.fits", overwrite=True)
    Table.from_pandas(spec_one).write(OUT / "expanded_spectroscopy_unique.fits", overwrite=True)
    # Compatibility aliases for earlier review-bundle scripts.
    Table.from_pandas(spec).write(OUT / "expanded_cached_spectroscopy_rows.fits", overwrite=True)
    Table.from_pandas(spec_one).write(OUT / "expanded_cached_spectroscopy_unique.fits", overwrite=True)

    base = gsp.merge(
        spec_one[["source_id", "FeH", "e_FeH", "alpha_proxy", "MgFe", "AlFe", "MnFe", "survey"]],
        on="source_id", how="left",
    )
    base["best_FeH"] = base["FeH"].where(base["FeH"].notna(), base["MH_gspphot"])
    base["best_alpha"] = base["alpha_proxy"]
    base["chem_source"] = np.where(base["FeH"].notna(), base["survey"], "GSP-Phot")
    base["population_context"] = [
        classify(float(feh) if pd.notna(feh) else np.nan,
                 float(alpha) if pd.notna(alpha) else np.nan)
        for feh, alpha in zip(base["best_FeH"], base["best_alpha"])
    ]
    Table.from_pandas(base).write(OUT / "expanded_population_context.fits", overwrite=True)

    tier_abc = base["tier"].isin(["A", "B", "C"])
    tier_ab = base["tier"].isin(["A", "B"])
    alpha_mask = tier_abc & base["FeH"].notna() & base["alpha_proxy"].notna()
    pops = ["Splash", "GSE", "Aurora", "disk", "unclassified"]
    alpha_counts = {p: int(((base["population_context"] == p) & alpha_mask).sum()) for p in pops}
    alpha_medians = {
        p: finite_stat(base.loc[(base["population_context"] == p) & alpha_mask, "FeH"])["p50"]
        for p in pops
    }
    summary = {
        "n_expanded_master": int(len(base)),
        "tier_counts": {k: int((base["tier"] == k).sum()) for k in ["A", "B", "C", "D", "X"]},
        "tier_AB": int(tier_ab.sum()),
        "tier_ABC": int(tier_abc.sum()),
        "gsp_phot": {
            "tier_AB_finite_MH": int((tier_ab & base["MH_gspphot"].notna()).sum()),
            "tier_ABC_finite_MH": int((tier_abc & base["MH_gspphot"].notna()).sum()),
            "tier_AB_MH": finite_stat(base.loc[tier_ab, "MH_gspphot"]),
            "tier_ABC_MH": finite_stat(base.loc[tier_abc, "MH_gspphot"]),
        },
        "spectroscopy": {
            "rows_total_in_expanded_master": int(len(spec)),
            "unique_total_in_expanded_master": int(len(spec_one)),
            "tier_ABC_rows": int(spec["tier"].isin(["A", "B", "C"]).sum()),
            "tier_ABC_unique": int(spec_one["tier"].isin(["A", "B", "C"]).sum()),
            "tier_ABC_with_finite_alpha": int(alpha_mask.sum()),
            "tier_ABC_APOGEE_rows": int(((spec["survey"] == "APOGEE") & spec["tier"].isin(["A", "B", "C"])).sum()),
            "tier_ABC_GALAH_rows": int(((spec["survey"] == "GALAH") & spec["tier"].isin(["A", "B", "C"])).sum()),
            "note": "APOGEE/GALAH rows come from the full expanded Phase 14W exact-source-id crossmatch." if (APOGEE.exists() and GALAH.exists()) else "APOGEE/GALAH rows are cached legacy crossmatches; they are synchronized to expanded tiers but are not a complete expanded-catalogue crossmatch.",
        },
        "alpha_subsample": {
            "n": int(alpha_mask.sum()),
            "counts": alpha_counts,
            "fractions": {p: float(alpha_counts[p] / max(int(alpha_mask.sum()), 1)) for p in pops},
            "median_FeH": finite_stat(base.loc[alpha_mask, "FeH"])["p50"],
            "median_FeH_by_population": alpha_medians,
        },
    }
    (OUT / "expanded_chemistry_summary.json").write_text(json.dumps(summary, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    bins = np.linspace(-3.0, 0.8, 39)
    for tiers, label, color in [
        (["A", "B"], "Tier A+B", "tab:red"),
        (["C"], "Tier C", "tab:green"),
        (["X"], "Tier X", "0.65"),
    ]:
        vals = base.loc[base["tier"].isin(tiers), "MH_gspphot"].dropna().to_numpy(dtype=float)
        if len(vals):
            ax.hist(vals, bins=bins, histtype="step", density=True, linewidth=1.7,
                    label=f"{label} (n={len(vals)})", color=color)
    ax.set_xlabel("[M/H] (Gaia DR3 GSP-Phot)")
    ax.set_ylabel("density")
    ax.set_title("Expanded catalogue GSP-Phot metallicity context")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "expanded_gspphot_metallicity.png", dpi=180)
    fig.savefig(REPO / "release/figures/fig04_metallicity_histogram.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    pal = {"Splash": "tab:orange", "GSE": "tab:blue", "Aurora": "tab:purple",
           "disk": "tab:cyan", "unclassified": "0.55"}
    for pop, color in pal.items():
        m = alpha_mask & (base["population_context"] == pop)
        if int(m.sum()):
            ax.scatter(base.loc[m, "FeH"], base.loc[m, "alpha_proxy"], s=28,
                       color=color, edgecolor="white", linewidth=0.3,
                       label=f"{pop} (n={int(m.sum())})")
    ax.axvline(-1.0, ls=":", color="black", alpha=0.45)
    ax.axvline(-0.7, ls=":", color="black", alpha=0.45)
    ax.axhline(0.15, ls=":", color="black", alpha=0.45)
    ax.axhline(0.25, ls=":", color="black", alpha=0.45)
    ax.set_xlabel("[Fe/H]")
    ax.set_ylabel(r"[$\alpha$/Fe] proxy")
    ax.set_title("APOGEE/GALAH alpha subset in expanded Tier A+B+C")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "expanded_alpha_population_context.png", dpi=180)
    fig.savefig(REPO / "release/figures/fig14_alpha_fe_by_population.pdf")
    plt.close(fig)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
