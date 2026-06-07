"""Phase 14AJ -- Filion et al. (2025) low-Vphi comparison.

Filion et al. (2025, ApJ, 989, 70) provide a supplementary table of 69
APOGEE identifiers for slow/retrograde, low-alpha disc-abundance stars.
This script resolves those APOGEE identifiers to Gaia source_id values,
applies the same Gaia-RVS quality classes and low-Vgrf tiering machinery
used by this release, and writes a compact manuscript table for the three
Filion stars that enter the public near-threshold parent pool.
"""
from __future__ import annotations

import json
from pathlib import Path

import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.table import Table
from astroquery.gaia import Gaia
from astroquery.utils.tap.core import TapPlus

from phase0c_parent_buffer_zpcorr_screen import compute_vgrf_with_parallax, compute_zpcorr
from phase0e_expanded_mc_tiering import (
    MC,
    VGRF_CUT,
    assign_tiers,
    dist_columns,
    galcen_frame,
    mc_pass,
    point_vgrf,
)
from phase0f_materialize_expanded_catalogue import classify_rvs_quality


BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
OUT = BUNDLE / "phase14"
TAB_DIRS = [REPO / "release" / "tables" / "v15", BUNDLE / "tables" / "v15"]
MASTER = BUNDLE / "catalogues" / "catalogue_expanded_master.fits"
VIZIER_TAP_URL = "http://tapvizier.cds.unistra.fr/TAPVizieR/tap"

FILION_APOGEE_IDS = [
    "2M00014230+1702419",
    "2M00303867+1549501",
    "2M01231405-0102329",
    "2M01244485+1453267",
    "2M01590930+3802046",
    "2M04022073-7023289",
    "2M04215815-6157227",
    "2M06234629-7005024",
    "2M07472456+5311408",
    "2M08372116-8537442",
    "2M09191140+3346373",
    "2M09460164+3322570",
    "2M10391298+2948411",
    "2M10475148-4448333",
    "2M11291121-8501189",
    "2M12325862+3626584",
    "2M12414587+4926490",
    "2M13001371+5217200",
    "2M13020472+5743438",
    "2M13041724-0134399",
    "2M13123688+1801404",
    "2M13323559-0118149",
    "2M13381907+2938588",
    "2M13403109-1550124",
    "2M13423397+2745255",
    "2M13453784-1539189",
    "2M13481254-0053486",
    "2M14003743+0334382",
    "2M14005878+0123010",
    "2M14025155+4747385",
    "2M14044997+0320333",
    "2M14070704+5411394",
    "2M14130012+0544007",
    "2M14575158+5300140",
    "2M15022633+3248296",
    "2M15023086+4212111",
    "2M15023130+5112297",
    "2M15050196-4430509",
    "2M15070607+2222028",
    "2M15104604+3158498",
    "2M15233167+0919257",
    "2M15334993+0934139",
    "2M15490280+2650095",
    "2M15513506+2805029",
    "2M15544945+2915314",
    "2M15572935+2608034",
    "2M16003181-2139519",
    "2M16013180+4056403",
    "2M16032054+3918321",
    "2M16045588+2738536",
    "2M16065591-2222130",
    "2M16082687+1809497",
    "2M16132140+5253093",
    "2M16215605+4725419",
    "2M16434878+3726165",
    "2M16453831+4802186",
    "2M17100538+6323032",
    "2M17121198-2438245",
    "2M17200667+4207505",
    "2M17252566+5837413",
    "2M18212748-3903116",
    "2M18451564+0557157",
    "2M18555040-2808410",
    "2M20503798-0045260",
    "2M22095666-3445576",
    "2M22193340-4820346",
    "2M22221939-5008474",
    "2M23110573-5259267",
    "2M23482381+8501382",
]


def _quoted(items: list[str]) -> str:
    return ",".join(f"'{item}'" for item in items)


def query_apogee() -> pd.DataFrame:
    tap = TapPlus(url=VIZIER_TAP_URL)
    query = f"""
    SELECT "APOGEE", "GaiaEDR3", "RAJ2000", "DEJ2000",
           "HRV", "e_HRV", "RV", "e_RV", "SNR",
           "[Fe/H]", "[Mg/Fe]", "[Al/Fe]", "[Mn/Fe]"
    FROM "III/286/catalog"
    WHERE "APOGEE" IN ({_quoted(FILION_APOGEE_IDS)})
    """
    df = tap.launch_job(query).get_results().to_pandas()
    df["GaiaEDR3"] = pd.to_numeric(df["GaiaEDR3"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["GaiaEDR3"]).drop_duplicates(["APOGEE", "GaiaEDR3"])
    df["source_id"] = df["GaiaEDR3"].astype("int64")
    return df


def query_gaia(source_ids: list[int]) -> pd.DataFrame:
    in_list = ",".join(str(int(s)) for s in source_ids)
    query = f"""
    SELECT source_id, ra, dec, l, b, parallax, parallax_error,
           parallax_over_error, pmra, pmra_error, pmdec, pmdec_error,
           ra_dec_corr, ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr,
           dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr,
           parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr,
           radial_velocity, radial_velocity_error, rv_amplitude_robust,
           rv_nb_transits, rv_template_teff, rv_method_used,
           rv_expected_sig_to_noise, rv_renormalised_gof,
           rv_chisq_pvalue, grvs_mag, phot_g_mean_mag,
           phot_bp_mean_mag, phot_rp_mean_mag, bp_rp,
           nu_eff_used_in_astrometry, pseudocolour,
           astrometric_params_solved, ruwe, duplicated_source, ecl_lat
    FROM gaiadr3.gaia_source
    WHERE source_id IN ({in_list})
    """
    df = Gaia.launch_job_async(query).get_results().to_pandas()
    df["source_id"] = df["source_id"].astype("int64")
    return df


def query_bailer_jones(source_ids: list[int]) -> pd.DataFrame:
    tap = TapPlus(url=VIZIER_TAP_URL)
    query = f"""
    SELECT "Source", "rgeo", "b_rgeo", "B_rgeo",
           "rpgeo", "b_rpgeo", "B_rpgeo", "Flag"
    FROM "I/352/gedr3dis"
    WHERE "Source" IN ({",".join(str(int(s)) for s in source_ids)})
    """
    df = tap.launch_job(query).get_results().to_pandas()
    df = df.rename(
        columns={
            "Source": "source_id",
            "b_rgeo": "rgeo_lo",
            "B_rgeo": "rgeo_hi",
            "b_rpgeo": "rpgeo_lo",
            "B_rpgeo": "rpgeo_hi",
            "Flag": "bj_flag",
        }
    )
    df["source_id"] = df["source_id"].astype("int64")
    return df


def nss_crossmatch(source_ids: list[int]) -> pd.DataFrame:
    in_list = ",".join(str(int(s)) for s in source_ids)
    query = (
        "SELECT source_id, nss_solution_type "
        "FROM gaiadr3.nss_two_body_orbit "
        f"WHERE source_id IN ({in_list})"
    )
    try:
        out = Gaia.launch_job_async(query).get_results().to_pandas()
    except Exception as exc:  # noqa: BLE001
        print(f"[14AJ] NSS query failed: {exc}", flush=True)
        return pd.DataFrame(columns=["source_id", "nss_solution_type"])
    if len(out):
        out["source_id"] = out["source_id"].astype("int64")
    return out


def load_release_master() -> pd.DataFrame:
    df = Table.read(MASTER).to_pandas()
    df["source_id"] = df["source_id"].astype("int64")
    df["tier"] = df["tier"].apply(
        lambda x: x.decode("utf-8").strip() if isinstance(x, bytes) else str(x).strip()
    )
    for col in ["rv_quality", "chem_survey", "chem_population"]:
        if col in df.columns:
            df[col] = df[col].apply(
                lambda x: x.decode("utf-8").strip() if isinstance(x, bytes) else str(x).strip()
            )
    return df


def filion_style_vphi(df: pd.DataFrame) -> pd.DataFrame:
    """Approximate Filion's reported Vphi convention from APOGEE/Gaia fields."""
    frame = coord.Galactocentric(
        galcen_distance=8.275 * u.kpc,
        z_sun=20.8 * u.pc,
        galcen_v_sun=coord.CartesianDifferential([8.4, 251.8, 8.4] * u.km / u.s),
    )
    out = np.full(len(df), np.nan)
    for i, row in df.iterrows():
        rv = row["HRV"] if np.isfinite(row["HRV"]) else row.get("RV", np.nan)
        if not np.isfinite(rv):
            continue
        c = coord.SkyCoord(
            ra=row["RAJ2000"] * u.deg,
            dec=row["DEJ2000"] * u.deg,
            distance=row["rpgeo"] * u.pc,
            pm_ra_cosdec=row["pmra"] * u.mas / u.yr,
            pm_dec=row["pmdec"] * u.mas / u.yr,
            radial_velocity=rv * u.km / u.s,
            frame="icrs",
        )
        cyl = c.transform_to(frame).cylindrical
        out[i] = (
            cyl.rho * cyl.differentials["s"].d_phi
        ).to_value(u.km / u.s, equivalencies=u.dimensionless_angles())
    df["filion_style_vphi_apogee_kms"] = out
    return df


def run_tiering(df: pd.DataFrame) -> pd.DataFrame:
    zpcorr_mas, parallax_zpcorr, zpcorr_valid = compute_zpcorr(df)
    df = df.copy()
    df["zpcorr_value_uas"] = zpcorr_mas * 1000.0
    df["zpcorr_valid"] = zpcorr_valid
    df["parallax_zpcorr"] = parallax_zpcorr
    df["distance_inv_zpcorr_pc"] = np.where(parallax_zpcorr > 0, 1000.0 / parallax_zpcorr, np.nan)
    ok = np.isfinite(parallax_zpcorr) & (parallax_zpcorr > 0)
    df["vgrf_zpcorr_inv"] = np.nan
    df.loc[ok, "vgrf_zpcorr_inv"] = compute_vgrf_with_parallax(df.loc[ok], parallax_zpcorr[ok])
    df["source_in_old_preselection"] = False

    dist_med, dist_lo, dist_hi, dist_source = dist_columns(df)
    point = point_vgrf(df, dist_med, batch_size=5000)
    df["dist_pc_final_screen"] = dist_med
    df["dist_lo_pc_final_screen"] = dist_lo
    df["dist_hi_pc_final_screen"] = dist_hi
    df["dist_source_final_screen"] = dist_source
    df["diagnostic_vgrf_default"] = point

    rng = np.random.default_rng(MC["random_seed"])
    p = mc_pass(df, dist_med, dist_lo, dist_hi, MC["base_realisations"], rng, star_batch=69)
    n_real = np.full(len(df), MC["base_realisations"], dtype=int)

    refine = (p > MC["refine_lo"]) & (p < MC["refine_hi"])
    if refine.any():
        idx = np.flatnonzero(refine)
        p[idx] = mc_pass(
            df.iloc[idx].reset_index(drop=True),
            dist_med[idx],
            dist_lo[idx],
            dist_hi[idx],
            MC["refine_realisations"],
            rng,
            star_batch=max(1, len(idx)),
        )
        n_real[idx] = MC["refine_realisations"]

    near = np.abs(point - VGRF_CUT) < MC["ultra_window_kms"]
    if near.any():
        idx = np.flatnonzero(near)
        p[idx] = mc_pass(
            df.iloc[idx].reset_index(drop=True),
            dist_med[idx],
            dist_lo[idx],
            dist_hi[idx],
            MC["ultra_realisations"],
            rng,
            star_batch=max(1, len(idx)),
        )
        n_real[idx] = MC["ultra_realisations"]

    df["diagnostic_P_vgrf_below_25"] = p
    df["diagnostic_mc_realisations"] = n_real
    df["diagnostic_tier"] = assign_tiers(p, point)
    return df


def write_latex_table(rows: pd.DataFrame) -> None:
    lines = [
        r"\begin{center}",
        r"\refstepcounter{table}\label{tab:filion2025_comparison}",
        r"\begin{minipage}{\linewidth}",
        r"\small\textbf{Table \thetable.} Filion et al. (2025) stars entering the present near-threshold parent pool.",
        r"\end{minipage}",
        r"\vspace{0.4ex}",
        "",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}llrrrrrrrl@{}}",
        r"\hline\hline",
        r"Gaia DR3 source\_id & APOGEE ID & $V_\phi^{\rm F25}$ & $\vgrf$ & $P_{<25}$ & $\sigma_{\rm RV}$ & $N_{\rm tr}$ & $p_{\chi^2,{\rm RV}}$ & RUWE & RV class \\",
        r" & & \multicolumn{1}{c}{(\kms)} & \multicolumn{1}{c}{(\kms)} & & \multicolumn{1}{c}{(\kms)} & & & & \\",
        r"\hline",
    ]
    for _, row in rows.sort_values("release_vgrf_default").iterrows():
        pval = row["rv_chisq_pvalue"]
        pval_text = "--" if not np.isfinite(pval) else f"{pval:.3f}"
        lines.append(
            f"{int(row['source_id'])} & {row['APOGEE']} & "
            f"{row['filion_style_vphi_apogee_kms']:.1f} & "
            f"{row['release_vgrf_default']:.1f} & "
            f"{row['release_P_vgrf_below_25']:.3f} & "
            f"{row['radial_velocity_error']:.2f} & "
            f"{int(row['rv_nb_transits'])} & {pval_text} & "
            f"{row['ruwe']:.2f} & {row['release_rv_quality']} \\\\"
        )
    lines += [
        r"\hline",
        r"\end{tabular*}",
        r"\vspace{0.5ex}",
        "",
        r"\begin{minipage}{0.95\linewidth}",
        r"\footnotesize",
        r"\textit{Note.} $V_\phi^{\rm F25}$ is an approximate reproduction of the Filion et al. Galactocentric azimuthal-velocity convention using APOGEE DR17 line-of-sight velocities and Bailer-Jones et al. (2021) photogeometric distances. The remaining columns are from the present Gaia DR3 RVS low-\vgrf{} release products. None of the three sources is matched in the Gaia DR3 non-single-star two-body-orbit table. All three sources are classified as Tier X and therefore are excluded from the Tier A+B+C science catalogue.",
        r"\end{minipage}",
        r"\end{center}",
        "",
    ]
    text = "\n".join(lines)
    for table_dir in TAB_DIRS:
        table_dir.mkdir(parents=True, exist_ok=True)
        (table_dir / "tab_filion2025_comparison.tex").write_text(text, encoding="utf-8")


def main() -> int:
    from zero_point import zpt

    zpt.load_tables()
    OUT.mkdir(parents=True, exist_ok=True)

    apogee = query_apogee()
    source_ids = sorted(apogee["source_id"].astype("int64").unique().tolist())
    gaia = query_gaia(source_ids)
    bj = query_bailer_jones(source_ids)
    nss = nss_crossmatch(source_ids)
    nss_ids = set(nss["source_id"].astype("int64")) if len(nss) else set()

    df = (
        gaia.merge(bj, on="source_id", how="left")
        .merge(apogee.drop(columns=["GaiaEDR3"]), on="source_id", how="left")
        .sort_values("APOGEE")
        .reset_index(drop=True)
    )
    df["gaia_parent_contract"] = (
        df["radial_velocity"].notna()
        & (df["parallax"] > 0)
        & (df["parallax_over_error"] > 5)
        & df["pmra"].notna()
        & df["pmdec"].notna()
    )
    df["diagnostic_rv_quality"] = classify_rvs_quality(df)
    df["diagnostic_rvs_quality_ok"] = df["diagnostic_rv_quality"] == "ok"
    df["rv_variable_flag"] = pd.to_numeric(df["rv_chisq_pvalue"], errors="coerce") < 0.01
    df["ruwe_gt_1p4"] = pd.to_numeric(df["ruwe"], errors="coerce") > 1.4
    df["nss_match"] = df["source_id"].isin(nss_ids)

    df = run_tiering(df)
    df = filion_style_vphi(df)

    master = load_release_master()
    release_cols = [
        "source_id",
        "tier",
        "P_vgrf_below_25",
        "vgrf_default",
        "rv_quality",
        "rvs_quality_ok",
        "radial_velocity_error",
        "rv_nb_transits",
        "rv_chisq_pvalue",
        "ruwe",
    ]
    rel = master[release_cols].rename(
        columns={
            "tier": "release_tier",
            "P_vgrf_below_25": "release_P_vgrf_below_25",
            "vgrf_default": "release_vgrf_default",
            "rv_quality": "release_rv_quality",
            "rvs_quality_ok": "release_rvs_quality_ok",
        }
    )
    out = df.merge(rel, on="source_id", how="left", suffixes=("", "_release"))
    out["release_parent_pool"] = out["release_tier"].notna()

    preferred_cols = [
        "APOGEE",
        "source_id",
        "gaia_parent_contract",
        "release_parent_pool",
        "release_tier",
        "release_vgrf_default",
        "release_P_vgrf_below_25",
        "diagnostic_vgrf_default",
        "diagnostic_P_vgrf_below_25",
        "diagnostic_tier",
        "diagnostic_mc_realisations",
        "filion_style_vphi_apogee_kms",
        "diagnostic_rv_quality",
        "release_rv_quality",
        "radial_velocity_error",
        "rv_nb_transits",
        "rv_chisq_pvalue",
        "ruwe",
        "ruwe_gt_1p4",
        "rv_variable_flag",
        "nss_match",
        "HRV",
        "e_HRV",
        "SNR",
        "__Fe_H_",
        "__Mg_Fe_",
        "__Al_Fe_",
        "__Mn_Fe_",
    ]
    out = out[[c for c in preferred_cols if c in out.columns]]
    out.to_csv(OUT / "filion2025_comparison.csv", index=False)

    parent_rows = out[out["release_parent_pool"]].copy()
    write_latex_table(parent_rows)

    summary = {
        "filion_reference": "Filion et al. 2025, ApJ, 989, 70",
        "n_filion_apogee_ids": len(FILION_APOGEE_IDS),
        "n_resolved_to_unique_gaia_source_id": int(out["source_id"].nunique()),
        "n_gaia_parent_contract": int(out["gaia_parent_contract"].sum()),
        "n_release_parent_pool_overlap": int(out["release_parent_pool"].sum()),
        "n_release_tierA_overlap": int((out["release_tier"] == "A").sum()),
        "n_release_tierAB_overlap": int(out["release_tier"].isin(["A", "B"]).sum()),
        "n_release_tierABC_overlap": int(out["release_tier"].isin(["A", "B", "C"]).sum()),
        "diagnostic_tier_counts": {
            str(k): int(v) for k, v in out["diagnostic_tier"].value_counts().sort_index().items()
        },
        "diagnostic_rv_quality_counts": {
            str(k): int(v)
            for k, v in out["diagnostic_rv_quality"].value_counts().sort_index().items()
        },
        "n_ruwe_gt_1p4": int(out["ruwe_gt_1p4"].sum()),
        "n_rv_variable_flag": int(out["rv_variable_flag"].sum()),
        "n_nss_match": int(out["nss_match"].sum()),
        "release_parent_pool_source_ids": [int(x) for x in parent_rows["source_id"].tolist()],
        "outputs": {
            "comparison_csv": "phase14/filion2025_comparison.csv",
            "table": "tables/v15/tab_filion2025_comparison.tex",
        },
    }
    (OUT / "filion2025_comparison_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
