"""Phase 14AF -- radial-velocity systematics and unresolved-multiplicity audit.

Referee response (deep-review Issue 9). A 25 km/s GRF-speed threshold is
unusually vulnerable to RV pathologies and unresolved binaries. This script
documents, for each tier:

  * the RV-quality distribution (sigma_RV, RVS transits, rvs_quality_ok,
    rv_chisq_pvalue as an RV-variability/binary indicator),
  * an explicit cross-match against the Gaia DR3 non-single-star catalogue
    (gaiadr3.nss_two_body_orbit) to flag known spectroscopic/astrometric/
    eclipsing binaries,
  * an external RV consistency check against APOGEE DR17 and GALAH DR3,
  * a top-N list of the slowest, highest-confidence sources with their
    RV-specific warning indicators.

Outputs
-------
phase14/expanded_rv_audit_summary.json
phase14/expanded_rv_audit_top.csv
tables/v15/tab_rv_quality.tex
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
MASTER = BUNDLE / "catalogues" / "catalogue_expanded_master.fits"
XMATCH = BUNDLE / "phase14" / "expanded_spectroscopic_crossmatch"
OUT = BUNDLE / "phase14"
TAB_DIRS = [REPO / "release" / "tables" / "v15", BUNDLE / "tables" / "v15"]

CHISQ_VAR = 0.01   # rv_chisq_pvalue below this => RV-variable (binary candidate)


def nss_crossmatch(source_ids: np.ndarray) -> pd.DataFrame:
    """Cross-match source_ids against gaiadr3.nss_two_body_orbit (chunked)."""
    from astroquery.gaia import Gaia
    import warnings
    warnings.filterwarnings("ignore")
    rows = []
    ids = [int(s) for s in source_ids]
    for i in range(0, len(ids), 500):
        chunk = ids[i:i + 500]
        in_list = ",".join(str(s) for s in chunk)
        q = (f"SELECT source_id, nss_solution_type FROM gaiadr3.nss_two_body_orbit "
             f"WHERE source_id IN ({in_list})")
        try:
            res = Gaia.launch_job_async(q).get_results().to_pandas()
            rows.append(res)
        except Exception as exc:  # noqa: BLE001
            print(f"  NSS chunk {i} failed: {exc}", flush=True)
    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(columns=["source_id", "nss_solution_type"])


def tier_rv_stats(s: pd.DataFrame) -> dict:
    return {
        "n": int(len(s)),
        "median_sigma_rv": round(float(np.nanmedian(s["radial_velocity_error"])), 2),
        "frac_sigma_rv_gt2_pct": round(float(np.mean(s["radial_velocity_error"] > 2.0) * 100), 1),
        "median_rv_nb_transits": int(np.nanmedian(s["rv_nb_transits"])),
        "frac_rvs_quality_ok_pct": round(float(np.mean(s["rvs_quality_ok"].astype(bool)) * 100), 1),
        "frac_ruwe_gt1p4_pct": round(float(np.mean(s["ruwe"] > 1.4) * 100), 1),
        "frac_rv_variable_pct": round(float(np.mean(s["rv_chisq_pvalue"] < CHISQ_VAR) * 100), 1),
    }


def external_rv_check(master: pd.DataFrame) -> dict:
    out = {}
    for survey, fname, rvcol in [("APOGEE", "apogee_dr17_expanded_crossmatch.csv", None),
                                 ("GALAH", "galah_dr3_expanded_crossmatch.csv", None)]:
        p = XMATCH / fname
        if not p.exists():
            continue
        x = pd.read_csv(p)
        # discover an external RV column
        cands = [c for c in x.columns if c.lower() in
                 ("vhelio_avg", "rv", "rv_galah", "radial_velocity", "vrad", "rv_value")]
        if "source_id" not in x.columns or not cands:
            out[survey] = {"n_overlap": int(len(x)), "note": f"no external RV column found among {list(x.columns)[:8]}"}
            continue
        rvcol = cands[0]
        j = x[["source_id", rvcol]].merge(
            master[["source_id", "radial_velocity", "tier"]], on="source_id", how="inner")
        j = j[j["tier"].isin(["A", "B", "C"]) & np.isfinite(j[rvcol]) & np.isfinite(j["radial_velocity"])]
        d = j["radial_velocity"] - j[rvcol]
        if len(d):
            out[survey] = {
                "rv_column": rvcol, "n_overlap_tierABC": int(len(d)),
                "median_delta_rv_kms": round(float(np.median(d)), 2),
                "robust_scatter_kms": round(float(1.4826 * np.median(np.abs(d - np.median(d)))), 2),
                "frac_abs_delta_gt5_pct": round(float(np.mean(np.abs(d) > 5) * 100), 1),
            }
    return out


def main() -> int:
    m = Table.read(MASTER).to_pandas()
    m["tier"] = m["tier"].str.decode("utf-8").str.strip() if m["tier"].dtype == object and isinstance(m["tier"].iloc[0], bytes) else m["tier"].astype(str).str.strip()
    abc = m[m["tier"].isin(["A", "B", "C"])].copy()

    print(f"[14AF] NSS cross-match for {len(abc)} Tier A+B+C sources ...", flush=True)
    nss = nss_crossmatch(abc["source_id"].to_numpy())
    nss_ids = set(nss["source_id"].astype("int64")) if len(nss) else set()
    abc["nss_match"] = abc["source_id"].astype("int64").isin(nss_ids)

    per_tier = {lab: tier_rv_stats(m[sel]) for lab, sel in [
        ("A", m["tier"] == "A"), ("B", m["tier"] == "B"),
        ("A+B", m["tier"].isin(["A", "B"])), ("A+B+C", m["tier"].isin(["A", "B", "C"]))]}
    for lab, sel in [("A", abc["tier"] == "A"), ("B", abc["tier"] == "B"),
                     ("A+B", abc["tier"].isin(["A", "B"])), ("A+B+C", abc["tier"].isin(["A", "B", "C"]))]:
        per_tier[lab]["n_nss_match"] = int(abc.loc[sel, "nss_match"].sum())

    # top-N slowest highest-confidence with RV warning indicators
    ab = abc[abc["tier"].isin(["A", "B"])].copy()
    ab = ab.sort_values("vgrf_default").head(50)
    top = ab[["source_id", "tier", "vgrf_default", "P_vgrf_below_25", "radial_velocity",
              "radial_velocity_error", "rv_nb_transits", "rv_chisq_pvalue", "ruwe", "nss_match"]]
    top.to_csv(OUT / "expanded_rv_audit_top.csv", index=False)

    summ = {
        "rv_variability_threshold_chisq_p": CHISQ_VAR,
        "per_tier": per_tier,
        "nss_total_matches_tierABC": int(abc["nss_match"].sum()),
        "nss_matches_tierAB": int(abc.loc[abc["tier"].isin(["A", "B"]), "nss_match"].sum()),
        "nss_solution_types": (nss["nss_solution_type"].value_counts().to_dict() if len(nss) else {}),
        "external_rv": external_rv_check(m),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "expanded_rv_audit_summary.json").write_text(json.dumps(summ, indent=2))
    _latex_table(per_tier)
    print(json.dumps(summ, indent=2))
    return 0


def _latex_table(per_tier: dict) -> None:
    order = ["A", "B", "A+B", "A+B+C"]
    lines = [
        r"% Auto-generated by scripts/phase14af_rv_audit.py",
        r"\begin{deluxetable}{lcccccc}",
        r"\tablecaption{Radial-velocity quality and unresolved-multiplicity audit by tier."
        r" $\sigma_{\rm RV}$ is the DR3 RV uncertainty; $f(\sigma_{\rm RV}{>}2)$ the fraction"
        r" above 2\kms; $f_{\rm var}$ the fraction with $p(\chi^2_{\rm RV}){<}0.01$ (an RV-variability"
        r" / binary indicator); $f_{\rm RUWE>1.4}$ the astrometric-binary fraction; and $N_{\rm NSS}$"
        r" the number matched in the Gaia DR3 non-single-star catalogue"
        r" (\texttt{nss\_two\_body\_orbit}).\label{tab:rv_quality}}",
        r"\tablehead{\colhead{tier} & \colhead{$\tilde\sigma_{\rm RV}$} & \colhead{$f(\sigma_{\rm RV}{>}2)$}"
        r" & \colhead{$\tilde{N}_{\rm tr}$} & \colhead{$f_{\rm var}$} & \colhead{$f_{\rm RUWE>1.4}$}"
        r" & \colhead{$N_{\rm NSS}$}}",
        r"\startdata",
    ]
    for lab in order:
        d = per_tier[lab]
        lines.append(f"{lab} & {d['median_sigma_rv']:.2f} & {d['frac_sigma_rv_gt2_pct']:.0f}\\% & "
                     f"{d['median_rv_nb_transits']} & {d['frac_rv_variable_pct']:.1f}\\% & "
                     f"{d['frac_ruwe_gt1p4_pct']:.1f}\\% & {d.get('n_nss_match', 0)} \\\\")
    lines += [r"\enddata", r"\end{deluxetable}"]
    text = "\n".join(lines) + "\n"
    for dd in TAB_DIRS:
        dd.mkdir(parents=True, exist_ok=True)
        (dd / "tab_rv_quality.tex").write_text(text)


if __name__ == "__main__":
    raise SystemExit(main())
