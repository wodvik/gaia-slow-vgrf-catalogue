"""Phase 14AE -- expanded distance-covariance stress test.

Referee response (deep-review Issue 6). The catalogue Monte Carlo samples
the Gaia astrometric covariance but draws the Bailer-Jones photogeometric
distance independently of the parallax. The manuscript quantified this
approximation with a 200-star transition-band test. This script expands
that test to the ENTIRE tier-sensitive band (adopted P>=0.30, ~3,200
sources) and reports:

  * the distribution of |Delta P| under a maximal parallax-distance
    Gaussian-copula coupling (z_distance = -z_parallax),
  * a tier-label confusion matrix before vs after coupled sampling,
  * the stability of the primary Tier A+B membership.

It recomputes BOTH the independent and the coupled probability fresh at
the adopted Vc=229 frame (so Delta P isolates the coupling, not the
realisation count), using the bundle's adopted candidate pool
(private_inputs/expanded_candidates_mc_tiered.csv; tiers 289/541/1952).

Outputs
-------
phase14/expanded_covariance_stress_per_star.csv
phase14/expanded_covariance_stress_summary.json
tables/v15/tab_covariance_stress.tex
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
CAND = BUNDLE / "private_inputs" / "expanded_candidates_mc_tiered.csv"
OUT = BUNDLE / "phase14"
TAB_DIRS = [REPO / "release" / "tables" / "v15", BUNDLE / "tables" / "v15"]
FIG_DIRS = [REPO / "release" / "figures", BUNDLE / "figures"]

VGRF_CUT = 25.0
TIERS = {"A": 0.95, "B": 0.84, "C": 0.50}
P_FLOOR = 0.30          # tier-sensitive band: everything that can reach Tier C
N_SAMP = 2000
SEED = 20260605
SOLAR = {"R0_kpc": 8.178, "z_sun_pc": 25.0, "Vc_kms": 229.0,
         "U_kms": 11.1, "V_kms": 12.24, "W_kms": 7.25}


def galcen_frame() -> coord.Galactocentric:
    return coord.Galactocentric(
        galcen_distance=SOLAR["R0_kpc"] * u.kpc, z_sun=SOLAR["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            SOLAR["U_kms"] * u.km / u.s,
            (SOLAR["Vc_kms"] + SOLAR["V_kms"]) * u.km / u.s,
            SOLAR["W_kms"] * u.km / u.s))


def cov3(df: pd.DataFrame) -> np.ndarray:
    sig = np.stack([df["parallax_error"].to_numpy(float),
                    df["pmra_error"].to_numpy(float),
                    df["pmdec_error"].to_numpy(float)], axis=-1)
    n = len(df)
    rho = np.zeros((n, 3, 3))
    for i in range(3):
        rho[:, i, i] = 1.0
    rho[:, 0, 1] = rho[:, 1, 0] = df["parallax_pmra_corr"].to_numpy(float)
    rho[:, 0, 2] = rho[:, 2, 0] = df["parallax_pmdec_corr"].to_numpy(float)
    rho[:, 1, 2] = rho[:, 2, 1] = df["pmra_pmdec_corr"].to_numpy(float)
    return rho * sig[:, :, None] * sig[:, None, :]


def cholesky_stack(cov: np.ndarray) -> np.ndarray:
    out = np.empty_like(cov)
    for i in range(cov.shape[0]):
        try:
            out[i] = np.linalg.cholesky(cov[i])
        except np.linalg.LinAlgError:
            out[i] = np.linalg.cholesky(cov[i] + 1e-9 * np.eye(3))
    return out


def split_normal_from_z(med, lo, hi, z):
    sig_lo = np.maximum(med - lo, 1.0)
    sig_hi = np.maximum(hi - med, 1.0)
    sig = np.where(z < 0, sig_lo[:, None], sig_hi[:, None])
    return np.maximum(med[:, None] + z * sig, 1.0)


def p_vgrf(df, mode, rng, n_samp):
    frame = galcen_frame()
    pvals = np.empty(len(df))
    chunk = 120
    for start in range(0, len(df), chunk):
        sub = df.iloc[start:start + chunk].reset_index(drop=True)
        n = len(sub)
        mu = np.stack([sub["parallax_zpcorr"].to_numpy(float),
                       sub["pmra"].to_numpy(float),
                       sub["pmdec"].to_numpy(float)], axis=-1)
        L = cholesky_stack(cov3(sub))
        z_ast = rng.standard_normal((n, n_samp, 3))
        ast = mu[:, None, :] + np.einsum("nij,nsj->nsi", L, z_ast)
        pmra_s, pmdec_s = ast[:, :, 1], ast[:, :, 2]
        if mode == "independent":
            z_dist = rng.standard_normal((n, n_samp))
        else:  # parallax_anticorrelated copula
            z_dist = -z_ast[:, :, 0]
        dist_s = split_normal_from_z(sub["dist_pc_final_screen"].to_numpy(float),
                                     sub["dist_lo_pc_final_screen"].to_numpy(float),
                                     sub["dist_hi_pc_final_screen"].to_numpy(float), z_dist)
        rv_s = (sub["radial_velocity"].to_numpy(float)[:, None]
                + rng.standard_normal((n, n_samp)) * sub["radial_velocity_error"].to_numpy(float)[:, None])
        ra = np.broadcast_to(sub["ra"].to_numpy(float)[:, None], (n, n_samp)).ravel()
        dec = np.broadcast_to(sub["dec"].to_numpy(float)[:, None], (n, n_samp)).ravel()
        icrs = coord.SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=dist_s.ravel() * u.pc,
                              pm_ra_cosdec=pmra_s.ravel() * u.mas / u.yr, pm_dec=pmdec_s.ravel() * u.mas / u.yr,
                              radial_velocity=rv_s.ravel() * u.km / u.s, frame="icrs")
        g = icrs.transform_to(frame)
        vgrf = np.sqrt(g.v_x.to_value(u.km / u.s) ** 2 + g.v_y.to_value(u.km / u.s) ** 2
                       + g.v_z.to_value(u.km / u.s) ** 2).reshape(n, n_samp)
        pvals[start:start + n] = np.mean(vgrf < VGRF_CUT, axis=1)
    return pvals


def tier_cat(p: np.ndarray) -> np.ndarray:
    out = np.full(len(p), "subC", dtype=object)
    out[p > TIERS["C"]] = "C"
    out[p > TIERS["B"]] = "B"
    out[p > TIERS["A"]] = "A"
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-samp", type=int, default=N_SAMP)
    args = ap.parse_args()

    cols = ["source_id", "ra", "dec", "parallax_zpcorr", "pmra", "pmdec",
            "parallax_error", "pmra_error", "pmdec_error", "parallax_pmra_corr",
            "parallax_pmdec_corr", "pmra_pmdec_corr", "radial_velocity",
            "radial_velocity_error", "dist_pc_final_screen", "dist_lo_pc_final_screen",
            "dist_hi_pc_final_screen", "P_vgrf_below_25", "tier"]
    df = pd.read_csv(CAND, usecols=cols)
    df = df[df["P_vgrf_below_25"] >= P_FLOOR].sort_values("source_id").reset_index(drop=True)
    n = len(df)
    print(f"[14AE] tier-sensitive band (adopted P>={P_FLOOR}): {n} sources; N_SAMP={args.n_samp}")

    p_ind = p_vgrf(df, "independent", np.random.default_rng(SEED), args.n_samp)
    p_cop = p_vgrf(df, "copula", np.random.default_rng(SEED), args.n_samp)
    dP = p_cop - p_ind
    t_ind, t_cop = tier_cat(p_ind), tier_cat(p_cop)

    cats = ["A", "B", "C", "subC"]
    confusion = {ri: {ci: int(np.sum((t_ind == ri) & (t_cop == ci))) for ci in cats} for ri in cats}

    # Primary-catalogue stability: adopted Tier A+B members (P_adopted>0.84) that fall
    # below 0.84 under the coupled draw.
    adopted_AB = df["P_vgrf_below_25"].to_numpy(float) > TIERS["B"]
    coupled_AB_lost = int(np.sum(adopted_AB & (p_cop <= TIERS["B"])))
    # ... isolating coupling: independent-fresh A+B that the copula flips
    indep_AB = p_ind > TIERS["B"]
    flip_AB = int(np.sum(indep_AB != (p_cop > TIERS["B"])))

    pd.DataFrame({
        "source_id": df["source_id"].astype("int64"),
        "adopted_P": df["P_vgrf_below_25"].astype(float), "adopted_tier": df["tier"].astype(str),
        "P_independent": p_ind, "P_copula": p_cop, "delta_P": dP,
        "tier_independent": t_ind, "tier_copula": t_cop,
    }).to_csv(OUT / "expanded_covariance_stress_per_star.csv", index=False)
    _figure(p_ind, p_cop, t_ind, t_cop)

    summ = {
        "n_sources_band": n, "p_floor": P_FLOOR, "n_samp": args.n_samp,
        "coupling": "Gaussian copula z_distance=-z_parallax (maximal anti-correlation)",
        "abs_delta_P": {
            "median": float(np.median(np.abs(dP))), "p84": float(np.percentile(np.abs(dP), 84)),
            "p95": float(np.percentile(np.abs(dP), 95)), "max": float(np.max(np.abs(dP))),
        },
        "tier_switches_independent_vs_copula": int(np.sum(t_ind != t_cop)),
        "confusion_matrix_indep_rows_copula_cols": confusion,
        "headline_AB_adopted": int(np.sum(adopted_AB)),
        "headline_AB_lost_under_copula": coupled_AB_lost,
        "AB_boundary_flips_indep_vs_copula": flip_AB,
        "crossings": {
            "p_gt_0p95": int(np.sum((p_ind > 0.95) != (p_cop > 0.95))),
            "p_gt_0p84": int(np.sum((p_ind > 0.84) != (p_cop > 0.84))),
            "p_gt_0p50": int(np.sum((p_ind > 0.50) != (p_cop > 0.50))),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "expanded_covariance_stress_summary.json").write_text(json.dumps(summ, indent=2))
    _latex_table(summ)
    print(json.dumps(summ, indent=2))
    return 0


def _figure(p_ind, p_cop, t_ind, t_cop) -> None:
    sw = t_ind != t_cop
    dP = p_cop - p_ind
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.scatter(p_ind[~sw], dP[~sw], s=5, c="#9aa0a6", alpha=0.5, edgecolors="none",
               label="tier unchanged", rasterized=True)
    ax.scatter(p_ind[sw], dP[sw], s=10, c="#d62728", alpha=0.8, edgecolors="none",
               label=f"tier changed (N={int(sw.sum())})")
    for x in (0.50, 0.84, 0.95):
        ax.axvline(x, color="k", ls=":", lw=0.7)
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlabel(r"$P(V_{\rm grf}<25\,{\rm km\,s^{-1}})$ (independent draw)")
    ax.set_ylabel(r"$\Delta P$ (copula $-$ independent)")
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    for d in FIG_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / "fig_joint_distance_test.pdf", bbox_inches="tight", dpi=250)
    plt.close(fig)


def _latex_table(summ: dict) -> None:
    cats = ["A", "B", "C", "subC"]
    disp = {"A": "A", "B": "B", "C": "C", "subC": r"$<$C"}
    cm = summ["confusion_matrix_indep_rows_copula_cols"]
    lines = [
        r"% Auto-generated by scripts/phase14ae_covariance_stress.py",
        r"\begin{deluxetable}{lcccc}",
        rf"\tablecaption{{Tier-label confusion matrix for the full tier-sensitive band"
        rf" ($N={summ['n_sources_band']:,}$ sources with adopted $P\geq{summ['p_floor']}$)"
        r" under a maximal parallax--distance Gaussian-copula coupling, relative to an"
        r" independent-distance recomputation at the same realisations. Rows are the"
        r" independent-draw tier, columns the coupled-draw tier; off-diagonal entries are"
        rf" tier switches. Only {summ['tier_switches_independent_vs_copula']} of"
        rf" {summ['n_sources_band']:,} sources switch. Relative to the adopted"
        rf" primary catalogue, {summ['headline_AB_adopted'] - summ['headline_AB_lost_under_copula']}"
        rf" of {summ['headline_AB_adopted']} Tier~A+B members remain above"
        r" $P=0.84$ under the coupled draw (median $|\Delta P|="
        rf"{summ['abs_delta_P']['median']:.4f}$, max ${summ['abs_delta_P']['max']:.4f}$)."
        r"\label{tab:covariance_stress}}",
        r"\tablehead{\colhead{indep.$\backslash$copula} & \colhead{A} & \colhead{B}"
        r" & \colhead{C} & \colhead{$<$C}}",
        r"\startdata",
    ]
    for ri in cats:
        cells = " & ".join(str(cm[ri][ci]) for ci in cats)
        lines.append(f"{disp[ri]} & {cells} " + r"\\")
    lines += [r"\enddata", r"\end{deluxetable}"]
    text = "\n".join(lines) + "\n"
    for d in TAB_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        (d / "tab_covariance_stress.tex").write_text(text)


if __name__ == "__main__":
    raise SystemExit(main())
