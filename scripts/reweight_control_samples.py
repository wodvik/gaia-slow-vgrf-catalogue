"""Expanded matched-control reweighting and comparison summary.

Uses the already-integrated control orbit products, but recomputes
observability weights against the expanded Tier A+B+C slow catalogue.
This avoids mixing the legacy 632-star slow marginal with the expanded
1,835-star catalogue.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table


REPO = Path(__file__).resolve().parents[1]
MASTER = REPO / "catalogues/catalogue_expanded_master.fits"
ORBITS = REPO / "catalogues/catalogue_expanded_orbits_tierABC.fits"
CONTROLS = REPO / "phase5/control_orbits.fits"
OUT = REPO / "analysis_products"
FIG = REPO / "figures"

BANDS = [
    ("slow", "Slow <25", "#d62728"),
    ("vgrf_25_50", "25-50", "#ff7f0e"),
    ("vgrf_50_100", "50-100", "#2ca02c"),
    ("vgrf_100_200", "100-200", "#1f77b4"),
    ("vgrf_200_260", "200-260", "#6a3d9a"),
]


def clean_strings(s: pd.Series) -> pd.Series:
    if s.dtype == object and len(s) and isinstance(s.iloc[0], (bytes, bytearray)):
        return s.str.decode("utf-8").str.strip()
    return s.astype(str).str.strip()


def weighted_quantile(values: np.ndarray, weights: np.ndarray, qs: list[float]) -> list[float]:
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[mask]
    weights = weights[mask]
    if len(values) == 0:
        return [float("nan") for _ in qs]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    return [float(np.interp(q, cdf, values)) for q in qs]


def kernel_weights(ctrl: pd.DataFrame, slow: pd.DataFrame) -> np.ndarray:
    # Product of 1D Gaussian-kernel density ratios against expanded slow
    # observability marginals. Normalize within each band to mean 1.
    slow_G = slow["phot_g_mean_mag"].to_numpy(float)
    slow_logd = np.log10(np.maximum(slow["dist_pc"].to_numpy(float), 1.0))
    slow_sb = np.abs(np.sin(np.deg2rad(slow["b"].to_numpy(float))))
    G = ctrl["G"].to_numpy(float)
    logd = np.log10(np.maximum(ctrl["dist_pc_v1"].to_numpy(float), 1.0))
    sb = np.abs(np.sin(np.deg2rad(ctrl["b_deg"].to_numpy(float))))
    h_G, h_logd, h_sb = 0.3, 0.15, 0.1

    def kde_ratio(x: np.ndarray, anchors: np.ndarray, h: float) -> np.ndarray:
        # Chunk anchors so this stays memory-light.
        out = np.zeros(len(x), dtype=float)
        for start in range(0, len(anchors), 500):
            a = anchors[start:start + 500]
            out += np.exp(-0.5 * ((x[None, :] - a[:, None]) / h) ** 2).sum(axis=0)
        return out / max(len(anchors), 1)

    w = kde_ratio(G, slow_G, h_G) * kde_ratio(logd, slow_logd, h_logd) * kde_ratio(sb, slow_sb, h_sb)
    med = np.nanmedian(w[np.isfinite(w) & (w > 0)])
    if np.isfinite(med) and med > 0:
        w = w / med
    w = np.clip(w, 0, np.nanpercentile(w[np.isfinite(w)], 99.5))
    mean = np.nanmean(w[np.isfinite(w) & (w > 0)])
    return w / mean if np.isfinite(mean) and mean > 0 else np.ones(len(ctrl))


def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    master = Table.read(MASTER).to_pandas()
    master["tier"] = clean_strings(master["tier"])
    slow = master[master["tier"].isin(["A", "B", "C"])].copy()
    orbits = Table.read(ORBITS).to_pandas()
    orbits["tier"] = clean_strings(orbits["tier"])
    controls = Table.read(CONTROLS).to_pandas()
    controls["band"] = clean_strings(controls["band"])
    return slow, orbits, controls


def summarize(slow: pd.DataFrame, orbits: pd.DataFrame, controls: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    slow_res = orbits["res_ratio_OmegaR_over_dPhi"].to_numpy(float)
    rows.append({
        "band": "slow",
        "label": "slow Tier A+B+C",
        "n": int(len(orbits)),
        "median_R_peri_pc": float(np.nanmedian(orbits["static_R_peri_kpc"]) * 1000.0),
        "median_R_apo_kpc": float(np.nanmedian(orbits["static_R_apo_kpc"])),
        "median_ecc": float(np.nanmedian(orbits["static_ecc"])),
        "olr_frac_pct": float(np.mean(np.abs(slow_res - 2.0) < 0.3) * 100.0),
        "four_to_one_frac_pct": float(np.mean(np.abs(slow_res - 4.0) < 0.3) * 100.0),
        "weight_note": "unweighted expanded slow sample",
    })

    weight_rows = []
    for band, label, _color in BANDS[1:]:
        sub = controls[controls["band"] == band].copy()
        w = kernel_weights(sub, slow)
        sub["expanded_kw_weight"] = w
        weight_rows.append(sub[["band", "source_id", "expanded_kw_weight"]])
        rperi = sub["R_peri_kpc"].to_numpy(float) * 1000.0
        rapo = sub["R_apo_kpc"].to_numpy(float)
        ecc = sub["ecc"].to_numpy(float)
        res = sub["res_ratio"].to_numpy(float)
        rows.append({
            "band": band,
            "label": f"control {label}",
            "n": int(len(sub)),
            "median_R_peri_pc": weighted_quantile(rperi, w, [0.5])[0],
            "median_R_apo_kpc": weighted_quantile(rapo, w, [0.5])[0],
            "median_ecc": weighted_quantile(ecc, w, [0.5])[0],
            "olr_frac_pct": float(np.sum(w[np.abs(res - 2.0) < 0.3]) / np.sum(w) * 100.0),
            "four_to_one_frac_pct": float(np.sum(w[np.abs(res - 4.0) < 0.3]) / np.sum(w) * 100.0),
            "weight_note": "expanded slow-sample kernel weights",
        })
    return pd.DataFrame(rows), pd.concat(weight_rows, ignore_index=True)


def plot_figures(summary: pd.DataFrame, orbits: pd.DataFrame, controls: pd.DataFrame, weights: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    c = controls.merge(weights, on=["band", "source_id"], how="left")

    # Eccentricity vs velocity band medians.
    fig, ax = plt.subplots(figsize=(4.2, 3.25))
    x = np.arange(len(summary))
    colors = [b[2] for b in BANDS]
    ax.plot(x, summary["median_ecc"], color="black", lw=0.8, zorder=1)
    ax.scatter(x, summary["median_ecc"], s=36, c=colors, edgecolors="black", linewidths=0.4, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(["slow", "25-50", "50-100", "100-200", "200-260"], rotation=25, ha="right")
    ax.set_ylabel("Median eccentricity")
    ax.set_xlabel(r"$V_\mathrm{GRF}$ band (km s$^{-1}$)")
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig09_eccentricity_vs_vgrf.pdf", bbox_inches="tight", dpi=250)
    fig.savefig(OUT / "expanded_control_eccentricity_trend.png", bbox_inches="tight", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.2, 3.25))
    ax.plot(x, summary["median_R_peri_pc"], color="black", lw=0.8, zorder=1)
    ax.scatter(x, summary["median_R_peri_pc"], s=36, c=colors, edgecolors="black", linewidths=0.4, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(["slow", "25-50", "50-100", "100-200", "200-260"], rotation=25, ha="right")
    ax.set_ylabel(r"Median $R_\mathrm{peri}$ (pc)")
    ax.set_xlabel(r"$V_\mathrm{GRF}$ band (km s$^{-1}$)")
    ax.set_yscale("log")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(FIG / "fig12_rperi_rapo_by_band.pdf", bbox_inches="tight", dpi=250)
    fig.savefig(OUT / "expanded_control_rperi_trend.png", bbox_inches="tight", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(5, 1, figsize=(3.7, 5.4), sharex=True)
    bins = np.linspace(0, 1, 36)
    axes[0].hist(orbits["static_ecc"], bins=bins, color=colors[0], alpha=0.75, density=True)
    axes[0].axvline(summary.loc[0, "median_ecc"], color="black", ls="--", lw=0.8)
    axes[0].text(0.02, 0.78, f"slow (N={len(orbits):,})", transform=axes[0].transAxes, fontsize=6.5)
    for ax, (band, label, color) in zip(axes[1:], BANDS[1:]):
        sub = c[c["band"] == band]
        ax.hist(sub["ecc"], bins=bins, color=color, alpha=0.75, density=True, weights=sub["expanded_kw_weight"])
        med = summary.loc[summary["band"] == band, "median_ecc"].iloc[0]
        ax.axvline(med, color="black", ls="--", lw=0.8)
        ax.text(0.02, 0.78, f"{label} (N={len(sub):,})", transform=ax.transAxes, fontsize=6.5)
    for ax in axes:
        ax.set_yticks([])
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("Eccentricity")
    fig.tight_layout()
    fig.savefig(FIG / "fig10_eccentricity_by_band.pdf", bbox_inches="tight", dpi=250)
    fig.savefig(OUT / "expanded_control_eccentricity_histograms.png", bbox_inches="tight", dpi=180)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    slow, orbits, controls = load()
    summary, weights = summarize(slow, orbits, controls)
    summary.to_csv(OUT / "expanded_control_comparison_summary.csv", index=False)
    weights.to_csv(OUT / "expanded_control_weights.csv", index=False)
    payload = {
        "n_slow_tierABC": int(len(orbits)),
        "summary": summary.to_dict(orient="records"),
        "outputs": {
            "summary_csv": str(OUT / "expanded_control_comparison_summary.csv"),
            "weights_csv": str(OUT / "expanded_control_weights.csv"),
        },
    }
    (OUT / "expanded_control_comparison_summary.json").write_text(json.dumps(payload, indent=2))
    plot_figures(summary, orbits, controls, weights)
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
