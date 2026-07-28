"""Regenerate core manuscript figures from expanded catalogue products."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from astropy.table import Table


BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
MASTER = BUNDLE / "catalogues" / "catalogue_expanded_master.fits"
ORBITS = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.fits"
ENERGY = BUNDLE / "phase14" / "expanded_static_energy.csv"

# --- population-prior retier switch (Phase 16F) -------------------------------
# With GAIA_RETIER=1 these scripts read the retiered catalogue products, whose
# `tier` column is defined on the population-prior probability rather than the
# forward Monte Carlo score. Default behaviour is unchanged.
import os as _os
_RETIER = _os.environ.get("GAIA_RETIER", "").lower() in ("1", "true", "yes")
if _RETIER:
    MASTER = BUNDLE / "catalogues" / "catalogue_retier_master.fits"
    ORBITS = BUNDLE / "catalogues" / "catalogue_retier_orbits_tierABC.fits"
# -----------------------------------------------------------------------------

FIG = BUNDLE / "figures"
OUT = BUNDLE / "phase14"

TIER_STYLE = {
    "A": {"color": "#b2182b", "marker": "o", "size": 12, "alpha": 0.85, "label": "Tier A"},
    "B": {"color": "#ef8a62", "marker": "^", "size": 11, "alpha": 0.78, "label": "Tier B"},
    "C": {"color": "#2166ac", "marker": "s", "size": 7, "alpha": 0.45, "label": "Tier C"},
    "D": {"color": "#7f7f7f", "marker": ".", "size": 5, "alpha": 0.25, "label": "Tier D"},
    "X": {"color": "#c7c7c7", "marker": ".", "size": 4, "alpha": 0.12, "label": "Tier X"},
}


def clean_tier(s: pd.Series) -> pd.Series:
    if s.dtype == object and len(s) and isinstance(s.iloc[0], (bytes, bytearray)):
        return s.str.decode("utf-8").str.strip()
    return s.astype(str).str.strip()


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    master = Table.read(MASTER).to_pandas()
    master["tier"] = clean_tier(master["tier"])
    orbits = Table.read(ORBITS).to_pandas()
    orbits["tier"] = clean_tier(orbits["tier"])
    return master, orbits


def save(fig: plt.Figure, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / name, bbox_inches="tight", dpi=300)
    fig.savefig(OUT / name.replace(".pdf", ".png"), bbox_inches="tight", dpi=180)
    plt.close(fig)
    print(f"wrote {FIG / name}")


def tier_marker_legend(
    ax: plt.Axes,
    df: pd.DataFrame,
    loc: str = "upper right",
    *,
    bbox_to_anchor: tuple[float, float] | None = None,
    ncol: int = 1,
    short_labels: bool = False,
) -> None:
    handles = []
    for tier in ("C", "B", "A"):
        st = TIER_STYLE[tier]
        label = tier if short_labels else st["label"]
        handles.append(Line2D(
            [0], [0],
            marker=st["marker"], linestyle="None",
            markersize=max(4.0, np.sqrt(st["size"]) * 1.15),
            markerfacecolor="0.25", markeredgecolor="black",
            markeredgewidth=0.45, color="black",
            label=f"{label} (N={int((df['tier'] == tier).sum()):,})",
        ))
    ax.legend(
        handles=handles,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor,
        ncol=ncol,
        fontsize=5.6 if ncol > 1 else 6,
        frameon=True,
        borderaxespad=0.2,
        handletextpad=0.45,
        columnspacing=0.9,
    )


def plot_pvgrf(master: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 3.3))
    for tier in ("X", "D", "C", "B", "A"):
        sub = master[master["tier"] == tier]
        st = TIER_STYLE[tier]
        ax.scatter(
            sub["vgrf_default"], sub["P_vgrf_below_25"],
            s=st["size"], c=st["color"], marker=st["marker"],
            alpha=st["alpha"], edgecolors="none", rasterized=True,
            label=f"{st['label']} (N={len(sub):,})",
        )
    ax.axvline(25, color="black", lw=0.8)
    for y, txt in [(0.50, "C"), (0.84, "B"), (0.95, "A")]:
        ax.axhline(y, color="black", ls="--", lw=0.55, alpha=0.45)
        ax.text(52, y + 0.012, txt, ha="right", va="bottom", fontsize=7)
    ax.set_xlim(0, 55)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel(r"Point-estimate $V_\mathrm{GRF}$ (km s$^{-1}$)")
    ax.set_ylabel(r"$P(V_\mathrm{GRF}<25\,\mathrm{km\,s}^{-1})$")
    ax.legend(loc="center right", fontsize=5.8, frameon=True)
    fig.tight_layout()
    save(fig, "fig_phase14_pvgrf_vs_vgrf.pdf")


def plot_sky(master: pd.DataFrame) -> None:
    abc = master[master["tier"].isin(["A", "B", "C"])].copy()
    l = np.deg2rad(((abc["l"].to_numpy() + 180.0) % 360.0) - 180.0)
    b = np.deg2rad(abc["b"].to_numpy())
    fig = plt.figure(figsize=(7.2, 3.8))
    ax = fig.add_subplot(111, projection="mollweide")
    for tier in ("C", "B", "A"):
        sub = abc[abc["tier"] == tier]
        ll = np.deg2rad(((sub["l"].to_numpy() + 180.0) % 360.0) - 180.0)
        bb = np.deg2rad(sub["b"].to_numpy())
        st = TIER_STYLE[tier]
        sc = ax.scatter(
            -ll, bb, c=sub["vgrf_default"], cmap="inferno", vmin=0, vmax=25,
            s=st["size"], marker=st["marker"], alpha=st["alpha"],
            edgecolors="none", rasterized=True, label=f"{st['label']} (N={len(sub):,})",
        )
    ax.grid(True, alpha=0.35)
    ax.set_xlabel("Galactic longitude")
    ax.set_ylabel("Galactic latitude")
    cb = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.08, fraction=0.055)
    cb.set_label(r"$V_\mathrm{GRF}$ (km s$^{-1}$)")
    fig.tight_layout()
    save(fig, "fig05_sky_aitoff.pdf")


def plot_rperi_rapo(orbits: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(3.6, 3.25))
    for tier in ("C", "B", "A"):
        sub = orbits[orbits["tier"] == tier]
        st = TIER_STYLE[tier]
        sc = ax.scatter(
            sub["static_R_peri_kpc"], sub["static_R_apo_kpc"],
            c=sub["static_ecc"], cmap="inferno", vmin=0.8, vmax=1.0,
            s=st["size"], marker=st["marker"], alpha=max(st["alpha"], 0.58),
            edgecolors="0.15", linewidths=0.18, rasterized=True,
        )
    ax.set_xlim(0, 0.65)
    ax.set_ylim(0, 19)
    ax.set_xlabel(r"$R_\mathrm{peri}$ (kpc)")
    ax.set_ylabel(r"$R_\mathrm{apo}$ (kpc)")
    tier_marker_legend(ax, orbits, loc="upper right")
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Eccentricity")
    fig.tight_layout()
    save(fig, "fig07_rperi_rapo.pdf")


def plot_toomre(orbits: pd.DataFrame) -> None:
    df = orbits.copy()
    R = np.hypot(df["x_kpc"], df["y_kpc"])
    df["V_R"] = (df["x_kpc"] * df["vx_kms"] + df["y_kpc"] * df["vy_kms"]) / R
    df["V_phi"] = (df["x_kpc"] * df["vy_kms"] - df["y_kpc"] * df["vx_kms"]) / R
    df["V_perp"] = np.sqrt(df["V_R"] ** 2 + df["vz_kms"] ** 2)
    fig, ax = plt.subplots(figsize=(3.6, 3.45))
    for tier in ("C", "B", "A"):
        sub = df[df["tier"] == tier]
        st = TIER_STYLE[tier]
        sc = ax.scatter(
            sub["V_phi"], sub["V_perp"], c=sub["vgrf_default_exact"],
            cmap="inferno", vmin=0, vmax=25, s=st["size"],
            marker=st["marker"], alpha=st["alpha"], edgecolors="none",
            rasterized=True, label=f"{st['label']} (N={len(sub):,})",
        )
    theta = np.linspace(0, np.pi, 240)
    for v in [50, 100, 150, 200, 250]:
        ax.plot(v * np.cos(theta), v * np.sin(theta), color="0.7", lw=0.45, alpha=0.6)
    ax.set_xlim(-300, 350)
    ax.set_ylim(0, 300)
    ax.set_xlabel(r"$V_\phi$ (km s$^{-1}$)")
    ax.set_ylabel(r"$\sqrt{V_R^2 + V_z^2}$ (km s$^{-1}$)")
    ax.legend(loc="upper right", fontsize=6, frameon=True)
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label(r"$V_\mathrm{GRF}$ (km s$^{-1}$)")
    fig.tight_layout()
    save(fig, "fig08_toomre_diagram.pdf")


def plot_elz(orbits: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(3.6, 3.25))
    if ENERGY.exists():
        energy = pd.read_csv(ENERGY)
        orbits = orbits.merge(energy, on="source_id", how="left", validate="one_to_one")
    if "static_E" not in orbits.columns or orbits["static_E"].isna().any():
        raise RuntimeError(f"Missing static energy cache: {ENERGY}")
    y = orbits["static_E"]
    ylabel = r"$E$ (km$^2$ s$^{-2}$)"
    for tier in ("C", "B", "A"):
        sub = orbits[orbits["tier"] == tier]
        yy = y.loc[sub.index] if hasattr(y, "loc") else y[sub.index]
        st = TIER_STYLE[tier]
        sc = ax.scatter(
            sub["static_Lz_kpc_kms"], yy, c=sub["static_ecc"], cmap="inferno",
            vmin=0.8, vmax=1.0, s=st["size"], marker=st["marker"],
            alpha=max(st["alpha"], 0.58), edgecolors="0.15",
            linewidths=0.18, rasterized=True,
        )
    ax.axvline(0, color="black", lw=0.7, alpha=0.7)
    ax.set_xlabel(r"$L_z$ (kpc km s$^{-1}$)")
    ax.set_ylabel(ylabel)
    tier_marker_legend(
        ax,
        orbits,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.015),
        ncol=3,
        short_labels=True,
    )
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Eccentricity")
    fig.tight_layout()
    save(fig, "fig15_e_lz.pdf")


def main() -> int:
    master, orbits = load()
    plot_pvgrf(master)
    plot_sky(master)
    plot_rperi_rapo(orbits)
    plot_toomre(orbits)
    plot_elz(orbits)
    summary = {
        "n_master": int(len(master)),
        "n_tier_ab": int(master["tier"].isin(["A", "B"]).sum()),
        "n_tier_abc": int(master["tier"].isin(["A", "B", "C"]).sum()),
        "figures": [
            "fig_phase14_pvgrf_vs_vgrf.pdf",
            "fig05_sky_aitoff.pdf",
            "fig07_rperi_rapo.pdf",
            "fig08_toomre_diagram.pdf",
            "fig15_e_lz.pdf",
        ],
    }
    (OUT / "expanded_core_figures_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
