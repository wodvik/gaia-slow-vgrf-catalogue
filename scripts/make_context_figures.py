"""Expanded edge-on and face-on context figures.

Computes orbit-event coordinates for the expanded Tier A+B+C catalogue
and regenerates:
  - figures/fig01_edge_on_Rgc_Zgc.pdf
  - figures/fig02_face_on_Xgc_Ygc.pdf

Run from WSL when the coordinate cache must be built because AGAMA is
installed in the WSL Python environment.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
from astropy.table import Table


REPO = Path(__file__).resolve().parents[1]
ORBITS = REPO / "catalogues/catalogue_expanded_orbits_tierABC.fits"
WORK = REPO / "phase3_agama/_hunter24_workdir"
OUT = REPO / "analysis_products"
FIG = REPO / "figures"
CACHE = OUT / "expanded_orbit_event_positions.csv"
SUMMARY = OUT / "expanded_context_figures_summary.json"

EDGE_PDF = FIG / "fig01_edge_on_Rgc_Zgc.pdf"
FACE_PDF = FIG / "fig02_face_on_Xgc_Ygc.pdf"
EDGE_PNG = OUT / "fig01_edge_on_Rgc_Zgc_expanded.png"
FACE_PNG = OUT / "fig02_face_on_Xgc_Ygc_expanded.png"

GYR = 1.0 / 0.9778
TRAJSIZE = 40001
CHUNK = 150
TIERS = ("A", "B", "C")
PLOT_ORDER = ("C", "B", "A")
STYLE = {
    "A": {"marker": "o", "size": 14, "alpha": 0.90},
    "B": {"marker": "^", "size": 13, "alpha": 0.82},
    "C": {"marker": "s", "size": 6, "alpha": 0.42},
}


def log(message: str) -> None:
    print(f"[14S-context t={time.time() - T0:7.1f}s] {message}", flush=True)


def clean_tier(s: pd.Series) -> pd.Series:
    if s.dtype == object and len(s) and isinstance(s.iloc[0], (bytes, bytearray)):
        return s.str.decode("utf-8").str.strip()
    return s.astype(str).str.strip()


def load_orbits() -> pd.DataFrame:
    df = Table.read(ORBITS).to_pandas()
    df["tier"] = clean_tier(df["tier"])
    df = df[df["tier"].isin(TIERS)].copy().reset_index(drop=True)
    return df


def compute_cache() -> pd.DataFrame:
    import agama

    agama.setUnits(length=1, mass=1, velocity=1)
    df = load_orbits()
    ic = df[["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]].to_numpy(float)
    pot = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))

    rows: list[dict] = []
    log(f"computing orbit-event cache for {len(df)} expanded Tier A+B+C stars")
    for start in range(0, len(df), CHUNK):
        stop = min(start + CHUNK, len(df))
        log(f"chunk {start}-{stop} trajsize={TRAJSIZE}")
        res = agama.orbit(potential=pot, ic=ic[start:stop], time=4.0 * GYR, trajsize=TRAJSIZE)
        for k in range(stop - start):
            row = df.iloc[start + k]
            traj = np.asarray(res[k, 1])
            x = traj[:, 0]
            y = traj[:, 1]
            z = traj[:, 2]
            R = np.hypot(x, y)
            j_peri = int(np.argmin(R))
            j_apo = int(np.argmax(R))
            j_zmax = int(np.argmax(np.abs(z)))
            rows.append({
                "source_id": int(row["source_id"]),
                "tier": str(row["tier"]),
                "x_current_kpc": float(row["x_kpc"]),
                "y_current_kpc": float(row["y_kpc"]),
                "z_current_kpc": float(row["z_kpc"]),
                "R_current_kpc": float(np.hypot(row["x_kpc"], row["y_kpc"])),
                "x_peri_kpc": float(x[j_peri]),
                "y_peri_kpc": float(y[j_peri]),
                "z_peri_kpc": float(z[j_peri]),
                "R_peri_kpc": float(R[j_peri]),
                "x_apo_kpc": float(x[j_apo]),
                "y_apo_kpc": float(y[j_apo]),
                "z_apo_kpc": float(z[j_apo]),
                "R_apo_kpc": float(R[j_apo]),
                "x_zmax_kpc": float(x[j_zmax]),
                "y_zmax_kpc": float(y[j_zmax]),
                "z_zmax_kpc": float(z[j_zmax]),
                "R_zmax_kpc": float(R[j_zmax]),
                "abs_zmax_kpc": float(abs(z[j_zmax])),
            })
    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(CACHE, index=False)
    log(f"wrote {CACHE}")
    return out


def load_or_compute_cache() -> pd.DataFrame:
    if CACHE.exists():
        cached = pd.read_csv(CACHE)
        if len(cached) == len(load_orbits()):
            log(f"loaded {CACHE} ({len(cached)} rows)")
            return cached
        log("cache row count does not match expanded catalogue; rebuilding")
    return compute_cache()


def merge_context(coords: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "source_id", "tier", "vgrf_default_exact", "static_ecc",
        "static_Lz_kpc_kms", "static_R_peri_kpc", "static_R_apo_kpc",
    ]
    orb = load_orbits()[cols].copy()
    df = coords.drop(columns=["tier"], errors="ignore").merge(orb, on="source_id", how="left", validate="one_to_one")
    df["tier"] = clean_tier(df["tier"])
    return df.sort_values("source_id").reset_index(drop=True)


def tier_legend(ax: plt.Axes, df: pd.DataFrame, loc: str = "upper right") -> None:
    handles = []
    labels = []
    for tier in TIERS:
        st = STYLE[tier]
        handles.append(Line2D([0], [0], marker=st["marker"], color="black", linestyle="None",
                              markersize=np.sqrt(st["size"] * 2.8)))
        labels.append(f"Tier {tier} (N={(df['tier'] == tier).sum():,})")
    ax.legend(handles, labels, loc=loc, fontsize=6.5, frameon=True, framealpha=0.9)


def scatter_tiers(ax: plt.Axes, df: pd.DataFrame, x: str, y: str, c: str,
                  vmin: float, vmax: float, scale: float = 1.0):
    last = None
    for tier in PLOT_ORDER:
        sub = df[df["tier"] == tier]
        st = STYLE[tier]
        last = ax.scatter(
            sub[x], sub[y], c=sub[c], cmap="inferno", vmin=vmin, vmax=vmax,
            marker=st["marker"], s=st["size"] * scale, alpha=st["alpha"],
            edgecolors="none",
        )
    return last


def draw_edge(df: pd.DataFrame) -> None:
    panels = [
        ("Current positions", "R_current_kpc", "z_current_kpc", "vgrf_default_exact", 0, 25),
        ("Modelled pericentres", "R_peri_kpc", "z_peri_kpc", "static_ecc", 0.8, 1.0),
        (r"Modelled maximum $|Z_{\rm GC}|$", "R_zmax_kpc", "z_zmax_kpc", "static_ecc", 0.8, 1.0),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 4.55), sharey=True)
    rmax = max(8.178, *(float(np.nanmax(df[p[1]])) for p in panels)) * 1.05
    zmax = max(4.0, *(float(np.nanmax(np.abs(df[p[2]]))) for p in panels)) * 1.10
    scatters = []
    for ax, (title, x, y, c, vmin, vmax) in zip(axes, panels):
        scatters.append(scatter_tiers(ax, df, x, y, c, vmin, vmax))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r"$R_{\rm GC}$ (kpc)")
        ax.axhline(0, color="0.65", ls=":", lw=0.5)
        ax.axvline(0, color="0.65", ls=":", lw=0.5)
        ax.axvline(8.178, color="0.55", ls="-.", lw=0.7, alpha=0.65)
        ax.set_xlim(0, rmax)
        ax.set_ylim(-zmax, zmax)
    axes[0].set_ylabel(r"$Z_{\rm GC}$ (kpc)")
    tier_legend(axes[0], df)
    axins = inset_axes(axes[1], width="38%", height="46%", loc="upper right", borderpad=0.9)
    scatter_tiers(axins, df, "R_peri_kpc", "z_peri_kpc", "static_ecc", 0.8, 1.0, scale=0.5)
    axins.set_xlim(0, 1.0)
    axins.set_ylim(-1.0, 1.0)
    axins.set_title("inner 1 kpc", fontsize=6, pad=1.5)
    axins.tick_params(labelsize=5, length=2, pad=1)

    fig.subplots_adjust(left=0.07, right=0.985, top=0.86, bottom=0.27, wspace=0.18)
    boxes = [ax.get_position() for ax in axes]
    cax_v = fig.add_axes([boxes[0].x0, 0.13, boxes[0].width, 0.023])
    sm_v = ScalarMappable(norm=Normalize(0, 25), cmap="inferno")
    cb_v = fig.colorbar(sm_v, cax=cax_v, orientation="horizontal")
    cb_v.set_label(r"$V_{\rm GRF}$ (km s$^{-1}$)", fontsize=8)
    cax_e = fig.add_axes([boxes[1].x0, 0.13, boxes[-1].x0 + boxes[-1].width - boxes[1].x0, 0.023])
    sm_e = ScalarMappable(norm=Normalize(0.8, 1.0), cmap="inferno")
    cb_e = fig.colorbar(sm_e, cax=cax_e, orientation="horizontal")
    cb_e.set_label("Eccentricity", fontsize=8)
    fig.savefig(EDGE_PDF, dpi=250)
    fig.savefig(EDGE_PNG, dpi=180)
    plt.close(fig)
    log(f"wrote {EDGE_PDF}")


def draw_face(df: pd.DataFrame) -> None:
    plot = df.copy()
    plot["x_current_plot_kpc"] = -plot["x_current_kpc"]
    plot["x_peri_plot_kpc"] = -plot["x_peri_kpc"]

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.15))
    sc0 = scatter_tiers(axes[0], plot, "x_current_plot_kpc", "y_current_kpc", "vgrf_default_exact", 0, 25)
    sc1 = scatter_tiers(axes[1], plot, "x_peri_plot_kpc", "y_peri_kpc", "static_ecc", 0.8, 1.0)

    x_all = np.concatenate([
        plot["x_current_plot_kpc"].to_numpy(float),
        plot["x_peri_plot_kpc"].to_numpy(float),
        np.array([0.0, 8.178]),
    ])
    y_all = np.concatenate([
        plot["y_current_kpc"].to_numpy(float),
        plot["y_peri_kpc"].to_numpy(float),
        np.array([0.0]),
    ])
    x_min = min(-1.4, float(np.nanmin(x_all)) - 0.4)
    x_max = max(9.6, float(np.nanmax(x_all)) + 0.4)
    y_half = max(abs(float(np.nanmin(y_all))), abs(float(np.nanmax(y_all))), (x_max - x_min) / 2.0)
    y_pad = 0.25

    for ax, title in zip(axes, ["Current positions", "Modelled projected pericentres"]):
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r"$-X_{\rm GC}$ (kpc)")
        ax.set_ylabel(r"$Y_{\rm GC}$ (kpc)")
        ax.axhline(0, color="0.65", ls=":", lw=0.5)
        ax.axvline(0, color="0.65", ls=":", lw=0.5)
        ax.scatter([0], [0], marker="+", s=45, c="black", linewidths=0.9)
        ax.scatter([8.178], [0], marker="*", s=36, c="0.2")
        ax.text(0, 0.45, "GC", ha="center", fontsize=7)
        ax.text(8.178, 0.45, "Sun", ha="center", fontsize=7)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-(y_half + y_pad), y_half + y_pad)
        ax.set_aspect("equal", adjustable="box")
    tier_legend(axes[0], df, loc="upper right")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.26, wspace=0.22)
    b0, b1 = [ax.get_position() for ax in axes]
    cax0 = fig.add_axes([b0.x0, 0.11, b0.width, 0.025])
    cb0 = fig.colorbar(ScalarMappable(norm=Normalize(0, 25), cmap="inferno"), cax=cax0, orientation="horizontal")
    cb0.set_label(r"$V_{\rm GRF}$ (km s$^{-1}$)", fontsize=8)
    cax1 = fig.add_axes([b1.x0, 0.11, b1.width, 0.025])
    cb1 = fig.colorbar(ScalarMappable(norm=Normalize(0.8, 1.0), cmap="inferno"), cax=cax1, orientation="horizontal")
    cb1.set_label("Eccentricity", fontsize=8)
    fig.savefig(FACE_PDF, dpi=250)
    fig.savefig(FACE_PNG, dpi=180)
    plt.close(fig)
    log(f"wrote {FACE_PDF}")


def write_summary(df: pd.DataFrame) -> None:
    payload = {
        "n": int(len(df)),
        "tier_counts": {tier: int((df["tier"] == tier).sum()) for tier in TIERS},
        "median_R_peri_pc": float(np.median(df["R_peri_kpc"]) * 1000.0),
        "max_R_peri_pc": float(np.max(df["R_peri_kpc"]) * 1000.0),
        "median_R_apo_kpc": float(np.median(df["R_apo_kpc"])),
        "median_abs_zmax_kpc": float(np.median(df["abs_zmax_kpc"])),
        "figures": [str(EDGE_PDF), str(FACE_PDF)],
        "cache": str(CACHE),
    }
    SUMMARY.write_text(json.dumps(payload, indent=2))
    log(json.dumps(payload, indent=2))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    coords = load_or_compute_cache()
    df = merge_context(coords)
    draw_edge(df)
    draw_face(df)
    write_summary(df)
    return 0


if __name__ == "__main__":
    T0 = time.time()
    sys.exit(main())
