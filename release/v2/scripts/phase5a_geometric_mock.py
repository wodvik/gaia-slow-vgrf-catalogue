"""Phase 5A — geometric-bias mock for the e(Vgrf) trend.

Question
========
The peer review wanted us to demonstrate that the eccentricity-vs-Vgrf
trend in the slow-Vgrf catalogue has a *physical* component beyond the
geometric apocentre-selection baseline (the trivial expectation that
stars selected at solar position with low Vgrf are biased toward orbits
with apocentres near R0).

Mock construction
=================
1. Sample N_mock stars from a thick-disk-like Gaussian velocity DF at
   solar position (R0 +/- 0.5 kpc, |z|<0.5 kpc):
     sigma_R = 60, sigma_phi = 40, sigma_z = 40 km/s
     v_phi mean = 200 km/s
   These dispersions span thin-thick blend.
2. Integrate each in static_full (Hunter+2024 axi) for 4 Gyr.
3. Compute v_grf at the solar position (the natural observable),
   eccentricity e = (R_apo - R_peri) / (R_apo + R_peri), and R_peri.
4. Compare e(Vgrf) trend against v2 observed.

Output
======
release/v2/phase5/geometric_mock.fits
release/v2/phase5/gate5A_e_vs_vgrf.png
release/v2/phase5/gate5A_geometric.json
release/v2/phase5/PHASE5A_NOTE.md
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import agama
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase5"
OUT.mkdir(parents=True, exist_ok=True)
WORK = REPO / "release/v2/phase3_agama/_hunter24_workdir"

R0 = float(CONFIG["solar_variants"]["default"]["R0_kpc"])
V0 = float(CONFIG["solar_variants"]["default"]["Vc_kms"])

agama.setUnits(length=1, mass=1, velocity=1)


def log(msg):
    print(f"[5A t={time.time()-T0:5.1f}s] {msg}", flush=True)


def make_mock(n_mock: int, rng: np.random.Generator,
              sig_iso_kms: float = 150.0):
    """Sample N stars at solar neighbourhood with isotropic halo-like DF.

    Isotropic Gaussian velocity (sigma=150 km/s) with NO net rotation —
    this is the geometric baseline for the slow-Vgrf catalogue, which is
    dominated by halo / Splash / Aurora stars whose underlying DF is
    closer to isotropic than to disk rotation. The 150 km/s dispersion
    populates the full Vgrf range from 0 to ~600 km/s with a fat tail
    in the slow-Vgrf bins (Vgrf < 50 km/s) where the observed catalogue
    lives.
    """
    R   = R0 + rng.uniform(-0.5, 0.5, size=n_mock)
    phi = rng.uniform(0, 2 * np.pi, size=n_mock)
    z   = rng.uniform(-0.5, 0.5, size=n_mock)
    vx  = rng.normal(0, sig_iso_kms, size=n_mock)
    vy  = rng.normal(0, sig_iso_kms, size=n_mock)
    vz  = rng.normal(0, sig_iso_kms, size=n_mock)
    x = R * np.cos(phi); y = R * np.sin(phi)
    ic = np.column_stack([x, y, z, vx, vy, vz])
    vgrf = np.sqrt(vx*vx + vy*vy + vz*vz)
    return ic, vgrf


def main():
    global T0; T0 = time.time()
    log(f"agama {agama.__version__}")

    # Load static_full Hunter+2024 axi
    log("loading static_full potential...")
    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))

    n_mock = 50000
    rng = np.random.default_rng(20260503)
    log(f"sampling {n_mock} mock stars (thick-disk DF) at solar nbd...")
    ic, vgrf_mock = make_mock(n_mock, rng)

    Gyr = 1.0 / 0.9778
    log("integrating mock orbits 4 Gyr...")
    t0 = time.time()
    res = agama.orbit(potential=pot_axi, ic=ic, time=4.0 * Gyr, trajsize=1001)
    log(f"  done in {time.time()-t0:.1f}s")

    log("computing R_peri / R_apo / e for mock orbits...")
    R_peri = np.zeros(n_mock); R_apo = np.zeros(n_mock); zmax = np.zeros(n_mock)
    for i in range(n_mock):
        traj = np.asarray(res[i, 1])
        Rcyl = np.sqrt(traj[:, 0]**2 + traj[:, 1]**2)
        R_peri[i] = float(Rcyl.min()); R_apo[i] = float(Rcyl.max())
        zmax[i] = float(np.abs(traj[:, 2]).max())
    ecc_mock = (R_apo - R_peri) / (R_apo + R_peri)

    # Save mock catalogue
    mock = Table()
    mock["x_kpc"] = ic[:, 0]; mock["y_kpc"] = ic[:, 1]; mock["z_kpc"] = ic[:, 2]
    mock["vx_kms"] = ic[:, 3]; mock["vy_kms"] = ic[:, 4]; mock["vz_kms"] = ic[:, 5]
    mock["vgrf_kms"] = vgrf_mock
    mock["R_peri_kpc"] = R_peri
    mock["R_apo_kpc"]  = R_apo
    mock["z_max_kpc"]  = zmax
    mock["ecc"] = ecc_mock
    mock.write(OUT / "geometric_mock.fits", overwrite=True)
    log(f"wrote geometric_mock.fits ({len(mock)} mock stars)")

    # Load v2 orbit catalogue
    log("loading v2 orbit catalogue for comparison...")
    v2 = Table.read(REPO / "release/v2/phase4/catalogue_v2_orbits.fits").to_pandas()
    if v2["tier"].dtype == object and isinstance(v2["tier"].iloc[0],
                                                  (bytes, bytearray)):
        v2["tier"] = v2["tier"].str.decode("utf-8")

    # ---- e(Vgrf) trend in 5 km/s bins ----
    bins = np.arange(0, 105, 5)
    centers = 0.5 * (bins[:-1] + bins[1:])

    def trend(vgrf, ecc):
        med = np.full(len(centers), np.nan)
        n = np.zeros(len(centers), dtype=int)
        for i in range(len(centers)):
            m = (vgrf >= bins[i]) & (vgrf < bins[i+1])
            if m.sum() > 5:
                med[i] = np.median(ecc[m]); n[i] = int(m.sum())
        return med, n

    mock_med, mock_n = trend(vgrf_mock, ecc_mock)
    tier_ABC = v2["tier"].isin(["A", "B", "C"]).to_numpy()
    obs_v = v2.loc[tier_ABC, "vgrf_default"].to_numpy()
    obs_e = v2.loc[tier_ABC, "static_ecc"].to_numpy()
    obs_med, obs_n = trend(obs_v, obs_e)

    # Difference at low-Vgrf bins (the slow-Vgrf interest range)
    excess = []
    for i, c in enumerate(centers):
        if c < 25 and obs_n[i] > 5 and not np.isnan(mock_med[i]):
            excess.append({"vgrf_bin_center": float(c),
                           "obs_med_e": float(obs_med[i]), "obs_n": int(obs_n[i]),
                           "mock_med_e": float(mock_med[i]), "mock_n": int(mock_n[i]),
                           "delta": float(obs_med[i] - mock_med[i])})

    summary = {
        "n_mock": int(n_mock),
        "v2_tierABC_n": int(tier_ABC.sum()),
        "mock_DF": {"sig_iso_kms": 150, "vphi_mean_kms": 0,
                    "rationale": "isotropic halo-like; populates slow-Vgrf tail"},
        "e_vs_vgrf_5kms_bins": {
            "centers": centers.tolist(),
            "obs_median_ecc": obs_med.tolist(),
            "obs_n": obs_n.tolist(),
            "mock_median_ecc": mock_med.tolist(),
            "mock_n": mock_n.tolist(),
        },
        "excess_at_low_vgrf_bins": excess,
        "median_excess_below_25kms": float(np.nanmedian([e["delta"] for e in excess])) if excess else None,
    }
    (OUT / "gate5A_geometric.json").write_text(json.dumps(summary, indent=2))
    log("wrote gate5A_geometric.json")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(centers, mock_med, "o-", color="C0", label=f"geometric mock (n={n_mock})")
    ax.plot(centers, obs_med,  "s-", color="C3", label=f"v2 Tier A+B+C (n={int(tier_ABC.sum())})")
    ax.set_xlim(0, 100); ax.set_ylim(0, 1)
    ax.set_xlabel("Vgrf (km/s)"); ax.set_ylabel("median eccentricity")
    ax.set_title("Phase 5A — e(Vgrf): mock vs observed")
    ax.legend(); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(OUT / "gate5A_e_vs_vgrf.png", dpi=140); plt.close(fig)
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
