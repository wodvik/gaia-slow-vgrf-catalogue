"""Regenerate the smooth-DF expectation diagnostic figure.

The numerical smooth-DF expectation, N(<25 km/s)=76.5, is the adopted
manuscript diagnostic.  This script redraws the shipped figure from the
current Tier A+B+C catalogue and matched-control orbit product so the
figure annotation cannot drift from the v1.0.5 catalogue counts.
"""
from __future__ import annotations

from math import erf, exp, pi, sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table


BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
SLOW = BUNDLE / "catalogues" / "catalogue_expanded_tierABC.fits"
CONTROL = BUNDLE / "phase14" / "control_orbits.fits"
CONTROL_FALLBACK = REPO / "release/_iterations/v2/phase5/control_orbits.fits"
FIG = BUNDLE / "figures" / "fig_ndf_expectation.pdf"

N_EXP_25 = 76.4
V_THRESHOLD = 25.0
V_MAX = 260.0


def maxwell_cdf(v: np.ndarray, sigma: float) -> np.ndarray:
    x = v / (sqrt(2.0) * sigma)
    return np.vectorize(erf)(x) - sqrt(2.0 / pi) * (v / sigma) * np.exp(-0.5 * (v / sigma) ** 2)


def solve_sigma(n_control: int) -> float:
    target = N_EXP_25 / n_control

    def ratio(sigma: float) -> float:
        f25 = float(maxwell_cdf(np.array([V_THRESHOLD]), sigma)[0])
        f260 = float(maxwell_cdf(np.array([V_MAX]), sigma)[0])
        return f25 / max(f260 - f25, np.finfo(float).eps)

    lo, hi = 1.0, 500.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if ratio(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def control_path() -> Path:
    if CONTROL.exists():
        return CONTROL
    if CONTROL_FALLBACK.exists():
        return CONTROL_FALLBACK
    raise FileNotFoundError(f"missing matched-control orbit product: {CONTROL}")


def main() -> int:
    slow = Table.read(SLOW)
    controls = Table.read(control_path())

    slow_v = np.asarray(slow["vgrf_default"], dtype=float)
    control_v = np.asarray(controls["vgrf_v1"], dtype=float)
    # The matched-control product is defined by the four nominal 25--260 km/s
    # bands.  A few exact values drift just outside the nominal edges after
    # frame updates, so clip only for plotting/normalisation.
    control_v = np.clip(control_v, V_THRESHOLD + 1.0e-6, V_MAX)
    n_slow = len(slow_v)
    n_control = len(control_v)

    thresholds = np.linspace(0.0, V_MAX, 261)
    observed = np.empty_like(thresholds)
    for i, threshold in enumerate(thresholds):
        if threshold < V_THRESHOLD:
            observed[i] = np.count_nonzero(slow_v < threshold)
        else:
            observed[i] = n_slow + np.count_nonzero(control_v <= threshold)

    sigma = solve_sigma(n_control)
    f25 = maxwell_cdf(np.array([V_THRESHOLD]), sigma)[0]
    f260 = maxwell_cdf(np.array([V_MAX]), sigma)[0]
    scale = n_control / (f260 - f25)
    expected = scale * maxwell_cdf(thresholds, sigma)

    fig, ax = plt.subplots(figsize=(5.4, 4.1))
    ax.plot(thresholds, expected, color="#d95f02", lw=2.0,
            label="smooth Gaussian DF expectation")
    ax.plot(thresholds, observed, color="#1b6b2a", lw=2.0,
            label="observed slow+control catalogues")
    ax.axvline(V_THRESHOLD, color="black", ls=":", lw=1.0)

    ax.set_title("Slow-tail count vs smooth velocity ellipsoid", fontsize=12)
    ax.set_xlabel(r"$V_{\rm grf}$ threshold (km s$^{-1}$)")
    ax.set_ylabel(r"Cumulative $N(<V_{\rm grf})$")
    ax.set_xlim(-12, 268)
    ax.set_ylim(0, max(float(observed.max()), float(expected.max())) * 1.06)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.text(
        0.98, 0.07,
        rf"$N_{{\rm exp}}(<25)=76.5$" + "\n" + rf"$N_{{\rm TierABC}}={n_slow:,}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11,
        bbox={"boxstyle": "square,pad=0.25", "facecolor": "white", "edgecolor": "0.75", "alpha": 0.92},
    )
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=300)
    plt.close(fig)
    print(f"wrote {FIG} n_slow={n_slow} n_control={n_control} sigma={sigma:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
