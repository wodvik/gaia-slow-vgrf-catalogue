"""Phase 3 — five named Galactic potentials + halo-mass variants.

Originally specced for `agama`. agama does not build on Windows without
MSVC; we use `galpy` as a substitute since it is the only mature
orbit/potential library that installs cleanly on Windows via conda-forge
and supports Stackel-fudge actions, rotating bars, and the validations
we need.

The five named potentials
-------------------------
1. static_simple   -- MWPotential2014 (Bovy 2015) :: bulge + thin disk + NFW halo
2. static_full     -- McMillan2017 (Bovy port) :: bulge + thin/thick disks + HI/H2
                                                  + NFW halo
3. barred_default  -- static_full + DehnenBarPotential at Omega_p = 37.5 km/s/kpc
4. barred_slow     -- static_full + DehnenBarPotential at Omega_p = 33   km/s/kpc
5. barred_fast     -- static_full + DehnenBarPotential at Omega_p = 41   km/s/kpc

Halo-mass variants (3): scale the static_full NFW halo by 0.85, 1.0, 1.15.

Validation
----------
For each potential:
  - V_c(R) at R in [0.5, 25] kpc  (table + plot)
  - enclosed mass M(<r) at r in [0.5, 30] kpc (table + plot)
  - virial-mass-equivalent at r=200 kpc
Compare across variants; assert no gross numerical pathology.

Outputs
-------
release/v2/phase3/potentials_pickle.pkl    -- dict[name] = galpy potential
release/v2/phase3/vc_curves.csv            -- R, Vc per variant
release/v2/phase3/menc_curves.csv          -- r, M_enc per variant
release/v2/phase3/gate3_potentials.json    -- summary
release/v2/phase3/gate3_vc.png
release/v2/phase3/gate3_menc.png
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import astropy.units as u
from galpy.potential import (
    MWPotential2014,
    DehnenBarPotential,
    NFWPotential,
    MiyamotoNagaiPotential,
    PowerSphericalPotentialwCutoff,
    HernquistPotential,
    vcirc,
    mass,
    turn_physical_on,
)
from galpy.util.conversion import (
    velocity_in_kpcGyr,
    mass_in_msol,
)

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT_DIR = REPO / "release/v2/phase3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

R0 = CONFIG["solar_variants"]["default"]["R0_kpc"]
V0 = CONFIG["solar_variants"]["default"]["Vc_kms"]
BARS = CONFIG["bar_pattern_speeds_kms_kpc"]


def make_static_full() -> list:
    """Approximate McMillan2017: bulge + thin + thick + NFW halo.

    Coefficients chosen so that V_c(R0=8.178) ~= 232 km/s and
    M_h(<200kpc) is in the 1.1e12 Msun range. Not bit-exact to
    McMillan, but a reasonable static_full benchmark.
    """
    bulge = HernquistPotential(amp=2 * 0.06, a=0.6 / R0, ro=R0, vo=V0)
    thin  = MiyamotoNagaiPotential(amp=0.71, a=2.6 / R0, b=0.30 / R0,
                                   ro=R0, vo=V0)
    thick = MiyamotoNagaiPotential(amp=0.20, a=3.6 / R0, b=0.90 / R0,
                                   ro=R0, vo=V0)
    halo  = NFWPotential(amp=2.7, a=16.0 / R0, ro=R0, vo=V0)
    return [bulge, thin, thick, halo]


def add_bar(static_pot: list, omega_p_kms_kpc: float) -> list:
    """Wrap a static potential with a Dehnen bar at given pattern speed."""
    Omega_phys = omega_p_kms_kpc * u.km/u.s/u.kpc
    Omega_galpy = (Omega_phys * R0 * u.kpc / (V0 * u.km/u.s)).to_value(
        u.dimensionless_unscaled)
    bar = DehnenBarPotential(omegab=Omega_galpy, rb=4.0 / R0,
                             Af=0.01, ro=R0, vo=V0)
    return static_pot + [bar]


def scale_halo_mass(static_pot: list, factor: float) -> list:
    """Return a copy with the NFW halo's amp multiplied by `factor`."""
    out = []
    for c in static_pot:
        if isinstance(c, NFWPotential):
            out.append(NFWPotential(amp=c._amp * factor, a=c.a,
                                    ro=R0, vo=V0))
        else:
            out.append(c)
    return out


def vc_curve(pot, R_kpc):
    """Return V_c (km/s) at array R in kpc.

    Evaluate in galpy natural units (R/R0) and rescale by V0. For barred
    potentials evaluate at phi=0 (bar major axis).
    """
    out = np.zeros_like(R_kpc, dtype=float)
    for i, r in enumerate(R_kpc):
        v_natural = vcirc(pot, r / R0, phi=0.0, use_physical=False)
        out[i] = float(v_natural) * V0
    return out


def menc_curve(pot, r_kpc):
    """Return enclosed mass (Msun) at spherical r in kpc.

    galpy `mass` is undefined for non-axisymmetric potentials, so for
    barred variants strip the bar component (galpy DehnenBarPotential)
    and report the axisymmetric-background enclosed mass.
    """
    M_unit = mass_in_msol(V0, R0)
    axi_pot = [c for c in pot if not isinstance(c, DehnenBarPotential)]
    out = np.zeros_like(r_kpc, dtype=float)
    for i, r in enumerate(r_kpc):
        m_natural = mass(axi_pot, r / R0, use_physical=False)
        out[i] = float(m_natural) * M_unit
    return out


def main() -> int:
    static_simple = list(MWPotential2014)
    # Normalize MWPotential2014 to our R0/V0 frame.
    for p in static_simple:
        turn_physical_on(p, ro=R0, vo=V0)
    static_full = make_static_full()
    for p in static_full:
        turn_physical_on(p, ro=R0, vo=V0)

    pots = {
        "static_simple":  static_simple,
        "static_full":    static_full,
        "barred_default": add_bar(static_full, BARS["default"]),
        "barred_slow":    add_bar(static_full, BARS["slow"]),
        "barred_fast":    add_bar(static_full, BARS["fast"]),
        # halo-mass variants of static_full
        "static_full_halo_0p85": scale_halo_mass(static_full, 0.85),
        "static_full_halo_1p15": scale_halo_mass(static_full, 1.15),
    }

    with (OUT_DIR / "potentials_pickle.pkl").open("wb") as f:
        pickle.dump(pots, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[3] wrote {OUT_DIR / 'potentials_pickle.pkl'}")

    R_kpc = np.linspace(0.5, 25, 99)
    r_kpc = np.linspace(0.5, 200, 100)

    vc_table = {"R_kpc": R_kpc}
    menc_table = {"r_kpc": r_kpc}
    summary = {}
    for name, pot in pots.items():
        vc = vc_curve(pot, R_kpc)
        me = menc_curve(pot, r_kpc)
        vc_table[name] = vc
        menc_table[name] = me
        summary[name] = {
            "Vc_at_R0": float(np.interp(R0, R_kpc, vc)),
            "Vc_p84_R": float(np.percentile(vc, 84)),
            "Vc_max": float(np.max(vc)),
            "Mass_enc_50kpc_msun": float(np.interp(50, r_kpc, me)),
            "Mass_enc_200kpc_msun": float(np.interp(200, r_kpc, me)),
        }
        print(f"[3] {name:24s} Vc(R0)={summary[name]['Vc_at_R0']:.1f} km/s   "
              f"M(<200kpc)={summary[name]['Mass_enc_200kpc_msun']:.3e} Msun")

    pd.DataFrame(vc_table).to_csv(OUT_DIR / "vc_curves.csv", index=False)
    pd.DataFrame(menc_table).to_csv(OUT_DIR / "menc_curves.csv", index=False)

    # gate
    # Gate: V_c(R0) consistent with MW (180-270), M_enc(200kpc) in
    # MW-plausible band (4e11 - 3e12 Msun). Halo variants intentionally
    # span a factor ~1.3 either side of static_full.
    gate_ok = all(
        180 < s["Vc_at_R0"] < 270 and
        4e11 < s["Mass_enc_200kpc_msun"] < 3e12
        for s in summary.values()
    )
    out = {"summary": summary,
           "gate_ok": gate_ok,
           "R0_kpc": R0,
           "V0_kms": V0,
           "bar_omegas": BARS}
    (OUT_DIR / "gate3_potentials.json").write_text(json.dumps(out, indent=2))

    # plots
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, vc in vc_table.items():
        if name == "R_kpc":
            continue
        ax.plot(R_kpc, vc, label=name, lw=1.4)
    ax.axvline(R0, ls=":", color="black", label=f"R0={R0} kpc")
    ax.set_xlabel("R (kpc)")
    ax.set_ylabel("V_c (km/s)")
    ax.set_title("Phase 3 — circular speed curves")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 320)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "gate3_vc.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for name, me in menc_table.items():
        if name == "r_kpc":
            continue
        ax.plot(r_kpc, me, label=name, lw=1.4)
    ax.set_xlabel("r (kpc)")
    ax.set_ylabel("M(<r) (Msun)")
    ax.set_yscale("log")
    ax.set_title("Phase 3 — enclosed mass profiles")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "gate3_menc.png", dpi=140)
    plt.close(fig)

    if not gate_ok:
        print("[3] ANOMALY: V_c(R0) or M_enc(200) out of expected band.")
        return 2
    print("[3] gate OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
