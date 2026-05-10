"""Phase 3 (agama) - MW potential via Hunter+2024.

Builds the Hunter+2024 agama potential configuration used by the release.


Approach:
  Use agama's bundled `example_mw_potential_hunter24.py` to generate the
  Hunter+2024 MW model — this is the literature-standard combination of
  the Sormani+2022 bar + Chatzopoulos NSC + Sormani+2020 NSD + thin/thick
  stellar disks (with a 2.7 kpc inner cutoff so the bar fills the
  central region) + HI/H2 gas + dark halo, with NO McMillan bulge to
  double-count.

Variants:
  static_simple       MWPotential2014 (axi)
  static_full_axi     Hunter+2024 axisymmetrised (Sormani bar -> axi-CylSpline)
  barred_default      Hunter+2024 with rotating bar @ Omega_p = -37.5 km/s/kpc
  barred_slow         Hunter+2024 with rotating bar @ Omega_p = -33   km/s/kpc
  barred_fast         Hunter+2024 with rotating bar @ Omega_p = -41   km/s/kpc
  static_full_halo_0p85   Hunter+2024 axi but with NFW halo mass × 0.85
  static_full_halo_1p15   Hunter+2024 axi but with NFW halo mass × 1.15
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml
import agama

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase3_agama"
WORK = OUT / "_hunter24_workdir"
WORK.mkdir(parents=True, exist_ok=True)

R0 = float(CONFIG["solar_variants"]["default"]["R0_kpc"])
V0 = float(CONFIG["solar_variants"]["default"]["Vc_kms"])
BARS = CONFIG["bar_pattern_speeds_kms_kpc"]

agama.setUnits(length=1, mass=1, velocity=1)
AGAMA_DATA = Path(agama.__file__).parent / "data"
AGAMA_PY = Path(agama.__file__).parent / "py"


def run_hunter24():
    """Run agama's Hunter+2024 example in WORK dir to generate .ini files.

    Skips re-run if MWPotentialHunter24_axi.ini already exists.
    """
    needed = ["MWPotentialHunter24_axi.ini",
              "MWPotentialHunter24_full.ini",
              "MWPotentialHunter24_rotating.ini"]
    if all((WORK / f).exists() for f in needed):
        print("[3v2] Hunter24 .ini files already present; skipping rebuild")
        return
    # Copy example + bar dependency into WORK (the example uses relative
    # imports to example_mw_bar_potential).
    shutil.copy(AGAMA_PY / "example_mw_potential_hunter24.py",
                WORK / "example_mw_potential_hunter24.py")
    shutil.copy(AGAMA_PY / "example_mw_bar_potential.py",
                WORK / "example_mw_bar_potential.py")
    print("[3v2] running Hunter+2024 example -> generates .ini files...")
    # Use --noplot to skip matplotlib; if not supported, just suppress show().
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    cmd = [sys.executable, "example_mw_potential_hunter24.py"]
    proc = subprocess.run(cmd, cwd=str(WORK), env=env,
                          capture_output=True, text=True, timeout=900)
    print("STDOUT:", proc.stdout[-2000:])
    if proc.returncode != 0:
        print("STDERR:", proc.stderr[-2000:])
        raise RuntimeError(f"Hunter24 example failed (rc={proc.returncode})")
    print("[3v2] Hunter24 .ini files generated")


def vc_at(pot, R, phi=0.0):
    pos = np.array([R * np.cos(phi), R * np.sin(phi), 0.0])
    f = pot.force(pos)
    Fr = -(f[0] * np.cos(phi) + f[1] * np.sin(phi))
    return float(np.sqrt(max(R * Fr, 0.0)))


def vc_phiavg(pot, R, n_phi=16):
    phis = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    return float(np.sqrt(np.mean([vc_at(pot, R, p)**2 for p in phis])))


def menc_kepler(pot, r, n_phi=16):
    G = 4.30091e-6  # kpc * (km/s)^2 / Msun
    return r * vc_phiavg(pot, r, n_phi)**2 / G


def make_halo_variant(axi_pot, halo_factor: float):
    """Add a delta-NFW to scale the halo by `halo_factor` (assumes Hunter+24
    uses an NFW halo with mass ~1.27e12 at scaleRadius~19.6 kpc)."""
    delta_halo = agama.Potential(type="NFW", scaleRadius=19.6,
                                 mass=1.27e12 * (halo_factor - 1.0))
    return agama.Potential(axi_pot, delta_halo)


def main():
    # Step 1: produce Hunter+2024 .ini files
    run_hunter24()

    # Step 2: load each into agama
    pot_axi = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    pot_full = agama.Potential(file=str(WORK / "MWPotentialHunter24_full.ini"))
    # Load the rotating .ini file (Omega encoded inside)
    pot_rot_default = agama.Potential(file=str(WORK / "MWPotentialHunter24_rotating.ini"))

    # Build slow/fast rotating variants by duplicating the rotating .ini and
    # patching Omega.
    rot_text = (WORK / "MWPotentialHunter24_rotating.ini").read_text()
    print("\n[3v2] rotating .ini contents (preview):")
    print(rot_text[:600])

    # The rotation block is `rotation=[[0,phi0],[1,phi0+Omega]]`. Replace
    # 37.5 with 33 / 41 numerically.
    def make_rotating(omega_p_kms_kpc: float, name: str) -> agama.Potential:
        # Re-open the full (non-rotating) .ini and write a custom rotating
        # wrapper to avoid string-rewriting brittleness.
        full_ini = WORK / "MWPotentialHunter24_full.ini"
        # Pattern speed sign convention: negative = clockwise rotation
        # (consistent with Hunter+ 2024).
        omega_signed = -abs(omega_p_kms_kpc)
        angle_bar = -0.44  # radians, Hunter+2024 default
        rot_ini = WORK / f"MWPotentialHunter24_rot_{name}.ini"
        rot_ini.write_text(
            "# Hunter+2024 with custom Omega_p = %.1f km/s/kpc\n"
            "[Potential]\n"
            "file=%s\n"
            "rotation=[[0,%.4f],[1,%.4f]]\n" % (
                omega_signed, str(full_ini), angle_bar, angle_bar + omega_signed))
        return agama.Potential(file=str(rot_ini))

    pots = {
        "static_simple":  agama.Potential(file=str(AGAMA_DATA / "MWPotential2014.ini")),
        "static_full":    pot_axi,
        "barred_default": make_rotating(BARS["default"], "default"),
        "barred_slow":    make_rotating(BARS["slow"],    "slow"),
        "barred_fast":    make_rotating(BARS["fast"],    "fast"),
        "static_full_halo_0p85": make_halo_variant(pot_axi, 0.85),
        "static_full_halo_1p15": make_halo_variant(pot_axi, 1.15),
    }

    # V_c phi-averaged + M_enc curves
    R_kpc = np.linspace(0.5, 25, 99)
    r_kpc = np.linspace(0.5, 200, 100)
    summary = {}
    vc = {"R_kpc": R_kpc.tolist()}
    me = {"r_kpc": r_kpc.tolist()}
    for name, pot in pots.items():
        vc_arr = np.array([vc_phiavg(pot, r) for r in R_kpc])
        me_arr = np.array([menc_kepler(pot, r) for r in r_kpc])
        vc[name] = vc_arr.tolist()
        me[name] = me_arr.tolist()
        # Also record phi=0 for barred variants (so we can see range)
        vc_phi0_R0 = vc_at(pot, R0, phi=0.0)
        vc_avg_R0 = vc_phiavg(pot, R0, n_phi=32)
        # Range
        phis = np.linspace(0, 2 * np.pi, 16, endpoint=False)
        vcs_R0 = [vc_at(pot, R0, p) for p in phis]
        summary[name] = {
            "Vc_phiavg_R0_kms": vc_avg_R0,
            "Vc_phi0_R0_kms":   vc_phi0_R0,
            "Vc_R0_min_kms":    float(min(vcs_R0)),
            "Vc_R0_max_kms":    float(max(vcs_R0)),
            "Vc_max_kms":       float(np.max(vc_arr)),
            "Mass_enc_50kpc_msun":  float(np.interp(50, r_kpc, me_arr)),
            "Mass_enc_200kpc_msun": float(np.interp(200, r_kpc, me_arr)),
        }
        print(f"  {name:24s} Vc_phi-avg(R0)={vc_avg_R0:6.1f} km/s "
              f"(phi=0: {vc_phi0_R0:6.1f}, range "
              f"[{min(vcs_R0):.1f},{max(vcs_R0):.1f}])  "
              f"M(<200)={summary[name]['Mass_enc_200kpc_msun']:.3e}")

    import csv
    with (OUT / "vc_curves_v2.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(list(vc.keys()))
        for i in range(len(R_kpc)):
            w.writerow([vc[k][i] for k in vc])
    with (OUT / "menc_curves_v2.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(list(me.keys()))
        for i in range(len(r_kpc)):
            w.writerow([me[k][i] for k in me])

    # Gate
    gate_ok = (
        all(180 < s["Vc_phiavg_R0_kms"] < 270 for s in summary.values()) and
        all(4e11 < s["Mass_enc_200kpc_msun"] < 3e12 for s in summary.values())
    )

    out = {
        "agama_version": agama.__version__,
        "approach": "Hunter+2024 example (peer-review-clean composition)",
        "no_bulge_double_counting": True,
        "summary_per_variant": summary,
        "R0_kpc": R0, "V0_kms": V0, "bar_omegas": BARS,
        "gate_physical_ok": gate_ok,
    }
    (OUT / "gate3_potentials_v2.json").write_text(json.dumps(out, indent=2))

    # Plots
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    for name in pots:
        ax.plot(R_kpc, vc[name], lw=1.4, label=name)
    ax.axvline(R0, ls=":", color="black", label=f"R0={R0} kpc")
    ax.set_xlabel("R (kpc)"); ax.set_ylabel("V_c phi-averaged (km/s)")
    ax.set_title("Phase 3 v2 (agama, Hunter+2024) — circular speed curves")
    ax.set_ylim(0, 320); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(OUT / "gate3_vc_v2.png", dpi=140); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for name in pots:
        ax.plot(r_kpc, me[name], lw=1.4, label=name)
    ax.set_xlabel("r (kpc)"); ax.set_ylabel("M(<r) (Msun)")
    ax.set_yscale("log")
    ax.set_title("Phase 3 v2 (agama, Hunter+2024) — enclosed mass profiles")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(OUT / "gate3_menc_v2.png", dpi=140); plt.close(fig)

    print(f"\n[3v2] gate_physical_ok = {gate_ok}")
    print(json.dumps(summary, indent=2))
    return 0 if gate_ok else 2


if __name__ == "__main__":
    sys.exit(main())
