"""Phase 14AI -- full-sample convergent action audit (referee Issue 1).

The review-stage catalogue ``J_R/J_z/J_phi`` are single-evaluation AGAMA
Staeckel-fudge actions (phase0g ``stackel_actions``).  For this highly radial
(e>0.9) population the Staeckel fudge is unreliable, and the earlier audit only
covered a 100-star Tier-A control subset (WP-5), leaving every other star with
an ``unsampled_global_caveat`` flag.

This script extends the convergent audit to ALL Tier A+B+C stars using the
IDENTICAL axisymmetric potential and integration settings as the catalogue and
WP-5: integrate each orbit for 4 Gyr in the static Hunter+2024 potential
(trajsize=1001, dprkn8, accuracy=1e-8), evaluate the interpolated ActionFinder
along the trajectory, and time-average to obtain a convergent action estimate
(``J_*_timeavg``).  Each star then gets a real ``action_max_fracdiff`` and
``action_reliability_flag`` instead of the global caveat.

Run from WSL (Agama installed):
    cd /mnt/c/Users/humbl/GAIA2026
    python3 release/gaia_slow_vgrf_catalogue_v1.0.5_review/scripts/phase14ai_full_action_audit.py

Updates ``catalogues/catalogue_expanded_orbits_tierABC.fits`` (+ .csv) and writes
``phase14/expanded_action_audit_per_star.csv`` and
``phase14/expanded_action_audit_summary.json``.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import agama
import numpy as np
import pandas as pd
from astropy.table import Table

T0 = time.time()


def log(msg: str) -> None:
    print(f"[14ai t={time.time() - T0:7.1f}s] {msg}", flush=True)


def find_repo() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / ".git").exists():
            return parent
    return Path("/mnt/c/Users/humbl/GAIA2026")


REPO = find_repo()
BUNDLE = Path(__file__).resolve().parents[1]
PRIVATE_WORK = REPO / "release" / "_iterations" / "v2" / "phase3_agama" / "_hunter24_workdir"
POTENTIALS = BUNDLE / "potentials"
WORK = POTENTIALS if (POTENTIALS / "MWPotentialHunter24_axi.ini").exists() else PRIVATE_WORK
AXI_INI = WORK / "MWPotentialHunter24_axi.ini"

ORBIT_FITS = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.fits"
ORBIT_CSV = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.csv"
AUDIT_CSV = BUNDLE / "phase14" / "expanded_action_audit_per_star.csv"
AUDIT_JSON = BUNDLE / "phase14" / "expanded_action_audit_summary.json"

GYR = 1.0 / 0.9778
T_INTEGRATE = 4.0 * GYR
TRAJSIZE = 1001
CHUNK = 400
DENOM_FLOOR = 1.0

agama.setUnits(length=1, mass=1, velocity=1)

IC_COLS = ["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]
CAT_J = ["J_R", "J_z", "J_phi"]


def timeavg_actions(ic: np.ndarray, pot_axi) -> np.ndarray:
    """Orbit-time-averaged convergent actions for each IC (N x 3)."""
    af = agama.ActionFinder(pot_axi, interp=True)
    out = np.full((len(ic), 3), np.nan, dtype=float)
    for start in range(0, len(ic), CHUNK):
        sl = slice(start, min(start + CHUNK, len(ic)))
        res = agama.orbit(potential=pot_axi, ic=ic[sl], time=T_INTEGRATE,
                          trajsize=TRAJSIZE, accuracy=1e-8, method="dprkn8")
        for k in range(sl.stop - sl.start):
            traj = np.asarray(res[k, 1], dtype=float)
            acts = np.asarray(af(traj), dtype=float)
            out[start + k] = np.nanmean(acts, axis=0)
        log(f"  time-averaged {sl.stop}/{len(ic)} orbits")
    return out


def flag_from_frac(fmax: np.ndarray) -> np.ndarray:
    return np.select([fmax <= 0.15, fmax <= 0.50],
                     ["sampled_ok", "sampled_caution"], default="sampled_poor")


def main() -> int:
    log(f"agama {agama.__version__}; repo {REPO}")
    if not AXI_INI.exists():
        raise FileNotFoundError(f"axi potential not found: {AXI_INI}")
    tab = Table.read(ORBIT_FITS)
    n = len(tab)
    log(f"loaded {n} Tier A+B+C orbit rows")
    ic = np.column_stack([np.asarray(tab[c], dtype=float) for c in IC_COLS])
    cat = {c: np.asarray(tab[c], dtype=float) for c in CAT_J}

    pot_axi = agama.Potential(file=str(AXI_INI))
    log("integrating + time-averaging actions for all stars")
    tavg = timeavg_actions(ic, pot_axi)

    denom = np.maximum(np.abs(tavg), DENOM_FLOOR)
    frac = np.abs(np.column_stack([cat["J_R"], cat["J_z"], cat["J_phi"]]) - tavg) / denom
    fmax = np.nanmax(frac, axis=1)
    flag = flag_from_frac(fmax)

    # --- update orbit catalogue columns in place ---
    tab["J_R_timeavg"] = tavg[:, 0]
    tab["J_z_timeavg"] = tavg[:, 1]
    tab["J_phi_timeavg"] = tavg[:, 2]
    tab["J_R_timeavg"].description = "Orbit-time-averaged (convergent) radial action over a 4 Gyr static integration."
    tab["J_z_timeavg"].description = "Orbit-time-averaged (convergent) vertical action over a 4 Gyr static integration."
    tab["J_phi_timeavg"].description = "Orbit-time-averaged (convergent) azimuthal action over a 4 Gyr static integration."
    tab["action_accuracy_sampled"] = np.ones(n, dtype=bool)
    tab["action_max_fracdiff"] = fmax
    if "action_reliability_flag" in tab.colnames:
        tab.remove_column("action_reliability_flag")
    tab["action_reliability_flag"] = np.asarray(flag, dtype="U24")
    tab.write(ORBIT_FITS, overwrite=True)
    log(f"updated {ORBIT_FITS.name}")

    # --- regenerate the CSV from the final table (chem + actions consistent) ---
    df = tab.to_pandas()
    for c in df.columns:
        if df[c].dtype == object and len(df) and isinstance(df[c].iloc[0], (bytes, bytearray)):
            df[c] = df[c].str.decode("utf-8")
    df.to_csv(ORBIT_CSV, index=False)
    log(f"regenerated {ORBIT_CSV.name}")

    # --- per-star audit product ---
    tier = np.asarray([t.decode() if isinstance(t, (bytes, bytearray)) else str(t) for t in tab["tier"]])
    audit = pd.DataFrame({
        "source_id": np.asarray(tab["source_id"]).astype("int64"),
        "tier": tier,
        "J_R_catalogue": cat["J_R"], "J_z_catalogue": cat["J_z"], "J_phi_catalogue": cat["J_phi"],
        "J_R_timeavg": tavg[:, 0], "J_z_timeavg": tavg[:, 1], "J_phi_timeavg": tavg[:, 2],
        "J_R_fracdiff": frac[:, 0], "J_z_fracdiff": frac[:, 1], "J_phi_fracdiff": frac[:, 2],
        "action_max_fracdiff": fmax,
        "action_reliability_flag": flag,
    })
    audit.to_csv(AUDIT_CSV, index=False)

    def pct(a):
        a = a[np.isfinite(a)]
        return {"n": int(len(a)), "median": float(np.median(a)),
                "p84": float(np.percentile(a, 84)), "max": float(np.max(a))}

    abc = np.isin(tier, ["A", "B", "C"])
    ab = np.isin(tier, ["A", "B"])
    summary = {
        "agama_version": str(agama.__version__),
        "potential": "Hunter+2024 axisymmetric (MWPotentialHunter24_axi.ini)",
        "integration": {"time_Gyr": 4.0, "trajsize": TRAJSIZE, "method": "dprkn8",
                        "accuracy": 1e-8, "denominator_floor_kpc_kms": DENOM_FLOOR},
        "n_total": int(n),
        "reliability_counts_all": {f: int((flag == f).sum())
                                   for f in ["sampled_ok", "sampled_caution", "sampled_poor"]},
        "reliability_counts_tierAB": {f: int(((flag == f) & ab).sum())
                                      for f in ["sampled_ok", "sampled_caution", "sampled_poor"]},
        "action_max_fracdiff": pct(fmax),
        "J_R_fracdiff": pct(frac[:, 0]),
        "J_z_fracdiff": pct(frac[:, 1]),
        "J_phi_fracdiff": pct(frac[:, 2]),
        "frac_poor_all": float((fmax > 0.50).mean()),
        "frac_reliable_le0p50_all": float((fmax <= 0.50).mean()),
        "note": "Convergent reference is the orbit-time-averaged interpolated ActionFinder; "
                "catalogue J_* are single-evaluation Staeckel-fudge actions at the IC.",
    }
    AUDIT_JSON.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    log("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
