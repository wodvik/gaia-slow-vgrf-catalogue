"""Phase 16M -- azimuth swept per radial period (apsidal geometry).

A small pericentre invites a comet-like reading in which the star reaches
inward, reverses, and retraces its path. For these orbits that picture is
wrong, and the orbit products already contain what is needed to say so.

The azimuth swept between successive apocentres is 2*pi*Omega_phi/Omega_R.
Reference values:
    360 deg  closed Kepler ellipse -- star returns along its inbound track
    254 deg  flat rotation curve / singular isothermal sphere
    180 deg  harmonic (constant-density) core -- successive apocentres lie
             on opposite sides of the centre

Caveat recorded with the result: these are Staeckel-fudge frequencies, and
the fudge is least reliable for the vertical action of plunging orbits. The
azimuthal quantities are the better-determined ones (J_phi = L_z is exact and
J_R is recovered to a median 8%), and the conclusion here rests on the
distribution being far from 360 deg rather than on any individual value.

Outputs: phase14/apsidal_geometry_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
ORBITS = BUNDLE / "catalogues" / "catalogue_retier_orbits_tierABC.fits"
OUT = BUNDLE / "phase14" / "apsidal_geometry_summary.json"


def main() -> None:
    t = Table.read(ORBITS)
    oR = np.abs(np.asarray(t["Omega_R"], dtype=float))
    op = np.abs(np.asarray(t["Omega_phi"], dtype=float))
    good = np.isfinite(oR) & np.isfinite(op) & (oR > 0)
    sweep = np.degrees(2 * np.pi * op[good] / oR[good])

    tier = np.asarray(t["tier"]).astype(str)[good]
    res = {
        "n_total": int(len(t)),
        "n_with_frequencies": int(good.sum()),
        "sweep_deg_percentiles": {p: float(np.percentile(sweep, p))
                                  for p in (5, 16, 50, 84, 95)},
        "n_above_300deg": int((sweep > 300).sum()),
        "n_below_200deg": int((sweep < 200).sum()),
        "frac_below_200deg": float(np.mean(sweep < 200)),
        "reference_values_deg": {"kepler_closed_ellipse": 360,
                                 "flat_rotation_curve": 254,
                                 "harmonic_core": 180},
        "by_tier": {k: {"n": int((tier == k).sum()),
                        "median_sweep_deg": float(np.median(sweep[tier == k]))}
                    for k in ("A", "B", "C") if (tier == k).any()},
        "caveat": ("Staeckel-fudge frequencies; the conclusion rests on the "
                   "distribution lying far from the 360 deg closed-ellipse "
                   "value, not on individual per-star values."),
    }
    OUT.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
