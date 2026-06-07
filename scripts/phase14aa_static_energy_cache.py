"""Regenerate the static-energy cache used by the E-Lz figure.

Run this helper from a Python environment with AGAMA installed, for example:

    python3 scripts/phase14aa_static_energy_cache.py
"""
from __future__ import annotations

from pathlib import Path

import agama
import numpy as np
import pandas as pd
from astropy.table import Table


BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
ORBITS = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.fits"
PRIVATE_WORK = REPO / "release/_iterations/v2/phase3_agama/_hunter24_workdir"
POTENTIALS = BUNDLE / "potentials"
WORK = POTENTIALS if (POTENTIALS / "MWPotentialHunter24_axi.ini").exists() else PRIVATE_WORK
OUT = BUNDLE / "phase14" / "expanded_static_energy.csv"


def main() -> int:
    agama.setUnits(length=1, mass=1, velocity=1)
    table = Table.read(ORBITS).to_pandas()
    pot = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))

    xyz = table[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(float)
    vel = table[["vx_kms", "vy_kms", "vz_kms"]].to_numpy(float)
    phi = np.asarray(pot.potential(xyz), dtype=float)
    kinetic = 0.5 * np.sum(vel * vel, axis=1)

    out = pd.DataFrame({
        "source_id": table["source_id"].astype("int64"),
        "static_E": phi + kinetic,
    })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT} rows={len(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
