"""
Phase 9 — generate convenience tier-subset FITS files for the Zenodo deposit.

Reads the master `catalogue_v2.fits` (2,859 rows; tiers A/B/C/D/X) and writes:
  - catalogue_tierA.fits   : Tier A only        (214 rows; P > 0.95)
  - catalogue_tierAB.fits  : Tier A+B headline  (334 rows; P > 0.84)

Both files preserve every column of the master so any downstream join keys
(`source_id`, `tier`, MC orbit posteriors, chemistry crossmatches) remain valid.

Outputs are written next to the master in release/v2/phase1/.
"""

from pathlib import Path

import numpy as np
from astropy.table import Table

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "release" / "v2" / "phase1" / "catalogue_v2.fits"
OUT_DIR = MASTER.parent

EXPECTED = {"A": 214, "AB": 334}


def main() -> None:
    master = Table.read(MASTER)
    tier = np.char.strip(np.array(master["tier"], dtype=str))

    subsets = {
        "A": tier == "A",
        "AB": np.isin(tier, ["A", "B"]),
    }

    for label, mask in subsets.items():
        n = int(mask.sum())
        assert n == EXPECTED[label], f"Tier {label} expected {EXPECTED[label]} rows, got {n}"
        sub = master[mask]
        out = OUT_DIR / f"catalogue_tier{label}.fits"
        sub.write(out, overwrite=True)
        print(f"wrote {out.relative_to(REPO_ROOT)}: {n} rows, {len(sub.colnames)} cols")


if __name__ == "__main__":
    main()
