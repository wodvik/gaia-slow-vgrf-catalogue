"""Phase 16I -- covariance-stress headline on the population-prior primary catalogue.

Phase 14AE stressed the parallax--distance coupling with a maximal Gaussian
copula and recorded, per source, the membership probability under independent
and coupled draws. Its headline was expressed against the forward-defined
Tier A+B (527 of 541 retained). With the catalogue now tiered on the
population-prior probability, the headline must be restated for the 276-star
primary sample.

The stress table itself is unchanged -- it is a property of the astrometric
covariance, not of the tier definition -- so this pass only re-derives the
headline and rewrites the caption. The confusion matrix, which is defined over
the full P >= 0.3 band, is likewise untouched.

Outputs: phase14/expanded_covariance_stress_summary.json  (headline updated)
         tables/v15/tab_covariance_stress.tex
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
STRESS = BUNDLE / "phase14" / "expanded_covariance_stress_per_star.csv"
SIDECAR = BUNDLE / "phase14" / "expanded_covariance_stress_summary.json"
RETIER = BUNDLE / "catalogues" / "catalogue_retier_tierAB.fits"
TAB_DIRS = [BUNDLE / "tables" / "v15", REPO / "release" / "tables" / "v15"]

THRESH_AB = 0.84


def main() -> None:
    cs = pd.read_csv(STRESS)
    ab_ids = {int(s) for s in Table.read(RETIER)["source_id"]}
    sub = cs[cs["source_id"].astype("int64").isin(ab_ids)]
    n_ab = len(sub)

    # A member is "retained" if it stays above the A+B threshold under the
    # coupled draw. The population-prior members sit well clear of the
    # boundary, so in practice none are lost.
    retained = int((sub["P_copula"] > THRESH_AB).sum())
    lost = n_ab - retained
    dp = np.abs(sub["delta_P"].to_numpy(dtype=float))

    d = json.loads(SIDECAR.read_text(encoding="utf-8"))
    d["headline_AB_adopted"] = n_ab
    d["headline_AB_lost_under_copula"] = lost
    d["headline_AB_retained_under_copula"] = retained
    d["headline_tier_definition"] = "population-prior (Phase 16F)"
    d["headline_AB_abs_delta_P"] = {
        "median": float(np.median(dp)), "p95": float(np.percentile(dp, 95)),
        "max": float(dp.max()),
    }
    d["headline_AB_forward_reference"] = {
        "n": 541, "retained": 527, "lost": 14,
        "note": "pre-retier headline, retained for provenance",
    }
    SIDECAR.write_text(json.dumps(d, indent=2), encoding="utf-8")

    med, mx = float(np.median(dp)), float(dp.max())
    sentence = (f"Relative to the adopted primary catalogue, {retained} of {n_ab} "
                f"Tier~A+B members remain above $P=0.84$ under the coupled draw "
                f"(median $|\\Delta P|={med:.4f}$, max ${mx:.4f}$).")

    for td in TAB_DIRS:
        f = td / "tab_covariance_stress.tex"
        if not f.exists():
            continue
        txt = f.read_text(encoding="utf-8")
        # lambda replacement: the sentence contains LaTeX backslashes, which
        # re.sub would otherwise interpret as group escapes.
        txt = re.sub(
            r"Relative to the adopted primary catalogue, \d+ of [\d,]+ Tier~A\+B members "
            r"remain above \$P=0\.84\$ under the coupled draw "
            r"\(median \$\|\\Delta P\|=[\d.]+\$, max \$[\d.]+\$\)\.",
            lambda _m: sentence, txt)
        f.write_text(txt, encoding="utf-8")
        print(f"updated {f}")

    print(f"population-prior Tier A+B: {retained}/{n_ab} retained, {lost} lost; "
          f"median |dP|={med:.4f}, max={mx:.4f}")


if __name__ == "__main__":
    main()
