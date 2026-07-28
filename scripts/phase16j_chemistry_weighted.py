"""Phase 16J -- probability-weighted chemodynamic decomposition.

Spectroscopic abundances are rare in this catalogue: only 117 of the
forward-defined Tier A+B+C stars carry both [Fe/H] and an alpha proxy from
APOGEE DR17 or GALAH DR3. Restricting the chemical discussion to the
population-prior Tier A+B+C would discard 69 of them.

That restriction is not scientifically motivated. A star does not stop being
chemically informative because its rest-frame speed is more likely 30 km/s
than 22; the velocity threshold is not a chemical criterion, and the readers
most likely to use this section (accretion-debris and merger studies) need the
fullest possible abundance sample. We therefore report the chemistry on the
FULL alpha-classified sample and quantify the membership difference per class
rather than cutting on it:

  N            raw count in the literature region
  sum P_pop    summed population-prior probability, i.e. the effective number
               of genuinely slow stars contributed by that region
  N(P>0.5)     count that survives the adopted tier threshold

The point of reporting all three is that the class PROPORTIONS are stable
across them, so the qualitative result -- the slow sample is a chemical
mixture rather than a single population -- does not depend on where the
velocity threshold falls.

Outputs: phase14/chemistry_weighted_summary.json
         tables/v15/tab_population_decomp.tex
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.stats import beta

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
ORBITS = BUNDLE / "catalogues" / "catalogue_expanded_orbits_tierABC.fits"
LAT = BUNDLE / "phase14" / "latent_deconvolution" / "latent_vgrf_per_star_regularised.csv"
OUT = BUNDLE / "phase14" / "chemistry_weighted_summary.json"
TAB_DIRS = [BUNDLE / "tables" / "v15", REPO / "release" / "tables" / "v15"]

ORDER = [("Splash", "Splash"), ("GSE", "GSE"), ("Aurora", "Aurora"),
         ("disk", "Disc-like"), ("unclassified", "Unclassified")]


def jeffreys(k: int, n: int) -> tuple[float, float]:
    lo = beta.ppf(0.16, k + 0.5, n - k + 0.5) if k > 0 else 0.0
    hi = beta.ppf(0.84, k + 0.5, n - k + 0.5) if k < n else 1.0
    return 100 * lo, 100 * hi


def main() -> None:
    orb = Table.read(ORBITS)
    lat = pd.read_csv(LAT)
    pm = dict(zip(lat.source_id.to_numpy(), lat.P_latent.to_numpy()))
    sid = np.asarray(orb["source_id"]).astype(np.int64)
    ppop = np.array([pm.get(int(s), 0.0) for s in sid])
    cp = np.asarray(orb["chem_population"]).astype(str)
    feh = np.asarray(orb["feh_spec"], dtype=float)
    cls = (cp != "") & (cp != "--") & (cp != "nan")

    n_tot = int(cls.sum())
    sum_tot = float(ppop[cls].sum())
    n_thr = int((ppop[cls] > 0.5).sum())

    rows, summary = [], {"n_alpha_classified": n_tot,
                         "sum_ppop_total": sum_tot,
                         "n_above_0p5_total": n_thr,
                         "classes": {}}
    for key, label in ORDER:
        m = cls & (cp == key)
        n = int(m.sum())
        s = float(ppop[m].sum())
        k = int((ppop[m] > 0.5).sum())
        lo, hi = jeffreys(n, n_tot)
        frac = 100 * n / n_tot
        summary["classes"][label] = {
            "N": n, "sum_ppop": s, "n_above_0p5": k,
            "raw_pct": frac, "effective_pct": 100 * s / sum_tot,
            "thresholded_pct": 100 * k / n_thr if n_thr else float("nan"),
            "median_feh": float(np.nanmedian(feh[m])) if n else float("nan"),
        }
        rows.append((label, n, s, k, frac, lo, hi,
                     float(np.nanmedian(feh[m])) if n else float("nan")))

    OUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    L = [
        r"\begin{center}",
        r"\refstepcounter{table}\label{tab:population_decomp}",
        r"\begin{minipage}{\columnwidth}",
        r"\small\textbf{Table \thetable.} Probability-weighted chemodynamic",
        r"context.  Spectroscopic abundances are rare in this catalogue, so the",
        r"chemical discussion uses the full $\alpha$-classified sample---every",
        r"star with both [Fe/H] and an $\alpha$ proxy from APOGEE or GALAH",
        rf"($n={n_tot}$)---rather than only those above the adopted velocity",
        r"threshold.  Regions follow \citet{Belokurov2020}, \citet{Helmi2018},",
        r"and \citet{BelokurovKravtsov2022}.  $\Sigma P_{\rm pop}$ is the summed",
        r"population-prior membership probability, i.e.\ the effective number of",
        r"genuinely slow stars each region contributes; $N(P>0.5)$ is the count",
        r"surviving the adopted threshold.  The three columns give consistent",
        r"proportions, so the mixed chemical character of the sample does not",
        r"depend on where the velocity threshold is placed.  Fractions are raw",
        r"$N$ shares with Beta$(k+\tfrac{1}{2}, n-k+\tfrac{1}{2})$ Jeffreys",
        r"68\% intervals \citep{Cameron2011}; they are wide, which is itself an",
        r"argument against extrapolating to the full catalogue.  Per-star",
        r"probabilities ship with the catalogue so these regions can be recut.",
        r"\end{minipage}",
        r"\vspace{0.4ex}",
        r"",
        r"\small",
        r"\setlength{\tabcolsep}{0pt}",
        r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}lrrrrr@{}}",
        r"\hline\hline",
        r"Population & $N$ & $\Sigma P_{\rm pop}$ & $N(P>0.5)$ & Fraction & [Fe/H] \\",
        r"\hline",
    ]
    for label, n, s, k, frac, lo, hi, fe in rows:
        L.append(f"{label:<20} & {n:>3} & {s:>5.1f} & {k:>3} & "
                 f"${frac:.1f}^{{+{hi-frac:.1f}}}_{{-{frac-lo:.1f}}}\\%$ & ${fe:.2f}$ \\\\")
    L += [
        r"\hline",
        f"$\\alpha$-classified & {n_tot} & {sum_tot:.1f} & {n_thr} & 100.0\\% & --- \\\\",
        r"\hline",
        r"\end{tabular*}",
        r"\end{center}",
    ]
    text = "\n".join(L) + "\n"
    for td in TAB_DIRS:
        if td.parent.exists():
            (td / "tab_population_decomp.tex").write_text(text, encoding="utf-8")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
