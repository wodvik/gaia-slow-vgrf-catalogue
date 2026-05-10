"""Generate reviewer-facing AAS-style machine-readable catalogue tables.

The published FITS files remain the authoritative data products. These ASCII
tables provide journal-friendly machine-readable counterparts and compact
sample-table rows for the manuscript.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[2].parent
OUT = REPO / "release/v2/mrt"
OUT.mkdir(parents=True, exist_ok=True)

COLUMNS = [
    ("source_id", "none", "Gaia DR3 source identifier"),
    ("ra", "deg", "ICRS right ascension"),
    ("dec", "deg", "ICRS declination"),
    ("parallax", "mas", "Gaia DR3 parallax"),
    ("parallax_error", "mas", "Gaia DR3 parallax uncertainty"),
    ("parallax_zpcorr", "mas", "Parallax after Lindegren et al. zero-point correction"),
    ("dist_pc", "pc", "Adopted distance"),
    ("dist_lo_pc", "pc", "Lower distance credible bound"),
    ("dist_hi_pc", "pc", "Upper distance credible bound"),
    ("dist_source", "none", "Distance estimator used"),
    ("pmra", "mas/yr", "Proper motion in right ascension"),
    ("pmdec", "mas/yr", "Proper motion in declination"),
    ("radial_velocity", "km/s", "Gaia DR3 radial velocity"),
    ("radial_velocity_error", "km/s", "Radial-velocity uncertainty"),
    ("rv_quality", "none", "Internal Gaia DR3 RVS-quality class"),
    ("vgrf_default", "km/s", "Default Galactic rest-frame speed"),
    ("vgrf_grav22", "km/s", "Galactic rest-frame speed using GRAVITY+2022 solar parameters"),
    ("vgrf_lsr6", "km/s", "Galactic rest-frame speed using six-dimensional LSR convention"),
    ("vgrf_rb20", "km/s", "Galactic rest-frame speed using Reid & Brunthaler 2020 convention"),
    ("P_vgrf_below_25", "none", "Monte Carlo probability that Vgrf is below 25 km/s"),
    ("mc_realisations", "none", "Number of Monte Carlo realisations used for threshold probability"),
    ("tier", "none", "Velocity-threshold membership tier"),
]

FORMATS = {
    "source_id": "{:19d}",
    "ra": "{:12.7f}",
    "dec": "{:12.7f}",
    "parallax": "{:10.5f}",
    "parallax_error": "{:10.5f}",
    "parallax_zpcorr": "{:10.5f}",
    "dist_pc": "{:10.2f}",
    "dist_lo_pc": "{:10.2f}",
    "dist_hi_pc": "{:10.2f}",
    "dist_source": "{:<28s}",
    "pmra": "{:11.5f}",
    "pmdec": "{:11.5f}",
    "radial_velocity": "{:10.3f}",
    "radial_velocity_error": "{:10.3f}",
    "rv_quality": "{:<8s}",
    "vgrf_default": "{:10.3f}",
    "vgrf_grav22": "{:10.3f}",
    "vgrf_lsr6": "{:10.3f}",
    "vgrf_rb20": "{:10.3f}",
    "P_vgrf_below_25": "{:8.5f}",
    "mc_realisations": "{:6d}",
    "tier": "{:<1s}",
}

SAMPLE_COLUMNS = [
    "source_id",
    "ra",
    "dec",
    "vgrf_default",
    "P_vgrf_below_25",
    "mc_realisations",
    "tier",
    "rv_quality",
]


def _value(row, name):
    value = row[name]
    if isinstance(value, (bytes, bytearray)):
        value = value.decode("utf-8")
    if isinstance(value, np.generic):
        value = value.item()
    return value


def _write_mrt(path: Path, title: str, table: Table) -> None:
    widths = {name: len(FORMATS[name].format(_value(table[0], name)))
              for name, _, _ in COLUMNS}
    for name, _, _ in COLUMNS:
        widths[name] = max(widths[name], len(name))
        for row in table:
            widths[name] = max(widths[name], len(FORMATS[name].format(_value(row, name))))

    lines = [
        f"# {title}",
        "# ASCII machine-readable companion to the FITS catalogue.",
        "# Units marked 'none' are dimensionless identifiers, labels, or probabilities.",
        "#",
        "# Columns:",
    ]
    for idx, (name, unit, desc) in enumerate(COLUMNS, start=1):
        lines.append(f"# {idx:02d} {name} [{unit}] - {desc}")
    lines.append("#")
    lines.append(" ".join(name.ljust(widths[name]) for name, _, _ in COLUMNS))
    lines.append(" ".join("-" * widths[name] for name, _, _ in COLUMNS))
    for row in table:
        parts = []
        for name, _, _ in COLUMNS:
            text = FORMATS[name].format(_value(row, name))
            parts.append(text.rjust(widths[name]) if name not in {"dist_source", "rv_quality", "tier"} else text.ljust(widths[name]))
        lines.append(" ".join(parts))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_sample_tex(path: Path, table: Table) -> None:
    sample = table[np.array([_value(row, "rv_quality") == "ok" for row in table])][:8]
    lines = [
        r"\begin{deluxetable*}{rrrrcrll}",
        r"\tablecaption{Sample rows from the Tier~A+B headline catalogue. The full table is available in machine-readable form in the review package as \path{release/v2/mrt/catalogue_tierAB_mrt.txt}; the FITS version is \path{release/v2/phase1/catalogue_tierAB.fits}.\label{tab:headline_catalogue_sample}}",
        r"\tablehead{",
        r"\colhead{Gaia DR3 source\_id} & \colhead{R.A.} & \colhead{Decl.} & \colhead{$\vgrf$} & \colhead{$P(\vgrf<25)$} & \colhead{$N_{\rm MC}$} & \colhead{Tier} & \colhead{RVS QC}\\",
        r"\colhead{} & \colhead{(deg)} & \colhead{(deg)} & \colhead{(\kms)} & \colhead{} & \colhead{} & \colhead{} & \colhead{}",
        r"}",
        r"\startdata",
    ]
    for row in sample:
        lines.append(
            f"{int(row['source_id'])} & "
            f"{float(row['ra']):.5f} & "
            f"{float(row['dec']):.5f} & "
            f"{float(row['vgrf_default']):.2f} & "
            f"{float(row['P_vgrf_below_25']):.3f} & "
            f"{int(row['mc_realisations'])} & "
            f"{str(_value(row, 'tier'))} & "
            f"{str(_value(row, 'rv_quality'))} \\\\"
        )
    lines.extend([
        r"\enddata",
        r"\tablecomments{Only eight rows are shown for form and content. The electronic table contains all 334 Tier~A+B sources.}",
        r"\end{deluxetable*}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    tier_a = Table.read(REPO / "release/v2/phase1/catalogue_tierA.fits")
    tier_ab = Table.read(REPO / "release/v2/phase1/catalogue_tierAB.fits")
    tier_a.sort("source_id")
    tier_ab.sort("source_id")
    _write_mrt(OUT / "catalogue_tierA_mrt.txt", "Tier A slow-Vgrf catalogue", tier_a)
    _write_mrt(OUT / "catalogue_tierAB_mrt.txt", "Tier A+B headline slow-Vgrf catalogue", tier_ab)
    _write_sample_tex(OUT / "tab_headline_catalogue_sample.tex", tier_ab)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
