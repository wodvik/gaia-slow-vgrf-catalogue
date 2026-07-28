"""Generate AAS/CDS Machine-Readable Tables (MRT) for the review bundle.

The published FITS files remain the authoritative data products. These ASCII
tables are the journal-/archive-standard machine-readable counterparts,
written in the AAS/CDS byte-by-byte MRT format via astropy so they are fully
regenerable on every release (no hand-authored byte ranges) and handle
missing values correctly (sparse chemistry/action columns appear as blank
fixed-width fields flagged "?" in the byte-by-byte description).

Products:
  mrt/catalogue_tierA_mrt.txt     kinematics + chemistry (Tier A)
  mrt/catalogue_tierAB_mrt.txt    kinematics + chemistry (Tier A+B, primary)
  mrt/catalogue_tierABC_mrt.txt   kinematics + chemistry (Tier A+B+C)
  mrt/catalogue_orbits_tierABC_mrt.txt  orbit summaries + actions + chemistry
  mrt/tab_headline_catalogue_sample.tex  manuscript 8-row sample (also -> tables/v15/)

Run from the bundle root:
    python scripts/make_mrt_tables.py
"""
from __future__ import annotations

import argparse
import io
from pathlib import Path

import numpy as np
from astropy.io import ascii
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
AUTHORS = "Humble R."

# (column, CDS unit string or None for dimensionless, description).
# dex abundances carry no CDS unit (CDS has no 'dex'); 'dex' is in the text.
TIER_COLUMNS = [
    ("source_id", None, "Gaia DR3 source identifier"),
    ("ra", "deg", "ICRS right ascension"),
    ("dec", "deg", "ICRS declination"),
    ("parallax", "mas", "Gaia DR3 parallax"),
    ("parallax_error", "mas", "Gaia DR3 parallax uncertainty"),
    ("parallax_zpcorr", "mas", "Parallax after Lindegren et al. (2021) zero-point correction"),
    ("dist_pc", "pc", "Adopted distance"),
    ("dist_lo_pc", "pc", "Lower distance credible bound"),
    ("dist_hi_pc", "pc", "Upper distance credible bound"),
    ("dist_source", None, "Distance estimator used"),
    ("pmra", "mas/yr", "Proper motion in right ascension (incl. cos dec)"),
    ("pmdec", "mas/yr", "Proper motion in declination"),
    ("radial_velocity", "km/s", "Gaia DR3 radial velocity"),
    ("radial_velocity_error", "km/s", "Radial-velocity uncertainty"),
    ("rv_quality", None, "Internal Gaia DR3 RVS-quality class (ok/marginal/poor)"),
    ("rvs_quality_ok", None, "True when rv_quality is ok"),
    ("vgrf_default", "km/s", "Adopted Galactic rest-frame speed"),
    ("P_vgrf_below_25_popprior", None,
     "Adopted probability that Vgrf < 25 km/s, after the population-prior correction; defines tier"),
    ("P_vgrf_below_25_forward", None,
     "Forward Monte Carlo probability that Vgrf < 25 km/s, before the population-prior correction"),
    ("mc_realisations", None, "Number of Monte Carlo realisations for the threshold probability"),
    ("tier", None, "Adopted membership tier from P_vgrf_below_25_popprior (A/B/C)"),
    ("tier_forward", None, "Membership tier from P_vgrf_below_25_forward (A/B/C/D/X)"),
    ("mh_gspphot", None, "Gaia GSP-Phot photometric metallicity [M/H] in dex (contextual; blank if none)"),
    ("feh_spec", None, "Spectroscopic [Fe/H] in dex from APOGEE DR17/GALAH DR3 (blank if no match)"),
    ("feh_spec_err", None, "Uncertainty on feh_spec in dex (blank if none)"),
    ("alpha_spec", None, "Spectroscopic alpha proxy in dex ([a/M] APOGEE, [a/Fe] GALAH; blank if none)"),
    ("chem_survey", None, "Spectroscopic source of feh_spec/alpha_spec (APOGEE/GALAH; blank if none)"),
    ("chem_population", None, "Chemodynamic class from [Fe/H]+alpha (Splash/GSE/Aurora/disk/unclassified; blank if none)"),
    ("nss_two_body", None, "Gaia DR3 nss_two_body_orbit non-single-star flag (True/False)"),
    ("nss_solution_type", None, "Gaia DR3 NSS solution type (SB1/Orbital/AstroSpectroSB1; blank if not matched)"),
]

ORBIT_COLUMNS = [
    ("source_id", None, "Gaia DR3 source identifier"),
    ("tier", None, "Adopted membership tier after the population-prior correction (A/B/C)"),
    ("tier_forward", None, "Membership tier from the uncorrected forward score (A/B/C)"),
    ("P_vgrf_below_25", None, "Forward Monte Carlo probability that Vgrf < 25 km/s"),
    ("static_R_peri_kpc", "kpc", "Pericentre radius in the static Hunter+2024 potential"),
    ("static_R_apo_kpc", "kpc", "Apocentre radius in the static potential"),
    ("static_z_max_kpc", "kpc", "Maximum absolute height in the static potential"),
    ("static_ecc", None, "Orbital eccentricity (R_apo-R_peri)/(R_apo+R_peri), static potential"),
    ("static_Lz_kpc_kms", "kpc km/s", "Mean z-angular momentum L_z in the static potential"),
    ("J_R", "kpc km/s", "Radial action (single-evaluation Staeckel fudge); caveated, see action_reliability_flag"),
    ("J_z", "kpc km/s", "Vertical action (single-evaluation Staeckel fudge); least reliable for these orbits"),
    ("J_phi", "kpc km/s", "Azimuthal action (single-evaluation Staeckel fudge)"),
    ("J_R_timeavg", "kpc km/s", "Convergent radial action (orbit-time-averaged over 4 Gyr)"),
    ("J_z_timeavg", "kpc km/s", "Convergent vertical action (orbit-time-averaged over 4 Gyr)"),
    ("J_phi_timeavg", "kpc km/s", "Convergent azimuthal action (orbit-time-averaged over 4 Gyr)"),
    ("action_max_fracdiff", None, "Max |single-eval - time-averaged|/max(|time-avg|,1) over J_R,J_z,J_phi"),
    ("action_reliability_flag", None, "Per-star action reliability (sampled_ok/caution/poor)"),
    ("mh_gspphot", None, "Gaia GSP-Phot photometric metallicity [M/H] in dex (blank if none)"),
    ("feh_spec", None, "Spectroscopic [Fe/H] in dex from APOGEE DR17/GALAH DR3 (blank if no match)"),
    ("alpha_spec", None, "Spectroscopic alpha proxy in dex (blank if none)"),
    ("chem_survey", None, "Spectroscopic source (APOGEE/GALAH; blank if none)"),
    ("chem_population", None, "Chemodynamic class (Splash/GSE/Aurora/disk/unclassified; blank if none)"),
    ("nss_two_body", None, "Gaia DR3 nss_two_body_orbit non-single-star flag (True/False)"),
    ("nss_solution_type", None, "Gaia DR3 NSS solution type (SB1/Orbital/AstroSpectroSB1; blank if not matched)"),
]


def write_mrt(fits_path: Path, out_path: Path, title: str, table_desc: str, colspec) -> int:
    t = Table.read(fits_path)
    t.sort("source_id")
    sub = t[[c for c, _, _ in colspec]]
    # The MRT writer cannot format boolean columns; render them as True/False.
    for name in list(sub.colnames):
        if sub[name].dtype == bool:
            sub[name] = np.array(["True" if bool(v) else "False" for v in sub[name]])
    for name, unit, desc in colspec:
        sub[name].unit = unit if unit else None
        sub[name].description = desc
    buf = io.StringIO()
    ascii.write(sub, buf, format="mrt")
    lines = buf.getvalue().splitlines()
    for i, line in enumerate(lines[:6]):
        if line == "Title:":
            lines[i] = "Title: " + title
        elif line == "Authors:":
            lines[i] = "Authors: " + AUTHORS
        elif line == "Table:":
            lines[i] = "Table: " + table_desc
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(sub)


def _value(row, name):
    value = row[name]
    if isinstance(value, (bytes, bytearray)):
        value = value.decode("utf-8")
    if isinstance(value, np.generic):
        value = value.item()
    return value


def write_sample_tex(path: Path, table: Table) -> None:
    table = table.copy()
    table.sort("source_id")
    sample = table[np.array([_value(row, "rv_quality") == "ok" for row in table])][:8]
    lines = [
        r"\begin{deluxetable*}{rrrrcrll}",
        r"\tablecaption{Sample rows from the Tier~A+B primary catalogue. The full table is available in machine-readable form in the review package as \path{mrt/catalogue_tierAB_mrt.txt}; the FITS version is \path{catalogues/catalogue_expanded_tierAB.fits}.\label{tab:headline_catalogue_sample}}",
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
        rf"\tablecomments{{Only eight rows are shown for form and content. The electronic table contains all {len(table)} Tier~A+B sources.}}",
        r"\end{deluxetable*}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-root", type=Path, default=BUNDLE)
    parser.add_argument("--forward", action="store_true",
                        help="build from the uncorrected forward-score tiers "
                             "instead of the adopted population-prior ones")
    args = parser.parse_args()
    bundle = args.bundle_root
    out = bundle / "mrt"
    out.mkdir(parents=True, exist_ok=True)
    table_out = bundle / "tables" / "v15"
    table_out.mkdir(parents=True, exist_ok=True)
    cat = bundle / "catalogues"
    # The population-prior tiers are the adopted catalogue (Phase 16F); the
    # forward-score products remain available behind --forward.
    stem = "catalogue_expanded" if args.forward else "catalogue_retier"
    tierdef = ("forward Monte Carlo score" if args.forward
               else "population-prior corrected score")

    jobs = [
        (f"{stem}_tierA.fits", "catalogue_tierA_mrt.txt",
         "Gaia DR3 slow-Vgrf catalogue, Tier A",
         "Tier A (P(Vgrf<25 km/s)>0.95) kinematics and chemistry", TIER_COLUMNS),
        (f"{stem}_tierAB.fits", "catalogue_tierAB_mrt.txt",
         "Gaia DR3 slow-Vgrf catalogue, Tier A+B (primary)",
         "Tier A+B (P>0.84) primary kinematics and chemistry", TIER_COLUMNS),
        (f"{stem}_tierABC.fits", "catalogue_tierABC_mrt.txt",
         "Gaia DR3 slow-Vgrf catalogue, Tier A+B+C",
         "Tier A+B+C (P>0.50) kinematics and chemistry", TIER_COLUMNS),
        (f"{stem}_orbits_tierABC.fits", "catalogue_orbits_tierABC_mrt.txt",
         "Gaia DR3 slow-Vgrf catalogue, Tier A+B+C orbit summaries and actions",
         "Tier A+B+C point-estimate orbit summaries, actions, and chemistry", ORBIT_COLUMNS),
    ]
    for fits_name, out_name, title, desc, spec in jobs:
        desc = f"{desc}. Tier membership uses the {tierdef}"
        n = write_mrt(cat / fits_name, out / out_name, title, desc, spec)
        print(f"wrote mrt/{out_name}  ({n} rows, {len(spec)} columns)")

    tier_ab = Table.read(cat / f"{stem}_tierAB.fits")
    write_sample_tex(out / "tab_headline_catalogue_sample.tex", tier_ab)
    write_sample_tex(table_out / "tab_headline_catalogue_sample.tex", tier_ab)
    print("wrote tab_headline_catalogue_sample.tex (mrt/ and tables/v15/)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
