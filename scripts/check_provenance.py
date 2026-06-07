"""Pre-submission provenance guard (deep-review Issue 2/13).

Fails (exit 1) if any tier count in the manuscript, tables, sidecar JSONs, or
the adopted candidate pool disagrees with the authoritative counts read live
from the released tier catalogues, or if a known stale total (the pre-Vc=229
517 / 1,835 working counts -> 1{,}835, 1{,}318) survives in the manuscript or
tables. Run before every freeze/packaging:

    python scripts/check_provenance.py

The authoritative counts are the row counts of
catalogues/catalogue_expanded_tier{A,AB,ABC}.fits, so this check cannot drift
from the catalogue itself.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
ADOPTED = (289, 541, 1952)          # expected (A, A+B, A+B+C); also derived live
STALE_TEX = ["1{,}835", "1{,}318"]  # LaTeX-formatted stale totals
STALE_BARE = [r"\b1835\b", r"\b1318\b"]
STALE_PATH_PATTERNS = [
    "v1.0.4_review",
    r"C:\\Users\\humbl\\GAIA2026",
    "/mnt/c/Users/humbl/GAIA2026",
    "/release/v2/phase3_agama/",
]
STALE_TEXT_PATTERNS = [
    "Tier A+B (N=517)",
    "N=517",
    "No Tier A+B pair lies within 25",
    "38.9\\,pc",
    "38.9 pc",
    "three close pairs",
]


def n_rows(p: Path) -> int:
    return len(Table.read(p))


def rel(p: Path) -> str:
    try:
        return str(p.relative_to(REPO))
    except ValueError:
        return str(p)


def pdf_text(path: Path) -> str | None:
    exe = shutil.which("pdftotext")
    if exe is None or not path.exists():
        return None
    proc = subprocess.run(
        [exe, "-layout", str(path), "-"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="ignore",
    )
    if proc.returncode != 0:
        return None
    return proc.stdout


def main() -> int:
    failures: list[str] = []
    warnings: list[str] = []

    A = n_rows(BUNDLE / "catalogues" / "catalogue_expanded_tierA.fits")
    AB = n_rows(BUNDLE / "catalogues" / "catalogue_expanded_tierAB.fits")
    ABC = n_rows(BUNDLE / "catalogues" / "catalogue_expanded_tierABC.fits")
    B_inc, C_inc = AB - A, ABC - AB
    print(f"Authoritative catalogue counts: A={A} B_inc={B_inc} C_inc={C_inc} "
          f"A+B={AB} A+B+C={ABC}")
    if (A, AB, ABC) != ADOPTED:
        failures.append(f"catalogue row counts {(A, AB, ABC)} != adopted {ADOPTED}")

    # --- sidecar JSON consistency ---
    sidecar = BUNDLE / "phase14" / "expanded_context_figures_summary.json"
    if sidecar.exists():
        d = json.loads(sidecar.read_text())
        if d.get("n") != ABC:
            failures.append(f"{sidecar.name}: n={d.get('n')} != {ABC}")
        tc = d.get("tier_counts", {})
        if (tc.get("A"), tc.get("B"), tc.get("C")) != (A, B_inc, C_inc):
            failures.append(f"{sidecar.name}: tier_counts={tc} != ({A},{B_inc},{C_inc})")

    # --- adopted candidate pool tiering ---
    cand = BUNDLE / "private_inputs" / "expanded_candidates_mc_tiered.csv"
    if cand.exists():
        t = pd.read_csv(cand, usecols=["tier"]).tier.value_counts().to_dict()
        if (t.get("A"), t.get("B"), t.get("C")) != (A, B_inc, C_inc):
            failures.append(f"bundle candidate pool tiers (A={t.get('A')},B={t.get('B')},"
                            f"C={t.get('C')}) != ({A},{B_inc},{C_inc})")

    # --- stale tokens in manuscript + tables (both trees) ---
    tex_files = [REPO / "release" / "main.tex", BUNDLE / "main.tex"]
    tex_files += list((REPO / "release" / "tables" / "v15").glob("*.tex"))
    tex_files += list((BUNDLE / "tables" / "v15").glob("*.tex"))
    for f in sorted(set(tex_files)):
        if not f.exists():
            continue
        txt = f.read_text(encoding="utf-8", errors="ignore")
        for lit in STALE_TEX:
            if lit in txt:
                # 1,835 is also the legitimate complement 1,952-117 (Tier A+B+C
                # minus the alpha-classified subset); accept only when written
                # with that explicit arithmetic, otherwise flag as stale.
                if lit == "1{,}835" and "1{,}952-117=1{,}835" in txt:
                    continue
                i = txt.find(lit)
                failures.append(f"stale '{lit}' in {f.relative_to(REPO)}: "
                                f"...{txt[max(0, i-30):i+30].strip()}...")
        for pat in STALE_BARE:
            m = re.search(pat, txt)
            if m:
                failures.append(f"stale {pat} in {f.relative_to(REPO)}: "
                                f"...{txt[max(0, m.start()-30):m.start()+30].strip()}...")

    # --- targeted scalar/caption guards for reviewer-visible drift ---
    threshold_tables = [REPO / "release" / "tables" / "v15" / "tab_threshold.tex",
                        BUNDLE / "tables" / "v15" / "tab_threshold.tex"]
    for f in threshold_tables:
        if not f.exists():
            continue
        txt = f.read_text(encoding="utf-8", errors="ignore")
        if "$<25$  & 2{,}755 & 1{,}952" not in txt:
            failures.append(f"{f.relative_to(REPO)}: <25 threshold row must carry "
                            f"N_ABC={ABC} alongside point-count 2,755")

    cov_tables = [REPO / "release" / "tables" / "v15" / "tab_covariance_stress.tex",
                  BUNDLE / "tables" / "v15" / "tab_covariance_stress.tex"]
    cov_sidecar = BUNDLE / "phase14" / "expanded_covariance_stress_summary.json"
    cov = json.loads(cov_sidecar.read_text()) if cov_sidecar.exists() else {}
    remain = cov.get("headline_AB_adopted", AB) - cov.get("headline_AB_lost_under_copula", 0)
    for f in cov_tables:
        if not f.exists():
            continue
        txt = f.read_text(encoding="utf-8", errors="ignore")
        if "no primary Tier~A+B member is lost" in txt:
            failures.append(f"{f.relative_to(REPO)}: covariance caption contradicts "
                            "adopted-catalogue A+B losses")
        expected = f"{remain} of {AB} Tier~A+B members remain above"
        if expected not in txt:
            failures.append(f"{f.relative_to(REPO)}: covariance caption must state "
                            f"'{expected}'")

    # --- stale review-bundle paths and stale prose in distributed sidecars ---
    sidecar_files = list((BUNDLE / "phase14").rglob("*.json"))
    sidecar_files += list((BUNDLE / "phase14").rglob("*.md"))
    sidecar_files += list((BUNDLE / "potentials").glob("*.ini"))
    sidecar_files += [BUNDLE / "README.md", BUNDLE / "README-submission.md"]
    for f in sorted(set(sidecar_files)):
        if not f.exists():
            continue
        txt = f.read_text(encoding="utf-8", errors="ignore")
        for pat in STALE_PATH_PATTERNS:
            if pat in txt:
                i = txt.find(pat)
                failures.append(f"stale path '{pat}' in {rel(f)}: "
                                f"...{txt[max(0, i-40):i+50].strip()}...")
        for pat in STALE_TEXT_PATTERNS:
            if pat in txt:
                i = txt.find(pat)
                failures.append(f"stale prose '{pat}' in {rel(f)}: "
                                f"...{txt[max(0, i-40):i+60].strip()}...")

    # --- rendered-PDF stale text check ---
    pdf_checks = [
        (BUNDLE / "main.pdf", ["N=517", "38.9 pc", "three close pairs"]),
        (BUNDLE / "figures" / "fig_sf_map.pdf", ["Tier A+B (N=517)", "N=517"]),
    ]
    for pdf, pats in pdf_checks:
        txt = pdf_text(pdf)
        if txt is None:
            warnings.append(f"could not extract text from {rel(pdf)}; skipped rendered stale-text check")
            continue
        for pat in pats:
            if pat in txt:
                failures.append(f"stale rendered text '{pat}' in {rel(pdf)}")

    # --- informational: stale scratch candidate pool on D: ---
    dpath = Path("D:/GAIA/parent_scan/expanded_candidates_mc_tiered.csv")
    if dpath.exists():
        warnings.append("D:/GAIA/parent_scan/expanded_candidates_mc_tiered.csv carries the "
                        "pre-Vc=229 tiering (517/1835); config.yml production path points here. "
                        "The adopted pool is bundle private_inputs/. Do not feed the D: copy to reruns.")

    print("\n=== PROVENANCE CHECK ===")
    for w in warnings:
        print("WARN:", w)
    if failures:
        for fl in failures:
            print("FAIL:", fl)
        print(f"\n{len(failures)} failure(s) -- provenance NOT clean.")
        return 1
    print("PASS: tier counts consistent across catalogues, sidecars, candidate pool, "
          "manuscript, tables, path sidecars, and rendered PDF text; no stale tokens.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
