"""Phase 5E — XP metallicity recompute via VizieR (Andrae+2023, Zhang+2023).

Cross-matches v2 Tier A+B+C (632 stars) against:
  - Andrae+2023 GSP-Phot stellar parameters from XP spectra (VizieR I/360)
  - Zhang+2023 stellar parameters from XP spectra (VizieR I/362)

Outputs the pooled (Andrae OR Zhang) [Fe/H] estimate per Tier A+B+C star
and tier-resolved [Fe/H] distributions.

Outputs
-------
release/v2/phase5/xp_metallicity_v2.fits
release/v2/phase5/gate5E_xp.json
release/v2/phase5/gate5E_xp_FeH.png
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import pyvo
from astropy.table import Table, vstack

REPO = Path(__file__).resolve().parents[2].parent
CONFIG = yaml.safe_load((REPO / "release/v2/config.yml").read_text())
OUT = REPO / "release/v2/phase5"
OUT.mkdir(parents=True, exist_ok=True)
CACHE = REPO / CONFIG["compute"]["cache_dir"]
CACHE.mkdir(parents=True, exist_ok=True)


def log(m): print(f"[5E t={time.time()-T0:5.1f}s] {m}", flush=True)


def vizier_xmatch(table_name: str, columns: list, source_ids: np.ndarray,
                  cache_path: Path, source_col: str = "Source",
                  chunk: int = 200) -> pd.DataFrame:
    if cache_path.exists():
        log(f"  using cache {cache_path.name}")
        return Table.read(cache_path).to_pandas()
    tap = pyvo.dal.TAPService("https://tapvizier.cds.unistra.fr/TAPVizieR/tap")
    chunks = []
    cols = ", ".join(f"t.{c}" for c in columns)
    for i in range(0, len(source_ids), chunk):
        ids = source_ids[i:i + chunk]
        in_clause = ",".join(str(int(x)) for x in ids)
        adql = f'SELECT {cols} FROM "{table_name}" AS t WHERE t.{source_col} IN ({in_clause})'
        log(f"  {table_name} chunk {i}-{i+len(ids)}")
        try:
            res = tap.search(adql).to_table()
        except Exception as e:
            log(f"    chunk failed: {e}")
            continue
        chunks.append(res)
    out = vstack(chunks) if chunks else Table()
    if len(out):
        out.write(cache_path, format="fits", overwrite=True)
    log(f"  cached {len(out)} rows -> {cache_path.name}")
    return out.to_pandas() if len(out) else pd.DataFrame()


def main():
    global T0; T0 = time.time()

    # Load Tier A+B+C source_ids
    v2 = Table.read(REPO / "release/v2/phase1/catalogue_v2.fits").to_pandas()
    if v2["tier"].dtype == object and isinstance(v2["tier"].iloc[0], (bytes, bytearray)):
        v2["tier"] = v2["tier"].str.decode("utf-8")
    tABC = v2["tier"].isin(["A", "B", "C"])
    ids = v2.loc[tABC, "source_id"].astype(int).to_numpy()
    log(f"querying XP catalogues for {len(ids)} Tier A+B+C stars")

    # Andrae+2023 GSP-Phot XP table on Vizier: I/360/gspphotxp
    # Available cols: Source, MH (metallicity), Teff, logg, etc.
    andrae = vizier_xmatch(
        "I/360/gspphotxp",
        ["Source", "MH", "Teff", "logg"],
        ids, CACHE / "andrae23_gspphotxp.fits",
        source_col="Source",
    )
    if not andrae.empty:
        andrae = andrae.rename(columns={"Source": "source_id",
                                         "MH": "MH_andrae",
                                         "Teff": "Teff_andrae",
                                         "logg": "logg_andrae"})
        log(f"Andrae+2023: {len(andrae)} matches")
    else:
        log("Andrae+2023: no matches (table may not exist or query failed)")
        andrae = pd.DataFrame(columns=["source_id"])

    # Zhang+2023 XP table: I/362/csv  (table name uncertain across mirrors)
    zhang = vizier_xmatch(
        "I/362/params",
        ["Source", "MH", "Teff", "logg"],
        ids, CACHE / "zhang23_params.fits",
        source_col="Source",
    )
    if not zhang.empty:
        zhang = zhang.rename(columns={"Source": "source_id",
                                       "MH": "MH_zhang",
                                       "Teff": "Teff_zhang",
                                       "logg": "logg_zhang"})
        log(f"Zhang+2023: {len(zhang)} matches")
    else:
        log("Zhang+2023: no matches (table may not exist or query failed)")
        zhang = pd.DataFrame(columns=["source_id"])

    # Merge to v2 base
    v2_subset = v2.loc[tABC, ["source_id", "tier", "P_vgrf_below_25",
                              "vgrf_default"]].copy()
    out = v2_subset.merge(andrae, on="source_id", how="left").merge(
        zhang, on="source_id", how="left")

    # Pool: prefer Andrae if both present (smaller uncertainties typically)
    out["MH_xp"] = out["MH_andrae"].where(out["MH_andrae"].notna(),
                                           out["MH_zhang"])
    out["MH_source"] = np.where(out["MH_andrae"].notna(), "andrae",
                       np.where(out["MH_zhang"].notna(), "zhang", "none"))

    Table.from_pandas(out).write(OUT / "xp_metallicity_v2.fits", overwrite=True)
    log(f"wrote xp_metallicity_v2.fits ({len(out)} rows; "
        f"{out['MH_xp'].notna().sum()} have XP MH)")

    # Tier-resolved [Fe/H] (XP) summary
    def stat(arr):
        a = np.asarray(arr, dtype=float); a = a[np.isfinite(a)]
        return ({"n": int(len(a)), "p16": float(np.percentile(a, 16)),
                 "p50": float(np.percentile(a, 50)),
                 "p84": float(np.percentile(a, 84))} if len(a) else None)

    summary = {"n_total_TierABC": int(tABC.sum()),
               "n_with_MH_xp":     int(out["MH_xp"].notna().sum()),
               "n_andrae":         int(out["MH_andrae"].notna().sum()) if "MH_andrae" in out else 0,
               "n_zhang":          int(out["MH_zhang"].notna().sum()) if "MH_zhang" in out else 0,
               "by_tier_xp_MH": {}}
    for tier_set, lab in [(["A", "B"], "tier_AB"),
                          (["A", "B", "C"], "tier_ABC")]:
        m = out["tier"].isin(tier_set)
        summary["by_tier_xp_MH"][lab] = {"n": int(m.sum()),
                                          "MH": stat(out.loc[m, "MH_xp"])}
    (OUT / "gate5E_xp.json").write_text(json.dumps(summary, indent=2))
    log("wrote gate5E_xp.json")
    print(json.dumps(summary, indent=2))

    # Plot histogram of XP MH per tier
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for tier, color in [("A","C3"), ("B","C1"), ("C","C2")]:
        m = out["tier"] == tier
        v = out.loc[m, "MH_xp"].dropna().to_numpy()
        if len(v) > 5:
            ax.hist(v, bins=np.linspace(-3, 0.5, 36), alpha=0.55,
                    label=f"Tier {tier} (n={len(v)})", color=color)
    ax.set_xlabel("[Fe/H] from XP (Andrae+2023 / Zhang+2023 pooled)")
    ax.set_ylabel("count")
    ax.set_title(f"Phase 5E — XP [Fe/H] for {out['MH_xp'].notna().sum()} Tier A+B+C stars")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "gate5E_xp_FeH.png", dpi=140); plt.close(fig)
    log(f"DONE in {time.time()-T0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
