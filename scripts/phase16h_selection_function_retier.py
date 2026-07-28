"""Phase 16H -- selection-function diagnostics and sky map on the retiered tiers.

The GaiaUnlimited DR3-RVS query itself does not need repeating: Phase 14Z already
stored the per-source selection-function value, weight, and parent count for all
20,829 candidates in phase14/expanded_selection_function.fits. Because the
population-prior tiers are strict nested subsets of the forward-defined ones,
every retiered star already has its selection-function entry, so this pass only
re-aggregates that table over the new tier definition and redraws the map.

Outputs: phase14/expanded_selection_function_summary_retier.json
         tables/v15/tab_selection_function.tex
         figures/fig_sf_map.pdf
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from astropy.table import Table

BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
SF = BUNDLE / "phase14" / "expanded_selection_function.fits"
RETIER = BUNDLE / "catalogues" / "catalogue_retier_master.fits"
OUT = BUNDLE / "phase14"
FIG_DIRS = [BUNDLE / "figures", REPO / "release" / "figures"]
TAB_DIRS = [BUNDLE / "tables" / "v15", REPO / "release" / "tables" / "v15"]

LOW_PARENT = 10
DL = 5.0   # degrees, sky-map pixel size
DB = 5.0

TABLE_HEAD = r"""\begin{center}
\refstepcounter{table}\label{tab:selection_function}
\begin{minipage}{\columnwidth}
\small\textbf{Table \thetable.} Gaia DR3 RVS selection-function
diagnostics for the catalogue.  The catalogue is selected from the
observed Gaia DR3 6D/RVS subset and is therefore not a volume-complete
Milky Way census.  $N_{\rm valid}$ is the number of sources for which
the GaiaUnlimited DR3-RVS implementation
\citep{CastroGinard2023} returned a finite selection-function weight at the source's $(G, BP{-}RP, l, b)$ cell;
the $N-N_{\rm valid}$ deficit is dominated by cells with zero parent
count.  $f_{n<10}$ is the fraction of $N_{\rm valid}$ sources whose
underlying GaiaUnlimited cell carries a parent count below 10 (the
regime in which the Beta$(k{+}1,\,n{-}k{+}1)$ posterior remains
competitive with its Beta$(1,1)$ prior: posterior standard error
$\gtrsim 0.16$ for $n<10$); $f_{\sum,n<10}$ is the share of
$\sum 1/S_{\rm RVS}$ contributed by those prior-dominated cells.  Tiers
are the adopted population-prior ones.  The
lower block restricts each tier to the inner-Galaxy quadrant
$l \in [330\degr,\,30\degr]$ where the catalogue concentrates
(Section~\ref{sec:sky}).  Inverse-selection weighted sums are reported
only as contextual observability diagnostics and are not used as
primary population counts.
\end{minipage}
\vspace{0.4ex}

\footnotesize
\setlength{\tabcolsep}{1pt}
\begin{tabular}{@{}lrrrrrr@{}}
\hline\hline
Sample & $N$ & $N_{\rm val}$ & $\Sigma 1/S$ & $\tilde S$ & $f_{<10}$ & $f_{\Sigma,<10}$ \\
\hline
"""


def clean(a) -> np.ndarray:
    a = np.asarray(a)
    if a.dtype.kind in "SO":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a]).astype(str)
    return a.astype(str)


def stats(sel: np.ndarray, valid: np.ndarray, w: np.ndarray, low: np.ndarray,
          sf_value: np.ndarray) -> dict:
    n = int(sel.sum())
    v = sel & valid
    nv = int(v.sum())
    sum_w = float(np.nansum(w[v]))
    lowv = v & low
    sum_low = float(np.nansum(w[lowv]))
    return {
        "N": n, "N_valid": nv,
        "sum_inv_S": sum_w,
        "median_S": float(np.nanmedian(sf_value[v])) if nv else float("nan"),
        "f_n_lt10_pct": 100.0 * int(lowv.sum()) / nv if nv else float("nan"),
        "n_lt10_count": int(lowv.sum()),
        "f_sum_n_lt10_pct": 100.0 * sum_low / sum_w if sum_w else float("nan"),
    }


def fmt_n(n: int) -> str:
    return f"{n:,}".replace(",", "{,}")


def main() -> None:
    t = Table.read(SF)
    sid = np.asarray(t["source_id"]).astype(np.int64)
    l = np.asarray(t["l"], dtype=float)
    b = np.asarray(t["b"], dtype=float)
    sf_value = np.asarray(t["sf_value"], dtype=float)
    w = np.asarray(t["sf_weight"], dtype=float)
    low = np.asarray(t["sf_prior_dominated_n_lt10"]).astype(bool)
    invalid = np.asarray(t["sf_invalid"]).astype(bool)
    valid = ~invalid & np.isfinite(w)

    # Population-prior tier labels, joined by source_id.
    r = Table.read(RETIER)
    rt = dict(zip(np.asarray(r["source_id"]).astype(np.int64).tolist(),
                  clean(r["tier"]).tolist()))
    tier = np.array([rt.get(int(s), "X") for s in sid])

    isA, isB, isC = tier == "A", tier == "B", tier == "C"
    isAB, isABC = isA | isB, isA | isB | isC
    pool = np.ones(len(t), bool)

    # Inner-Galaxy quadrant l in [330, 30]
    inner = (l >= 330.0) | (l <= 30.0)

    samples = [
        ("Candidate pool", pool, False), ("Tier A", isA, False), ("Tier B", isB, False),
        ("Tier C", isC, False), ("Tier A+B", isAB, False), ("Tier A+B+C", isABC, False),
        ("Tier A (inner)", isA & inner, True), ("Tier A+B (inner)", isAB & inner, True),
        ("Tier A+B+C (inner)", isABC & inner, True),
    ]
    summary = {"source_table": SF.name, "tier_definition": "population-prior (Phase 16F)",
               "low_parent_count_threshold": LOW_PARENT, "samples": {}}
    rows = []
    for name, sel, is_inner in samples:
        s = stats(sel, valid, w, low, sf_value)
        summary["samples"][name] = s
        rows.append((name, s, is_inner))

    (OUT / "expanded_selection_function_summary_retier.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    # ---- table ----
    # The header is a literal template. It must NOT be recovered by re-reading
    # the output file: doing so makes the script append a fresh body block on
    # every run, silently duplicating the table.
    lines = [TABLE_HEAD.rstrip("\n"),
             r"\multicolumn{7}{l}{\textit{Full sample}} \\"]
    for name, s, is_inner in rows:
        if is_inner:
            continue
        f10 = "--" if name == "Candidate pool" else f"{s['f_n_lt10_pct']:.1f}\\%"
        fs10 = "--" if name == "Candidate pool" else f"{s['f_sum_n_lt10_pct']:.1f}\\%"
        lines.append(f"{name:<23} & {fmt_n(s['N']):<8} & {fmt_n(s['N_valid']):<8} & "
                     f"{s['sum_inv_S']:,.1f}".replace(",", "{,}") +
                     f" & {s['median_S']:.2f} & {f10:<7} & {fs10:<7} \\\\")
    lines.append(r"\hline")
    lines.append(r"\multicolumn{7}{l}{\textit{Inner-Galaxy quadrant, $l\in[330\degr,\,30\degr]$}} \\")
    for name, s, is_inner in rows:
        if not is_inner:
            continue
        lines.append(f"{name:<23} & {fmt_n(s['N']):<8} & {fmt_n(s['N_valid']):<8} & "
                     f"{s['sum_inv_S']:,.1f}".replace(",", "{,}") +
                     f" & {s['median_S']:.2f} & {s['f_n_lt10_pct']:.1f}\\% & "
                     f"{s['f_sum_n_lt10_pct']:.1f}\\% \\\\")
    lines += [r"\hline", r"\end{tabular}", r"\end{center}"]
    text = "\n".join(lines) + "\n"
    for td in TAB_DIRS:
        if td.parent.exists():
            td.mkdir(parents=True, exist_ok=True)
            (td / "tab_selection_function.tex").write_text(text, encoding="utf-8")

    # ---- sky map ----
    # l is plotted increasing to the LEFT with l=0 centred, matching the released
    # figure: wrap to [-180, 180] and invert the axis.
    lw = np.where(l > 180.0, l - 360.0, l)
    le = np.arange(-180.0, 180.0 + DL, DL)
    be = np.arange(-90.0, 90.0 + DB, DB)

    vpool = pool & valid
    li = np.digitize(lw[vpool], le) - 1
    bi = np.digitize(b[vpool], be) - 1
    logw = np.log10(np.maximum(w[vpool], 1e-9))
    grid = np.full((len(be) - 1, len(le) - 1), np.nan)
    acc = np.zeros_like(grid)
    cnt = np.zeros_like(grid)
    ok = (li >= 0) & (li < len(le) - 1) & (bi >= 0) & (bi < len(be) - 1)
    np.add.at(acc, (bi[ok], li[ok]), logw[ok])
    np.add.at(cnt, (bi[ok], li[ok]), 1.0)
    grid = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)

    # Hatch pixels where >= 50% of catalogue (Tier A+B+C) stars are prior-dominated.
    cat = isABC & valid
    cli = np.digitize(lw[cat], le) - 1
    cbi = np.digitize(b[cat], be) - 1
    cok = (cli >= 0) & (cli < len(le) - 1) & (cbi >= 0) & (cbi < len(be) - 1)
    tot = np.zeros_like(grid)
    lown = np.zeros_like(grid)
    np.add.at(tot, (cbi[cok], cli[cok]), 1.0)
    np.add.at(lown, (cbi[cok], cli[cok]), low[cat][cok].astype(float))
    hatch = (tot > 0) & (lown / np.maximum(tot, 1) >= 0.5)

    fig, ax = plt.subplots(figsize=(5.77, 3.43))
    im = ax.pcolormesh(le, be, grid, cmap="viridis", shading="auto")
    # Hatch prior-dominated pixels individually. contourf would interpolate
    # between grid centres and smear the hatching across unflagged cells.
    matplotlib.rcParams["hatch.linewidth"] = 0.45
    for j, i in zip(*np.nonzero(hatch)):
        ax.add_patch(Rectangle((le[i], be[j]), DL, DB, fill=False,
                               hatch="/////", edgecolor="white",
                               linewidth=0.0, zorder=3))
    nab = int(isAB.sum())
    ax.scatter(lw[isAB], b[isAB], s=4, c="#d62728", alpha=0.85,
               linewidths=0, label=f"Tier A+B (N={nab})", zorder=5)
    ax.set_xlim(180, -180)
    ax.set_ylim(-90, 90)
    ax.set_xticks([-150, -100, -50, 0, 50, 100, 150])
    ax.set_xticklabels(["150", "100", "50", "0", "50", "100", "150"])
    ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80])
    ax.set_yticklabels(["80", "60", "40", "20", "0", "20", "40", "60", "80"])
    ax.set_xlabel("Galactic longitude $l$ (deg)")
    ax.set_ylabel("Galactic latitude $b$ (deg)")
    ax.set_title("GaiaUnlimited DR3-RVS selection-function weights", fontsize=9)
    ax.legend(loc="upper right", fontsize=6.5, framealpha=0.9, markerscale=2)
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label(r"$\log_{10}\,\langle 1/S_{\rm RVS}\rangle$", fontsize=8)
    fig.tight_layout()
    for fd in FIG_DIRS:
        if fd.parent.exists():
            fd.mkdir(parents=True, exist_ok=True)
            fig.savefig(fd / "fig_sf_map.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)

    for name, s, _ in rows:
        print(f"{name:<22} N={s['N']:>6}  Nval={s['N_valid']:>6}  "
              f"sum1/S={s['sum_inv_S']:>9.1f}  f<10={s['f_n_lt10_pct']:>5.1f}%  "
              f"fsum<10={s['f_sum_n_lt10_pct']:>5.1f}%")


if __name__ == "__main__":
    main()
