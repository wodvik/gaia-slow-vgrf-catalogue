# Star counts — quick reference

This page is the plain-English version of `release/v2/numbers.yml`.
Every number a reader has to track for the Humble (2026) slow-Vgrf
catalogue is here, with the FITS file that provides it.

## Headline

> **The headline catalogue is Tier A+B: 334 stars with
> P(Vgrf<25 km/s) > 0.84.**

It is the empirical, conservative, paper-defining sample. Use
`catalogue_tierAB.fits`.

## The full ladder

| Sample | N | Selection | FITS file |
|---|---:|---|---|
| Propagated Gaia DR3 6D candidate pool | 2,859 | Vgrf-eligible 6D sources, pre-MC | (see master) |
| Corrected point-estimate Vgrf<25 km/s | 709 | After ZP+BJ corrections | (filter master) |
| Tier A | 214 | P(Vgrf<25) > 0.95 | `catalogue_tierA.fits` |
| **Tier A+B (HEADLINE)** | **334** | P > 0.84 | **`catalogue_tierAB.fits`** |
| Tier A+B+C | 632 | P > 0.50 | `catalogue_v2.fits` (master) |
| Tier D | 77 | Point-estimate <25 km/s, P ≤ 0.50 | (filter master) |

The master `catalogue_v2.fits` carries every row (all 2,859) plus a
`tier` column (A/B/C/D/X), so the convenience subsets are simple
filters; they are shipped pre-cut for casual reuse.

## Where each count is used in the paper

| Used for | Sample | N |
|---|---|---:|
| Headline catalogue, Gold subsets, abstract count | **Tier A+B** | **334** |
| Highest-confidence empirical anchor, robustness checks | Tier A | 214 |
| Orbit medians (e, R_peri, R_apo), MC propagation, population distributions | Tier A+B+C | 632 |
| Splash/GSE/Aurora/disc fractions (chemistry) | alpha-classified subset within Tier A+B+C | 53 |
| Selection-function context (full inverse weighting) | Tier A+B effective | 539 |
| Selection-function context (conservative bracket) | Tier A+B effective | 398 |
| Sgr A* candidate evaluation | Tier A+B+C | 4 |

## What the convenience subsets contain

`catalogue_tierA.fits` and `catalogue_tierAB.fits` carry the **same 26
columns** as the master (`source_id`, `ra`, `dec`, `parallax`,
`parallax_zpcorr`, `dist_pc`, `dist_lo_pc`, `dist_hi_pc`, `dist_source`,
`pmra`, `pmdec`, `radial_velocity`, four Vgrf flavours,
`P_vgrf_below_25`, `mc_realisations`, `tier`, etc.). Joining them back to
`catalogue_mc_orbits.fits`, `catalogue_v2_sf.fits`, `chemistry_v2.fits`,
or `populations_v2.fits` on `source_id` works without modification.

## What is NOT a per-star quantity

- The **539 / 398 selection-function effective counts** are
  population-scale weights, not per-star membership flags.
- The **53 alpha-classified subset** is bounded by spectroscopic coverage,
  not by tier; chemistry-region fractions are quoted only inside it.
- The **0–10 bridger range** is a sensitivity bracket across barred
  pattern speeds, not a robust population fraction.
