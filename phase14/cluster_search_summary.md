# Phase 14P compact-group search

Purpose: check whether the slow-Vgrf release catalogue contains an obvious
present-day open/globular-cluster-like aggregate.

Method: DBSCAN/FoF searches in 3D Galactocentric position, with optional
3D velocity coherence. Inputs are the master catalogue and
point-estimate orbit table; no Gaia archive rescan is required.

## Result

No candidate group was found under any configured compact-cluster search.

## Sample summary

| sample | N | min pair sep pc | pairs <25 pc | pairs <25 pc and <5 km/s | DBSCAN groups |
|---|---:|---:|---:|---:|---:|
| tierAB_headline | 517 | 38.9 | 0 | 0 | 0 |
| tierABC_statistical | 1835 | 20.2 | 3 | 0 | 0 |

Interpretation: an actual compact bound cluster in this catalogue would be
expected to produce multiple stars within tens of parsecs, usually also
with small relative velocities. The closest headline Tier A+B pair is
already about 39 pc apart, and no three-star compact clump appears when
the sample is broadened to Tier A+B+C.

Caveat: this is a present-day compactness check, not a full cluster
membership paper. A complete literature-grade test would also crossmatch
source IDs against Gaia DR3 open-cluster and globular-cluster membership
catalogues, then inspect CMD/isochrone consistency for any matches.
