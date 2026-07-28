# Phase 14P compact-group search

Purpose: check whether the slow-Vgrf release catalogue contains an obvious
present-day open/globular-cluster-like aggregate.

Method: DBSCAN/FoF searches in 3D Galactocentric position, with optional
3D velocity coherence. Inputs are the expanded master catalogue and
expanded point-estimate orbit table; no Gaia archive rescan is required.

## Result

No candidate group was found under any configured compact-cluster search.

## Sample summary

| sample | N | min pair sep pc | pairs <25 pc | pairs <25 pc and <5 km/s | DBSCAN groups |
|---|---:|---:|---:|---:|---:|
| tierAB_primary | 541 | 20.2 | 1 | 0 | 0 |
| tierABC_statistical | 1952 | 20.2 | 2 | 0 | 0 |

Interpretation: an actual compact bound cluster in this catalogue would be
expected to produce multiple stars within tens of parsecs, usually also
with small relative velocities. The closest primary Tier A+B pair is 20.2 pc apart; Tier A+B contains 1 pair below 25 pc, and Tier A+B+C contains 2 pairs below 25 pc.
No DBSCAN compact group appears in either sample.

Caveat: this is a present-day compactness check, not a full cluster
membership paper. A complete literature-grade test would also crossmatch
source IDs against Gaia DR3 open-cluster and globular-cluster membership
catalogues, then inspect CMD/isochrone consistency for any matches.
