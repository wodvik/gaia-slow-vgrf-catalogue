# Gaia DR3 Online Verification Report

This audit checked the expanded v2 catalogue products against the public ESA
Gaia Archive TAP service, querying `gaiadr3.gaia_source` by `source_id`.

Endpoint used:

- <https://gea.esac.esa.int/tap-server/tap/sync>

Columns compared:

- `ra`
- `dec`
- `parallax`
- `pmra`
- `pmdec`
- `radial_velocity`
- `phot_g_mean_mag`

Tolerances:

- `ra`, `dec`, `parallax`, `pmra`, `pmdec`: `1e-9`
- `radial_velocity`, `phot_g_mean_mag`: `1e-6`

## Results

| Catalogue | Local rows | Unique local source IDs | Gaia Archive rows recovered | Missing IDs | Discrepancies |
|---|---:|---:|---:|---:|---:|
| `catalogue_expanded_tierAB.fits` | 517 | 517 | 517 | 0 | 0 |
| `catalogue_expanded_tierABC.fits` | 1,835 | 1,835 | 1,835 | 0 | 0 |
| `catalogue_expanded_master.fits` | 20,829 | 20,829 | 20,829 | 0 | 0 |

Maximum absolute deltas in the full master audit:

| Column | Maximum absolute delta |
|---|---:|
| `ra` | `1.1368683772161603e-13` |
| `dec` | `7.105427357601002e-15` |
| `parallax` | `7.105427357601002e-15` |
| `pmra` | `8.526512829121202e-14` |
| `pmdec` | `5.684341886080802e-14` |
| `radial_velocity` | `0.0` |
| `phot_g_mean_mag` | `0.0` |

## Output Files

- `verify_gaia_online.py`: reproducible verifier script.
- `catalogue_expanded_tierAB_online_verification_summary.json`
- `catalogue_expanded_tierABC_online_verification_summary.json`
- `catalogue_expanded_master_online_verification_summary.json`
- `*_missing_source_ids.csv`: header-only when no source IDs are missing.
- `*_discrepancies.csv`: header-only when no compared values disagree.

## Conclusion

All checked Gaia DR3 `source_id` values exist in the public ESA Gaia Archive,
and the local identifying astrometric, radial-velocity, and G-band photometric
values match the Gaia DR3 online records within the configured tolerances.
