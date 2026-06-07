# Reproducibility / release-validation targets (deep-review Issue 10).
PYTHON ?= python

.PHONY: validate-release checksums smoke provenance

## Run the full pre-submission validation suite.
validate-release:
	$(PYTHON) scripts/validate_release.py

## (Re)generate SHA-256 checksums for the released products.
checksums:
	$(PYTHON) scripts/make_release_checksums.py

## Product presence / schema / primary-count smoke test.
smoke:
	$(PYTHON) tests/smoke_regression.py --bundle-root .

## Tier-count consistency + stale-token provenance guard.
provenance:
	$(PYTHON) scripts/check_provenance.py
