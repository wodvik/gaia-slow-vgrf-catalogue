"""Run the full pre-submission release validation (deep-review Issue 10).

Chains the bundle's integrity checks into one command:
  1. tests/smoke_regression.py   -- product presence, schema, primary counts
  2. scripts/check_provenance.py -- tier-count consistency + stale-token scan
  3. scripts/make_release_checksums.py --verify -- release-product integrity

    python scripts/validate_release.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

BUNDLE = Path(__file__).resolve().parents[1]
STEPS = [
    [sys.executable, "tests/smoke_regression.py", "--bundle-root", "."],
    [sys.executable, "scripts/check_provenance.py"],
    [sys.executable, "scripts/make_release_checksums.py", "--verify"],
]


def main() -> int:
    failures = 0
    for step in STEPS:
        print(f"\n==> {' '.join(step[1:])}", flush=True)
        rc = subprocess.run(step, cwd=BUNDLE).returncode
        if rc != 0:
            failures += 1
            print(f"  STEP FAILED (rc={rc})")
    print("\n" + ("=" * 40))
    if failures:
        print(f"validate-release: {failures} step(s) FAILED")
        return 1
    print("validate-release: ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
