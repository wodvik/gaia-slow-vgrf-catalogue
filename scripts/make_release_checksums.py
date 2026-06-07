"""Generate or verify SHA-256 checksums for the released products.

Referee response (deep-review Issue 10). Writes release_checksums.sha256 over
the figures, tables, catalogues, and MRT products so the frozen release can be
integrity-checked.

    python scripts/make_release_checksums.py            # write checksums
    python scripts/make_release_checksums.py --verify    # verify against file
    python scripts/make_release_checksums.py --file PATH  # hash one file
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

BUNDLE = Path(__file__).resolve().parents[1]
CHECKSUMS = BUNDLE / "release_checksums.sha256"
DIRS = ["figures", "tables/v15", "catalogues", "mrt"]
EXTS = {".pdf", ".tex", ".fits", ".txt", ".csv", ".json"}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def iter_products():
    for d in DIRS:
        for p in sorted((BUNDLE / d).rglob("*")):
            if p.is_file() and p.suffix.lower() in EXTS:
                yield p


def main(argv: list[str]) -> int:
    if "--file" in argv:
        path = Path(argv[argv.index("--file") + 1])
        print(f"{sha256(path)}  {path}")
        return 0
    if "--verify" in argv:
        if not CHECKSUMS.exists():
            print("no release_checksums.sha256 to verify against")
            return 1
        expected = {}
        for line in CHECKSUMS.read_text().splitlines():
            if line.strip():
                digest, rel = line.split(None, 1)
                expected[rel.strip()] = digest
        bad = []
        for p in iter_products():
            rel = p.relative_to(BUNDLE).as_posix()
            if rel not in expected:
                bad.append(f"UNLISTED {rel}")
            elif sha256(p) != expected[rel]:
                bad.append(f"MISMATCH {rel}")
        for rel in expected:
            if not (BUNDLE / rel).exists():
                bad.append(f"MISSING  {rel}")
        if bad:
            print("\n".join(bad))
            print(f"\n{len(bad)} checksum problem(s).")
            return 1
        print(f"OK: {len(expected)} released products match release_checksums.sha256")
        return 0
    lines = [f"{sha256(p)}  {p.relative_to(BUNDLE).as_posix()}" for p in iter_products()]
    CHECKSUMS.write_text("\n".join(lines) + "\n")
    print(f"wrote {CHECKSUMS.name} ({len(lines)} products)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
