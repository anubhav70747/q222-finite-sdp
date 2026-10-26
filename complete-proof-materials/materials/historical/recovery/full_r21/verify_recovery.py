#!/usr/bin/env python3
"""Verify delivered file identity and nested ZIP CRCs, not mathematics."""
from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "MANIFEST_SHA256.json").read_text(encoding="utf-8"))
    errors: list[str] = []
    checked = 0
    for entry in manifest["files"]:
        relative = Path(entry["path"])
        path = (root / relative).resolve()
        if relative.is_absolute() or root not in path.parents:
            errors.append(f"Unsafe manifest path: {relative}")
            continue
        if not path.is_file():
            errors.append(f"Missing: {relative}")
            continue
        checked += 1
        if path.stat().st_size != entry["size_bytes"]:
            errors.append(f"Size mismatch: {relative}")
        if digest(path) != entry["sha256"]:
            errors.append(f"SHA-256 mismatch: {relative}")

    archives = 0
    for entry in manifest["files"]:
        if not entry["path"].endswith(".zip"):
            continue
        path = root / entry["path"]
        try:
            with zipfile.ZipFile(path) as archive:
                bad = archive.testzip()
            archives += 1
            if bad is not None:
                errors.append(f"ZIP CRC mismatch: {entry['path']} / {bad}")
        except (OSError, zipfile.BadZipFile, RuntimeError) as exc:
            errors.append(f"ZIP unreadable: {entry['path']}: {exc}")

    print(f"Manifest-bound files checked: {checked}")
    print(f"Nested ZIP archives checked: {archives}")
    if errors:
        print("RESULT: FAIL")
        for error in errors:
            print(error)
        return 1
    print("RESULT: PASS — file identity and archive integrity only.")
    print("No mathematical claim was tested or certified by this program.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
