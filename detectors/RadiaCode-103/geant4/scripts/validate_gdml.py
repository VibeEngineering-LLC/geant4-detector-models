#!/usr/bin/env python3
"""Well-formed XML check for gdml/**/*.gdml under this geant4 tree."""
from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

GEANT4_ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    gdml_root = GEANT4_ROOT / "gdml"
    files = sorted(gdml_root.rglob("*.gdml"))
    if not files:
        print(f"No *.gdml under {gdml_root}", file=sys.stderr)
        return 1
    failed = False
    for path in files:
        rel = path.relative_to(GEANT4_ROOT)
        try:
            ET.parse(path)
        except ET.ParseError as exc:
            print(f"FAIL {rel}: {exc}", file=sys.stderr)
            failed = True
        else:
            print(f"OK   {rel}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
