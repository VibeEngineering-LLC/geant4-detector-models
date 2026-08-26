#!/usr/bin/env python3
"""Binary STL -> GDML tessellated. See ../docs/GEANT4-MODEL.md."""
from __future__ import annotations

import re
import struct
import sys
from pathlib import Path

GEANT4_ROOT = Path(__file__).resolve().parents[1]
RC103_ROOT = GEANT4_ROOT.parent

def resolve_stl(root: Path | None = None) -> Path:
    """Pick Rc-103.stl or Rc-103.stl.bak under RC103 root."""
    base = root if root is not None else RC103_ROOT
    for name in ("Rc-103.stl", "Rc-103.stl.bak"):
        candidate = base / name
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"No STL found under {base} (Rc-103.stl or Rc-103.stl.bak)")


STL_PATH = resolve_stl()
DETECTOR_GDML = GEANT4_ROOT / "gdml" / "detector" / "RC103_detector.gdml"
OUT_MESH = GEANT4_ROOT / "gdml" / "case_mesh" / "RC103_case_mesh.gdml"
OUT_FULL = GEANT4_ROOT / "gdml" / "full" / "RC103_full.gdml"




def remap(x: float, y: float, z: float) -> tuple[float, float, float]:
    return (x, y, z)  # Rc-103.stl already in device mm frame


def read_stl_triangles(path: Path) -> list[tuple[tuple[float, float, float], ...]]:
    data = path.read_bytes()
    n = struct.unpack_from("<I", data, 80)[0]
    tris: list[tuple[tuple[float, float, float], ...]] = []
    off = 84
    for _ in range(n):
        vals = struct.unpack_from("<12fH", data, off)
        verts = [remap(vals[3 + i * 3], vals[4 + i * 3], vals[5 + i * 3]) for i in range(3)]
        tris.append((verts[0], verts[1], verts[2]))
        off += 50
    return tris


def write_mesh_gdml(tris: list, out_path: Path) -> None:
    n = len(tris)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n')
        f.write(f"<!-- RC103 case mesh; {n} triangles from {STL_PATH.name} -->\n")
        f.write(
            '<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n'
            '      xsi:noNamespaceSchemaLocation="'
            "http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd"
            '">\n\n'
        )
        f.write("  <define><position name=\"p_origin\" x=\"0\" y=\"0\" z=\"0\" unit=\"mm\"/></define>\n\n")
        f.write("  <materials>\n")
        f.write('    <material Z="6" name="Case_ABS" state="solid">\n')
        f.write('      <D unit="g/cm3" value="1.05"/>\n')
        f.write('      <fraction n="0.85" ref="G4_C"/>\n')
        f.write('      <fraction n="0.08" ref="G4_H"/>\n')
        f.write('      <fraction n="0.07" ref="G4_N"/>\n')
        f.write("    </material>\n")
        f.write('    <material Z="1" name="G4_AIR" state="gas">\n')
        f.write('      <D unit="g/cm3" value="0.00120479"/>\n')
        f.write('      <fraction n="0.75527" ref="G4_N"/>\n')
        f.write('      <fraction n="0.231781" ref="G4_O"/>\n')
        f.write('      <fraction n="0.012827" ref="G4_Ar"/>\n')
        f.write("    </material>\n  </materials>\n\n")
        f.write("  <solids>\n")
        f.write('    <box lunit="mm" name="WorldBox" x="400" y="400" z="400"/>\n')
        f.write('    <tessellated aunit="deg" lunit="mm" name="Case_STL_mesh">\n')
        for v1, v2, v3 in tris:
            f.write('      <triangular type="ABSOLUTE">\n')
            for vx, vy, vz in (v1, v2, v3):
                f.write(f'        <point x="{vx:.6f}" y="{vy:.6f}" z="{vz:.6f}" unit="mm"/>\n')
            f.write("      </triangular>\n")
        f.write("    </tessellated>\n  </solids>\n\n")
        f.write("  <structure>\n")
        f.write('    <volume name="Case_STL_log"><materialref ref="Case_ABS"/><solidref ref="Case_STL_mesh"/></volume>\n')
        f.write('    <volume name="World"><materialref ref="G4_AIR"/><solidref ref="WorldBox"/>\n')
        f.write('      <physvol name="pv_case_stl"><volumeref ref="Case_STL_log"/><positionref ref="p_origin"/></physvol>\n')
        f.write("    </volume>\n  </structure>\n\n")
        f.write('  <setup name="Default" version="1.0"><world ref="World"/></setup>\n</gdml>\n')


def extract_gdml_section(text: str, tag: str) -> str:
    m = re.search(rf"<{tag}>(.*?)</{tag}>", text, flags=re.DOTALL)
    if not m:
        raise ValueError(f"GDML section <{tag}> not found")
    return m.group(1)


def filter_solid_primitives(solids: str) -> str:
    solids = re.sub(r'\s*<box lunit="mm" name="Case_inner".*?/>\n', "", solids, count=1)
    solids = re.sub(
        r'\s*<subtraction name="Case_shell">.*?</subtraction>\n', "", solids, flags=re.DOTALL, count=1
    )
    return solids


def filter_structure_shell(structure: str) -> str:
    structure = re.sub(
        r'\s*<volume name="Case_shell_log">.*?</volume>\n', "", structure, flags=re.DOTALL, count=1
    )
    return structure.replace(
        '      <physvol name="pv_case_shell">\n'
        '        <volumeref ref="Case_shell_log"/>\n'
        '        <positionref ref="p_origin"/>\n'
        "      </physvol>\n",
        '      <physvol name="pv_case_stl">\n'
        '        <volumeref ref="Case_STL_log"/>\n'
        '        <positionref ref="p_origin"/>\n'
        "      </physvol>\n",
    )


def write_tessellated_block(tris: list) -> str:
    lines = ['    <tessellated aunit="deg" lunit="mm" name="Case_STL_mesh">']
    for v1, v2, v3 in tris:
        lines.append('      <triangular type="ABSOLUTE">')
        for vx, vy, vz in (v1, v2, v3):
            lines.append(f'        <point x="{vx:.6f}" y="{vy:.6f}" z="{vz:.6f}" unit="mm"/>')
        lines.append("      </triangular>")
    lines.append("    </tessellated>")
    return "\n".join(lines) + "\n"


def write_full_gdml(tris: list, out_path: Path) -> None:
    n = len(tris)
    det = DETECTOR_GDML.read_text(encoding="utf-8")
    define = extract_gdml_section(det, "define")
    materials = extract_gdml_section(det, "materials")
    solids = filter_solid_primitives(extract_gdml_section(det, "solids"))
    structure = filter_structure_shell(extract_gdml_section(det, "structure"))
    structure = (
        '    <volume name="Case_STL_log">\n'
        '      <materialref ref="Case_ABS"/>\n'
        '      <solidref ref="Case_STL_mesh"/>\n'
        "    </volume>\n"
    ) + structure

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n')
        f.write("<!-- RC103 full: STL case + all internals -->\n")
        f.write(f"<!-- STL triangles: {n}; merged from {DETECTOR_GDML.name} -->\n")
        f.write(
            '<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n'
            '      xsi:noNamespaceSchemaLocation="'
            "http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd"
            '">\n\n'
        )
        f.write(f"  <define>{define}  </define>\n\n")
        f.write(f"  <materials>{materials}  </materials>\n\n")
        f.write("  <solids>\n")
        f.write(write_tessellated_block(tris))
        f.write(solids)
        f.write("  </solids>\n\n")
        f.write(f"  <structure>{structure}  </structure>\n\n")
        f.write('  <setup name="Default" version="1.0"><world ref="World"/></setup>\n</gdml>\n')


def main() -> int:
    tris = read_stl_triangles(STL_PATH)
    print(f"Triangles: {len(tris)}", file=sys.stderr)
    write_mesh_gdml(tris, OUT_MESH)
    write_full_gdml(tris, OUT_FULL)
    for p in (OUT_MESH, OUT_FULL):
        print(f"{p}: {p.stat().st_size / 1e6:.1f} MB", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
