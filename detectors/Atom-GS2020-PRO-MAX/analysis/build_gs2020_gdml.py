import json
import struct
import sys
import os
import re
import math
from pathlib import Path

# Constants
Z_SHIFT = 41.5
Y_TO_Z = 37.75
WORLD_HALF_SIZE = 100.0

def main():
    if len(sys.argv) != 5:
        print("Usage: python build_gs2020_gdml.py <profile.json> <parts_dir> <materials.json> <out.gdml>")
        sys.exit(1)

    profile_path = sys.argv[1]
    parts_dir = sys.argv[2]
    materials_path = sys.argv[3]
    out_gdml_path = sys.argv[4]

    # Load inputs
    with open(profile_path, 'r', encoding='utf-8') as f:
        profile_data = json.load(f)
    
    manifest_path = os.path.join(parts_dir, "_manifest.json")
    with open(manifest_path, 'r', encoding='utf-8') as f:
        manifest = json.load(f)
    
    with open(materials_path, 'r', encoding='utf-8') as f:
        materials_data = json.load(f)

    # Process Materials
    custom_materials = materials_data.get("custom", [])
    profile_parts_map = materials_data.get("profile_parts", {})
    stl_rules = materials_data.get("stl_rules", [])
    stl_default = materials_data.get("stl_default", "G4_AIR")
    world_material = materials_data.get("world_material", "G4_AIR")

    # Helper to get material for STL part
    def get_stl_material(name):
        for pattern, mat in stl_rules:
            if re.search(pattern, name):
                return mat
        return stl_default

    # Prepare GDML sections
    define_section = []
    materials_section = []
    solids_section = []
    structure_volumes = []
    report_lines = []
    
    total_facets = 0
    total_skipped = 0
    parts_dropped = 0
    
    # Process Profile Parts (genericPolycone)
    for part_name, part_info in profile_data.get("parts", {}).items():
        mat_name = profile_parts_map.get(part_name)
        if not mat_name:
            print(f"Error: No material defined for profile part '{part_name}'")
            sys.exit(2)
        
        solid_name = f"S_{part_name}"
        profile_r_y = part_info["profile_r_y"]
        
        # Build rzpoints
        rzpoints_xml = ""
        for r, y in profile_r_y:
            z = Y_TO_Z - y
            r_str = repr(round(r, 6))
            z_str = repr(round(z, 6))
            rzpoints_xml += f'    <rzpoint r="{r_str}" z="{z_str}"/>\n'
        
        solid_xml = f'<genericPolycone name="{solid_name}" startphi="0" deltaphi="360" aunit="deg" lunit="mm">\n{rzpoints_xml}</genericPolycone>'
        solids_section.append(solid_xml)
        
        vol_name = "Crystal_log" if (os.environ.get("GS_DEVICE_MODE") == "1" and part_name == "NaI_crystal") else f"LV_{part_name}"
        vol_xml = f'<volume name="{vol_name}"><materialref ref="{mat_name}"/><solidref ref="{solid_name}"/></volume>'
        structure_volumes.append(vol_xml)
        
        n_points = len(profile_r_y)
        report_lines.append(f"VOL {part_name} genericPolycone {mat_name} {n_points} skipped=0")

    # Process STL Parts (tessellated)
    for idx, item in enumerate(manifest):
        original_name = item["name"]
        file_name = item["file"]
        safe_name = Path(file_name).stem
        solid_name = f"S_{safe_name}"
        
        mat_name = get_stl_material(original_name)
        
        stl_path = os.path.join(parts_dir, file_name)
        if not os.path.exists(stl_path):
            print(f"Warning: STL file {stl_path} not found. Skipping.")
            parts_dropped += 1
            continue
            
        try:
            with open(stl_path, 'rb') as f:
                data = f.read()
        except Exception as e:
            print(f"Error reading {stl_path}: {e}")
            parts_dropped += 1
            continue

        # Parse Binary STL
        # Header 80 bytes
        # Count uint32
        # Triangles: normal(3f), v1(3f), v2(3f), v3(3f), attr(1h) -> 12*4 + 2 = 50 bytes
        
        if len(data) < 84:
            print(f"Warning: STL file {file_name} too small. Skipping.")
            parts_dropped += 1
            continue
            
        count = struct.unpack_from('<I', data, 80)[0]
        
        vertices_map = {} # key: rounded tuple -> index
        vertex_list = [] # list of (x,y,z)
        facets_xml = ""
        skipped_triangles = 0
        valid_facets_count = 0
        tri_list = []  # fix: collect facets, flip all if the mesh is inside-out (negative signed volume)
        signed_vol = 0.0
        
        offset = 84
        for _ in range(count):
            if offset + 50 > len(data):
                break
                
            # Unpack triangle data
            # Normal (3 floats)
            nx, ny, nz = struct.unpack_from('<fff', data, offset); offset += 12
            # V1 (3 floats)
            x1, y1, z1 = struct.unpack_from('<fff', data, offset); offset += 12
            # V2 (3 floats)
            x2, y2, z2 = struct.unpack_from('<fff', data, offset); offset += 12
            # V3 (3 floats)
            x3, y3, z3 = struct.unpack_from('<fff', data, offset); offset += 12
            # Attribute byte count (uint16)
            attr = struct.unpack_from('<H', data, offset)[0]; offset += 2
            
            # Apply Z_SHIFT: z -> z - Z_SHIFT
            z1 -= Z_SHIFT
            z2 -= Z_SHIFT
            z3 -= Z_SHIFT
            
            # Round to 1e-5 mm for deduplication key
            def round_key(x, y, z):
                return (round(x, 5), round(y, 5), round(z, 5))
            
            k1 = round_key(x1, y1, z1)
            k2 = round_key(x2, y2, z2)
            k3 = round_key(x3, y3, z3)
            
            # Map keys to indices
            if k1 not in vertices_map:
                vertices_map[k1] = len(vertex_list)
                vertex_list.append((x1, y1, z1))
            if k2 not in vertices_map:
                vertices_map[k2] = len(vertex_list)
                vertex_list.append((x2, y2, z2))
            if k3 not in vertices_map:
                vertices_map[k3] = len(vertex_list)
                vertex_list.append((x3, y3, z3))
                
            i1 = vertices_map[k1]
            i2 = vertices_map[k2]
            i3 = vertices_map[k3]
            
            # Skip if duplicate vertices in triangle
            if i1 == i2 or i2 == i3 or i1 == i3:
                skipped_triangles += 1
                continue
                
            # Calculate area to skip degenerate triangles
            # Area = 0.5 * |cross(b-a, c-a)|
            ax, ay, az = vertex_list[i1]
            bx, by, bz = vertex_list[i2]
            cx, cy, cz = vertex_list[i3]
            
            abx = bx - ax
            aby = by - ay
            abz = bz - az
            
            acx = cx - ax
            acy = cy - ay
            acz = cz - az
            
            cross_x = aby * acz - abz * acy
            cross_y = abz * acx - abx * acz
            cross_z = abx * acy - aby * acx
            
            area = 0.5 * math.sqrt(cross_x**2 + cross_y**2 + cross_z**2)
            
            if area < 1e-8:
                skipped_triangles += 1
                continue
                
            # Add facet
            v1_name = f"v{idx}_{i1}"
            v2_name = f"v{idx}_{i2}"
            v3_name = f"v{idx}_{i3}"
            
            tri_list.append((v1_name, v2_name, v3_name))
            signed_vol += (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6.0
            valid_facets_count += 1
            
        flipped = signed_vol < 0
        for a_, b_, c_ in tri_list:
            if flipped: b_, c_ = c_, b_
            facets_xml += f'    <triangular vertex1="{a_}" vertex2="{b_}" vertex3="{c_}" type="ABSOLUTE"/>\n'
        if flipped: print(f"FLIPPED {safe_name} signed_volume_mm3={signed_vol:.3f}")
        if valid_facets_count < 4:
            print(f"Warning: Part {original_name} has fewer than 4 valid triangles ({valid_facets_count}). Skipping.")
            parts_dropped += 1
            continue
            
        # Add vertices to define section
        for i, (vx, vy, vz) in enumerate(vertex_list):
            vx_str = repr(round(vx, 6))
            vy_str = repr(round(vy, 6))
            vz_str = repr(round(vz, 6))
            v_name = f"v{idx}_{i}"
            define_section.append(f'  <position name="{v_name}" x="{vx_str}" y="{vy_str}" z="{vz_str}" unit="mm"/>')
            
        # Build solid XML
        solid_xml = f'<tessellated name="{solid_name}">\n{facets_xml}</tessellated>'
        solids_section.append(solid_xml)
        
        vol_name = f"LV_{safe_name}"  # fix: ASCII-safe names (original names contain spaces/colons)
        vol_xml = f'<volume name="{vol_name}"><materialref ref="{mat_name}"/><solidref ref="{solid_name}"/></volume>'
        structure_volumes.append(vol_xml)
        
        total_facets += valid_facets_count
        total_skipped += skipped_triangles
        
        report_lines.append(f"VOL {safe_name} tessellated {mat_name} {valid_facets_count} skipped={skipped_triangles}")

    # Build Materials Section XML
    materials_xml = ""
    for mat in custom_materials:
        m_name = mat["name"]
        m_state = mat.get("state", "solid")
        m_density = mat["density_g_cm3"]
        
        mat_str = f'<material name="{m_name}" state="{m_state}"><D value="{m_density}" unit="g/cm3"/>'
        for elem, frac in mat["fractions"]:
            mat_str += f'<fraction n="{frac}" ref="{elem}"/>'
        mat_str += '</material>'
        materials_xml += mat_str + "\n"

    # Build Define Section XML
    define_xml = ""
    if define_section:
        define_xml = "<define>\n" + "\n".join(define_section) + "\n</define>"

    # Build Solids Section XML
    # World Box Solid (fix: was built but never added to <solids>)
    world_box_size = 2 * WORLD_HALF_SIZE
    world_solid = f'<box name="WorldBox" x="{world_box_size}" y="{world_box_size}" z="{world_box_size}" lunit="mm"/>'
    # DEVICE mode (env GS_DEVICE_MODE=1, 2026-09-24): parts go into an air box "RC103_device_log" (name expected by
    # rc103_field, run_field/Rc103FieldDetectorConstruction.cc:260), NaI LV renamed "Crystal_log" (same file, :268).
    DEVICE = os.environ.get("GS_DEVICE_MODE") == "1"
    if DEVICE:
        world_solid += '\n<box name="DeviceBox" x="64.0" y="64.0" z="84.0" lunit="mm"/>'
    solids_xml = "<solids>\n" + world_solid + "\n" + "\n".join(solids_section) + "\n</solids>"

    # Build Structure Section XML
    
    structure_xml = "<structure>\n"
    # Add all part volumes first
    for vol in structure_volumes:
        structure_xml += f"  {vol}\n"
        
    # World Volume (in DEVICE mode the parts are placed into RC103_device_log, which is placed into World)
    world_vol_name = "RC103_device_log" if DEVICE else "World"
    structure_xml += f'  <volume name="{world_vol_name}"><materialref ref="{world_material}"/><solidref ref="{"DeviceBox" if DEVICE else "WorldBox"}"/>'
    
    # Add physvols for each part volume
    for vol in structure_volumes:
        # Extract volume name from XML string
        match = re.search(r'name="([^"]+)"', vol)
        if match:
            v_name = match.group(1)
            structure_xml += f'<physvol name="PV_{v_name}"><volumeref ref="{v_name}"/></physvol>'
            
    structure_xml += "</volume>\n"
    if DEVICE:
        structure_xml += (f'  <volume name="World"><materialref ref="{world_material}"/><solidref ref="WorldBox"/>'
                          '<physvol name="PV_RC103_device"><volumeref ref="RC103_device_log"/></physvol></volume>\n')
    structure_xml += "</structure>"

    # Setup Section
    setup_xml = '<setup name="Default" version="1.0"><world ref="World"/></setup>'

    # Assemble GDML
    gdml_header = '<?xml version="1.0" encoding="UTF-8"?>\n'
    gdml_root_start = '<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">\n'
    
    full_gdml = gdml_header + gdml_root_start
    if define_xml:
        full_gdml += define_xml + "\n"
    if materials_xml:
        full_gdml += "<materials>\n" + materials_xml + "</materials>\n"
    full_gdml += solids_xml + "\n"
    full_gdml += structure_xml + "\n"
    full_gdml += setup_xml + "\n</gdml>"

    # Write GDML
    with open(out_gdml_path, 'w', encoding='utf-8') as f:
        f.write(full_gdml)

    # Report
    print("\n--- Build Report ---")
    for line in report_lines:
        print(line)
        
    total_volumes = len(report_lines)
    final_line = f"TOTAL volumes={total_volumes} facets={total_facets} skipped={total_skipped} parts_dropped={parts_dropped}"
    print(final_line)

    # Write Summary JSON
    summary_data = []
    for line in report_lines:
        parts = line.split()
        name = parts[1]
        solid_type = parts[2]
        material = parts[3]
        count_str = parts[4]
        skipped_str = parts[5].split('=')[1]
        
        summary_data.append({
            "name": name,
            "solid_type": solid_type,
            "material": material,
            "count": int(count_str),
            "skipped": int(skipped_str)
        })
        
    summary_data.append({
        "summary": True,
        "total_volumes": total_volumes,
        "total_facets": total_facets,
        "total_skipped": total_skipped,
        "parts_dropped": parts_dropped
    })

    summary_path = out_gdml_path + ".summary.json"
    with open(summary_path, 'w', encoding='utf-8') as f:
        json.dump(summary_data, f, indent=2)

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
