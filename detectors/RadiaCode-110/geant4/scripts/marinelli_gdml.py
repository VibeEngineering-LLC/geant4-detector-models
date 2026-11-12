import sys
sys.stdout.reconfigure(encoding="utf-8")
import json
import math
import numpy as np
import trimesh
from matplotlib.path import Path

if len(sys.argv) != 5:
    sys.stderr.write("Usage: python marinelli_gdml.py <can.stl> <cap.stl> <name> <out_dir>\n")
    sys.exit(2)

can_path, cap_path, name, out_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

# Step 1: Load and center in fit frame F
def load_centered(path):
    m = trimesh.load(path, process=False)
    c = (m.bounds[0] + m.bounds[1]) / 2
    m.apply_translation([-c[0], -c[1], -m.bounds[0][2]])
    return m

can = load_centered(can_path)
cap = load_centered(cap_path)
z_neck = float(can.bounds[1][2])

# Step 2: Cap geometry parameters
v = cap.vertices
r = np.hypot(v[:, 0], v[:, 1])
mask_bore = v[:, 2] > 0.5 * cap.bounds[1][2]
r_bore = float(np.min(r[mask_bore]))
# plate inner face has vertices only at its rim (bore corner): lowest z of bore-radius vertices above the outer face
mask_plate = (r < r_bore + 0.05) & (v[:, 2] > 0.05)
t_plate = round(float(np.min(v[mask_plate, 2])), 4)

# Step 3: Seat the cap (180 deg about X, then translate)
cap.apply_transform(trimesh.transformations.rotation_matrix(math.pi, [1, 0, 0]))
cap.apply_translation([0, 0, z_neck + t_plate])
print(f"SEAT t_plate={t_plate} z_neck={z_neck} cap_z={cap.bounds[0][2]}..{cap.bounds[1][2]}")

# Step 4: Choose cap angle theta about Z by minimal overlap volume
DZ = 0.1
NA = 72

def inside(mesh, z, pts):
    s = mesh.section(plane_origin=[0, 0, z], plane_normal=[0, 0, 1])
    result = np.zeros(len(pts), dtype=bool)
    if s is None:
        return result
    for loop in s.discrete:
        result ^= Path(loop[:, :2]).contains_points(pts)
    return result

rad = np.linspace(r_bore - 1.5, cap.bounds[1][0] + 0.2, 40)
ang = np.linspace(0, 2 * math.pi, 720, endpoint=False)
RR, AA = np.meshgrid(rad, ang)
P = np.c_[(RR * np.cos(AA)).ravel(), (RR * np.sin(AA)).ravel()]
cell = (rad[1] - rad[0]) * (2 * math.pi / 720) * RR.ravel() * DZ
zs = np.arange(cap.bounds[0][2] + DZ / 2, z_neck - 1e-3, DZ)

# Precompute can cross-sections once
can_in = [inside(can, float(z), P) for z in zs]

overlap = np.zeros(NA)
for k in range(NA):
    theta = 2 * math.pi * k / NA
    c_val, s_val = math.cos(theta), math.sin(theta)
    Pr = np.c_[c_val * P[:, 0] + s_val * P[:, 1], -s_val * P[:, 0] + c_val * P[:, 1]]
    ov = 0.0
    for i in range(len(zs)):
        ov += cell[can_in[i] & inside(cap, float(zs[i]), Pr)].sum()
    overlap[k] = ov

k_best = int(np.argmin(overlap))
theta_best = 2 * math.pi * k_best / NA
cap.apply_transform(trimesh.transformations.rotation_matrix(theta_best, [0, 0, 1]))
print(f"THETA best_deg={math.degrees(theta_best)} overlap_min_mm3={overlap[k_best]} overlap_max_mm3={overlap.max()}")

# Step 5: Well floor on the axis (section of can at z=0, normal x)
s = can.section(plane_origin=[0, 0, 0], plane_normal=[1, 0, 0])
crossings = []
if s is not None:
    for loop in s.discrete:
        for j in range(len(loop) - 1):
            a = loop[j]
            b = loop[j + 1]
            if (a[1] <= 0 < b[1]) or (b[1] <= 0 < a[1]):
                z_cross = a[2] + (b[2] - a[2]) * (0 - a[1]) / (b[1] - a[1])
                crossings.append(float(z_cross))
crossings.sort()
z_floor = crossings[0] if crossings else 0.0
print(f"AXIS crossings={[round(c, 4) for c in crossings]}")

# Step 6: Working orientation (frame W): 180 deg about X then translate [0,0,z_floor]
M_rot = trimesh.transformations.rotation_matrix(math.pi, [1, 0, 0])
M_full = np.eye(4)
M_full[0, 0] = 1.0
M_full[1, 1] = -1.0
M_full[2, 2] = -1.0
M_full[2, 3] = z_floor
can.apply_transform(M_full)
cap.apply_transform(M_full)
# 1 um clearance (cap moved away from can, -Z in W): touching thread/end faces otherwise give
# a 2 nm G4 overlap from %.6f rounding (checked 2026-11-12: -1 um -> OK, +1 um -> 1 um overlap)
CAP_CLEARANCE_MM = 0.001
cap.apply_translation([0, 0, -CAP_CLEARANCE_MM])

# Step 7: GDML output
def dedup_vertices(vertices):
    seen = {}
    order = []
    for idx in range(len(vertices)):
        v = vertices[idx]
        key = f"{v[0]:.6f},{v[1]:.6f},{v[2]:.6f}"
        if key not in seen:
            seen[key] = len(order)
            order.append(idx)
    return seen, order

can_map, can_order = dedup_vertices(can.vertices)
cap_map, cap_order = dedup_vertices(cap.vertices)
n_vert_can = len(can_order)
n_vert_cap = len(cap_order)

# World box dimensions
all_coords = np.vstack([can.vertices, cap.vertices])
max_abs = np.max(np.abs(all_coords), axis=0)
box_x = float(2 * (max_abs[0] + 20))
box_y = float(2 * (max_abs[1] + 20))
box_z = float(2 * (max_abs[2] + 20))

lines = []
lines.append('<?xml version="1.0" encoding="UTF-8"?>')
lines.append('<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">')

# define: vertices
lines.append('  <define>')
for i in range(len(can_order)):
    orig = can_order[i]
    v = can.vertices[orig]
    lines.append(f'    <position name="can_v{i}" unit="mm" x="{v[0]:.6f}" y="{v[1]:.6f}" z="{v[2]:.6f}"/>')
for i in range(len(cap_order)):
    orig = cap_order[i]
    v = cap.vertices[orig]
    lines.append(f'    <position name="cap_v{i}" unit="mm" x="{v[0]:.6f}" y="{v[1]:.6f}" z="{v[2]:.6f}"/>')
lines.append('  </define>')

# materials
lines.append('  <materials>')
lines.append('    <element name="H" formula="H" Z="1"><atom value="1.008"/></element>')
lines.append('    <element name="C" formula="C" Z="6"><atom value="12.011"/></element>')
lines.append('    <element name="O" formula="O" Z="8"><atom value="15.999"/></element>')
lines.append('    <element name="N" formula="N" Z="7"><atom value="14.007"/></element>')
lines.append('    <element name="Ar" formula="Ar" Z="18"><atom value="39.948"/></element>')
lines.append('    <material name="PETG" state="solid"><D value="1.27" unit="g/cm3"/><composite n="10" ref="C"/><composite n="8" ref="H"/><composite n="4" ref="O"/></material>')
lines.append('    <material name="AIR" state="gas"><D value="0.00120479" unit="g/cm3"/><fraction n="0.755267" ref="N"/><fraction n="0.231781" ref="O"/><fraction n="0.012827" ref="Ar"/><fraction n="0.000124" ref="C"/></material>')
lines.append('  </materials>')

# solids
lines.append('  <solids>')
lines.append(f'    <tessellated name="Marinelli_{name}_Can_solid" aunit="deg" lunit="mm">')
for face in can.faces:
    a, b, c = int(face[0]), int(face[1]), int(face[2])
    ka = f"{can.vertices[a][0]:.6f},{can.vertices[a][1]:.6f},{can.vertices[a][2]:.6f}"
    kb = f"{can.vertices[b][0]:.6f},{can.vertices[b][1]:.6f},{can.vertices[b][2]:.6f}"
    kc = f"{can.vertices[c][0]:.6f},{can.vertices[c][1]:.6f},{can.vertices[c][2]:.6f}"
    lines.append(f'      <triangular vertex1="can_v{can_map[ka]}" vertex2="can_v{can_map[kb]}" vertex3="can_v{can_map[kc]}" type="ABSOLUTE"/>')
lines.append('    </tessellated>')

lines.append(f'    <tessellated name="Marinelli_{name}_Cap_solid" aunit="deg" lunit="mm">')
for face in cap.faces:
    a, b, c = int(face[0]), int(face[1]), int(face[2])
    ka = f"{cap.vertices[a][0]:.6f},{cap.vertices[a][1]:.6f},{cap.vertices[a][2]:.6f}"
    kb = f"{cap.vertices[b][0]:.6f},{cap.vertices[b][1]:.6f},{cap.vertices[b][2]:.6f}"
    kc = f"{cap.vertices[c][0]:.6f},{cap.vertices[c][1]:.6f},{cap.vertices[c][2]:.6f}"
    lines.append(f'      <triangular vertex1="cap_v{cap_map[ka]}" vertex2="cap_v{cap_map[kb]}" vertex3="cap_v{cap_map[kc]}" type="ABSOLUTE"/>')
lines.append('    </tessellated>')

lines.append(f'    <box name="World_solid" x="{box_x}" y="{box_y}" z="{box_z}" lunit="mm"/>')
lines.append('  </solids>')

# structure
lines.append('  <structure>')
# GDML: materialref/solidref are child elements (attributes are silently ignored by G4GDMLParser)
lines.append(f'    <volume name="Marinelli_{name}_Can"><materialref ref="PETG"/><solidref ref="Marinelli_{name}_Can_solid"/></volume>')
lines.append(f'    <volume name="Marinelli_{name}_Cap"><materialref ref="PETG"/><solidref ref="Marinelli_{name}_Cap_solid"/></volume>')
lines.append(f'    <volume name="World"><materialref ref="AIR"/><solidref ref="World_solid"/>')
lines.append(f'      <physvol name="Marinelli_{name}_Can_pv"><volumeref ref="Marinelli_{name}_Can"/></physvol>')
lines.append(f'      <physvol name="Marinelli_{name}_Cap_pv"><volumeref ref="Marinelli_{name}_Cap"/></physvol>')
lines.append('    </volume>')
lines.append('  </structure>')

# setup
lines.append('  <setup name="Default" version="1.0"><world ref="World"/></setup>')
lines.append('</gdml>')

gdml_path = f"{out_dir}/RC110_marinelli_{name}.gdml"
with open(gdml_path, "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(lines) + "\n")

# Step 8: JSON report
bounds_W = {
    "min": [float(min(can.bounds[0][i], cap.bounds[0][i])) for i in range(3)],
    "max": [float(max(can.bounds[1][i], cap.bounds[1][i])) for i in range(3)]
}

report = {
    "can": can_path,
    "cap": cap_path,
    "name": name,
    "n_tri_can": len(can.faces),
    "n_tri_cap": len(cap.faces),
    "n_vert_can": n_vert_can,
    "n_vert_cap": n_vert_cap,
    "t_plate_mm": t_plate,
    "z_neck_F_mm": float(z_neck),
    "r_bore_mm": r_bore,
    "theta_deg": math.degrees(theta_best),
    "overlap_min_mm3": float(overlap[k_best]),
    "overlap_max_mm3": float(overlap.max()),
    "overlap_scan": [round(float(x), 3) for x in overlap],
    "z_floor_F_mm": float(z_floor),
    "axis_crossings_F_mm": [round(c, 4) for c in crossings],
    "can_volume_mm3": float(can.volume),
    "cap_volume_mm3": float(cap.volume),
    "bounds_W_mm": bounds_W,
    "gdml": gdml_path
}

json_path = f"{out_dir}/RC110_marinelli_{name}.json"
with open(json_path, "w", encoding="utf-8", newline="\n") as f:
    json.dump(report, f, ensure_ascii=False, indent=1)
    f.write("\n")

# Print JSON without overlap_scan
print_report = {k: v for k, v in report.items() if k != "overlap_scan"}
print(json.dumps(print_report, ensure_ascii=False, indent=1))
