"""Blender (CLI): STEP-mesh of AtomSpectra PRO / GS2020 + Capescint NaI-51x51 in place of the STEP crystal placeholder.
Run: blender -b -P blender_gs2020_capescint.py -- <device.glb> <capescint_stl_dir> <out.blend>"""
import bpy, glob, os, math, sys
G, C, OUT = sys.argv[sys.argv.index("--") + 1:][:3]
bpy.ops.wm.read_factory_settings(use_empty=True); sc = bpy.context.scene; sc.unit_settings.length_unit = "MILLIMETERS"
bpy.ops.import_scene.gltf(filepath=G)
zr = lambda o: [f((o.matrix_world @ v.co).z for v in o.data.vertices) * 1000 for f in (min, max)]
ph = [o for o in sc.objects if o.type == "MESH" and "crystal" in (o.name + o.data.name).lower()]; print("PH", [o.name for o in ph]); ph = ph[0]; z0 = zr(ph)[0] / 1000; ph.name = "STEP_crystal_placeholder"; ph.hide_set(True); ph.hide_render = True
col = bpy.data.collections.new("Capescint_NaI_51x51"); sc.collection.children.link(col)
RGBA = {"NaI_crystal": (0.2, 0.75, 1, 1), "Reflector": (0.95, 0.95, 0.95, 1), "Silica_window": (0.75, 0.9, 1, 0.35),
        "Al_housing_1": (0.62, 0.62, 0.68, 1), "Al_housing_2": (0.55, 0.55, 0.62, 1), "Epoxy_seal_mid": (0.9, 0.55, 0.1, 1), "Epoxy_seal_top": (0.9, 0.55, 0.1, 1)}
for f in sorted(glob.glob(os.path.join(C, "*.stl"))):
    n = os.path.basename(f)[:-4]; bpy.ops.wm.stl_import(filepath=f); o = bpy.context.selected_objects[0]; o.name = n
    for c in list(o.users_collection): c.objects.unlink(o)
    col.objects.link(o); o.rotation_euler = (math.pi / 2, 0, 0); o.scale = (1e-3,) * 3; o.location = (0, 0, z0)
    m = bpy.data.materials.new(n); m.diffuse_color = RGBA[n]; o.data.materials.append(m)
bpy.context.view_layer.update()
for o in sorted([o for o in sc.objects if o.type == "MESH"], key=lambda o: zr(o)[0]):
    xs = [(o.matrix_world @ v.co).x * 1000 for v in o.data.vertices]
    print(f"ZR {o.name[:40]:40s} z {zr(o)[0]:7.2f}..{zr(o)[1]:7.2f}  x {min(xs):7.2f}..{max(xs):7.2f}")
from mathutils import Matrix  # 1 BU = 1 mm, иначе прибор 0,083 BU не виден в стандартном виде
FLIP = Matrix.Translation((0, 0, 83)) @ Matrix.Rotation(3.141592653589793, 4, "X")  # оператор 24.09: переворот 180°, кристалл сверху
for o in [o for o in sc.objects if o.parent is None]: o.matrix_world = FLIP @ Matrix.Scale(1000, 4) @ o.matrix_world
sc.unit_settings.scale_length = 0.001
for s in [a.spaces[0] for scr in bpy.data.screens for a in scr.areas if a.type == "VIEW_3D"]:
    s.clip_start, s.clip_end, s.shading.show_xray = 0.1, 10000, True; s.region_3d.view_location, s.region_3d.view_distance = (0, 0, 41.5), 220
bpy.ops.wm.save_as_mainfile(filepath=OUT); print("SAVED", OUT)
