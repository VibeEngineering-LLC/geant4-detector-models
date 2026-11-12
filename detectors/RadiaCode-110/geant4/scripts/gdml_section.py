"""XZ section (y=0) of marinelli GDML tessellated solids read back from the GDML. Usage: gdml_section.py out.png a.gdml [b.gdml ...]"""
import re, sys, numpy as np, trimesh, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, len(sys.argv) - 2, figsize=(6 * (len(sys.argv) - 2), 6))
for ax, g in zip(np.atleast_1d(axs), sys.argv[2:]):
    t = open(g, encoding="utf-8").read()
    pos = {m[0]: (float(m[1]), float(m[2]), float(m[3])) for m in re.findall(r'<position name="(\w+)" unit="mm" x="([^"]+)" y="([^"]+)" z="([^"]+)"', t)}
    for part, col in (("can", "tab:blue"), ("cap", "tab:red")):
        tri = re.findall(rf'<triangular vertex1="({part}_v\d+)" vertex2="({part}_v\d+)" vertex3="({part}_v\d+)"', t)
        m = trimesh.Trimesh(np.array([[pos[a], pos[b], pos[c]] for a, b, c in tri]).reshape(-1, 3), np.arange(3 * len(tri)).reshape(-1, 3), process=False)
        s = m.section(plane_origin=[0, 0, 0], plane_normal=[0, 1, 0])
        for e in s.discrete: ax.plot(e[:, 0], e[:, 2], col, lw=0.7)
    ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_title(g.replace("\\", "/").split("/")[-1] + "  (blue can, red cap; mm)")
plt.tight_layout(); plt.savefig(sys.argv[1], dpi=90); print("saved", sys.argv[1])
