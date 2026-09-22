import sys, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
dt = np.dtype([("evt","<i4"),("pdg","<i4")]+[(k,"<f4") for k in ("E","phi","y","ux","uy","uz","t","Ep","wp")]+[(k,"<i2") for k in ("proc","vol","ip","pad")])
assert dt.itemsize == 52
a = np.fromfile(sys.argv[1], dtype=dt)
print("records", len(a), "events", len(np.unique(a["evt"])))
for p in np.unique(a["pdg"]):
    s = a[a["pdg"] == p]
    print("pdg", p, "n", len(s), "Emed", np.median(s["E"]), "vols", dict(zip(*np.unique(s["vol"], return_counts=True))), "procs", dict(zip(*np.unique(s["proc"], return_counts=True))))
