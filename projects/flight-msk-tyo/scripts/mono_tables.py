"""Таблицы моноэнергетического изотропного источника для CABIN_DIRECT (поток 1 см-2 с-1): tables/mono/<name>_<E>.tab + список заданий results/eff/jobs.lst."""
import numpy as np
SETS = {"gamma": (22, 0.02, 10), "neutron": (2112, 1e-8, 1000), "em": (11, 0.05, 1000), "ep": (-11, 0.05, 1000),
        "mum": (13, 1, 10000), "mup": (-13, 1, 10000), "proton": (2212, 1, 10000)}
NPTS = {"gamma": 16, "neutron": 16}
jobs = []
for name, (pdg, lo, hi) in SETS.items():
    for E in np.geomspace(lo, hi, NPTS.get(name, 12)):
        fn = f"tables/mono/{name}_{E:.4g}.tab"
        with open(fn, "w", newline="\n") as f:
            f.write(f"# mono {name} E={E:.6g} MeV isotropic\n1 100 0 0 0 1.0\n{E*0.9995:.9e} {E*1.0005:.9e}\n")
            f.write(" ".join(["1.0"] * 100) + "\n")
        jobs.append(f"{name} {pdg} {E:.6g} {fn}")
open("results/eff/jobs.lst", "w", newline="\n").write("\n".join(jobs) + "\n")
print(len(jobs), "заданий")
