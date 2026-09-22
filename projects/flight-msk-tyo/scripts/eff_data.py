"""A_eff(E) = скорость счёта 20-10000 кэВ при единичной плотности потока (изотропно, R=10 см) -> results/eff/eff_curves.csv"""
import glob, os, re, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
D = "results/eff"; S = 11.25  # средняя проекция кристалла 18x15x60 мм = S/4, см2
rows = []
for f in glob.glob(f"{D}/*_total.csv"):
    m = re.match(r"(gamma|neutron|proton|mum|mup|em|ep)_([0-9.e+-]+)_total\.csv$", os.path.basename(f))
    if not m:
        continue
    d = np.genfromtxt(f, delimiter=",", skip_header=2)
    a = d[20:10000, 1].sum(); e = np.sqrt(d[20:10000, 2].sum())
    rows.append((m.group(1), float(m.group(2)) * 1000.0, a, e, a / S))
rows.sort()
out = ["comp,E_keV,A_eff_cm2,sigma_cm2,eps"] + [f"{c},{E:.6g},{a:.6g},{e:.3g},{a / S:.5g}" for c, E, a, e, _ in rows]
open(f"{D}/eff_curves.csv", "w", encoding="utf-8").write("\n".join(out) + "\n")
print(len(rows), "точек")
