"""Сравнение r900base и r900fuel по полосам -> results/final/fuel.md (отношение с ошибкой по Σw²)."""
import math, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
B = [(20, 100), (100, 300), (300, 1000), (1000, 3000), (3000, 10000)]
N = [("neutron", "нейтроны"), ("gamma", "γ-кванты"), ("em", "электроны"), ("ep", "позитроны"), ("all", "все считавшиеся компоненты")]
rd = lambda h, c: np.genfromtxt(f"results/final/r900{h}_total_{c}.csv", delimiter=",", skip_header=2)
f = lambda a, lo, hi: (a[lo:hi, 1].sum(), math.sqrt(a[lo:hi, 2].sum()))
out = ["| Составляющая | " + " | ".join(f"{lo}–{hi} кэВ" for lo, hi in B) + " | 20–10 000 кэВ |", "|---|" + "---|" * (len(B) + 1)]
for c, n in N:
    a, b = rd("base", c), rd("fuel", c)
    cells = []
    for lo, hi in B + [(20, 10000)]:
        s0, e0 = f(a, lo, hi); s1, e1 = f(b, lo, hi)
        r = s1 / s0; er = r * math.hypot(e0 / s0, e1 / s1)
        cells.append(f"{r:.2f} ± {er:.2f}")
    out.append(f"| {n} | " + " | ".join(cells) + " |")
open("results/final/fuel.md", "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")
print("\n".join(out))
