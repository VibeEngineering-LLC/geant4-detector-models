"""Таблица A_eff и ε=A_eff/11,25 см² по сортам излучения на опорных энергиях -> results/final/eff.md"""
import csv, sys
sys.stdout.reconfigure(encoding="utf-8")
N = {"gamma": "γ-кванты", "neutron": "нейтроны", "proton": "протоны", "mum": "μ⁻", "mup": "μ⁺", "em": "e⁻", "ep": "e⁺"}
G = [(0.025, "25 кэВ"), (0.1, "100 кэВ"), (0.5, "500 кэВ"), (1, "1 МэВ"), (3, "3 МэВ"), (10, "10 МэВ"), (100, "100 МэВ"), (1000, "1 ГэВ")]
d = {}
for r in csv.DictReader(open("results/eff/eff_curves.csv", encoding="utf-8")):
    d.setdefault(r["comp"], []).append((float(r["E_keV"]) / 1000, float(r["A_eff_cm2"])))
out = ["| Излучение | " + " | ".join(g[1] for g in G) + " |", "|---|" + "---|" * len(G)]
for c, n in N.items():
    cells = []
    for e, _ in G:
        p = min(d[c], key=lambda x: abs(x[0] / e - 1) if x[0] > e else abs(e / x[0] - 1))
        ok = max(p[0] / e, e / p[0]) < 1.6
        cells.append(f"{p[1]:.2f}" if ok else "—")
    out.append(f"| {n} | " + " | ".join(cells) + " |")
open("results/final/eff.md", "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")
print("\n".join(out))
