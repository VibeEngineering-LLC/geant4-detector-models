"""График кривых эффективной площади из results/eff/eff_curves.csv -> results/eff/eff_curves.svg"""
import csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding="utf-8")
N = {"gamma": "γ-кванты", "neutron": "нейтроны", "proton": "протоны", "mum": "μ⁻", "mup": "μ⁺", "em": "e⁻", "ep": "e⁺"}
C = {"gamma": "#c0392b", "neutron": "#27ae60", "proton": "#2980b9", "mum": "#8e44ad", "mup": "#d35400", "em": "#16a085", "ep": "#7f8c8d"}
d = {}
for r in csv.DictReader(open("results/eff/eff_curves.csv", encoding="utf-8")):
    d.setdefault(r["comp"], []).append((float(r["E_keV"]), max(float(r["A_eff_cm2"]), 1e-4)))
plt.rcParams["svg.fonttype"] = "none"
fig, ax = plt.subplots(figsize=(9, 5))
for c, v in d.items():
    ax.loglog([x[0] for x in v], [x[1] for x in v], "o-", ms=3, lw=1.5, color=C[c], label=N[c])
ax.axhline(11.25, ls=":", color="k", lw=1)
ax.set_ylim(1e-3, 1e2); ax.set_xlabel("энергия частицы, кэВ"); ax.set_ylabel("эффективная площадь, см²")
ax.grid(True, which="both", alpha=.3); ax.legend(ncol=4, fontsize=9)
ax.set_title("Окно 20–10 000 кэВ; пунктир — средняя проекция кристалла 11,25 см²", fontsize=10)
fig.tight_layout(); fig.savefig("results/eff/eff_curves.svg"); print("svg ok")
