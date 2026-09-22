"""Таблицы вкладов по полосам: компоненты и объёмы рождения (из слитых файлов). Использование: band_tables.py <prefix> <out.md>"""
import sys, csv, numpy as np
pre, out = sys.argv[1], sys.argv[2]
B = [(20, 60), (60, 100), (100, 300), (300, 500), (500, 522), (522, 1000), (1000, 1500), (1500, 3000), (3000, 10000)]
comps = ["neutron", "proton", "mup", "mum", "em", "ep", "gamma", "k40"]
def rate(fn, lo, hi):
    d = np.genfromtxt(fn, delimiter=",", skip_header=2)
    m = (d[:, 0] >= lo) & (d[:, 0] < hi)
    return d[m, 1].sum() * 3600
RU = {"neutron": "нейтроны", "proton": "протоны", "mup": "μ⁺", "mum": "μ⁻", "em": "электроны", "ep": "позитроны", "gamma": "γ-кванты", "k40": "K-40 в людях", "Blanket": "Изоляция", "Cargo": "Багаж", "Floor": "Пол", "Pax": "Люди", "Skin": "Обшивка", "Trim": "Отделка", "Tube": "Воздух у прибора", "World": "Воздух вне конструкции", "other": "Прочее"}
rows = ["| полоса, кэВ | всего, отсч./ч | " + " | ".join(RU[c] for c in comps) + " |", "|" + "---|" * (len(comps) + 2)]
for lo, hi in B:
    t = rate(pre + "_total_all.csv", lo, hi)
    rows.append(f"| {lo}–{hi} | {t:.0f} | " + " | ".join(f"{100*rate(pre+'_total_'+c+'.csv', lo, hi)/t:.1f} %" for c in comps) + " |")
cats = list(csv.DictReader(open(pre + "_cat_all.csv", encoding="utf-8")))
vols = sorted({c["volname"] for c in cats})
rows += ["", "| полоса, кэВ | " + " | ".join(RU.get(v, v) for v in vols) + " |", "|" + "---|" * (len(vols) + 1)]
for lo, hi in B:
    v = {x: sum(sum(float(c[f"bin_{i}"]) for i in range(lo, hi)) for c in cats if c["volname"] == x) for x in vols}
    tot = sum(v.values())
    rows.append(f"| {lo}–{hi} | " + " | ".join(f"{100*v[x]/tot:.1f} %" for x in vols) + " |")
open(out, "w", encoding="utf-8", newline="\n").write("\n".join(rows) + "\n")
