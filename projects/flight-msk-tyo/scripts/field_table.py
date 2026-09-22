"""Состав поля PARMA на эшелоне: поток и доли по энергиям для каждой частицы (markdown). Использование: field_table.py <префикс таблиц> <out.md>"""
import sys, numpy as np
P, out = sys.argv[1], sys.argv[2]
C = [("нейтроны", 0), ("протоны", 1), ("μ⁺", 29), ("μ⁻", 30), ("электроны", 31), ("позитроны", 32), ("γ-кванты", 33)]
B = [(0, 1e-6, "< 1 эВ"), (1e-6, 0.1, "1 эВ – 100 кэВ"), (0.1, 10, "0,1–10 МэВ"), (10, 1000, "10 МэВ – 1 ГэВ"), (1000, 1e9, "> 1 ГэВ")]
f = lambda v: f"{v:.3g}".replace(".", ",")
rows, tot = [], 0
for name, ip in C:
    L = open(f"{P}{ip}.tab").read().split("\n")
    h = L[1].split(); nb = int(h[0]); ie = int(h[2]); line = float(h[4]); T = float(h[5])
    e = np.array(L[2].split(), float); D = np.array([np.array(L[3 + k].split(), float).sum() for k in range(nb)]) * np.diff(e)
    if ie > 0: D[ie - 1] += line
    D *= T / D.sum(); em = np.sqrt(e[:-1] * e[1:])
    rows.append((name, T, [D[(em >= a) & (em < b)].sum() / T for a, b, _ in B], line / T if ie > 0 else 0)); tot += T
o = ["| частица | поток, см⁻²·с⁻¹ | доля в сумме | " + " | ".join(x[2] for x in B) + " |", "|" + "---|" * (len(B) + 3)]
for name, T, fr, l5 in rows:
    o.append(f"| {name} | {f(T)} | {f"{100*T/tot:.1f}".replace(".", ",")} % | " + " | ".join(f"{100*x:.1f} %".replace(".", ",") for x in fr) + " |")
o.append(f"| **всего** | **{f(tot)}** | 100 % | " + " | ".join("" for _ in B) + " |")
o.append(""); o.append("Линия аннигиляции 511 кэВ в поле γ-квантов: " + f"{100*rows[-1][3]:.1f} %".replace(".", ",") + " потока γ.")
open(out, "w", encoding="utf-8", newline="\n").write("\n".join(o) + "\n"); print("\n".join(o))
