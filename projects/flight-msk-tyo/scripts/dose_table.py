"""Таблицы дозы фантома из results/final/dose_all.csv -> results/final/dose.md (мкГр/ч и мкЗв/ч по слоям и составляющим поля; по классу частицы у тела)."""
import csv, math, sys
sys.stdout.reconfigure(encoding="utf-8")
H = 3600e6   # Гр/с -> мкГр/ч (1 Гр = 1e6 мкГр)
LAY = ["кожа (0,07 мм)", "глубина 10 мм", "глубина 30 мм", "весь шар (среднее)"]
GR = [("neutron", ["neutron"], "нейтроны"), ("proton", ["proton"], "протоны"), ("mu", ["mup", "mum"], "мюоны"), ("e", ["em", "ep"], "e⁻ и e⁺"), ("gamma", ["gamma"], "γ-кванты")]
# K-40 в дозу пассажира не включён: это γ-излучение калия-40 в теле ДРУГИХ людей — фоновый источник,
# не зависящий от высоты полёта, а не часть дозы от космического излучения (решение оператора 22.09.2026)
R = list(csv.DictReader(open("results/final/dose_all.csv", encoding="utf-8")))
def s(fn, keys, lay, cls=None):
    v = [r for r in R if r["comp"] in keys and int(r["layer"]) == lay and (cls is None or int(r["class"]) == cls)]
    return sum(float(r[fn]) for r in v) * H, math.sqrt(sum(float(r["sig_eq"]) ** 2 for r in v)) * H
def tab(fn, name, note):
    o = [f"**{name}**", "", note, "", "| Слой | " + " | ".join(g[2] for g in GR) + " | Итого |", "|---|" + "---|" * (len(GR) + 1)]
    for l, ln in enumerate(LAY):
        c = [s(fn, g[1], l) for g in GR]; tot = sum(x[0] for x in c); et = math.sqrt(sum(x[1] ** 2 for x in c))
        fmt = lambda v: f"{v:.4f}".rstrip("0").rstrip(".") if v < 0.01 else f"{v:.3g}"
        o.append(f"| {ln} | " + " | ".join(fmt(x[0]) for x in c) + f" | **{tot:.3g}** ± {et:.2g} |")
    return o
n1 = "Поглощённая доза — энергия, выделенная в ткани, без поправки на вид излучения (1 мкГр/ч = 0,001 мГр/ч)."
n2 = "Эквивалентная доза — поглощённая доза, умноженная на коэффициент wR (ICRP 103): у нейтронов и протонов он больше 1, поэтому их доля в эквивалентной дозе выше, чем в поглощённой."
o = tab("abs_Gy_s", "Поглощённая доза, мкГр/ч", n1) + [""] + tab("eq_Sv_s", "Эквивалентная доза (wR по ICRP 103), мкЗв/ч", n2)
open("results/final/dose.md", "w", encoding="utf-8", newline="\n").write("\n".join(o) + "\n"); print("\n".join(o))
