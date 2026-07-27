# -*- coding: utf-8 -*-
"""Проникающая способность бета-излучения пробы.

Путь от пробы до кристалла складывается из гильзы колодца (1.25 мм PLA),
воздушного зазора, стенки корпуса (1.50 мм ABS) и отражающей чашки (1.25 мм),
что даёт около 0.49 г/см². По практическому пробегу R = 0.412*E^1.27 (г/см²,
Кац-Пенфолд) это соответствует Emax порядка 1.15 МэВ: электроны слабее этого
до кристалла не доходят вовсе.

Скрипт строит зависимость доли дошедших от энергии и от плотности пробы
(самопоглощение беты в самой пробе) и показывает, где лежит порог.
"""
import os

import numpy as np

import sys
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = rcspec.RESULTS
THRESH = 20.0        # кэВ, порог счёта прибора

# Бета-излучатели, встречающиеся в пробах ЕРН и техногенных
EMITTERS = [
    ("Cs-137", 514, "не доходит"),
    ("Pb-214", 1020, "на пороге"),
    ("K-40", 1311, "доходит краем"),
    ("Tl-208", 1803, "доходит"),
    ("Ac-228", 2076, "доходит"),
    ("Pa-234m", 2269, "доходит"),
    ("Y-90", 2280, "доходит"),
    ("Bi-214", 3272, "доходит уверенно"),
]


def katz_penfold(mass_thickness=0.49):
    """Обратная задача: при какой Emax практический пробег равен заданной
    массовой толщине. R = 0.412 * E^1.27 г/см² для 0.01 < E < 2.5 МэВ."""
    return (mass_thickness / 0.412) ** (1.0 / 1.27)


def main():
    base = rcspec.rdir("beta")
    if not os.path.isdir(base):
        raise SystemExit("нет результатов в " + base)

    print("порог по массовой толщине 0.49 г/см² => Emax = %.2f МэВ (Кац-Пенфолд)"
          % katz_penfold())

    data = {}
    for cfg in sorted(os.listdir(base)):
        d = os.path.join(base, cfg)
        if not os.path.isdir(d):
            continue
        rows = []
        for fn in sorted(os.listdir(d)):
            if not fn.startswith("E") or not fn.endswith(".csv"):
                continue
            meta, hist = rcspec.read_spec(os.path.join(d, fn))
            n = float(meta["N_primaries"])
            cnt = hist[int(THRESH):].sum()
            edep = (hist * (np.arange(len(hist)) + 0.5)).sum()
            rows.append(dict(E=float(meta["E_prim_keV"]), n=n,
                             frac=cnt / n, d=np.sqrt(max(cnt, 1)) / n,
                             edep=edep / n))
        rows.sort(key=lambda r: r["E"])
        if rows:
            data[cfg] = rows

    cfgs = sorted(data)
    print("\nДоля электронов пробы, давших сигнал в кристалле")
    print("%9s" % "Emax,кэВ", end="")
    for c in cfgs:
        print(" %16s" % c.replace("full_", ""), end="")
    print()
    energies = [r["E"] for r in data[cfgs[0]]]
    for i, E in enumerate(energies):
        print("%9.0f" % E, end="")
        for c in cfgs:
            r = data[c][i]
            print("  %.3e±%.0e" % (r["frac"], r["d"]), end="")
        print()

    out = rcspec.rdir("beta_transmission.csv")
    with open(out, "w", encoding="utf-8") as f:
        f.write("config,E_keV,frac,d_frac,mean_edep_keV\n")
        for c in cfgs:
            for r in data[c]:
                f.write("%s,%.1f,%.6e,%.2e,%.4f\n"
                        % (c, r["E"], r["frac"], r["d"], r["edep"]))
    print("\nтаблица:", out)

    print("\nБета-излучатели проб (Emax, доходит ли до кристалла)")
    ref = data.get("full_water_1.00") or data[cfgs[0]]
    xs = np.array([r["E"] for r in ref])
    ys = np.array([r["frac"] for r in ref])
    for nuc, emax, note in EMITTERS:
        v = np.interp(emax, xs, ys)
        print("  %-9s %5d кэВ  доля %.2e  %s" % (nuc, emax, v, note))


if __name__ == "__main__":
    main()
