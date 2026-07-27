# -*- coding: utf-8 -*-
"""Поправка на самопоглощение — одной формулой вместо таблицы.

Пик полного поглощения набирают только НЕрассеянные фотоны, поэтому падение
пика относительно случая без вещества есть средняя вероятность пройти пробу
без взаимодействия. Для равномерного распределения длин пути это

    f = (1 - exp(-mu*d)) / (mu*d),    mu = rho * (mu/rho)(E),

где d — эффективная толщина пробы, чисто геометрическая величина: она не должна
зависеть ни от энергии, ни от матрицы, ни от плотности. Проверка этого и есть
смысл подгонки: если один d описывает все 125 точек, поправка сводится к
формуле, а шум отдельных точек на высоких энергиях перестаёт мешать.
"""
import sys
import csv
import os

import numpy as np

# Корни путей — из переменных окружения (common/py/paths.py).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec
RESULTS = rcspec.RESULTS


def load_mu():
    # Таблица массовых коэффициентов ослабления — это РЕЗУЛЬТАТ (её считает
    # mucalc), поэтому она коммитится в results/. В каталоге расчётов её
    # ищем только как запасной вариант, для свежепосчитанной.
    p = os.path.join(RESULTS, "mu.csv")
    if not os.path.exists(p):
        p = os.path.join(str(paths.build("RadiaCode-103")), "mu.csv")
    mu = {}
    for r in csv.DictReader(l for l in open(p, encoding="utf-8")
                            if not l.startswith("#")):
        mu[(r["matrix"], round(float(r["rho_gcm3"]), 4),
            round(float(r["E_keV"]), 1))] = float(r["mu_over_rho_cm2_g"])
    return mu


def load_eff():
    p = rcspec.rdir("efficiency.csv")
    rows = list(csv.DictReader(open(p, encoding="utf-8")))
    for r in rows:
        r["rho"] = float(r["rho_gcm3"])
        r["E"] = float(r["E_keV"])
        r["eps"] = float(r["eps_p"])
        r["d_eps"] = float(r["d_eps_p"])
    return rows


def model(mu_d):
    """(1-exp(-x))/x, устойчиво при малых x."""
    out = np.where(mu_d > 1e-6, (1.0 - np.exp(-mu_d)) / np.maximum(mu_d, 1e-12),
                   1.0 - 0.5 * mu_d)
    return out


def main():
    mu_tab = load_mu()
    rows = load_eff()

    ref = {r["E"]: (r["eps"], r["d_eps"]) for r in rows if r["matrix"] == "air"}
    E, MU, F, W = [], [], [], []
    for r in rows:
        if r["matrix"] == "air" or r["E"] not in ref:
            continue
        e0, d0 = ref[r["E"]]
        if e0 <= 0 or r["eps"] <= 0:
            continue
        f = r["eps"] / e0
        # относительная ошибка отношения
        rel = np.hypot(r["d_eps"] / r["eps"], d0 / e0)
        key = (r["matrix"], round(r["rho"], 4), round(r["E"], 1))
        if key not in mu_tab:
            continue
        E.append(r["E"])
        MU.append(mu_tab[key] * r["rho"])       # 1/см
        F.append(f)
        W.append(1.0 / max(rel * f, 1e-6) ** 2)
    E, MU, F, W = map(np.array, (E, MU, F, W))
    print("точек для подгонки: %d" % len(E))

    # одномерный поиск d по взвешенному хи-квадрат
    ds = np.linspace(0.2, 8.0, 2000)
    chi = [(W * (model(MU * d) - F) ** 2).sum() for d in ds]
    d_best = ds[int(np.argmin(chi))]
    pred = model(MU * d_best)
    resid = (pred - F) / F
    print("эффективная толщина d = %.3f см" % d_best)
    print("невязка: среднее %+.1f %%, ско %.1f %%, максимум %.1f %%"
          % (100 * resid.mean(), 100 * resid.std(), 100 * np.abs(resid).max()))

    # проверка независимости d от матрицы и энергии
    print("\n%-22s %8s %8s" % ("группа", "d, см", "ско, %"))
    for name, m in [("E < 150 кэВ", E < 150), ("E 150..600", (E >= 150) & (E < 600)),
                    ("E > 600 кэВ", E >= 600)]:
        if m.sum() < 3:
            continue
        c = [(W[m] * (model(MU[m] * d) - F[m]) ** 2).sum() for d in ds]
        db = ds[int(np.argmin(c))]
        rr = (model(MU[m] * db) - F[m]) / F[m]
        print("%-22s %8.3f %8.1f" % (name, db, 100 * rr.std()))

    out = rcspec.rdir("selfabsorption_fit.txt")
    with open(out, "w", encoding="utf-8") as f:
        f.write("Поправка на самопоглощение в сосуде Маринелли 200 мл\n\n")
        f.write("  eps_p(матрица, rho, E) = eps_p_предел(E) * f,\n")
        f.write("  f = (1 - exp(-mu*d)) / (mu*d),  mu = rho * (mu/rho)(E), 1/см\n")
        f.write("  d = %.3f см   (эффективная толщина пробы)\n\n" % d_best)
        f.write("Подгонка по %d точкам: ско %.1f %%, максимум %.1f %%\n"
                % (len(E), 100 * resid.std(), 100 * np.abs(resid).max()))
        f.write("Проверено на матрицах органика/грунт/вода при плотности "
                "0.5..1.6 г/см³ и энергиях 30..3000 кэВ.\n")
    print("\n", out)


if __name__ == "__main__":
    main()
