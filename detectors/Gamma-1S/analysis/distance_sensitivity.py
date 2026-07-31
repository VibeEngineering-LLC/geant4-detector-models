# -*- coding: utf-8 -*-
"""Чувствительность эффективности к расстоянию источник-кристалл — АНАЛИТИКА.

ЗАЧЕМ. §5.4 отчёта отбраковывает гипотезу «постоянное смещение расстояния»,
сравнивая наблюдаемое отношение остатков (1,30) с ТРЕБУЕМЫМ — отношением
чувствительностей на 5 и 25 см. Использованные там числа 3,05/0,89 %/мм
унаследованы из более раннего расчёта БЕЗ производящего скрипта в дереве
(найдено вторым независимым аудитом 31.07.2026, task 95). Этот скрипт даёт
НЕЗАВИСИМУЮ, воспроизводимую оценку той же величины аналитической формулой —
не замену унаследованных чисел (они получены другим путём, возможно, полным
Монте-Карло с реальной физикой поглощения, а не голым телесным углом), а
независимую проверку порядка величины.

ФОРМУЛА. Телесный угол диска радиуса a на оси на расстоянии d от источника:
    Omega/4pi(d) = (1 - d/sqrt(d^2+a^2)) / 2
Это ВЕРХНЯЯ ГРАНКА геометрической части эффективности (без поглощения в торце
и в самом кристалле) — но производная по d в ОТНОСИТЕЛЬНЫХ единицах близка к
производной полной эффективности, потому что поглощение в первом приближении
зависит от угла падения, а не от d напрямую, при фиксированной геометрии торца.
Аналитическая производная:
    d(Omega/4pi)/dd = -a^2 / (2*(d^2+a^2)^1.5)
Радиус кристалла `a` берётся НЕ числом здесь, а вычисляется из cryDia,
прочитанного прямо из G1SDetector.hh — та же дисциплина, что стамп
провенанса: константа геометрии живёт в ОДНОМ месте (чертёж), скрипт её не
дублирует.

ПЛОСКОСТЬ ОТСЧЁТА — подтверждено оператором 31.07.2026: «5 см»/«25 см»
относятся к НАРУЖНОЙ плоскости корпуса, не к кристаллу. Это ровно то, что уже
разыгрывает геометрия (`run_all_grids.py`: `ZFACE + dist`, `ZFACE = 41,0 мм`
от центра кристалла — сумма стека торца). Расстояние ДО ПЛОСКОСТИ ДИСКА
(кристалл, полутолщина 31,5 мм) поэтому на `ZFACE − cryLen/2 = 9,5 мм` больше
заявленного: 59,5 мм и 259,5 мм, а не 50/250. Прежняя редакция считала по
50/250 — ошибка, исправлена.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

RESULTS = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results"))

# ZFACE — та же константа, что в run_all_grids.py (наружная плоскость торца
# = 41,0 мм от центра кристалла, сумма стека Al 2,0/воздух 1,0/банка Al 0,5/
# MgO 6,0 + полутолщина кристалла 31,5). Расстояние до ПЛОСКОСТИ ДИСКА
# (кристалл) = ZFACE + dist - cryLen/2, не dist само по себе.
ZFACE = 41.0
D_STATED = {"точечная 5 см": 50.0, "точечная 25 см": 250.0}


def _hh_field(name):
    """Число поля HeadGeom из G1SDetector.hh — не копия, а чтение чертежа."""
    hh = os.path.join(str(paths.geometry("Gamma-1S")), "G1SDetector.hh")
    text = open(hh, encoding="utf-8").read()
    m = re.search(name + r"\s*=\s*([\d.]+)", text)
    if not m:
        raise SystemExit("distance_sensitivity: не нашёл %s в %s" % (name, hh))
    return float(m.group(1))


def crystal_radius_mm():
    return _hh_field("cryDia") / 2.0


def crystal_half_len_mm():
    return _hh_field("cryLen") / 2.0


def omega_over_4pi(d, a):
    return (1.0 - d / math.sqrt(d * d + a * a)) / 2.0


def sensitivity_pct_per_mm(d, a):
    """|d(Omega/4pi)/dd| / (Omega/4pi) * 100 — аналитическая производная."""
    deriv = -a * a / (2.0 * (d * d + a * a) ** 1.5)
    return -deriv / omega_over_4pi(d, a) * 100.0


if __name__ == "__main__":
    a = crystal_radius_mm()
    half_len = crystal_half_len_mm()
    offset = ZFACE - half_len
    D_POINTS = {lbl: dist + offset for lbl, dist in D_STATED.items()}
    print("Радиус кристалла (из G1SDetector.hh): a = %.2f мм" % a)
    print("Плоскость отсчёта -> плоскость диска: ZFACE(%.1f) - cryLen/2(%.1f)"
          " = %.1f мм\n" % (ZFACE, half_len, offset))
    print("%-16s %10s %8s %14s %16s"
          % ("геометрия", "заявлено", "d, мм", "Omega/4pi", "%/мм"))
    rows = []
    for label, d in D_POINTS.items():
        f = omega_over_4pi(d, a)
        s = sensitivity_pct_per_mm(d, a)
        print("%-16s %10.1f %8.1f %14.6f %16.3f"
              % (label, D_STATED[label], d, f, s))
        rows.append((label, D_STATED[label], d, a, f, s))

    s5 = sensitivity_pct_per_mm(D_POINTS["точечная 5 см"], a)
    s25 = sensitivity_pct_per_mm(D_POINTS["точечная 25 см"], a)
    print("\nОтношение чувствительностей 5/25 см: %.2f" % (s5 / s25))

    OBS = {
        "quantity": "относительная чувствительность телесного угла диска"
                    " к расстоянию источник-кристалл; не эффективность"
                    " регистрации и не результат счёта",
        "area": "н/п — аналитическая формула; не съём площади пика",
        "window": "н/п",
        "shelf": "н/п — подложка не вычитается; расчёт не по спектру",
        "blurred": "н/п",
        "formula": "Omega/4pi(d)=(1-d/sqrt(d^2+a^2))/2; производная по d"
                   " аналитическая; не конечная разность",
        "caveat": "верхняя граница геометрической части эффективности;"
                  " поглощение в торце и кристалле не учтено",
    }
    csvio.write(
        os.path.join(RESULTS, "distance_sensitivity.csv"),
        ["geometry", "d_stated_mm", "d_to_disc_mm", "crystal_radius_mm",
         "omega_over_4pi", "sensitivity_pct_per_mm"],
        [(lbl, "%.1f" % ds, "%.1f" % d, "%.2f" % a, "%.6f" % f, "%.3f" % s)
         for lbl, ds, d, a, f, s in rows],
        comments=[
            "Аналитическая чувствительность телесного угла диска-кристалла к"
            " расстоянию источник-кристалл; НЕ замена унаследованным числам"
            " 3;05/0;89 %%/мм в report.md — независимая проверка порядка"
            " величины (см. докстринг скрипта).",
            "d_stated_mm — заявленное расстояние (5/25 см от НАРУЖНОЙ"
            " плоскости корпуса; подтверждено оператором 31.07.2026).",
            "d_to_disc_mm = d_stated + ZFACE - cryLen/2 — фактическое"
            " расстояние до плоскости кристалла; та же величина, что"
            " разыгрывает геометрия (run_all_grids.py: ZFACE + dist).",
            "sensitivity_pct_per_mm = |d(Omega/4pi)/dd| / (Omega/4pi) * 100;"
            " формула и вывод — в докстринге distance_sensitivity.py.",
            "Отношение чувствительностей 5/25 см = %.2f." % (s5 / s25),
        ],
        stamp=stamp.lines(
            "detectors/Gamma-1S/analysis/distance_sensitivity.py", OBS,
            geometry_dir=str(paths.geometry("Gamma-1S")),
            names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO)))
    print("\nтаблица: %s" % os.path.join(RESULTS, "distance_sensitivity.csv"))
