# -*- coding: utf-8 -*-
"""Сверка расчёта с реальным измерением: черника в маринелли 500 мл, RC-103.

Определяем активность Cs-137 по площади пика 662 и посчитанной эффективности,
и сравниваем с тем, что выдало приложение. Это проверка всей цепочки: промер
STL, геометрия прибора, эффективность, самопоглощение.
"""
import csv
import os
import sys

import numpy as np

# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec
import read_rcxml

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'common', 'py'))
import paths  # noqa: E402


BASE = str(paths.measured("RadiaCode-103"))  # измеренные спектры прибора
SAMPLE = os.path.join(BASE, "RC103 черника маринелли авторская домик 246 гр.xml")
BG_SHIELD = os.path.join(BASE, "Фон домик 23 дня.xml")

E0 = 661.657          # кэВ, Cs-137
YIELD = 0.851
FWHM_REL = 0.084      # паспорт RC-103 на 662 кэВ
APP_BQ_PER_KG = 3340.0   # что показало приложение

# окна: пик и две подложки по бокам (кэВ)
PEAK = 2.5            # полуширина окна в сигмах
SIDE_LO = (520.0, 575.0)
SIDE_HI = (745.0, 800.0)


def rebin_to(spec, other):
    """Приведение спектра other к энергетической шкале spec по каналам."""
    e_t = spec.energy
    edges_t = np.concatenate(([e_t[0] - (e_t[1] - e_t[0]) / 2],
                              (e_t[:-1] + e_t[1:]) / 2,
                              [e_t[-1] + (e_t[-1] - e_t[-2]) / 2]))
    e_o = other.energy
    edges_o = np.concatenate(([e_o[0] - (e_o[1] - e_o[0]) / 2],
                              (e_o[:-1] + e_o[1:]) / 2,
                              [e_o[-1] + (e_o[-1] - e_o[-2]) / 2]))
    # накопленная сумма позволяет переложить в другие границы без предположений
    cum = np.concatenate(([0.0], np.cumsum(other.counts)))
    return np.diff(np.interp(edges_t, edges_o, cum))


def window(spec, lo, hi):
    e = spec.energy
    return (e >= lo) & (e <= hi)


def eps_p_662(vessel="m500", cfg="full_organic_0.50"):
    for r in csv.DictReader(open(rcspec.rdir("efficiency.csv", v=vessel),
                                 encoding="utf-8")):
        if r["config"] == cfg and abs(float(r["E_keV"]) - 661.7) < 1:
            return float(r["eps_p"]), float(r["d_eps_p"])
    raise SystemExit("нет eps_p(662) для " + cfg)


def main():
    smp = read_rcxml.read(SAMPLE)[0]
    bg = read_rcxml.read(BG_SHIELD)[0]
    print("проба: живое %d с, %d отсчётов, %.3f имп/с"
          % (smp.live, smp.counts.sum(), smp.counts.sum() / smp.live))
    print("фон:   живое %d с, %d отсчётов, %.3f имп/с"
          % (bg.live, bg.counts.sum(), bg.counts.sum() / bg.live))

    # фон приводим к шкале пробы и к её времени
    bg_on_smp = rebin_to(smp, bg) * (smp.live / bg.live)
    net = smp.counts - bg_on_smp
    print("вычтено фона: %.0f отсчётов, остаток пробы %.0f (%.3f имп/с)"
          % (bg_on_smp.sum(), net.sum(), net.sum() / smp.live))

    sigma = FWHM_REL * E0 / 2.35482
    lo, hi = E0 - PEAK * sigma, E0 + PEAK * sigma
    mp = window(smp, lo, hi)
    ml = window(smp, *SIDE_LO)
    mh = window(smp, *SIDE_HI)
    print("\nокно пика %.0f..%.0f кэВ (%d каналов), sigma %.1f кэВ"
          % (lo, hi, mp.sum(), sigma))

    e = smp.energy
    # линейная подложка по средним уровням боковых окон (на канал)
    yl, yh = net[ml].mean(), net[mh].mean()
    xl, xh = e[ml].mean(), e[mh].mean()
    base = yl + (yh - yl) * (e[mp] - xl) / (xh - xl)
    gross = net[mp].sum()
    cont = base.sum()
    area = gross - cont
    # доля гауссианы вне окна +-2.5 sigma
    frac = 0.9876
    area_full = area / frac
    d_area = np.sqrt(smp.counts[mp].sum() + bg_on_smp[mp].sum()
                     + cont) / frac

    print("подложка слева %.2f, справа %.2f отсч/канал" % (yl, yh))
    print("в окне: всего %.0f, подложка %.0f, ЧИСТАЯ ПЛОЩАДЬ %.0f ± %.0f"
          % (gross, cont, area, d_area))
    print("полная площадь пика (с поправкой на хвосты): %.0f ± %.0f"
          % (area_full, d_area))

    eps, d_eps = eps_p_662()
    A = area_full / (eps * YIELD * smp.live)
    dA_stat = d_area / (eps * YIELD * smp.live)
    dA_eps = A * d_eps / eps
    a = A / smp.weight
    print("\neps_p(662, органика 0.50, m500) = %.4e ± %.1f %%"
          % (eps, 100 * d_eps / eps))
    print("АКТИВНОСТЬ ПО РАСЧЁТУ: %.0f ± %.0f Бк  (стат. %.0f, эффект. %.0f)"
          % (A, np.hypot(dA_stat, dA_eps), dA_stat, dA_eps))
    print("удельная: %.0f Бк/кг при массе %.3f кг" % (a, smp.weight))
    print("\nПРИЛОЖЕНИЕ: %.0f Бк/кг" % APP_BQ_PER_KG)
    print("расхождение: %+.1f %%" % (100 * (a / APP_BQ_PER_KG - 1)))


if __name__ == "__main__":
    main()
