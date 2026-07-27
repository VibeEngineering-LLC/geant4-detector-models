# -*- coding: utf-8 -*-
"""Площадь пика подгонкой «гауссиана + линейная подложка».

Оконные суммы зависят от выбора боковых окон (у Cs-137 слева комптоновский
континуум, справа почти пусто, и линейная интерполяция между ними завышает
подложку). Подгонка снимает этот произвол: подложка и пик определяются
одновременно из одних данных.
"""
import sys
import csv
import os

import numpy as np
from scipy.optimize import curve_fit

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
BG = os.path.join(BASE, "Фон домик 23 дня.xml")

E0_CS, YIELD_CS = 661.657, 0.851
APP_BQ_PER_KG = 3340.0     # приложение; его калибровка сверена с АТОМТЕХ (0.5 %)


def rebin_to(spec, other):
    def edges(s):
        e = s.energy
        return np.concatenate(([e[0] - (e[1] - e[0]) / 2],
                               (e[:-1] + e[1:]) / 2,
                               [e[-1] + (e[-1] - e[-2]) / 2]))
    cum = np.concatenate(([0.0], np.cumsum(other.counts)))
    return np.diff(np.interp(edges(spec), edges(other), cum))


def peak_model(x, area, mu, sigma, b0, b1):
    """Пик задан ПЛОЩАДЬЮ, а не высотой: тогда её ошибка выходит из подгонки."""
    g = area / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    return g + b0 + b1 * (x - mu)


def eps_at(E0, vessel, cfg):
    E, ep = [], []
    for r in csv.DictReader(open(rcspec.rdir("efficiency.csv", v=vessel),
                                 encoding="utf-8")):
        if r["config"] == cfg:
            E.append(float(r["E_keV"]))
            ep.append(float(r["eps_p"]))
    o = np.argsort(E)
    E, ep = np.array(E)[o], np.array(ep)[o]
    return float(np.exp(np.interp(np.log(E0), np.log(E), np.log(ep))))


def main():
    smp = read_rcxml.read(SAMPLE)[0]
    bg = read_rcxml.read(BG)[0]
    bgs = rebin_to(smp, bg) * (smp.live / bg.live)
    net = smp.counts - bgs
    e = smp.energy
    # на канал, чтобы подложка была линейной по энергии
    dE = np.gradient(e)

    fit = (e > 560) & (e < 790)
    x, y = e[fit], net[fit] / dE[fit]
    err = np.sqrt(np.maximum(smp.counts[fit] + bgs[fit], 1)) / dE[fit]

    p0 = [net[fit].sum(), E0_CS, 0.084 * 662 / 2.355, y.min(), 0.0]
    p, cov = curve_fit(peak_model, x, y, p0=p0, sigma=err, absolute_sigma=True)
    d = np.sqrt(np.diag(cov))
    area, mu, sigma = p[0], p[1], abs(p[2])
    resid = (y - peak_model(x, *p)) / err
    print("подгонка в окне 560..790 кэВ, %d точек" % fit.sum())
    print("  площадь пика   %.0f ± %.0f имп" % (area, d[0]))
    print("  центр          %.2f ± %.2f кэВ  (табл. %.2f)" % (mu, d[1], E0_CS))
    print("  сигма          %.2f ± %.2f кэВ => FWHM %.1f %% на 662"
          % (sigma, d[2], 100 * sigma * 2.35482 / mu))
    print("  хи-квадрат/степень свободы %.2f"
          % ((resid ** 2).sum() / (len(x) - len(p))))

    eps = eps_at(E0_CS, "m500", "full_organic_0.50")
    A = area / (eps * YIELD_CS * smp.live)
    a = A / smp.weight
    print("\nпо расчётной кривой eps_p = %.4e:" % eps)
    print("  активность %.0f Бк, удельная %.0f Бк/кг" % (A, a))
    print("  приложение (сверено с АТОМТЕХ): %.0f Бк/кг" % APP_BQ_PER_KG)
    print("  МОЯ МОДЕЛЬ ЗАВЫШАЕТ eps_p в %.3f раза" % (APP_BQ_PER_KG / a))
    print("  => нормировочный множитель к расчёту: %.3f" % (a / APP_BQ_PER_KG))
    print("\nэффективность, восстановленная из измерения: %.4e" % (eps * a / APP_BQ_PER_KG))


if __name__ == "__main__":
    main()
