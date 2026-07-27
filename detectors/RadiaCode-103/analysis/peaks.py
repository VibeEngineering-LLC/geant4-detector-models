# -*- coding: utf-8 -*-
"""Площади пиков в измеренном спектре и активности по расчётным кривым.

Независимая проверка кривых внутри одного измерения: у пробы известно
содержание K-40, поэтому активность, восстановленная по линии 1461, проверяет
кривую там, где она никак не связана с определением цезия.
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

# линия, нуклид, выход на распад
LINES = [(661.657, "Cs-137", 0.851), (1460.82, "K-40", 0.1055)]
FWHM_REL = 0.084


def curve(vessel, cfg):
    E, ep = [], []
    for r in csv.DictReader(open(rcspec.rdir("efficiency.csv", v=vessel),
                                 encoding="utf-8")):
        if r["config"] == cfg:
            E.append(float(r["E_keV"]))
            ep.append(float(r["eps_p"]))
    o = np.argsort(E)
    return np.array(E)[o], np.array(ep)[o]


def eps_at(E0, Eg, epg):
    return float(np.exp(np.interp(np.log(E0), np.log(Eg), np.log(epg))))


def rebin_to(spec, other):
    def edges(s):
        e = s.energy
        return np.concatenate(([e[0] - (e[1] - e[0]) / 2],
                               (e[:-1] + e[1:]) / 2,
                               [e[-1] + (e[-1] - e[-2]) / 2]))
    cum = np.concatenate(([0.0], np.cumsum(other.counts)))
    return np.diff(np.interp(edges(spec), edges(other), cum))


def peak_area(spec, net, gross_raw, bg_scaled, E0, nsig=2.5, gap=1.6):
    """Площадь за вычетом линейной подложки по двум боковым окнам."""
    sigma = FWHM_REL * E0 / 2.35482 * np.sqrt(E0 / 662.0) / np.sqrt(E0 / 662.0)
    sigma = FWHM_REL * 662.0 / 2.35482 * np.sqrt(E0 / 662.0)   # FWHM ~ sqrt(E)
    e = spec.energy
    mp = (e >= E0 - nsig * sigma) & (e <= E0 + nsig * sigma)
    ml = (e >= E0 - (gap + 2.2) * sigma) & (e <= E0 - gap * sigma * 1.6)
    mh = (e >= E0 + gap * sigma * 1.6) & (e <= E0 + (gap + 2.2) * sigma)
    if ml.sum() < 2 or mh.sum() < 2:
        return None
    yl, yh = net[ml].mean(), net[mh].mean()
    xl, xh = e[ml].mean(), e[mh].mean()
    base = yl + (yh - yl) * (e[mp] - xl) / (xh - xl)
    gross = net[mp].sum()
    cont = base.sum()
    area = (gross - cont) / 0.9876
    var = gross_raw[mp].sum() + bg_scaled[mp].sum() + abs(cont)
    return area, np.sqrt(var) / 0.9876, sigma, e[mp][0], e[mp][-1], cont


def main():
    smp = read_rcxml.read(os.path.join(
        BASE, "RC103 черника маринелли авторская домик 246 гр.xml"))[0]
    bg = read_rcxml.read(os.path.join(BASE, "Фон домик 23 дня.xml"))[0]
    bgs = rebin_to(smp, bg) * (smp.live / bg.live)
    net = smp.counts - bgs

    Eg, epg = curve("m500", "full_organic_0.50")
    print("проба %.0f с, фон приведён; чистый счёт пробы %.3f имп/с\n"
          % (smp.live, net.sum() / smp.live))
    print("%-8s %-8s %9s %8s %10s %12s %10s"
          % ("линия", "нуклид", "окно, кэВ", "подложка", "площадь", "A, Бк",
             "a, Бк/кг"))
    for E0, nuc, y in LINES:
        r = peak_area(smp, net, smp.counts, bgs, E0)
        if r is None:
            print("%-8.1f %-8s окно вне спектра" % (E0, nuc))
            continue
        area, d, sigma, lo, hi, cont = r
        eps = eps_at(E0, Eg, epg)
        A = area / (eps * y * smp.live)
        dA = d / (eps * y * smp.live)
        print("%-8.1f %-8s %4.0f..%-4.0f %8.0f %6.0f±%-4.0f %6.0f±%-4.0f %8.0f"
              % (E0, nuc, lo, hi, cont, area, d, A, dA, A / smp.weight))

    print("\nдля справки: eps_p(662) = %.3e, eps_p(1461) = %.3e"
          % (eps_at(661.7, Eg, epg), eps_at(1460.8, Eg, epg)))
    print("ожидаемое содержание калия в сушёной чернике: 230..250 Бк/кг")


if __name__ == "__main__":
    main()
