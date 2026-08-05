# -*- coding: utf-8 -*-
"""Опорная сверка: модель Geant4 и собственная запись Cs-137, 10 см от торца.

ЧТО СРАВНИВАЕТСЯ. Наблюдаемая объявляется явно, потому что «эффективность» без
определения окна и без указания грани здесь не задана (прибор анизотропен, а
запись снята с торца). Считаются ДВЕ величины, обе по одинаковой конвенции
для модели и для измерения:

  A) счёт в широком окне 500–870 кэВ. Именно так набрана площадь в записи:
     весь пик целиком плюс подложка многократного рассеяния над комптоновским
     краем (477 кэВ). Величина конвенционная, но воспроизводимая;
  B) площадь пика подгонкой «гауссиана + линейная подложка» — той же, что в
     `RadiaCode-103/analysis/fit_peak.py`. Модельный спектр для этого сперва
     размывается приборным ПШПВ, снятым с самой записи, и подгоняется ТЕМ ЖЕ
     кодом в ТОМ ЖЕ окне: иначе сравниваются разные вещи — дельта-функция и
     гауссиана.

     Оконная сумма по ±1,25 ПШПВ здесь не годится: запись идёт по 8192
     каналам (0,38 кэВ на канал, ~60 отсчётов в канале у вершины), и поиск
     полувысоты по точкам спотыкается на пуассоновском шуме — первая же
     флуктуация вниз обрывает поиск и даёт ПШПВ 10 кэВ вместо ~40.

Обе величины делятся на A * p_gamma * t_live и сравниваются с модельной
eps(4pi) = counts / N_primaries * solid_angle_frac.

ЧЕГО ЭТА СВЕРКА НЕ ДАЁТ. Активность источника не метрологическая: паспортная
дата пересчитана распадом, погрешность паспорта неизвестна. Число годится как
анкор порядка величины, а не как поверочная точка.

    python analysis/compare_cs137.py <спектр.xml> [<модель.csv>]
"""
import os
import sys
import math

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(_HERE, "..", "..", "..",
                                                 "common", "py")))
import becqmoni as bm   # noqa: E402

# --- паспорт источника -------------------------------------------------------
# Реквизиты со слов оператора (паспорт Эксплораниум): 9,25 кБк на 02.01.2002.
# Период полураспада и выход линии — IAEA Live Chart of Nuclides (ENSDF),
# nds.iaea.org/relnsd/v1/data?fields=decay_rads&nuclides=137cs&rad_types=g:
# half_life = 30,08(9) года; 661,657 кэВ, интенсивность 85,1(2) %.
A0_BQ = 9250.0
T_HALF_A = 30.08
P_GAMMA = 0.851
E0_KEV = 661.657
DT_YEARS = 20.925          # 02.01.2002 -> 03.12.2022, дата записи из XML
WIN = (500.0, 870.0)       # широкое окно, как в записи

FIT = (560.0, 780.0)       # окно подгонки пика 662 кэВ

DEF_MODEL = "cs137_end10cm.csv"


def peak_model(x, area, mu, sigma, b0, b1):
    """Гауссиана площади `area` плюс линейная подложка, всё на кэВ."""
    g = area / (sigma * math.sqrt(2 * math.pi)) \
        * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    return g + b0 + b1 * (x - mu)


def fit_peak(x, y_per_keV, e0):
    """Подгонка одиночного пика. -> (площадь, центроида, ПШПВ)."""
    from scipy.optimize import curve_fit
    m = (x >= FIT[0]) & (x <= FIT[1])
    xx, yy = np.asarray(x)[m], np.asarray(y_per_keV)[m]
    if len(xx) < 8:
        raise SystemExit("В окне подгонки %s меньше 8 точек." % (FIT,))
    a0 = max(float(yy.sum() * (xx[-1] - xx[0]) / len(xx)), 1.0)
    p0 = [a0, e0, 0.06 * e0 / 2.3548, 0.0, 0.0]
    p, _ = curve_fit(peak_model, xx, yy, p0=p0, maxfev=20000)
    return abs(p[0]), p[1], abs(p[2]) * 2.3548


def read_model(path):
    """CSV прогона -> (массив по 1 кэВ, словарь {E: отсчёты}, N, доля, шапка).

    Словарь нужен `becqmoni.broaden`: она принимает именно {энергия: счёт},
    а не массив — так же, как её зовут скрипты Гамма-1С.
    """
    head, e, c = {}, [], []
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#"):
                if "=" in ln:
                    k, v = ln[1:].split("=", 1)
                    head[k.strip()] = v.strip()
            elif "," in ln and not ln.startswith("E_keV"):
                a, b = ln.split(",")
                e.append(float(a))
                c.append(float(b))
    n = int(head["N_primaries"])
    frac = float(head["solid_angle_frac"])
    emax = int(max(e)) + 2
    hist = np.zeros(emax)
    dic = {}
    for ee, cc in zip(e, c):
        hist[int(ee)] += cc
        dic[ee] = dic.get(ee, 0.0) + cc
    return hist, dic, n, frac, head


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("Нет пути к записи. Файл BecqMoni XML с образцом и встроенным "
              "фоном\nпередаётся первым аргументом; в репозиторий он не "
              "коммитится (личная запись).")
        return 2
    xml = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else DEF_MODEL
    if not os.path.exists(xml):
        print("Нет файла записи: %s" % xml)
        return 2
    if not os.path.exists(model):
        print("Нет файла модели: %s\nСначала прогон: asn16 "
              "point_cs137_end10cm_long.mac" % model)
        return 2

    # read_checked, а не read: `becqmoni` предписывает всем разбирающим
    # скриптам ходить через него, иначе правило «калибровка фона проверяется
    # независимо» держится на памяти автора каждого скрипта.
    sam, bg, bgdiag = bm.read_checked(xml)
    if bg is None:
        print("В записи нет встроенного фона — вычитать нечем.")
        return 2
    # Поканальный вычет законен только при совпадающих калибровках и длинах:
    # у спектров BecqMoni они бывают разными, и тогда фон вычтется не оттуда,
    # молча. В этой записи обе калибровки совпадают до последнего знака.
    if len(sam.n) != len(bg.n) or sam.cal != bg.cal:
        print("Калибровки или длины образца и фона не совпадают — "
              "поканальный вычет неприменим.")
        return 2
    # В XML этого прибора элемента LiveTime нет вовсе, и becqmoni подставляет
    # РЕАЛЬНОЕ время. Говорить «живое» было бы неверно: мёртвое время в
    # расчёт не входит.
    t_label = "живое" if sam.live != sam.real else "реальное (LiveTime в XML нет)"

    # --- активность на дату записи ------------------------------------------
    act = A0_BQ * math.exp(-math.log(2.0) / T_HALF_A * DT_YEARS)

    # --- измерение -----------------------------------------------------------
    k = sam.live / bg.live                      # масштаб фона по времени набора

    # Подгонка пика идёт ПЕРВОЙ: она даёт истинное положение линии, по
    # которому потом ставится широкое окно. Спектр приводится к отсчётам НА
    # КЭВ: шаг канала на этом приборе не постоянен (калибровка 4-й степени),
    # и без деления на ширину канала площадь под гауссианой не была бы числом
    # отсчётов.
    chx = np.arange(len(sam.n), dtype=float)
    en = sam.energy(chx)
    dE = np.gradient(en)
    net = sam.n - k * bg.n
    net_roi, cen, fwhm = fit_peak(en, net / dE, E0_KEV)
    # Погрешность площади пика — по отсчётам образца и фона в ±1,25 ПШПВ:
    # подгонка снимает произвол подложки, но дисперсию отсчётов не убирает.
    s_roi, _ = sam.counts_between(cen - 1.25 * fwhm, cen + 1.25 * fwhm)
    b_roi, _ = bg.counts_between(cen - 1.25 * fwhm, cen + 1.25 * fwhm)
    d_net_roi = math.sqrt(s_roi + k * k * b_roi)

    # Широкое окно приводится к ШКАЛЕ ЗАПИСИ. Калибровка прибора смещена:
    # подгонка даёт центроиду около 654,5 кэВ при табличных 661,657, то есть
    # −1,1 %. Если окно измерения ставить по номиналу, а модельное — по
    # истинной шкале, сравниваются два разных интервала; поправка даёт около
    # +1 % к нетто. Систематика ниже статистики, но она не усредняется.
    scale = cen / E0_KEV
    win_m = (WIN[0] * scale, WIN[1] * scale)
    s_win, nch = sam.counts_between(*win_m)
    b_win, _ = bg.counts_between(*win_m)
    net_win = s_win - k * b_win
    d_net_win = math.sqrt(s_win + k * k * b_win)

    denom = act * P_GAMMA * sam.live

    # --- модель --------------------------------------------------------------
    hist, dic, npri, frac, head = read_model(model)
    ch = np.arange(len(hist))
    m_win = hist[(ch >= WIN[0]) & (ch < WIN[1])].sum()
    # Размытие приборным ПШПВ и ТА ЖЕ подгонка в ТОМ ЖЕ окне.
    wide = bm.broaden(dic, fwhm_at_662=fwhm)
    gx = np.arange(len(wide), dtype=float)
    m_roi, m_cen, m_fwhm = fit_peak(gx, wide, E0_KEV)
    eps_m_win = m_win / npri * frac
    eps_m_roi = m_roi / npri * frac
    d_rel_m = 1.0 / math.sqrt(max(m_win, 1.0))

    # --- вывод ---------------------------------------------------------------
    print("ЗАПИСЬ %s" % os.path.basename(xml))
    print("  %s время образца %.0f с, фона %.0f с, масштаб фона %.5f"
          % (t_label, sam.live, bg.live, k))
    print("  окно A приведено к шкале записи: %.1f-%.1f кэВ (номинал %.0f-%.0f)"
          % (win_m[0], win_m[1], WIN[0], WIN[1]))
    print("  проверка калибровки фона: %s" % (bgdiag,))
    print("  подгонка пика: центроида %.2f кэВ, ПШПВ %.2f кэВ (%.1f %% на 662)"
          % (cen, fwhm, 100 * fwhm / cen))
    print("  активность на дату записи %.0f Бк (паспорт %.0f Бк, распад %.3f лет)"
          % (act, A0_BQ, DT_YEARS))
    print("МОДЕЛЬ %s" % os.path.basename(model))
    print("  N = %d, доля телесного угла %.6f, штамп %s"
          % (npri, frac, head.get("src_sha1", "?")))
    print()
    hdr = "%-34s %12s %12s %8s" % ("наблюдаемая", "измерение", "модель",
                                   "модель/изм")
    print(hdr)
    print("-" * len(hdr))
    for name, ne, dn, em, ee in (
            ("A) окно %.0f-%.0f кэВ" % WIN, net_win, d_net_win, m_win, eps_m_win),
            ("B) площадь пика подгонкой", net_roi, d_net_roi, m_roi,
             eps_m_roi)):
        eps_e = ne / denom
        print("%-34s %12.4e %12.4e %8.3f"
              % (name + ", eps(4pi)", eps_e, ee, ee / eps_e))
        print("%-34s %8.0f +- %-3.0f %12.0f" % ("   отсчётов", ne, dn, em))
    print()
    # Round-trip размытия: подгонка по размытой модели обязана вернуть ту
    # ПШПВ, которой размывали, и центроиду линии. Расхождение означало бы, что
    # размытие или подгонка сдвигают пик, и тогда сравнение площадей неверно.
    print("сверка размытия: подгонка по модели даёт центроиду %.2f кэВ "
          "(вход %.3f) и ПШПВ %.2f (вход %.2f)"
          % (m_cen, E0_KEV, m_fwhm, fwhm))
    print("статистика: измерение %.1f / %.1f %% (окно / пик), модель %.2f %%"
          % (100 * d_net_win / net_win, 100 * d_net_roi / net_roi,
             100 * d_rel_m))
    print("АКТИВНОСТЬ НЕ МЕТРОЛОГИЧЕСКАЯ: погрешность паспорта неизвестна, "
          "дата\nпересчитана распадом. Число — анкор порядка, не поверочная "
          "точка.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
