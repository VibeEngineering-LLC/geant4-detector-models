"""Зонная аппроксимация кривой эффективности — модуль БЕЗ побочных эффектов.

Вынесено из compare_point.py по замечанию аудитора: тот модуль при импорте
проверяет каталог расчётных спектров и падает в чистом клоне, а сшивка нужна
скриптам, которым спектры не нужны вовсе (compare_effcalcmc работает по двум
файлам из репозитория). Правило: модуль, из которого импортируют функции,
не имеет права требовать окружение, нужное лишь части его пользователей.

Метод — по образцу ЛСРМ (Efficiency/EffCalcMC): кривая строится не одним
полиномом по всему диапазону, а зонами С ПЕРЕКРЫТИЕМ; в перекрытии соседние
ветви смешиваются линейно по log E («сшивка»).

Число зон и степени подобраны перебором, критерий — скользящее исключение
узла (не χ² подгонки: его занижают лишние степени). Победили ДВЕ зоны с
широким перекрытием 338–662 кэВ, степени 5/2: RMS предсказания 0,72 % (5 см)
и 1,7 % (25 см) против 0,85 % и 5,2 % у трёхзонной ЛСРМ-разметки. Смысл:
ниже ~660 кэВ кривая горбатая (нужна 5-я степень), выше — почти прямая в
log-log (хватает 2-й). Их Efficiency на расчётной кривой EffCalcMC сам
выбрал два интервала (степени 4/4) — сходится.
"""
import math

import numpy as np

ZONES = [(None, 661.7, 5), (338.3, None, 2)]


def zoned_fit(Eg, yg, dyg, zones=None):
    """[(lo, hi, deg, coeffs, chi2/nu, n)] по зонам + интерполятор log-log."""
    zones = ZONES if zones is None else zones
    Eg = np.asarray(Eg, float)
    yg = np.asarray(yg, float)
    dyg = np.asarray(dyg, float)
    lE, ly = np.log(Eg), np.log(yg)
    w = yg / np.maximum(dyg, 1e-30)
    fits = []
    for lo, hi, deg in zones:
        lo = Eg[0] if lo is None else lo
        hi = Eg[-1] if hi is None else hi
        m = (Eg >= lo * 0.999) & (Eg <= hi * 1.001)
        if m.sum() < deg + 2:      # зоне нужен запас узлов над степенью
            deg = max(1, int(m.sum()) - 2)
        cf = np.polyfit(lE[m], ly[m], deg, w=w[m])
        rr = (ly[m] - np.polyval(cf, lE[m])) * w[m]
        chi2 = (rr ** 2).sum() / max(1, m.sum() - deg - 1)
        fits.append((lo, hi, deg, cf, chi2, int(m.sum())))

    def ev(E):
        x = math.log(E)
        hit = [(lo, hi, cf) for lo, hi, _d, cf, _c, _n in fits
               if lo * 0.999 <= E <= hi * 1.001]
        if not hit:
            lo, hi, _d, cf, _c, _n = fits[0] if E < fits[0][1] else fits[-1]
            return math.exp(np.polyval(cf, x))
        if len(hit) == 1:
            return math.exp(np.polyval(hit[0][2], x))
        (l1, h1, c1), (l2, h2, c2) = hit[0], hit[1]
        a, b = math.log(max(l1, l2)), math.log(min(h1, h2))
        t = 0.5 if b <= a else min(1.0, max(0.0, (x - a) / (b - a)))
        return math.exp((1 - t) * np.polyval(c1, x) + t * np.polyval(c2, x))

    return fits, ev


def local_quad(Ec, yc):
    """Интерполятор ГУСТОЙ кривой: локальная квадратика в log-log.

    Для сравнения двух кривых интерполировать надо ГУСТУЮ в узлы редкой, а
    не наоборот (замечание аудитора): ошибка локальной квадратики на сетке
    из 50 лог-точек ничтожна, тогда как глобальная аппроксимация редкой
    кривой сглаживает ровно тот край, где вся структура.
    """
    lE = np.log(np.asarray(Ec, float))
    ly = np.log(np.asarray(yc, float))

    def ev(E):
        x = math.log(E)
        i = int(np.clip(np.searchsorted(lE, x), 1, len(lE) - 2))
        sl = slice(i - 1, i + 2)
        cf = np.polyfit(lE[sl], ly[sl], 2)
        return math.exp(np.polyval(cf, x))

    return ev
