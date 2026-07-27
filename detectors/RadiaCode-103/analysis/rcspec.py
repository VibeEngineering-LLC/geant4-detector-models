# -*- coding: utf-8 -*-
"""Общие операции над спектрами RadiaCode: чтение, свёртка с разрешением,
площади в областях интереса.

Монте-Карло даёт спектр ЭНЕРГОВЫДЕЛЕНИЯ — без разрешения прибора пик полного
поглощения стоит в одном канале. Чтобы сравнивать с реальным набором и считать
площади ROI, спектр сворачивается с гауссианой прибора.
"""
import os
import sys

import numpy as np

RESULTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       "..", "results"))
# Консоль Windows работает в CP-1251/CP-866, где нет ни греческой «бета», ни
# знака «≈». Печать β-вклада падала с UnicodeEncodeError на полпути, теряя
# уже посчитанный результат. Пусть лучше подставит вопросительный знак.
try:
    sys.stdout.reconfigure(errors="replace")
except (AttributeError, ValueError):   # не консоль или старый python
    pass

VESSELS = ("m200", "m500")

NBINS = 3201          # каналы по 1 кэВ, как пишет rc_curves


def vessel():
    """Сосуд берётся из аргументов командной строки, по умолчанию m200."""
    for a in sys.argv[1:]:
        if a in VESSELS:
            return a
    return "m200"


def rdir(*parts, v=None):
    """Путь внутри results для выбранного сосуда."""
    return os.path.join(RESULTS, v or vessel(), *parts)
# Паспортное разрешение на 662 кэВ: RC-102 9.5 %, RC-103 8.4 %.
# Ход по энергии принят статистическим, FWHM ~ sqrt(E); у CsI(Tl) на низких
# энергиях реальное разрешение хуже из-за непропорциональности отклика.
R662 = {"102": 0.095, "103": 0.084, "101": 0.095}


def read_spec(path):
    """-> (метаданные, массив каналов 1 кэВ длиной NBINS)."""
    meta = {}
    hist = np.zeros(NBINS)
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                if "=" in line:
                    k, v = line[1:].split("=", 1)
                    meta[k.strip()] = v.strip()
            elif line and line[0].isdigit():
                a, b = line.split(",")
                i = int(float(a))
                if 0 <= i < NBINS:
                    hist[i] += float(b)
    return meta, hist


def fwhm(E, model="103"):
    """Полная ширина на полувысоте, кэВ."""
    r = R662[model] if isinstance(model, str) else float(model)
    return r * 662.0 * np.sqrt(np.maximum(E, 1.0) / 662.0)


def fold(hist, model="103"):
    """Свёртка с разрешением прибора. Ширина зависит от энергии, поэтому
    раскладываем каждый канал в свою гауссиану."""
    E = np.arange(len(hist)) + 0.5
    out = np.zeros_like(hist)
    sig = fwhm(E, model) / 2.35482
    nz = np.nonzero(hist)[0]
    for i in nz:
        s = sig[i]
        lo = max(0, int(i - 5 * s))
        hi = min(len(hist), int(i + 5 * s) + 1)
        x = E[lo:hi]
        g = np.exp(-0.5 * ((x - E[i]) / s) ** 2)
        gs = g.sum()
        if gs > 0:
            out[lo:hi] += hist[i] * g / gs
    return out


def roi(hist, E0, model="103", nsig=1.2):
    """Площадь в окне +-nsig*sigma вокруг линии (около 80 % площади пика при
    nsig=1.2) и границы окна."""
    s = fwhm(E0, model) / 2.35482
    lo, hi = int(E0 - nsig * s), int(E0 + nsig * s) + 1
    lo, hi = max(0, lo), min(len(hist), hi)
    return hist[lo:hi].sum(), lo, hi


def peak_fraction(model="103", nsig=1.2):
    """Доля площади гауссианы, попадающая в окно +-nsig*sigma."""
    from math import erf, sqrt
    return erf(nsig / sqrt(2.0))
