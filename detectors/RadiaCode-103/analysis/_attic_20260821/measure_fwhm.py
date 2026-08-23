# -*- coding: utf-8 -*-
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import numpy as np
import rcspec
import read_rcxml

sys.path.insert(0, os.path.join(_HERE, "..", "..", "..", "common", "py"))
import paths

BASE = str(paths.measured("RadiaCode-103"))
FIELD = os.path.join(BASE, "Фон 7 дней без домика.xml")
BERRY = os.path.join(BASE, "RC103 черника маринелли авторская домик 246 гр.xml")
SHIELD = os.path.join(BASE, "Фон домик 23 дня.xml")

def spec_arrays(path, subtract=None):
    s = read_rcxml.read(path)[0]
    y = s.counts.astype(float) / s.live
    if subtract is not None:
        s2 = read_rcxml.read(subtract)[0]
        y2 = s2.counts.astype(float) / s2.live
        y2 = np.interp(s.energy, s2.energy, y2)
        y -= y2
    return (s.energy, y)

def net_peak(e, y, E0, half=3.0, gap=1.0):
    s = rcspec.fwhm(E0) / 2.35482
    inner = abs(e - E0) <= half * s
    left_band = (e < E0 - half*s) & (e >= E0 - (half + gap + 1.0)*s)
    right_band = (e > E0 + half*s) & (e <= E0 + (half + gap + 1.0)*s)
    band = left_band | right_band
    if np.sum(band) < 4:
        return (None, None)
    p = np.polyfit(e[band], y[band], 1)
    line = p[0] * e[inner] + p[1]
    return (e[inner], y[inner] - line)

def smooth(net, npt):
    """Скользящее среднее. Слабые линии открытого фона зашумлены: без
    сглаживания argmax садится на одиночный статистический выброс, и ширина
    выходит в один канал (проверено — 583/911/1764/2614 давали 2-7 кэВ)."""
    if npt < 3:
        return net
    k = np.ones(npt) / npt
    return np.convolve(net, k, mode="same")


def half_max_width(e, net, E0=None, win=1.0):
    """ПШПВ по половине максимума. Вершина ищется не по всему окну, а рядом с
    табличной энергией (+-win сигма), иначе её уводит соседняя линия/выброс."""
    if E0 is not None:
        s = rcspec.fwhm(E0) / 2.35482
        near = np.abs(e - E0) <= win * s
        i = int(np.arange(len(net))[near][np.argmax(net[near])])
    else:
        i = int(np.argmax(net))
    top = net[i]
    if top <= 0 or i == 0 or i == len(net) - 1:
        return (None, None)
    half = top / 2.0

    def cross(step):
        """Уходим от вершины, пока отсчёты выше половины; точка пересечения —
        линейной интерполяцией между последним «выше» и первым «ниже»."""
        j = i
        while 0 <= j + step < len(net) and net[j + step] > half:
            j += step
        k = j + step
        if not (0 <= k < len(net)):
            return None            # склон обрывается на краю окна
        if net[j] == net[k]:
            return e[j]
        return e[j] + (e[k] - e[j]) * (net[j] - half) / (net[j] - net[k])

    left, right = cross(-1), cross(+1)
    if left is None or right is None:
        return (None, e[i])
    return (right - left, e[i])

def main():
    cases = [
        ("661.7 Cs-137 (черника минус фон домика)", BERRY, SHIELD, 661.657),
        ("583.2 Tl-208 (открытый фон)", FIELD, None, 583.19),
        ("911.2 Ac-228 (открытый фон)", FIELD, None, 911.20),
        ("1460.8 K-40 (открытый фон)", FIELD, None, 1460.82),
        ("1764.5 Bi-214 (открытый фон)", FIELD, None, 1764.49),
        ("2614.5 Tl-208 (открытый фон)", FIELD, None, 2614.51)
    ]
    
    print("%-42s %9s %9s %8s %7s"
          % ("линия", "изм,кэВ", "модель", "изм/мод", "изм,%"))
    print("-" * 80)
    
    rows = []
    for name, path, subtract, E0 in cases:
        e, y = spec_arrays(path, subtract)
        enet, ynet = net_peak(e, y, E0)
        if enet is None:
            print("%-42s %8s" % (name, "окно вне шкалы"))
            continue
        # окно сглаживания ~0,4 сигма: мельче шума, много крупнее канала
        dE = float(np.median(np.diff(enet))) if len(enet) > 2 else 1.0
        npt = int(round(0.4 * rcspec.fwhm(E0) / 2.35482 / max(dE, 1e-6)))
        width, _ = half_max_width(enet, smooth(ynet, npt), E0)
        if width is None:
            print("%-42s %8s" % (name, "не разр."))
            continue
        model_fwhm = rcspec.fwhm(E0)
        ratio = width / model_fwhm
        percent = 100.0 * width / E0
        print("%-42s %9.1f %9.1f %8.2f %7.2f" % (name, width, model_fwhm, ratio, percent))
        rows.append((E0, width))
    
    if len(rows) >= 2:
        E = np.array([r[0] for r in rows])
        W = np.array([r[1] for r in rows])
        A2 = np.sum(W**2 * E) / np.sum(E**2)
        print("\nМНК по измеренным ширинам в форме FWHM^2 = A2*E:")
        print("  A2 = %.3f кэВ   (в rcspec сейчас %.3f)" % (A2, rcspec.FWHM_A2))
        fwhm_662_model = rcspec.fwhm(661.657)
        print("  FWHM(662) = %.1f кэВ = %.2f %%  (в rcspec сейчас %.1f кэВ = %.2f %%)"
              % (np.sqrt(A2 * 661.657), 100.0 * np.sqrt(A2 / 661.657),
                 fwhm_662_model, 100.0 * fwhm_662_model / 661.657))
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
