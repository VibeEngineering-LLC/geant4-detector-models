# -*- coding: utf-8 -*-
"""Параметр хвоста T функции отклика RC-103 по измеренным линиям фона.

Модель формы — LSRM peak-image донора SpectraVibe (gauss + левый
экспоненциальный хвост + комптоновская ступенька). Ступенька при подгонке
ВКЛЮЧЕНА: она есть в измерении, и без неё её счёт впитал бы хвост.
Скан по окну обязателен (урок W-036): число без размаха по параметру метода
результатом не считается.
"""
import os, sys
_H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _H)
import numpy as np, read_rcxml, calibrate as C
sys.path.insert(0, C.DONOR)
from gamma.peaks.peak_image import fit_peak_image

# ПШПВ задаётся, а не подгоняется (донор, F-449): при свободной sigma подгонка
# вырождается — T гулял от 1,1 до 8,1, sigma уходила в минус при converged=True.
# T известен независимо: скан по согласию ВСЕГО разложения дал минимум 0,75
# (см. web-background, профиль по параметру хвоста). Поэтому здесь T
# фиксируется, а свободной остаётся ширина — задача обратная вчерашней и
# обусловлена лучше: вчера свободными были и T, и sigma, и подгонка вырождалась.
T_FIX = 0.75
# Окна — в долях ожидаемой ПШПВ, а не произвольные: слишком широкое окно
# захватывает соседние линии и континуум, и подгонка уходит в абсурд
# (при полуокне 260 кэВ на 2614 она дала отрицательную ширину при
# converged=True). Ближайшие соседи: у 1460,8 это 1764,5, у 2614,5 — край
# диапазона 2830.
LINES = [("K-40 1460,8", 1460.82, 75.0, (1.2, 1.4, 1.6, 1.8, 2.0)),
         ("Tl-208 2614,5", 2614.51, 110.0, (1.2, 1.4, 1.6, 1.8, 2.0))]

def main():
    sys.stdout.reconfigure(encoding="utf-8")
    c = read_rcxml.read(C.FIELD)[0].counts[:-1].astype(float)
    ch = np.arange(len(c))
    E = C.CAL0[0] + C.CAL0[1] * ch + C.CAL0[2] * ch * ch
    print("%-14s %7s %8s %8s %8s %8s" % ("линия", "полуокно", "ПШПВ с хв.", "ПШПВ гаусс", "центр", "сошлось"))
    for name, e0, w0, ks in LINES:
        for k in ks:
            hw = k * w0
            sel = (E >= e0 - hw) & (E <= e0 + hw)
            x, y = E[sel], c[sel]
            nb = max(4, len(x)//10)
            base = np.linspace(y[:nb].mean(), y[-nb:].mean(), len(x))
            g = fit_peak_image(x, y - base, mu0=e0, T0=T_FIX, fit_T=False,
                               fit_step=True)
            h = fit_peak_image(x, y - base, mu0=e0, T0=0.0, fit_T=False,
                               fit_step=True)
            print("%-14s %7d %9.2f %9.2f %8.2f %8s"
                  % (name, hw, 2.35482 * g.sigma, 2.35482 * h.sigma,
                     g.mu, g.converged))

if __name__ == "__main__":
    main()
