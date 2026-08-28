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
LINES = [("K-40", 1460.82, (150, 200, 250), (70.0, 74.2, 80.0)),
         ("комплекс ~600", 600.2, (90, 120, 150), (45.0, 52.5, 58.0))]

def main():
    sys.stdout.reconfigure(encoding="utf-8")
    c = read_rcxml.read(C.FIELD)[0].counts[:-1].astype(float)
    ch = np.arange(len(c))
    E = C.CAL0[0] + C.CAL0[1] * ch + C.CAL0[2] * ch * ch
    print("%-14s %7s %8s %8s %8s %8s" % ("линия", "полуокно", "ПШПВ", "T", "sigma", "сошлось"))
    for name, e0, halfs, fwhms in LINES:
        for hw in halfs:
            s = (E >= e0 - hw) & (E <= e0 + hw)
            x, y = E[s], c[s]
            base = np.linspace(y[:4].mean(), y[-4:].mean(), len(x))
            for fw in fwhms:
                # sigma_fixed — в единицах абсциссы (у нас кэВ); fwhm_channels
                # ждёт КАНАЛЫ и на кэВ-сетке молча не действует: sigma
                # выходила одинаковой при всех заданных ПШПВ.
                r = fit_peak_image(x, y - base, mu0=e0,
                                   sigma_fixed=fw / 2.35482,
                                   fit_sigma=False, fit_T=True, fit_step=True)
                print("%-14s %7d %8.1f %8.3f %8.2f %8s"
                      % (name, hw, fw, r.T, r.sigma, r.converged))

if __name__ == "__main__":
    main()
