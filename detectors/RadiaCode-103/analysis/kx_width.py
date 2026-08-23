# -*- coding: utf-8 -*-
"""
Спектр Cs-137 содержит пик около 32 кэВ, который не является одной линией,
а представляет собой комплекс K-X бария, возникающий при внутреннем преобразовании
в Ba-137m. Инструментальное программное обеспечение подогнало одну гауссову кривую
к целому комплексу и сообщило FWHM = 9.999 кэВ при центроиде 31.900 кэВ.

Это значение FWHM завышено, поскольку оно включает как реальное разрешение детектора,
так и распределение компонентов K-X (Kalpha около 32 кэВ, Kbeta около 36.5 кэВ).
Используя 9.999 кэВ напрямую как точку калибровки разрешения, мы получим завышенное
разрешение детектора на низких энергиях.

Этот скрипт восстанавливает ИСТИННОЕ разрешение детектора при ~32 кэВ методом прямого моделирования:
строим комплекс из табличных линий, распределяем его с пробной шириной детектора,
подгоняем одну гауссову кривую так же, как делает ПО инструмента, и находим пробную ширину,
которая воспроизводит наблюдаемое значение 9.999 кэВ.

Источник данных: IAEA NDS, ENSDF (E. Browne and J. K. Tuli,
cut-off 1-Oct-2006), retrieved 2026-08-21 via
`https://nds.iaea.org/relnsd/v1/data?fields=decay_rads&nuclides=137cs&rad_types=x`.
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import numpy as np
from scipy.optimize import curve_fit, brentq
import rcspec

# Комплекс K-X бария (в кэВ, проценты интенсивности)
BA_KX = [
    (31.816, 1.9905),   # Kalpha2
    (32.193, 3.6671),   # Kalpha1
    (36.482, 1.0786),   # K'beta1
    (37.255, 0.2718),   # K'beta2
]

# ВНИМАНИЕ: строка ENSDF с меткой `KB` (36.827, 1.3504) является суммой K'beta1 и K'beta2,
# и НЕ должна добавляться дополнительно — это приведет к двойному учету Kbeta.

# Референсные данные
REF_FWHM = 9.999
REF_CENTROID = 31.900
REF_CS_E = 661.556
REF_CS_FWHM = 58.068
REF_CS_FWHM_ALT = 56.87
REF_CS_E_ALT = 658.03

def complex_profile(x, fwhm_true):
    """Возвращает модельную профиль комплекса K-X на сетке x."""
    sigma = fwhm_true / 2.35482
    profile = np.zeros_like(x)
    for energy, intensity in BA_KX:
        profile += intensity * np.exp(-0.5 * ((x - energy) / sigma)**2)
    return profile

def fit_single_gaussian(x, y):
    """Подгоняет одну гауссову кривую и возвращает центроид и FWHM."""
    def gaussian(x, a, mu, s, c):
        return a * np.exp(-0.5 * ((x - mu) / s)**2) + c

    ymax = y.max()
    mu_guess = x[np.argmax(y)]
    s_guess = 4.0
    a_guess = ymax
    c_guess = 0.0

    try:
        popt, _ = curve_fit(gaussian, x, y, p0=[a_guess, mu_guess, s_guess, c_guess], maxfev=10000)
        amplitude, mu, s, c = popt
        return mu, abs(s) * 2.35482
    except:
        return mu_guess, s_guess * 2.35482

def observed_fwhm(fwhm_true):
    """Вычисляет наблюдаемую FWHM для заданной истинной ширины."""
    x = np.linspace(20.0, 50.0, 3001)
    profile = complex_profile(x, fwhm_true)
    _, fwhm_observed = fit_single_gaussian(x, profile)
    return fwhm_observed

def observed_centroid(fwhm_true):
    """Вычисляет наблюдаемый центроид для заданной истинной ширины."""
    x = np.linspace(20.0, 50.0, 3001)
    profile = complex_profile(x, fwhm_true)
    mu, _ = fit_single_gaussian(x, profile)
    return mu

def solve_true_fwhm():
    """Находит истинную ширину детектора методом бисекции."""
    try:
        result = brentq(lambda f: observed_fwhm(f) - REF_FWHM, 0.5, 15.0)
        return result
    except ValueError:
        print(u"Ошибка: не удалось найти корень в заданном интервале")
        return None

def curve_two_point(E1, W1, E2, W2):
    """Решает уравнение FWHM = k*sqrt(E + alpha*E^2) через две точки."""
    u1 = W1**2
    u2 = W2**2
    alpha = (u1 * E2 - u2 * E1) / (u2 * E1**2 - u1 * E2**2)
    k = W1 / np.sqrt(E1 + alpha * E1**2)
    return k, alpha

def main():
    print(u"Комплекс K-X бария:")
    total_weighted_energy = sum(e * i for e, i in BA_KX) / sum(i for _, i in BA_KX)
    alpha_weighted_energy = sum(e * i for e, i in BA_KX[:2]) / sum(i for _, i in BA_KX[:2])
    print(u"  Взвешенная энергия комплекса: %.3f кэВ" % total_weighted_energy)
    print(u"  Взвешенная энергия Kalpha: %.3f кэВ" % alpha_weighted_energy)
    for e, i in BA_KX:
        print(u"  %.3f кэВ (%.4f%%)" % (e, i))

    true_fwhm = solve_true_fwhm()
    if true_fwhm is None:
        return 1

    print(u"")
    print(u"Наблюдаемая ПШПВ комплекса: %.3f кэВ" % REF_FWHM)
    print(u"Восстановленная истинная ПШПВ прибора: %.3f кэВ" % true_fwhm)
    diff_kev = REF_FWHM - true_fwhm
    diff_percent = (diff_kev / true_fwhm) * 100
    print(u"Разница: %.3f кэВ (%.2f%%)" % (diff_kev, diff_percent))

    observed_cent = observed_centroid(true_fwhm)
    print(u"")
    print(u"Модельный центроид при восстановленной ширине: %.3f кэВ" % observed_cent)
    if abs(observed_cent - REF_CENTROID) < 0.01:
        print(u"Центры совпадают — модель корректна")
    else:
        print(u"Центры не совпадают — возможная ошибка модели")

    print(u"")
    print(u"Кривая A (по двум точкам):")
    kA, alphaA = curve_two_point(alpha_weighted_energy, true_fwhm, REF_CS_E, REF_CS_FWHM)
    print(u"  k = %.6f, alpha = %.8f" % (kA, alphaA))

    print(u"Кривая B:")
    kB, alphaB = curve_two_point(alpha_weighted_energy, true_fwhm, REF_CS_E_ALT, REF_CS_FWHM_ALT)
    print(u"  k = %.6f, alpha = %.8f" % (kB, alphaB))

    print(u"")
    print(u"Сравнение разрешений:")
    print(u"Энергия кэВ   Кривая A FWHM   A %%   Кривая B FWHM   Текущее FWHM   Отношение A/текущее")
    energies = [32, 60, 100, 200, 400, 662, 1000, 1461, 2000, 2614]
    for E in energies:
        A_FWHM = kA * np.sqrt(E + alphaA * E**2)
        B_FWHM = kB * np.sqrt(E + alphaB * E**2)
        current_FWHM = rcspec.fwhm(E)
        ratio = A_FWHM / current_FWHM if current_FWHM > 0 else 0
        print(u"  %6d     %8.2f      %5.1f    %8.2f      %8.2f       %6.3f" %
              (E, A_FWHM, (A_FWHM / current_FWHM - 1) * 100,
               B_FWHM, current_FWHM, ratio))

    print(u"")
    print(u"Наибольшие отклонения от кривых:")
    max_diff_A = 0
    max_diff_B = 0
    for E in energies:
        A_FWHM = kA * np.sqrt(E + alphaA * E**2)
        B_FWHM = kB * np.sqrt(E + alphaB * E**2)
        current_FWHM = rcspec.fwhm(E)
        diff_A = abs(A_FWHM - current_FWHM) / current_FWHM
        diff_B = abs(B_FWHM - current_FWHM) / current_FWHM
        if diff_A > max_diff_A:
            max_diff_A = diff_A
        if diff_B > max_diff_B:
            max_diff_B = diff_B

    print(u"Максимальное отклонение от кривой A: %.2f%%" % (max_diff_A * 100))
    print(u"Максимальное отклонение от кривой B: %.2f%%" % (max_diff_B * 100))

    return 0

if __name__ == "__main__":
    sys.exit(main())
