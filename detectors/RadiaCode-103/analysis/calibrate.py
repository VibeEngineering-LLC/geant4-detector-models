# -*- coding: utf-8 -*-
"""Перекалибровка RadiaCode-103 по фоновому спектру: E(канал) и FWHM(E).

Скрипт ИЗМЕРЯЕТ и ПРЕДЛАГАЕТ, но ничего не переписывает: смена разрешения
обнуляет все свёртки контура, поэтому правка rcspec.py делается человеком.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import numpy as np
import read_rcxml
import rcspec

sys.path.insert(0, os.path.join(_HERE, "..", "..", "..", "common", "py"))
import paths

DONOR = r"D:\GoogleDrive\Дозиметрия\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis\scripts"
sys.path.insert(0, DONOR)
from gamma.calibration.fwhm_measure import measure_fwhm
from gamma.calibration.fwhm_fit import fit_fwhm_scintillator
from gamma.calibration.energy_fit import polynomial_energy_fit

BASE = str(paths.measured("RadiaCode-103"))
FIELD = os.path.join(BASE, "Фон 7 дней без домика.xml")

# Действующая калибровка контура — стартовое приближение и база сравнения.
CAL0 = [-3.711311, 2.444318, 0.000321]

# Реперные линии фона: нуклид, энергия кэВ. Набор классический для фонового
# спектра — сильные линии рядов урана и тория плюс калий.
LINES = [
    ("Pb-212", 238.63), ("Pb-214", 351.93), ("Tl-208", 583.19),
    ("Bi-214", 609.31), ("Ac-228", 911.20), ("Bi-214", 1120.29),
    ("K-40", 1460.82), ("Bi-214", 1764.49), ("Tl-208", 2614.51),
]

MIN_SIGMA = 5.0        # ниже — линия не считается измеренной
MAX_ITER = 4           # итераций «измерить -> пересчитать шкалу -> измерить»

def measure_lines(counts, cal):
    accepted = []
    rejected = []
    for nuclide, energy_keV in LINES:
        try:
            res = measure_fwhm(counts, energy_keV=energy_keV, energy_cal=cal, window_factor=1.25)
        except Exception:
            res = type('obj', (object,), {'passed': False, 'reason': 'исключение в measure_fwhm'})()
        if res.passed and res.significance_sigma is not None and res.significance_sigma >= MIN_SIGMA and res.centroid_channel is not None:
            accepted.append((nuclide, energy_keV, res))
        else:
            rejected.append((nuclide, energy_keV, res.reason))
    return accepted, rejected

def main():
    sys.stdout.reconfigure(encoding="utf-8")
    smp = read_rcxml.read(FIELD)[0]
    counts = smp.counts[:-1].astype(float)

    cal = list(CAL0)
    fit = None
    prev_fit = None
    for i in range(MAX_ITER):
        accepted, rejected = measure_lines(counts, cal)
        if len(accepted) < 3:
            print(f"Итерация {i+1}: принято {len(accepted)} линий (<3), остановка")
            break

        channels = [res.centroid_channel for _, _, res in accepted]
        energies = [e for _, e, _ in accepted]

        fit = polynomial_energy_fit(channels, energies, max_degree=2, target_residual_keV=1.0, min_degree=1)
        cal = list(fit.coefficients)

        print(f"Итерация {i+1}: принято {len(accepted)} линий, полином степени {fit.degree}, max={fit.max_residual_keV:.2f} кэВ, rms={fit.rms_residual_keV:.2f} кэВ")
        # Сходимость сравнивает ПРЕДЫДУЩУЮ подгонку с текущей. Сравнение fit с
        # самим собой давало бы тождественный ноль и объявляло сходимость на
        # второй итерации всегда.
        if prev_fit is not None:
            old_energy = np.asarray(prev_fit.predict([0, len(counts)-1]))
            new_energy = np.asarray(fit.predict([0, len(counts)-1]))
            diff = float(np.max(np.abs(old_energy - new_energy)))
            if diff < 0.05:
                print(f"Сходимость достигнута: изменение шкалы на краях {diff:.3f} кэВ")
                prev_fit = fit
                break
        prev_fit = fit

    if fit is None:
        print("Не удалось сойтись: нет данных для подгонки энергетической шкалы")
        return

    # Финальная таблица линий
    accepted, rejected = measure_lines(counts, cal)
    print("\nТаблица измеренных линий:")
    print("Нуклид     Энергия   Центроид  Отклонение  Канал  FWHM      Сигма   Действующая  Отношение")
    print("           кэВ       кэВ       кэВ          -      кэВ ± погр. -       кэВ        -")
    for nuclide, energy_keV, res in accepted:
        centroid_keV = res.centroid_keV
        deviation = centroid_keV - energy_keV
        fwhm = res.fwhm_keV
        uncertainty = res.fwhm_uncertainty_keV
        sigma = res.significance_sigma
        rcspec_fwhm = rcspec.fwhm(energy_keV, "103")
        ratio = fwhm / rcspec_fwhm if rcspec_fwhm > 0 else np.nan
        print(f"{nuclide:8} {energy_keV:6.2f} {centroid_keV:6.2f} {deviation:6.2f}     {res.centroid_channel:6.1f} {fwhm:5.2f} ± {uncertainty:5.2f} {sigma:5.1f}   {rcspec_fwhm:6.2f}    {ratio:5.3f}")
    for nuclide, energy_keV, reason in rejected:
        print(f"{nuclide:8} {energy_keV:6.2f} отказ     {reason}")

    # Подгонка кривой разрешения
    if len(accepted) < 2:
        print("\nНедостаточно данных для подгонки FWHM")
        return

    energies = [e for _, e, _ in accepted]
    fwhms = [res.fwhm_keV for _, _, res in accepted]

    fit_fwhm = fit_fwhm_scintillator(energies, fwhms)
    k, alpha = fit_fwhm.coefficients
    print(f"\nПодгонка FWHM сцинтиллятора:")
    print(f"k = {k:.4f}, alpha = {alpha:.6f}, точек = {fit_fwhm.n_points}")
    if alpha < 0:
        print("ПРЕДУПРЕЖДЕНИЕ: alpha < 0, кривая изгибается вниз — неправдоподобно для сцинтиллятора")

    # Подгонка FWHM^2 = C + A2*E
    E_array = np.array(energies)
    fwhm_array = np.array(fwhms)
    y = fwhm_array**2
    A = np.column_stack([np.ones_like(E_array), E_array])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        C, A2 = coeffs
        print(f"\nПодгонка FWHM^2 = C + A2*E:")
        print(f"C = {C:.4f} кэВ², A2 = {A2:.6f} кэВ")
        noise_level = np.sqrt(max(0.0, C))
        print(f"Шумовой член (sqrt(C)) = {noise_level:.3f} кэВ")

        # Остатки
        model_fwhm_sq = C + A2 * E_array
        residuals_sq = y - model_fwhm_sq
        print("\nОстатки по точкам:")
        for i, (e, res_sq) in enumerate(zip(energies, residuals_sq)):
            print(f"  {e:5.0f} кэВ: {res_sq:.3f} кэВ²")
    except Exception as e:
        print(f"\nОшибка подгонки FWHM^2 = C + A2*E: {e}")

    # Сравнительная таблица
    test_energies = [32, 60, 100, 200, 400, 662, 1000, 1461, 2000, 2614]
    print("\nСравнение FWHM:")
    print("Энергия   Действующая   Новая (C+A2*E)   Отношение")
    print("кэВ       кэВ           кэВ              -")
    for E in test_energies:
        old_fwhm = rcspec.fwhm(E, "103")
        new_fwhm = np.sqrt(max(0.0, C + A2 * E))
        ratio = new_fwhm / old_fwhm if old_fwhm > 0 else np.nan
        print(f"{E:5d}     {old_fwhm:6.3f}       {new_fwhm:6.3f}         {ratio:5.3f}")

    # Требует толкования
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    issues = []
    if len(accepted) < 5:
        issues.append("принято меньше 5 линий")
    if alpha < 0:
        issues.append("alpha < 0 в модели сцинтиллятора")
    if C < 0:
        issues.append("C < 0 в форме C + A2*E (отрицательный шумовой член)")
    if fit.max_residual_keV > 3.0:
        issues.append("максимальный остаток энергетической подгонки > 3 кэВ")
    for nuclide, energy_keV, res in accepted:
        if abs(res.centroid_keV - energy_keV) > 5.0:
            issues.append(f"отклонение центроиды > 5 кэВ для {nuclide}")
    for nuclide, energy_keV, res in accepted:
        rcspec_fwhm = rcspec.fwhm(energy_keV, "103")
        ratio = res.fwhm_keV / rcspec_fwhm if rcspec_fwhm > 0 else np.nan
        if not (0.8 <= ratio <= 1.25):
            issues.append(f"отношение FWHM < 0.8 или > 1.25 для {nuclide}")
    if issues:
        for issue in issues:
            print(f"- {issue}")
    else:
        print("вопросов нет")

    # Предлагаемые константы
    print("\nПРЕДЛАГАЕМЫЕ КОНСТАНТЫ ДЛЯ rcspec.py:")
    print(f"FWHM_C = {C:.6f}")
    print(f"FWHM_A2 = {A2:.6f}")
    print(f"FWHM_B = 0.0")
    print("k = {:.6f}, alpha = {:.6f} (относятся к альтернативной форме)".format(k, alpha))
    print("\nПРЕДУПРЕЖДЕНИЕ: скрипт ничего не записал и правка делается человеком осознанно,")
    print("поскольку смена разрешения обнуляет все свёртки контура.")

if __name__ == "__main__":
    main()
