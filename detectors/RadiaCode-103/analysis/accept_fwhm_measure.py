# -*- coding: utf-8 -*-
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import numpy as np
import read_rcxml

sys.path.insert(0, os.path.join(_HERE, "..", "..", "..", "common", "py"))
import paths

DONOR = r"D:\GoogleDrive\Дозиметрия\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis\scripts"
sys.path.insert(0, DONOR)
from gamma.calibration.fwhm_measure import measure_fwhm

BASE = str(paths.measured("RadiaCode-103"))
FIELD = os.path.join(BASE, "Фон 7 дней без домика.xml")
BERRY = os.path.join(BASE, "RC103 черника маринелли авторская домик 246 гр.xml")
SHIELD = os.path.join(BASE, "Фон домик 23 дня.xml")

def load(path):
    s = read_rcxml.read(path)[0]
    return (np.asarray(s.counts, dtype=float), np.asarray(s.energy, dtype=float), float(s.live))

def ecal_from(energy):
    coeffs = np.polyfit(np.arange(len(energy)), energy, 2)
    return list(coeffs[::-1])

def call(counts, energy, E0, **kw):
    try:
        res = measure_fwhm(counts, energy_keV=E0, energy_cal=ecal_from(energy), **kw)
        return res
    except Exception as exc:
        return {"passed": False, "fwhm_keV": None, "reason": u"ИСКЛЮЧЕНИЕ: " + str(exc)}

def val(res, name, default=None):
    if hasattr(res, name):
        return getattr(res, name)
    elif isinstance(res, dict) and name in res:
        return res[name]
    else:
        return default

def synth(fwhm_keV, E0, n_ch, gain, area, base_level, base_slope, seed):
    np.random.seed(seed)
    energy = gain * np.arange(n_ch)
    # Гауссова пикировка
    sigma = fwhm_keV / (2 * np.sqrt(2 * np.log(2)))
    counts = area * np.exp(-0.5 * ((energy - E0) / sigma)**2) / (sigma * np.sqrt(2 * np.pi))
    # Фон
    continuum = base_level + base_slope * (E0 - energy)
    continuum = np.clip(continuum, 0, None)
    counts += continuum
    # Пуассоновский шум
    counts = np.random.poisson(counts)
    return (np.asarray(counts, dtype=float), energy)

def check_reference():
    counts, energy, live = load(BERRY)
    res = call(counts, energy, 661.657)
    fwhm = val(res, "fwhm_keV")
    unc = val(res, "fwhm_keV_unc", default=val(res, "uncertainty_keV", default="-"))
    sig = val(res, "significance", default="-")
    reason = val(res, "reason", default="")
    print(u"Cs-137 661.657:")
    print(u"  FWHM = %.2f кэВ, неопределенность = %s, значимость = %s" % (fwhm, unc, sig))
    print(u"  Причина: %s" % reason)
    ref1 = 56.87
    ref2 = 58.068
    dev1 = abs(fwhm - ref1)
    dev2 = abs(fwhm - ref2)
    print(u"  Сравнение с эталонами: отклонение от %.2f = %.2f кэВ, от %.3f = %.3f кэВ" % (ref1, dev1, ref2, dev2))
    passed = val(res, "passed", False) and 55.5 <= fwhm <= 59.0
    print(u"[ПРОШЛА]" if passed else u"[НЕ ПРОШЛА]")
    return passed

def check_k40():
    counts, energy, live = load(FIELD)
    res = call(counts, energy, 1460.82)
    fwhm = val(res, "fwhm_keV")
    reason = val(res, "reason", default="")
    print(u"K-40 1460.82:")
    print(u"  FWHM = %.2f кэВ" % fwhm)
    print(u"  Причина: %s" % reason)
    passed = val(res, "passed", False) and 70.0 <= fwhm <= 80.0
    print(u"[ПРОШЛА]" if passed else u"[НЕ ПРОШЛА]")
    return passed

def check_honest_refusal():
    counts, energy, live = load(FIELD)
    res1 = call(counts, energy, 2614.51)
    res2 = call(counts, energy, 911.20)
    reason1 = val(res1, "reason", default="")
    reason2 = val(res2, "reason", default="")
    print(u"Tl-208 2614.51:")
    print(u"  Причина: %s" % reason1)
    print(u"Ac-228 911.20:")
    print(u"  Причина: %s" % reason2)
    passed = (not val(res1, "passed", False)) and (not val(res2, "passed", False))
    print(u"[ПРОШЛА]" if passed else u"[НЕ ПРОШЛА]")
    return passed

def check_absent_line():
    counts, energy, live = load(SHIELD)
    res = call(counts, energy, 661.657)
    reason = val(res, "reason", default="")
    print(u"Cs-137 661.657 в фоне:")
    print(u"  Причина: %s" % reason)
    passed = not val(res, "passed", False)
    print(u"[ПРОШЛА]" if passed else u"[НЕ ПРОШЛА]")
    return passed

def check_window_dependence():
    # Свип идёт по ЭТАЛОННОЙ линии (Cs-137 в чернике), а не по 2614 в фоне:
    # 2614 заведомо неизмерима (значимость 2,7 сигма), на ней свип показал бы
    # не зависимость от окна, а шум. Генератор подставил её ошибочно.
    counts, energy, live = load(BERRY)
    factors = [1.0, 1.25, 1.5, 2.0, 3.0]
    results = []
    print(u"Зависимость FWHM от window_factor (эталонная линия Cs-137 661,657):")
    print(u"window_factor\tFWHM")
    for factor in factors:
        res = call(counts, energy, 661.657, window_factor=factor)
        fwhm = val(res, "fwhm_keV")
        results.append(fwhm)
        print(u"%.2f\t\t%s" % (factor,
                               "отказ: " + str(val(res, "reason", ""))
                               if fwhm is None else "%.2f" % fwhm))
    ok3 = [w for w in results[:3] if w is not None]
    if len(ok3) < 3:
        print(u"[НЕ ПРОШЛА] - при window_factor <= 1.5 измеритель отказал")
        return False
    spread = (max(ok3) - min(ok3)) / np.median(ok3) * 100
    print(u"Разброс при window_factor <= 1.5: %.2f%%" % spread)
    passed = all(55.5 <= w <= 59.0 for w in ok3)
    if not passed:
        print(u"[НЕ ПРОШЛА]")
        return False
    # Проверка на window_factor >= 2.0
    for i, factor in enumerate(factors):
        if factor >= 2.0 and results[i] is not None and not (55.5 <= results[i] <= 59.0):
            print(u"ВНИМАНИЕ: значение FWHM при window_factor=%.2f выходит за пределы допустимого диапазона!" % factor)
    print(u"[ПРОШЛА]")
    return True

def check_mutation():
    # Синтетический спектр с известным FWHM
    synth_args = (56.87, 661.657, 1024, 2.76, 30000, 60, 0.05, 1)
    counts1, energy1 = synth(*synth_args)
    res1 = call(counts1, energy1, 661.657)
    w1 = val(res1, "fwhm_keV")
    
    # Синтетический спектр с увеличенным FWHM
    synth_args2 = (56.87 * 1.20, 661.657, 1024, 2.76, 30000, 60, 0.05, 1)
    counts2, energy2 = synth(*synth_args2)
    res2 = call(counts2, energy2, 661.657)
    w2 = val(res2, "fwhm_keV")
    
    # Синтетический спектр с уменьшенным FWHM
    synth_args3 = (56.87 * 0.75, 661.657, 1024, 2.76, 30000, 60, 0.05, 1)
    counts3, energy3 = synth(*synth_args3)
    res3 = call(counts3, energy3, 661.657)
    w3 = val(res3, "fwhm_keV")
    
    print(u"Проверка чувствительности к изменению FWHM:")
    print(u"  w1 = %.2f кэВ (ожидаемое = %.2f)" % (w1, 56.87))
    print(u"  w2 = %.2f кэВ (ожидаемое = %.2f)" % (w2, 56.87 * 1.20))
    print(u"  w3 = %.2f кэВ (ожидаемое = %.2f)" % (w3, 56.87 * 0.75))
    
    # Печатаем САМИ отношения, а не отклонения от них: сгенерированная версия
    # печатала abs(отношение-1) под подписью "отношение", и строка читалась как
    # "измеритель вернул 0,02 вместо 1,2" — то есть как провал вместо успеха.
    print(u"  Отношения: w1/56,87 = %.3f (ждём 1,000); w2/w1 = %.3f (ждём 1,200); "
          u"w3/w1 = %.3f (ждём 0,750)" % (w1 / 56.87, w2 / w1, w3 / w1))
    
    # Проверка условий
    cond1 = abs(w1 / 56.87 - 1) < 0.05
    cond2 = abs(w2 / w1 - 1.20) < 0.06
    cond3 = abs(w3 / w1 - 0.75) < 0.06
    
    passed = cond1 and cond2 and cond3
    if not passed:
        print(u"[НЕ ПРОШЛА]")
        return False
    else:
        print(u"[ПРОШЛА] (важно: w2 и w3 должны изменяться, иначе измеритель не чувствителен к FWHM)")
        return True

def check_overflow_channel():
    counts, energy, live = load(FIELD)
    last_ch = counts[-1]
    med_ch = np.median(counts[1000:1023])
    ratio = last_ch / med_ch if med_ch > 0 else float('inf')
    print(u"Проверка переполнения канала:")
    print(u"  Значение последнего канала: %.2f" % last_ch)
    print(u"  Медиана каналов 1000-1022: %.2f" % med_ch)
    print(u"  Отношение: %.2f" % ratio)
    
    res_full = call(counts, energy, 2614.51)
    res_no_last = call(counts[:-1], energy[:-1], 2614.51)
    
    passed = (val(res_full, "passed", False) == val(res_no_last, "passed", False))
    if not passed:
        print(u"[НЕ ПРОШЛА] - несоответствие между результатами с и без последнего канала")
        return False
    else:
        fwhm1 = val(res_full, "fwhm_keV")
        fwhm2 = val(res_no_last, "fwhm_keV")
        if fwhm1 is not None and fwhm2 is not None:
            diff = abs(fwhm1 - fwhm2) / fwhm1 * 100
            print(u"  Разница FWHM: %.2f%%" % diff)
        print(u"[ПРОШЛА]")
        return True

def check_subtracted_input():
    berry_counts, berry_energy, live_berry = load(BERRY)
    shield_counts, shield_energy, live_shield = load(SHIELD)
    
    # Интерполируем фон
    shield_interp = np.interp(np.arange(len(berry_counts)), np.arange(len(shield_counts)), shield_counts)
    bg_subtracted = berry_counts - shield_interp * (live_berry / live_shield)
    bg_subtracted = np.clip(bg_subtracted, 0, None)
    
    # Измеряем на сыром и фоне-вычтенном
    res_raw = call(berry_counts, berry_energy, 661.657)
    res_bg = call(bg_subtracted, berry_energy, 661.657)
    
    fwhm_raw = val(res_raw, "fwhm_keV")
    fwhm_bg = val(res_bg, "fwhm_keV")
    
    print(u"Сравнение с сырыми и фоном-вычтенными данными:")
    print(u"  Сырой спектр: %.2f кэВ" % fwhm_raw)
    print(u"  Фон-вычтенный: %.2f кэВ" % fwhm_bg)
    
    if fwhm_raw is not None and fwhm_bg is not None:
        diff = abs(fwhm_raw - fwhm_bg) / fwhm_raw * 100
        print(u"  Разница: %.2f%%" % diff)
    
    print(u"[СПРАВОЧНО]")
    return True

def main():
    print(u"Приемочное тестирование модуля FWHM")
    print(u"=====================================")
    
    checks = [
        ("check_reference", check_reference),
        ("check_k40", check_k40),
        ("check_honest_refusal", check_honest_refusal),
        ("check_absent_line", check_absent_line),
        ("check_window_dependence", check_window_dependence),
        ("check_mutation", check_mutation),
        ("check_overflow_channel", check_overflow_channel),
        ("check_subtracted_input", check_subtracted_input)
    ]
    
    results = []
    for name, func in checks:
        print(u"\n--- %s ---" % name)
        try:
            passed = func()
            results.append((name, passed))
        except Exception as e:
            print(u"[НЕ ПРОШЛА] - исключение: %s" % str(e))
            results.append((name, False))
    
    print(u"\n--- Результаты ---")
    print(u"Проверка\t\t\tРезультат")
    print(u"----------------------------------------")
    all_passed = True
    for name, passed in results:
        status = u"[ПРОШЛА]" if passed else u"[НЕ ПРОШЛА]"
        if not passed and name != "check_subtracted_input":
            all_passed = False
        print(u"%s\t%s" % (name.ljust(25), status))
    
    verdict = u"ПРИЕМКА ПРОЙДЕНА" if all_passed else u"ПРИЕМКА НЕ ПРОЙДЕНА"
    print(u"\n--- ИТОГ ---")
    print(verdict)
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())
