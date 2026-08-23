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

# ДОНОР: SpectraVibe, только чтение (§12 - чужая зона)
DONOR = r"D:\GoogleDrive\Дозиметрия\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis\scripts"
sys.path.insert(0, DONOR)
from gamma.calibration.fwhm_provider import _measure_fwhm_channels
from gamma.calibration.fwhm_fit import fit_fwhm_scintillator

BASE = str(paths.measured("RadiaCode-103"))
FIELD = os.path.join(BASE, "Фон 7 дней без домика.xml")
BERRY = os.path.join(BASE, "RC103 черника маринелли авторская домик 246 гр.xml")
SHIELD = os.path.join(BASE, "Фон домик 23 дня.xml")

# Пороги донора (fwhm_provider.py:62-66), взяты как есть
MIN_HEIGHT = 30.0     # отсчётов в пике над подложкой
MIN_SIGMA = 8.0       # значимость
ISO_FAC = 1.5         # изоляция: сосед не ближе 1.5*ПШПВ
MIN_ANCHORS = 3

LINES = [
    (661.657, "Cs-137", "berry"),
    (238.63, "Pb-212", "field"),
    (295.22, "Pb-214", "field"),
    (351.93, "Pb-214", "field"),
    (583.19, "Tl-208", "field"),
    (609.31, "Bi-214", "field"),
    (911.20, "Ac-228", "field"),
    (968.97, "Ac-228", "field"),
    (1120.29, "Bi-214", "field"),
    (1460.82, "K-40", "field"),
    (1764.49, "Bi-214", "field"),
    (2614.51, "Tl-208", "field")
]

NEIGHBOURS = [338.32, 463.0, 510.77, 727.33, 794.95, 860.56, 934.06, 1155.19, 1238.11,
              1281.0, 1377.67, 1408.01, 1588.2, 1620.5, 1729.6, 1847.4, 2118.55, 2204.21,
              2447.86]

def counts_of(key):
    if key == "field":
        data = read_rcxml.read(FIELD)[0]
        return data.energy, data.counts.astype(float), data.live
    else:  # berry
        berry = read_rcxml.read(BERRY)[0]
        shield = read_rcxml.read(SHIELD)[0]
        # Порядок обязателен: фон СНАЧАЛА переносится на сетку пробы, ПОТОМ
        # масштабируется по живому времени, и только затем вычитается.
        # Генератор выдал (berry-shield)*k с интерполяцией уже разности —
        # это и не вычитание фона, и не та сетка.
        bg = np.interp(berry.energy, shield.energy, shield.counts.astype(float))
        net = berry.counts.astype(float) - bg * (berry.live / shield.live)
        return berry.energy, net, berry.live

def channel_of(e, E0):
    return int(np.argmin(np.abs(e - E0)))


def flatten(counts, ch, seed_ch):
    """Снять локальный НАКЛОН континуума прямой по боковым полосам.

    Донор берёт подложку как среднее внешних пятых долей окна
    (fwhm_provider.py:86-89) — верно для пика на РОВНОЙ подложке. На фоновом
    спектре континуум круто падает, левое крыло выше правого, их среднее
    завышает подложку: уровень половины поднимается, ширина занижается, а на
    слабых линиях net_height уходит в минус и метод отказывает (K-40 на сыром
    спектре дал 60,0 кэВ против 73 кэВ со снятым наклоном). Поэтому наклон
    снимается ДО вызова донора: это подготовка входа под тот случай, для
    которого донорский метод и написан, а не правка чужого метода.
    """
    n = len(counts)
    half = max(3, int(round(3.0 * seed_ch)))
    lo, hi = max(0, ch - half), min(n, ch + half + 1)
    x = np.arange(lo, hi)
    y = np.asarray(counts, dtype=float)[lo:hi]
    side = np.abs(x - ch) > 1.2 * seed_ch
    if side.sum() < 4:
        return np.asarray(counts, dtype=float)
    k, b = np.polyfit(x[side], y[side], 1)
    out = np.zeros(n, dtype=float)
    out[lo:hi] = y - (k * x + b)
    return out

def dE_at(e, ch):
    de = (e[min(ch+1, len(e)-1)] - e[max(ch-1, 0)]) / 2.0
    return de if de > 0 else 1.0

def isolation_ok(E0, fwhm_keV, key):
    """Изоляция проверяется по линиям, которые ЕСТЬ В ЭТОМ спектре.

    В чернике фон домика уже вычтен, естественных соседей (Bi-212 727 и др.)
    там нет — общий список соседей забраковал бы единственную сильную линию
    прибора без всякой физической причины.
    """
    if key == "berry":
        return True
    for E in NEIGHBOURS + [x[0] for x in LINES if x[2] == "field"]:
        if abs(E - E0) > 0.1 and abs(E - E0) < ISO_FAC * fwhm_keV:
            return False
    return True

def significance(counts, ch, fwhm_ch):
    w = max(1, int(round(fwhm_ch / 2.0)))
    lo = max(0, ch - w)
    hi = min(len(counts), ch + w + 1)
    baseline = (counts[lo] + counts[hi-1]) / 2.0
    gross = counts[lo:hi].sum()
    base = baseline * (hi - lo)
    net = gross - base
    if gross <= 0:
        return 0.0
    else:
        return net / np.sqrt(gross)

def measure_all():
    sources = {}
    rows = []
    for E0_keV, label, source_key in LINES:
        if source_key not in sources:
            sources[source_key] = counts_of(source_key)
        e, counts, live = sources[source_key]
        ch = channel_of(e, E0_keV)
        dE = dE_at(e, ch)
        seed_ch = 1.5 * rcspec.fwhm(E0_keV) / dE
        counts = flatten(counts, ch, seed_ch)
        w_ch = _measure_fwhm_channels(counts, ch, seed_ch)
        fwhm_keV = w_ch * dE if w_ch is not None else None
        # Вычисляем высоту пика (без учета подложки)
        w = max(1, int(round(w_ch / 2.0)) if w_ch is not None else int(round(seed_ch / 2.0)))
        lo = max(0, ch - w)
        hi = min(len(counts), ch + w + 1)
        baseline = (counts[lo] + counts[hi-1]) / 2.0
        height = counts[ch] - baseline
        sig = significance(counts, ch, w_ch if w_ch is not None else seed_ch)
        iso = isolation_ok(E0_keV,
                           fwhm_keV if fwhm_keV is not None else rcspec.fwhm(E0_keV),
                           source_key)
        ok = True
        reason = ""
        if w_ch is None:
            ok = False
            reason = "ширина не читается"
        elif height < MIN_HEIGHT:
            ok = False
            reason = "мало отсчётов"
        elif sig < MIN_SIGMA:
            ok = False
            reason = "низкая значимость"
        elif not iso:
            ok = False
            reason = "линия не изолирована"
        rows.append({
            "E0": E0_keV,
            "label": label,
            "fwhm_keV": fwhm_keV,
            "height": height,
            "sig": sig,
            "iso": iso,
            "ok": ok,
            "reason": reason
        })
    return rows

def main():
    rows = measure_all()
    print("%-26s %10s %10s %10s %10s %10s %s" % ("Энергия, кэВ", "ПШПВ изм", "Модель", "Отношение",
                                                   "Значимость", "Результат", "Причина"))
    for r in rows:
        fwhm = "%10.2f" % r["fwhm_keV"] if r["fwhm_keV"] is not None else "-".rjust(10)
        model = "%10.2f" % rcspec.fwhm(r["E0"])
        ratio = ("%10.2f" % (r["fwhm_keV"] / rcspec.fwhm(r["E0"]))
                 if r["fwhm_keV"] is not None else "-".rjust(10))
        sig = "%10.2f" % r["sig"]
        verdict = "OK" if r["ok"] else "FAIL"
        reason = r["reason"] if not r["ok"] else ""
        print("%-26s %s %s %s %s %10s %s" % (r["label"] + " (%.2f)" % r["E0"], fwhm, model,
                                             ratio, sig, verdict, reason))
    good = [r for r in rows if r["ok"]]
    print("\nКоличество опорных точек: %d из %d" % (len(good), MIN_ANCHORS))
    if len(good) < MIN_ANCHORS:
        print("ПРЕДУПРЕЖДЕНИЕ: Не выполнено условие донора (%d точек)" % MIN_ANCHORS)
        if len(good) < 2:
            return 1
    res = fit_fwhm_scintillator([r["E0"] for r in good], [r["fwhm_keV"] for r in good])
    print("\nМодель: FWHM = k*sqrt(E + alpha*E^2)")
    # У FwhmFitResult нет полей .k/.alpha — коэффициенты лежат кортежем
    # в .coefficients (донор, fwhm_fit.py:48). Генератор их выдумал.
    k_fit, alpha_fit = res.coefficients
    print("k = %.4f" % k_fit)
    print("alpha = %.6f" % alpha_fit)
    print("Точек: %d" % res.n_points)
    print("Макс. ошибка: %.4f кэВ" % res.max_residual_keV)
    print("RMS ошибка: %.4f кэВ" % res.rms_residual_keV)
    print("Сходится: %s" % str(res.converged))
    print("\nСравнение на 3 энергиях:")
    for E in [662, 1461, 2614]:
        new = res.fwhm_at(E)
        old = rcspec.fwhm(E)
        print("E = %d кэВ: новая = %.2f кэВ (%.1f%%), старая = %.2f кэВ" % (E, new, new/E*100, old))
    return 0

if __name__ == "__main__":
    sys.exit(main())
