# -*- coding: utf-8 -*-
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

BASE = str(paths.measured("RadiaCode-103"))
FIELD = os.path.join(BASE, "Фон 7 дней без домика.xml")
BERRY = os.path.join(BASE, "RC103 черника маринелли авторская домик 246 гр.xml")
SHIELD = os.path.join(BASE, "Фон домик 23 дня.xml")

try:
    from gamma.calibration.fwhm_sanity_gate import check_fwhm_curve_consistency
except Exception as _exc:
    check_fwhm_curve_consistency = None
    _GATE_ERR = str(_exc)

CANDIDATES = [
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
    (2614.51, "Tl-208", "field"),
    (2614.51, "Tl-208 (домик 23 сут)", "shield"),
    (1460.82, "K-40 (домик 23 сут)", "shield")
]

_CACHE = {}

def load(source):
    if source in _CACHE:
        return _CACHE[source]
    if source == "field":
        fname = FIELD
    elif source == "berry":
        fname = BERRY
    elif source == "shield":
        fname = SHIELD
    else:
        raise ValueError("Unknown source: %s" % source)
    # read() отдаёт СПИСОК спектров, нужен элемент [0]; поля называются
    # .counts / .energy / .live (генератор выдумал energy_keV и live_seconds).
    data = read_rcxml.read(fname)[0]
    _CACHE[source] = (np.asarray(data.counts, dtype=float),
                      np.asarray(data.energy, dtype=float), float(data.live))
    return _CACHE[source]

def ecal_from(energy):
    poly = np.polyfit(np.arange(len(energy)), energy, 2)
    return [poly[2], poly[1], poly[0]]

# Линии естественных рядов, реально присутствующие в спектре открытого фона.
# Нужны измерителю для проверки изоляции — иначе слитый комплекс принимается
# за одиночную линию и завышает ширину.
KNOWN_LINES = [238.63, 241.98, 295.22, 300.09, 338.32, 351.93, 463.00, 510.77,
               583.19, 609.31, 665.45, 727.33, 768.36, 794.95, 806.17, 860.56,
               911.20, 934.06, 964.77, 968.97, 1120.29, 1155.19, 1238.11,
               1281.00, 1377.67, 1401.50, 1408.01, 1460.82, 1509.23, 1588.20,
               1620.50, 1661.28, 1729.60, 1764.49, 1847.42, 2118.55, 2204.21,
               2447.86, 2614.51]


def val(res, name, default=None):
    if hasattr(res, name):
        return getattr(res, name)
    if isinstance(res, dict) and name in res:
        return res[name]
    return default

def measure_all():
    rows = []
    for E0, label, source in CANDIDATES:
        try:
            counts, energy, live = load(source)
            ecal = ecal_from(energy)
            # known_lines_keV — список линий, которые ЕСТЬ в этом спектре; без него
            # измеритель не может знать про соседей и принимает слитый комплекс за
            # линию (Bi-214 1120 дал 135 кэВ при модели 70 — это 1120+1155 вместе).
            # Для черники список НЕ передаём: там цезий доминирует над фоновыми
            # линиями примерно в 30 раз, а с полным списком 661,657 честно объявится
            # блендированной с Bi-214 609/665 — и эталон, сошедшийся с двумя
            # приложениями, был бы отброшен без физической причины.
            # Полный список соседей оказался слишком жёстким: при ПШПВ ~80 кэВ
            # даже K-40 1460,8 формально не изолирован (1377,7 и 1509,2 внутри
            # 1,5 ПШПВ), и отбраковывалось ВСЁ. Поэтому измеряем без списка, а
            # выбросы отдаём донорскому гейту — он ищет их по остаткам, то есть
            # по данным, а не по таблице.
            res = measure_fwhm(counts, energy_keV=E0, energy_cal=ecal,
                               window_factor=1.25)
            fwhm = val(res, "fwhm_keV")
            unc = val(res, "fwhm_keV_unc") or val(res, "uncertainty_keV")
            sig = val(res, "significance_sigma") or val(res, "significance")
            passed = True
            reason = None
        except Exception as exc:
            fwhm = None
            unc = None
            sig = None
            passed = False
            reason = str(exc)
        rows.append({
            "E0": E0,
            "label": label,
            "source": source,
            "fwhm": fwhm,
            "unc": unc,
            "sig": sig,
            "passed": passed,
            "reason": reason
        })
    return rows

def print_table(rows):
    print("Метка и источник".ljust(30), "E0", "FWHM".rjust(10), "unc".rjust(10),
          "sig".rjust(10), "rcspec FWHM".rjust(12), "ratio".rjust(10))
    print("-" * 100)
    for row in rows:
        label = "%s (%s)" % (row["label"], row["source"])
        e0 = "%.3f" % row["E0"]
        fwhm = "%.3f" % row["fwhm"] if row["fwhm"] is not None else "отказ"
        unc = "%.3f" % row["unc"] if row["unc"] is not None else "-"
        sig = "%.2f" % row["sig"] if row["sig"] is not None else "-"
        rcspec_fwhm = "%.3f" % rcspec.fwhm(row["E0"])
        ratio = "%.3f" % (row["fwhm"] / rcspec.fwhm(row["E0"])) if row["fwhm"] is not None else "-"
        print(label.ljust(30), e0, fwhm.rjust(10), unc.rjust(10),
              sig.rjust(10), rcspec_fwhm.rjust(12), ratio.rjust(10))
    print()
    refused = [r for r in rows if not r["passed"]]
    if refused:
        print("Отказы:")
        for r in refused:
            print("  %s (%s): %s" % (r["label"], r["source"], r["reason"]))

def drop_outliers(rows, tol=0.25, min_keep=3):
    """Итеративно снимает точки, выпадающие из подгонки сильнее tol.

    Отбраковка ВСЕГДА печатается: молча урезанный набор якорей читается потом
    как «померили всё и сошлось», хотя половина точек выброшена.
    """
    cur = list(rows)
    while len(cur) > min_keep:
        E = [r["E0"] for r in cur]
        W = [r["fwhm"] for r in cur]
        f = fit_fwhm_scintillator(E, W)
        dev = [abs(w / f.fwhm_at(e) - 1.0) for e, w in zip(E, W)]
        i = int(np.argmax(dev))
        if dev[i] <= tol:
            break
        print("  ОТБРАКОВАНО: %s %.1f кэВ, измерено %.1f при подогнанных %.1f "
              "(отклонение %.0f %%)"
              % (cur[i]["label"], cur[i]["E0"], W[i], f.fwhm_at(E[i]), 100 * dev[i]))
        cur.pop(i)
    return cur


def fit_curve(rows):
    good_rows = [r for r in rows if r["passed"] and r["fwhm"] is not None]
    if len(good_rows) < 2:
        print("Невозможно построить кривую: меньше двух точек")
        return None
    good_rows = drop_outliers(good_rows)
    if len(good_rows) < 2:
        print("Невозможно построить кривую: после отбраковки осталось <2 точек")
        return None
    energies = [r["E0"] for r in good_rows]
    fwhms = [r["fwhm"] for r in good_rows]
    res = fit_fwhm_scintillator(energies, fwhms)
    k, alpha = res.coefficients
    print("k = %.6f" % k)
    print("alpha = %.6f" % alpha)
    print("n_points = %d" % len(energies))
    print("max_residual_keV = %.6f" % res.max_residual_keV)
    print("rms_residual_keV = %.6f" % res.rms_residual_keV)
    print("converged = %s" % str(res.converged))
    # FWHM^2 = k^2*E + k^2*alpha*E^2
    C = 0.0
    A2 = k**2
    B = k * np.sqrt(alpha) if alpha > 0 else 0.0
    if alpha < 0:
        print("ПРЕДУПРЕЖДЕНИЕ: отрицательное alpha означает, что кривая изгибается ВНИЗ,")
        print("что неправдоподобно для скинтиллятора. Использовать формулу нельзя.")
    print("C = %.6f" % C)
    print("A2 = %.6f" % A2)
    print("B = %.6f" % B)
    return (k, alpha, C, A2, B, energies)

def run_gate(rows, fit):
    if check_fwhm_curve_consistency is None:
        print("Гейт недоступен: %s" % _GATE_ERR)
        return
    # ГЕЙТУ ПОДАЮТСЯ ИЗМЕРЕННЫЕ ПАРЫ (E, ПШПВ), а не пересчитанные по кривой.
    # Сгенерированная версия подставляла k*sqrt(E+alpha*E^2), то есть проверяла
    # кривую саму против себя — такой гейт зелёный всегда и не значит ничего.
    pairs = [(r["E0"], r["fwhm"]) for r in rows
             if r.get("passed") and r.get("fwhm") is not None]
    if len(pairs) < 2:
        print("Гейт не запускался: якорей меньше двух")
        return
    print("Гейту переданы ИЗМЕРЕННЫЕ якоря: %s"
          % ", ".join("%.1f->%.1f" % p for p in pairs))
    try:
        result = check_fwhm_curve_consistency(pairs)
    except Exception as exc:
        import inspect
        print("Гейт не отработал: %s" % exc)
        print("сигнатура: %s" % inspect.signature(check_fwhm_curve_consistency))
        return
    print("Результат гейта: %s" % repr(result))
    if hasattr(result, 'passed'):
        print("passed = %s" % str(result.passed))
    if hasattr(result, 'reasons'):
        print("reasons:")
        for r in result.reasons:
            print("  %s" % r)
    if hasattr(result, 'warnings'):
        print("warnings:")
        for w in result.warnings:
            print("  %s" % w)
    if isinstance(result, dict):
        if 'passed' in result:
            print("passed = %s" % str(result['passed']))
        if 'reasons' in result:
            print("reasons:")
            for r in result['reasons']:
                print("  %s" % r)
        if 'warnings' in result:
            print("warnings:")
            for w in result['warnings']:
                print("  %s" % w)

def compare(fit):
    k, alpha, C, A2, B, anchors = fit
    energies = [32, 60, 100, 200, 400, 662, 1000, 1461, 2000, 2614]
    print("E".rjust(8), "новая FWHM".rjust(12), "процент от E".rjust(12),
          "rcspec FWHM".rjust(12), "процент от E".rjust(12), "отношение".rjust(10))
    print("-" * 80)
    for E in energies:
        new_fwhm = k * np.sqrt(E + alpha * E**2)
        new_percent = (new_fwhm / E) * 100
        old_fwhm = rcspec.fwhm(E)
        old_percent = (old_fwhm / E) * 100
        ratio = new_fwhm / old_fwhm if old_fwhm > 0 else 0.0
        marker = "<-- эталон" if E == 662 else ""
        print("%8d %12.3f %12.3f %12.3f %12.3f %10.3f %s" %
              (E, new_fwhm, new_percent, old_fwhm, old_percent, ratio, marker))

def main():
    print("Сборка кривой разрешения для RadiaCode-103")
    print("=" * 50)
    rows = measure_all()
    print_table(rows)
    fit = fit_curve(rows)
    if fit is None:
        return 1
    run_gate(rows, fit)
    compare(fit)
    print("\nКонстанты не записаны в файл, правка rcspec.py делается человеком осознанно,")
    print("потому что смена разрешения обнуляет все свёртки контура.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
