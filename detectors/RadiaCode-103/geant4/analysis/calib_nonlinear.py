# -*- coding: utf-8 -*-
"""
CAL-1 (03.09.2026). Проверка, нужна ли нелинейная энергетическая шкала измеренному спектру RadiaCode-103:
по реперным линиям фона измеряются центроиды пиков в КАНАЛАХ, затем сравниваются модели шкалы «энергия от канала»:
полиномы степени 2, 3, 4 и модель с насыщением SiPM. Штатное ПО прибора ограничено степенью 2;
в работе допустима любая до 4 (решение оператора 03.09.2026). Прогонов Geant4 нет.
"""
import os
import sys
import json
import argparse
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "analysis"))  # read_rcxml
DONOR = "D:/GoogleDrive/Дозиметрия/ИИ/1 Скилы/0_Work/gamma-spectrum-analysis/scripts"  # forward slash: без \-ловушек
sys.path.insert(0, DONOR)

def measure_lines(counts, device_coef):
    # Реперные линии (кэВ, имя)
    known_lines_keV = [
        (74.2,  "Pb K-alpha (72.80+74.97, смесь)"),
        (85.5,  "Pb K-beta (84.94+87.3, смесь)"),
        (186.0, "Ra-226 186.2 + U-235 185.7 (смесь)"),
        (238.63,"Pb-212"), (295.22,"Pb-214"), (351.93,"Pb-214"),
        (511.0, "аннигиляция"), (583.19,"Tl-208"), (609.31,"Bi-214"),
        (911.20,"Ac-228"), (1120.29,"Bi-214"), (1460.82,"K-40"),
        (1764.49,"Bi-214"), (2614.51,"Tl-208"),
    ]
    
    # Импорт внутри main, чтобы --selftest работал без донора
    import read_rcxml
    import gamma.calibration.fwhm_measure as FM

    lines = []
    for energy_keV, name in known_lines_keV:
        try:
            result = FM.measure_fwhm(
                counts,
                energy_keV=energy_keV,
                energy_cal=device_coef,
                window_factor=1.25,
                methods=(FM.METHOD_HALF_MAX, FM.METHOD_GAUSSIAN),
                known_lines_keV=[energy_keV],
                reject_blended=False
            )
        except Exception as e:
            print(f"Ошибка при измерении линии {name}: {e}")
            result = None

        if result and result.passed and result.centroid_channel is not None:
            lines.append({
                "name": name,
                "energy_keV": energy_keV,
                "channel": float(result.centroid_channel),
                "fwhm_keV": getattr(result, "fwhm_keV", None),
                "significance_sigma": getattr(result, "significance_sigma", None),
                "status": "измерена"
            })
        else:
            lines.append({
                "name": name,
                "energy_keV": energy_keV,
                "channel": None,
                "fwhm_keV": None,
                "significance_sigma": None,
                "status": "не измерена"
            })

    return lines

def fit_models(lines, device_coef):
    # Фильтрация прошедших линий
    passed = [l for l in lines if l["channel"] is not None]
    if not passed:
        raise ValueError("Нет прошедших линий для подгонки")

    channels = np.array([l["channel"] for l in passed])
    energies = np.array([l["energy_keV"] for l in passed])

    # Подгонка полиномов
    models = {}
    for d in [2, 3, 4]:
        try:
            from gamma.calibration.energy_fit import polynomial_energy_fit
            res = polynomial_energy_fit(channels, energies, min_degree=d, max_degree=d)
            coef = getattr(res, "coefficients", None)
            if coef is not None:
                models[f"poly{d}"] = {
                    "coefficients": coef,
                    "degree": d
                }
        except Exception as e:
            print(f"Ошибка подгонки полинома степени {d}: {e}")
            models[f"poly{d}"] = {"error": str(e)}

    # Подгонка модели насыщения SiPM
    def sat_model(ch, a0, a1, a2, ch0):
        g = -ch0 * np.log(1 - ch / ch0)
        return a0 + a1 * g + a2 * g**2

    def residuals(params, ch, E_ref):
        a0, a1, a2, ch0 = params
        E_fit = sat_model(ch, a0, a1, a2, ch0)
        return E_ref - E_fit

    # Начальные значения
    a0, a1, a2 = [c for c in device_coef]
    ch0_start = 5 * max(channels)
    x0 = [a0, a1, a2, ch0_start]

    bounds = ([float('-inf'), float('-inf'), float('-inf'), 1.05 * max(channels)],
              [float('inf'), float('inf'), float('inf'), 1e6])

    try:
        res = least_squares(residuals, x0, args=(channels, energies), bounds=bounds)
        models["sat"] = {
            "params": list(res.x),
            "success": True
        }
    except Exception as e:
        print(f"Ошибка подгонки модели насыщения: {e}")
        models["sat"] = {"error": str(e)}

    # Добавление "референсных" моделей
    models["device"] = {
        "coefficients": device_coef,
        "degree": 2
    }
    models["cal_room"] = {
        "coefficients": [-3.711311, 2.444318, 0.000321],
        "degree": 2
    }

    return models

def residual_table(model_name, coef_or_params, lines):
    # Фильтрация прошедших линий
    passed = [l for l in lines if l["channel"] is not None]
    channels = np.array([l["channel"] for l in passed])
    energies_ref = np.array([l["energy_keV"] for l in passed])

    if model_name == "device" or model_name == "cal_room":
        coef = coef_or_params
        energies_fit = np.polynomial.polynomial.polyval(channels, coef)
    elif model_name == "sat":
        a0, a1, a2, ch0 = coef_or_params
        g = -ch0 * np.log(1 - channels / ch0)
        energies_fit = a0 + a1 * g + a2 * g**2
    else:
        # Полином
        coef = coef_or_params
        energies_fit = np.polynomial.polynomial.polyval(channels, coef)

    residuals_keV = energies_ref - energies_fit

    # Вычисление σ для χ²
    sigmas = []
    for l in passed:
        if l["fwhm_keV"] is not None and l["significance_sigma"] is not None:
            sigma = l["fwhm_keV"] / 2.355 / np.sqrt(l["significance_sigma"])
        else:
            sigma = 1.0
        sigmas.append(sigma)

    chi2 = np.sum((residuals_keV / np.array(sigmas))**2) if all(s != 0 for s in sigmas) else float('inf')

    # RMS и максимумы
    rms = np.sqrt(np.mean(residuals_keV**2))
    max_abs_res_300 = max(abs(r) for r, e in zip(residuals_keV, energies_ref) if e < 300)
    max_abs_res_1500 = max(abs(r) for r, e in zip(residuals_keV, energies_ref) if e > 1500)

    return {
        "model": model_name,
        "rms": rms,
        "max_abs_res_300": max_abs_res_300,
        "max_abs_res_1500": max_abs_res_1500,
        "chi2": chi2,
        "residuals": list(residuals_keV),
        "fwhm_ratios": [r / f if f else None for r, f in zip(residuals_keV, [l["fwhm_keV"] for l in passed])]
    }

def selftest():
    print("Запуск самопроверки...")
    
    # 1. Позитив: сгенерировать 12 точек
    ch = np.linspace(30, 900, 12)
    E_true = 2.0 + 2.4 * ch + 3e-4 * ch**2 + 1e-7 * ch**3

    # Подгонка полиномов
    poly3_coef = np.polynomial.polynomial.polyfit(ch, E_true, 3)
    poly4_coef = np.polynomial.polynomial.polyfit(ch, E_true, 4)
    poly2_coef = np.polynomial.polynomial.polyfit(ch, E_true, 2)

    E3 = np.polynomial.polynomial.polyval(ch, poly3_coef)
    E4 = np.polynomial.polynomial.polyval(ch, poly4_coef)
    E2 = np.polynomial.polynomial.polyval(ch, poly2_coef)

    max_res3 = float(np.max(np.abs(E3 - E_true)))
    max_res4 = float(np.max(np.abs(E4 - E_true)))
    max_res2 = float(np.max(np.abs(E2 - E_true)))

    if max_res3 >= 1e-6 or max_res4 >= 1e-6:
        print("SELFTEST FAIL: полиномы 3 и 4 не достигли требуемой точности")
        return 1

    if max_res2 <= 0.5:
        print("SELFTEST FAIL: полином 2 не показал существенности кубического члена")
        return 1

    # 2. Негатив-мутация
    E_true_mut = E_true.copy()
    E_true_mut[3] += 15  # Испортить одну точку
    poly3_mut = np.polynomial.polynomial.polyfit(ch, E_true_mut, 3)  # ПЕРЕподгонка на испорченных данных
    E3_mut = np.polynomial.polynomial.polyval(ch, poly3_mut)
    max_res3_mut = float(np.max(np.abs(E3_mut - E_true_mut)))
    if max_res3_mut <= 1:
        print("SELFTEST FAIL: полином 3 не обнаружил порчу")
        return 1

    # 3. Модель насыщения
    ch0_true = 2000
    a0, a1, a2 = 0, 2.4, 0
    E_sat_true = a0 + a1 * (-ch0_true * np.log(1 - ch / ch0_true)) + a2 * (-ch0_true * np.log(1 - ch / ch0_true))**2

    def sat_model(ch, a0, a1, a2, ch0):
        g = -ch0 * np.log(1 - ch / ch0)
        return a0 + a1 * g + a2 * g**2

    def residuals(params, ch, E_ref):
        a0, a1, a2, ch0 = params
        E_fit = sat_model(ch, a0, a1, a2, ch0)
        return E_ref - E_fit

    x0 = [a0, a1, a2, 5 * max(ch)]
    bounds = ([float('-inf'), float('-inf'), float('-inf'), 1.05 * max(ch)],
              [float('inf'), float('inf'), float('inf'), 1e6])

    res = least_squares(residuals, x0, args=(ch, E_sat_true), bounds=bounds)
    ch0_fit = res.x[3]
    if abs(ch0_fit - ch0_true) / ch0_true >= 0.05:
        print("SELFTEST FAIL: параметр ch0 модели насыщения не совпадает с ожидаемым")
        return 1

    E_sat_fit = sat_model(ch, *res.x)
    max_res_sat = float(np.max(np.abs(E_sat_fit - E_sat_true)))
    if max_res_sat >= 0.05:
        print("SELFTEST FAIL: модель насыщения не достигла требуемой точности")
        return 1

    print("SELFTEST OK")
    return 0

def main(argv):
    sys.stdout.reconfigure(encoding="utf-8")

    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--selftest", action="store_true")

    args = parser.parse_args(argv)

    if args.selftest:
        return selftest()

    # Пути по умолчанию
    if args.spectrum is None:
        measured_dir = os.environ.get("G4MODELS_MEASURED", r"C:\g4work\measured\RadiaCode-103")
        args.spectrum = os.path.join(measured_dir, "Фон 7 дней без домика.xml")

    if args.out is None:
        out_dir = os.path.join(HERE, "out")
        os.makedirs(out_dir, exist_ok=True)
        args.out = os.path.join(out_dir, "calib_nonlinear.json")

    # Импорт донора
    try:
        import read_rcxml
        import gamma.calibration.fwhm_measure as FM
        from gamma.calibration.energy_fit import polynomial_energy_fit
    except ImportError as exc:
        print(f"Ошибка: не удалось импортировать необходимые модули: {exc!r}")
        return 1

    # Чтение спектра
    try:
        smp = read_rcxml.read(args.spectrum)[0]
        counts = smp.counts[:-1]  # Исключаем переполнение
        device_coef = smp.coef
        live_time = smp.live
    except Exception as e:
        print(f"Ошибка чтения спектра: {e}")
        return 1

    lines = measure_lines(counts, device_coef)
    models = fit_models(lines, device_coef)

    # Вывод в stdout
    print("Линии:")
    print("Имя\t\t\t\tE_ref\tКанал\tFWHM\tЗначимость\tСтатус")
    for l in lines:
        fwhm_str = f"{l['fwhm_keV']:.2f}" if l["fwhm_keV"] is not None else "-"
        sig_str = f"{l['significance_sigma']:.1f}" if l["significance_sigma"] is not None else "-"
        ch_str = f"{l['channel']:.1f}" if l["channel"] is not None else "-"
        print(f"{l['name']}\t{float(l['energy_keV']):.2f}\t{ch_str}\t{fwhm_str}\t{sig_str}\t{l['status']}")

    print("\nМодели:")
    print("Модель\t\tПараметры\tRMS\t\tMax|res|<300\tMax|res|>1500\tχ²")
    model_results = {}
    for name, model in models.items():
        if "error" in model:
            continue
        try:
            res_table = residual_table(name, model.get("coefficients", model.get("params", [])), lines)
            model_results[name] = res_table
            print(f"{name}\t\t{model.get('degree', len(model.get('params', [])))}\t\t"
                  f"{res_table['rms']:.4f}\t\t{res_table['max_abs_res_300']:.4f}\t\t"
                  f"{res_table['max_abs_res_1500']:.4f}\t\t{res_table['chi2']:.2f}")
        except Exception as e:
            print(f"Ошибка при вычислении таблицы для модели {name}: {e}")

    # Вывод коэффициентов
    print("\nКоэффициенты:")
    for name, model in models.items():
        if "error" in model:
            continue
        coef = model.get("coefficients", model.get("params", []))
        print(f"{name}: {coef}")

    # Сохранение в JSON
    output_data = {
        "spectrum": args.spectrum,
        "device_coef": device_coef,
        "lines": lines,
        "models": model_results
    }

    try:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(output_data, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"Ошибка записи в файл: {e}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
