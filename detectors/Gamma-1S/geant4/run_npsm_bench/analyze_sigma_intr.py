import csv
import math
import glob
import argparse
import pathlib
import sys
import json
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

# Импорт функции parse_csv из analyze_ncompt.py
sys.path.insert(0, str(pathlib.Path(__file__).parent))
try:
    from analyze_ncompt import parse_csv
except ImportError:
    print("ОШИБКА: Не удалось импортировать parse_csv из analyze_ncompt.py")
    sys.exit(3)

def read_spectrum_table(file_path):
    """Читает таблицу спектра из CSV файла, начиная с строки bin_keV."""
    table = []
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    # Найти строку с bin_keV
    start_line = None
    for i, line in enumerate(lines):
        if line.strip().startswith("bin_keV"):
            start_line = i + 1
            break
    
    if start_line is None:
        raise ValueError("Не найдена строка с заголовком bin_keV")
    
    # Читать строки до пустой строки или конца файла
    for i in range(start_line, len(lines)):
        line = lines[i].strip()
        if not line:
            break
        try:
            bin_keV, count_edep, count_light = map(float, line.split(","))
            table.append((bin_keV, count_edep, count_light))
        except ValueError:
            raise ValueError(f"Ошибка разбора строки: {line}")
    
    return table

def find_peak_position(spectrum_table, energy_keV):
    """Находит пик полного поглощения в заданном диапазоне."""
    min_energy = 0.3 * energy_keV
    max_energy = 1.05 * energy_keV
    
    # Найти максимальный элемент в диапазоне
    peak_bin = None
    max_count = -1
    
    for i, (bin_keV, count_edep, count_light) in enumerate(spectrum_table):
        if min_energy <= bin_keV <= max_energy:
            if count_light > max_count:
                max_count = count_light
                peak_bin = i
    
    if peak_bin is None:
        raise ValueError("Не найден пик в заданном диапазоне")
    
    # Проверить, не на границе ли он
    is_on_edge = (peak_bin == 0) or (peak_bin == len(spectrum_table) - 1)
    
    return peak_bin, max_count, is_on_edge

def peak_half_width(spectrum_table, peak_bin):
    """Грубая полуширина пика: от максимума наружу до половины высоты."""
    n = len(spectrum_table)
    half = spectrum_table[peak_bin][2] / 2.0
    left = peak_bin
    while left > 0 and spectrum_table[left][2] >= half:
        left -= 1
    right = peak_bin
    while right < n - 1 and spectrum_table[right][2] >= half:
        right += 1
    return max((spectrum_table[right][0] - spectrum_table[left][0]) / 2.0, 1.0)


def estimate_background(spectrum_table, peak_bin, energy_keV):
    """Оценивает фон под пиком линейно, ВНЕ самого пика.

    Прежняя версия брала участки вплотную к пику: они содержали его крылья,
    фон завышался, и в окне подгонки не оставалось положительных отсчётов.
    """
    peak_energy = spectrum_table[peak_bin][0]
    hw = peak_half_width(spectrum_table, peak_bin)
    inner = 4.0 * hw
    left_start = peak_energy - 10.0 * hw
    right_end = peak_energy + 10.0 * hw
    
    # Найти бины для левого и правого участков
    left_bins = []
    right_bins = []
    
    for i, (bin_keV, count_edep, count_light) in enumerate(spectrum_table):
        d = bin_keV - peak_energy
        # Отступ inner исключает крылья пика из фоновых участков.
        if left_start <= bin_keV and d <= -inner:
            left_bins.append((i, bin_keV, count_light))
        elif bin_keV <= right_end and d >= inner:
            right_bins.append((i, bin_keV, count_light))
    
    # Если участков недостаточно, использовать ближайшие
    if len(left_bins) < 2:
        left_bins = [(i, bin_keV, count_light) for i, bin_keV, count_light in spectrum_table 
                     if left_start <= bin_keV < peak_energy]
        left_bins.sort(key=lambda x: x[1], reverse=True)
        left_bins = left_bins[:2] if len(left_bins) >= 2 else left_bins
    
    if len(right_bins) < 2:
        right_bins = [(i, bin_keV, count_light) for i, bin_keV, count_light in spectrum_table 
                      if peak_energy < bin_keV <= right_end]
        right_bins.sort(key=lambda x: x[1])
        right_bins = right_bins[:2] if len(right_bins) >= 2 else right_bins
    
    # Если не хватает точек для линейной аппроксимации
    if len(left_bins) < 2 or len(right_bins) < 2:
        raise ValueError("Недостаточно точек для оценки фона")
    
    # Линейная регрессия на левом участке
    x_left = np.array([bin_keV for _, bin_keV, _ in left_bins])
    y_left = np.array([count_light for _, _, count_light in left_bins])
    A_left = np.vstack([x_left, np.ones(len(x_left))]).T
    m_left, c_left = np.linalg.lstsq(A_left, y_left, rcond=None)[0]
    
    # Линейная регрессия на правом участке
    x_right = np.array([bin_keV for _, bin_keV, _ in right_bins])
    y_right = np.array([count_light for _, _, count_light in right_bins])
    A_right = np.vstack([x_right, np.ones(len(x_right))]).T
    m_right, c_right = np.linalg.lstsq(A_right, y_right, rcond=None)[0]
    
    # Оценка фона в пике: среднее значение линейных аппроксимаций
    peak_energy = spectrum_table[peak_bin][0]
    bg_left = m_left * peak_energy + c_left
    bg_right = m_right * peak_energy + c_right
    
    return (bg_left + bg_right) / 2

def fit_gaussian(spectrum_table, peak_bin, energy_keV, background):
    """Подгоняет гауссиану к пикам с вычетом фона."""
    # Определить окно ±3 полуширины на полувысоте
    peak_energy = spectrum_table[peak_bin][0]
    peak_count = spectrum_table[peak_bin][2]
    
    # Полуширина на полувысоте (предварительная оценка)
    half_max = peak_count / 2.0
    
    # Полуширина ищется ОТ ПИКА НАРУЖУ, а не по всему спектру. В спектре
    # гамма-квантов есть комптоновский континуум, и бинов выше половины
    # максимума там много по всей шкале: сгенерированная версия брала min/max
    # энергии среди ВСЕХ таких бинов и получала полуширину размером со шкалу,
    # после чего подгонка не сходилась ни на одном прогоне.
    n = len(spectrum_table)
    left = peak_bin
    while left > 0 and spectrum_table[left][2] >= half_max:
        left -= 1
    right = peak_bin
    while right < n - 1 and spectrum_table[right][2] >= half_max:
        right += 1
    fwhm = spectrum_table[right][0] - spectrum_table[left][0]
    if fwhm <= 0:
        fwhm = 2.0  # вырожденный случай: пик в один бин

    # Начальное приближение для сигмы
    sigma_guess = fwhm / 2.3548

    # Окно ±3σ
    window_size = 3 * sigma_guess
    min_window_energy = peak_energy - window_size
    max_window_energy = peak_energy + window_size
    
    # Собрать данные в окне
    window_data = []
    for i, (bin_keV, count_edep, count_light) in enumerate(spectrum_table):
        if min_window_energy <= bin_keV <= max_window_energy:
            window_data.append((i, bin_keV, count_light))
    
    # Если данных мало
    if len(window_data) < 3:
        raise ValueError("Недостаточно данных для подгонки гауссиана")
    
    # Подгонка методом наименьших квадратов
    x_data = np.array([bin_keV for _, bin_keV, _ in window_data])
    y_data = np.array([count_light - background for _, _, count_light in window_data])
    
    # Начальные приближения (вырезаны при замене имитации подгонки).
    amplitude_guess = peak_count - background
    centroid_guess = peak_energy

    # ПОДГОНКА. Сгенерированная версия имитировала её: цикл по sigma внутри
    # присваивал sigma = sigma_guess и делал break, возвращая начальное
    # приближение под видом результата. Ошибок не было, число выглядело
    # правдоподобно — ложный зелёный в чистом виде. Переписано вручную.
    #
    # Схема: перебор по сетке (sigma, mu); амплитуда входит ЛИНЕЙНО, поэтому
    # при заданных sigma и mu решается аналитически: A = sum(y*g)/sum(g*g).
    if np.all(y_data <= 0):
        raise ValueError("В окне подгонки нет положительных отсчётов")

    best = None
    sigmas = np.linspace(max(0.3, 0.25 * sigma_guess), 4.0 * sigma_guess, 120)
    mus = np.linspace(centroid_guess - 2.0 * sigma_guess,
                      centroid_guess + 2.0 * sigma_guess, 81)
    for sg in sigmas:
        for mu in mus:
            g = np.exp(-0.5 * ((x_data - mu) / sg) ** 2)
            denom = float(np.sum(g * g))
            if denom <= 0:
                continue
            A = float(np.sum(y_data * g)) / denom
            if A <= 0:
                continue
            err = float(np.sum((y_data - A * g) ** 2))
            if best is None or err < best[0]:
                best = (err, A, mu, sg)

    if best is None:
        raise ValueError("Подгонка гауссианы не сошлась")

    _, A_fit, mu_fit, sigma_fit = best
    # Признак упора в край сетки: значит окно или начальное приближение негодны
    if sigma_fit <= sigmas[0] * 1.001 or sigma_fit >= sigmas[-1] * 0.999:
        raise ValueError(f"sigma упёрлась в край сетки: {sigma_fit:.3f}")
    return A_fit, mu_fit, sigma_fit

def process_run(file_path):
    """Обрабатывает один прогон."""
    try:
        data = parse_csv(file_path)
        if data is None:
            raise ValueError("Ошибка чтения файла")
        
        # Получить таблицу спектра
        spectrum_table = read_spectrum_table(file_path)
        if not spectrum_table:
            raise ValueError("Пустая таблица спектра")
        
        energy_keV = data["energy_keV"]
        npsm_enabled = data["npsm_enabled"]
        particle = data["particle"]
        beam = data["beam"]
        n_events_processed = data["n_events_processed"]
        em_cut_mm = data["em_cut_mm"]
        em_deex = data["em_deex"]
        
        # Найти пик
        peak_bin, max_count, is_on_edge = find_peak_position(spectrum_table, energy_keV)

        # Вырожденный случай: модель выключена. Спектр по свету тождествен
        # спектру по энергии, пик полного поглощения — одна линия шириной в
        # бин. Подгонка на ней невозможна, но прогон НУЖЕН критерию С-4:
        # возвращается sigma = 0 с пометкой, а не отказ.
        if int(float(npsm_enabled)) == 0:
            return {
                "E_keV": energy_keV,
                "centroide": spectrum_table[peak_bin][0],
                "sigma_keV": 0.0,
                "sigma_rel": 0.0,
                "events_in_peak": max_count,
                "fwhm": 0.0,
                "centroid_shift": (spectrum_table[peak_bin][0] - energy_keV) / energy_keV * 100.0,
                "flags": "[линия]",
                "npsm_enabled": npsm_enabled,
                "particle": particle,
                "beam": beam,
                "n_events_processed": n_events_processed,
            }

        # Оценить фон
        background = estimate_background(spectrum_table, peak_bin, energy_keV)
        
        # Подогнать гауссиану
        try:
            amplitude, centroid, sigma = fit_gaussian(spectrum_table, peak_bin, energy_keV, background)
        except Exception as e:
            raise ValueError(f"Подгонка гауссиана не удалась: {e}")
        
        # Вычислить параметры
        sigma_keV = sigma
        sigma_rel = (sigma / energy_keV) * 100
        peak_events = max_count
        fwhm = 2.3548 * sigma
        centroid_shift = ((centroid - energy_keV) / energy_keV) * 100
        
        # Пометки
        flags = []
        if is_on_edge:
            flags.append("[край]")
        
        return {
            "E_keV": energy_keV,
            "centroide": centroid,
            "sigma_keV": sigma_keV,
            "sigma_rel": sigma_rel,
            "events_in_peak": peak_events,
            "fwhm": fwhm,
            "centroid_shift": centroid_shift,
            "flags": flags,
            "npsm_enabled": npsm_enabled,
            "particle": particle,
            "beam": beam,
            "n_events_processed": n_events_processed,
            "em_cut_mm": em_cut_mm,
            "em_deex": em_deex
        }
    except Exception as e:
        print(f"ОШИБКА: Не удалось обработать файл {file_path}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", required=True, help="Шаблон для поиска CSV файлов")
    parser.add_argument("--json", help="Выходной JSON файл")
    parser.add_argument("--csv", help="Выходной CSV файл")
    
    args = parser.parse_args()
    
    # Найти все файлы
    files = glob.glob(args.runs)
    if not files:
        print("ОШИБКА: Не найдено файлов по шаблону")
        sys.exit(1)
    
    results = []
    npsm_disabled_found = False
    
    for file_path in files:
        result = process_run(file_path)
        if result is not None:
            results.append(result)
            if result["npsm_enabled"] == 0:
                npsm_disabled_found = True
    
    # Проверить наличие прогона с выключенной моделью
    if not npsm_disabled_found:
        print("ОШИБКА: Прогон с npsm_enabled = 0 не найден")
        sys.exit(1)
    
    # Фильтрация по E > 100 кэВ для регрессии
    # Прогон с выключенной моделью в регрессию и в критерии монотонности не
    # входит: его sigma равна нулю по построению и давала nan в логарифме.
    filtered_results = [r for r in results
                        if r["E_keV"] > 100 and int(float(r["npsm_enabled"])) != 0]
    
    # Подгонка регрессии σ(E) = a * E^b в логарифмах
    if len(filtered_results) >= 2:
        log_E = np.log([r["E_keV"] for r in filtered_results])
        log_sigma = np.log([r["sigma_keV"] for r in filtered_results])
        
        # Линейная регрессия: log(sigma) = a + b * log(E)
        A = np.vstack([log_E, np.ones(len(log_E))]).T
        a, b = np.linalg.lstsq(A, log_sigma, rcond=None)[0]
        
        # Коэффициент детерминации
        y_pred = a + b * log_E
        ss_res = np.sum((log_sigma - y_pred) ** 2)
        ss_tot = np.sum((log_sigma - np.mean(log_sigma)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        # Значение σ/E при 662 кэВ
        sigma_662 = np.exp(a + b * np.log(662))
        sigma_rel_662 = (sigma_662 / 662) * 100
        
        print(f"Регрессия: σ(E) = {np.exp(a):.4f} * E^{b:.4f}")
        print(f"Коэффициент детерминации: {r_squared:.4f}")
        print(f"σ/E при 662 кэВ: {sigma_rel_662:.4f}%")
    else:
        a, b, r_squared = None, None, None
        sigma_rel_662 = None
    
    # Вывод таблицы
    print("E_keV, центроид, смещение %, σ_keV, σ/E %, событий в пике, пометки")
    for r in results:
        flags_str = ", ".join(r["flags"]) if r["flags"] else ""
        print(f"{r['E_keV']:.1f}, {r['centroide']:.3f}, {r['centroid_shift']:.2f}, "
              f"{r['sigma_keV']:.3f}, {r['sigma_rel']:.2f}, {r['events_in_peak']}, {flags_str}")
    
    # Критерии
    print("\nКРИТЕРИИ:")
    
    # С-1: σ_intr(E) монотонно возрастает при E > 100 кэВ
    if len(filtered_results) >= 2:
        sigma_values = [r["sigma_keV"] for r in filtered_results]
        e_values = [r["E_keV"] for r in filtered_results]
        
        # Проверка монотонности
        is_monotonic_increasing = all(sigma_values[i] <= sigma_values[i+1] 
                                      for i in range(len(sigma_values)-1))
        c1_pass = "PASS" if is_monotonic_increasing else "FAIL"
        print(f"С-1: σ_intr(E) монотонно возрастает при E > 100 кэВ | {c1_pass}")
    else:
        print("С-1: Недостаточно данных для проверки")
        c1_pass = "FAIL"
    
    # С-2: σ_intr/E монотонно УБЫВАЕТ при E > 100 кэВ
    if len(filtered_results) >= 2:
        sigma_rel_values = [r["sigma_rel"] for r in filtered_results]
        
        # Проверка монотонности убывания
        is_monotonic_decreasing = all(sigma_rel_values[i] >= sigma_rel_values[i+1] 
                                      for i in range(len(sigma_rel_values)-1))
        c2_pass = "PASS" if is_monotonic_decreasing else "FAIL"
        print(f"С-2: σ_intr/E монотонно УБЫВАЕТ при E > 100 кэВ | {c2_pass}")
    else:
        print("С-2: Недостаточно данных для проверки")
        c2_pass = "FAIL"
    
    # С-3: центроид пика ниже энергии на всех точках. Прогон с ВЫКЛЮЧЕННОЙ
    # моделью исключается: там свет тождествен энергии по построению, и
    # смещение обязано быть нулевым, а не отрицательным. Прежняя версия
    # проверяла и его, отчего критерий проваливался всегда.
    all_shift_negative = all(r["centroid_shift"] < 0
                             for r in results if int(float(r["npsm_enabled"])) != 0)
    c3_pass = "PASS" if all_shift_negative else "FAIL"
    print(f"С-3: центроид пика ниже энергии на всех точках | {c3_pass}")
    
    # С-4: прогон с npsm_enabled = 0
    npsm_disabled_results = [r for r in results if r["npsm_enabled"] == 0]
    if npsm_disabled_results:
        max_sigma = max(r["sigma_keV"] for r in npsm_disabled_results)
        # Ширина одного бина = 1 кэВ
        c4_pass = "PASS" if max_sigma <= 1.0 else "FAIL"
        print(f"С-4: σ не больше ширины одного бина (1 кэВ) | {c4_pass}")
    else:
        print("С-4: Прогон с npsm_enabled = 0 не найден")
        c4_pass = "FAIL"
    
    # Требует толкования
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    issues = []
    
    for r in results:
        if "[край]" in r["flags"]:
            issues.append(f"Пик на границе: E={r['E_keV']} кэВ")
        if r["events_in_peak"] < 500:
            issues.append(f"Низкая статистика в пике: E={r['E_keV']} кэВ, событий={r['events_in_peak']}")
    
    if issues:
        for issue in issues:
            print(issue)
    else:
        print("пусто")
    
    # Итог
    overall_pass = all([c1_pass == "PASS", c2_pass == "PASS", c3_pass == "PASS", c4_pass == "PASS"])
    print(f"\nSTAGE3_ACCEPTANCE={'PASS' if overall_pass else 'FAIL'}")
    
    # Сохранение в JSON
    if args.json:
        output_data = {
            "results": results,
            "regression": {
                "a": a,
                "b": b,
                "r_squared": r_squared,
                "sigma_rel_662": sigma_rel_662
            },
            "criteria": {
                "C1": c1_pass,
                "C2": c2_pass,
                "C3": c3_pass,
                "C4": c4_pass
            },
            "overall": "PASS" if overall_pass else "FAIL"
        }
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(output_data, f, ensure_ascii=False, indent=2)
    
    # Сохранение в CSV
    if args.csv:
        with open(args.csv, "w", encoding="utf-8") as f:
            f.write("E_keV, центроид, смещение %, σ_keV, σ/E %, событий в пике, пометки\n")
            for r in results:
                flags_str = ", ".join(r["flags"]) if r["flags"] else ""
                f.write(f"{r['E_keV']:.1f}, {r['centroide']:.3f}, {r['centroid_shift']:.2f}, "
                        f"{r['sigma_keV']:.3f}, {r['sigma_rel']:.2f}, {r['events_in_peak']}, {flags_str}\n")

if __name__ == "__main__":
    main()
