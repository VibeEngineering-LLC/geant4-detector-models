import os
import sys
import csv
import json
import math
import numpy as np
from collections import namedtuple

sys.stdout.reconfigure(encoding="utf-8")

def read_csv(filename):
    with open(filename, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = {}
        in_spectrum = False
        data = []
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            # В файле ДВЕ таблицы: сначала гистограмма числа рассеяний
            # (n_compt,count_absorbed,count_all), затем спектр. Обе трёхколоночные,
            # поэтому различать их по числу колонок нельзя — собираем строки
            # ТОЛЬКО после заголовка bin_keV.
            if row[0] == 'bin_keV':
                in_spectrum = True
                continue
            if row[0] == 'n_compt':
                in_spectrum = False
                continue
            if not in_spectrum:
                if len(row) == 2:
                    header[row[0]] = row[1]
                continue
            if len(row) >= 3:
                data.append({'bin_keV': float(row[0]),
                             'count_edep': float(row[1]),
                             'count_light': float(row[2])})
    if not data:
        print("Ошибка: таблица спектра не найдена в", filename, file=sys.stderr)
        sys.exit(3)
    return header, data

# Функция для поиска пика полного поглощения
def read_spe(filename):
    try:
        # Путь к скриптам гамма-анализа — из переменной окружения:
        # зашитый абсолютный путь работал только на машине автора и
        # раскрывал структуру его дисков в публичном репозитории
        # (убрано 07.09.2026 перед публикацией).
        _gsa = os.environ.get("GAMMA_SPECTRUM_ANALYSIS_SCRIPTS")
        if _gsa:
            sys.path.insert(0, _gsa)
        from gamma.io.lsrm_spe import read_lsrm_spe
        return read_lsrm_spe(filename)
    except ImportError:
        print("Ошибка: не удалось импортировать модуль read_lsrm_spe", file=sys.stderr)
        sys.exit(3)

def calculate_sigma(E):
    # ПШПВ(E) = -6.2914323508675 + 1.7442782325498*sqrt(E) + 0.0049483216398*E
    fwhm = -6.2914323508675 + 1.7442782325498 * math.sqrt(E) + 0.0049483216398 * E
    return fwhm / (2 * math.sqrt(2 * math.log(2)))

def find_peak_position(bins, counts, energy_keV):
    # ⚠ Таблица спектра в CSV РАЗРЕЖЕНА: записаны только непустые бины, и их
    # индекс в массиве НЕ равен энергии. Поэтому область поиска задаётся в
    # энергиях по массиву bins, а не по индексам.
    # Нижняя граница 0.3*E0 (не 3*E0: десятичная запятая уже подводила).
    # Бин ниже 1 кэВ отбрасывается всегда: там сидят события без депозита —
    # в пробном прогоне это 9602 отсчёта из 10000, они забивают максимум.
    lo, hi = 0.3 * energy_keV, 1.05 * energy_keV
    idx = [i for i, b in enumerate(bins) if lo <= b <= hi and b >= 1.0]
    if not idx:
        return None

    max_pos = max(idx, key=lambda i: counts[i])
    max_count = counts[max_pos]
    if max_count <= 0:
        return None
    
    # Уточняем положение центром тяжести ПО ЭНЕРГИИ (bins), а не по индексу.
    total = 0.0
    weighted_sum = 0.0
    for i in range(max(0, max_pos - 5), min(len(counts), max_pos + 6)):
        if counts[i] > max_count / 2:
            total += counts[i]
            weighted_sum += bins[i] * counts[i]

    if total == 0:
        return None

    return weighted_sum / total

def calibrate_light_to_energy(bins, count_light, energy_keV):
    # Найти пик полного поглощения В ШКАЛЕ СВЕТА (положение уже в кэВ света)
    peak_pos = find_peak_position(bins, count_light, energy_keV)
    if peak_pos is None:
        raise ValueError("Не удалось найти пик полного поглощения")

    # Калибровка умножает ШКАЛУ, счёт не трогается.
    k = energy_keV / peak_pos
    calibrated_centers = [b * k for b in bins]
    return calibrated_centers, k

def convolve_with_gaussian(counts, centers, sigma_func):
    # Свертка с гауссианой переменной ширины
    result = np.zeros_like(counts)
    
    for i in range(len(counts)):
        sigma = sigma_func(centers[i])
        # Определяем диапазон для свёртки (±4σ)
        start_bin = max(0, int(i - 4 * sigma))
        end_bin = min(len(counts), int(i + 4 * sigma) + 1)
        
        norm = 0.0
        for j in range(start_bin, end_bin):
            dx = centers[j] - centers[i]
            gaussian = math.exp(-0.5 * (dx / sigma)**2)
            result[i] += counts[j] * gaussian
            norm += gaussian
        
        if norm > 0:
            result[i] /= norm
    
    return result

def interpolate_spectrum(spectrum, energy_cal, grid):
    # Интерполируем спектр на общую сетку
    # energy_cal - коэффициенты полинома (по возрастанию степени)
    # grid - массив энергий (шаг 1 кэВ от 50 до 1000)
    
    def channel_to_energy(ch):
        return sum(c * (ch ** i) for i, c in enumerate(energy_cal))
    
    # ⚠ Полином переводит КАНАЛ в ЭНЕРГИЮ. Подставлять в него энергию, чтобы
    # получить канал, нельзя — так шкала уезжала, и в районе 662 кэВ выходили
    # нули. Строим энергию каждого канала и интерполируем по этой шкале.
    n = len(spectrum)
    energies = np.array([channel_to_energy(ch) for ch in range(n)], dtype=float)

    # Плотность отсчётов на кэВ: шаг канала переменный (полином 3-й степени),
    # без деления на ширину сравнение формы было бы искажено.
    widths = np.gradient(energies)
    widths[widths <= 0] = np.nan
    density = np.asarray(spectrum, dtype=float) / widths

    good = np.isfinite(density) & np.isfinite(energies)
    return np.interp(np.asarray(grid, dtype=float),
                     energies[good], density[good], left=0.0, right=0.0)

def calculate_fwhm(peak_spectrum, peak_energy, grid):
    # Находим пик в окне ±60 кэВ вокруг 661.658
    start_idx = max(0, int(peak_energy - 60) - 50)
    end_idx = min(len(grid), int(peak_energy + 60) - 50 + 1)
    
    if start_idx >= len(peak_spectrum) or end_idx <= 0:
        return None
    
    # Находим максимум и ЕГО ПОЛОЖЕНИЕ — спуск пойдёт от вершины
    window = list(peak_spectrum[start_idx:end_idx])
    if not window or max(window) <= 0:
        return None
    max_val = max(window)
    max_idx = start_idx + window.index(max_val)
    half_max = max_val / 2.0
    
    # Ищем левую и правую точки половины высоты
    left_ch = -1
    right_ch = -1
    
    # Спуск от ВЕРШИНЫ в обе стороны. Искать первое превышение половины от
    # края окна нельзя: слева от пика лежит комптоновский континуум, и точка
    # полуспада уехала бы на край окна.
    for i in range(max_idx, start_idx - 1, -1):
        if peak_spectrum[i] < half_max:
            left_ch = i
            break
    for i in range(max_idx, end_idx):
        if peak_spectrum[i] < half_max:
            right_ch = i
            break

    if left_ch == -1 or right_ch == -1:
        return None
    
    # Линейная интерполяция между соседними бинами
    def interpolate_val(ch):
        ch_floor = int(ch)
        ch_ceil = ch_floor + 1
        if ch_ceil >= len(peak_spectrum):
            return peak_spectrum[ch_floor]
        weight = ch - ch_floor
        return peak_spectrum[ch_floor] * (1 - weight) + peak_spectrum[ch_ceil] * weight
    
    left_energy = grid[left_ch]
    right_energy = grid[right_ch]
    
    # Используем линейную интерполяцию для точного определения
    # Проверим, что между соседними бинами есть пересечение
    if left_ch + 1 < len(peak_spectrum):
        val_left = interpolate_val(left_ch)
        val_right = interpolate_val(left_ch + 1)
        if val_left >= half_max and val_right < half_max:
            # Линейная интерполяция для левой точки
            t = (half_max - val_left) / (val_right - val_left)
            left_energy = grid[left_ch] + t * (grid[left_ch + 1] - grid[left_ch])
    
    if right_ch > 0 and right_ch < len(peak_spectrum):
        val_left = interpolate_val(right_ch - 1)
        val_right = interpolate_val(right_ch)
        if val_left >= half_max and val_right < half_max:
            # Линейная интерполяция для правой точки
            t = (half_max - val_left) / (val_right - val_left)
            right_energy = grid[right_ch - 1] + t * (grid[right_ch] - grid[right_ch - 1])
    
    return right_energy - left_energy

def calculate_chi2(spectrum1, spectrum2, sigma1, sigma2):
    # Вычисляем χ² по бинам
    chi2 = 0.0
    n = 0
    
    for i in range(len(spectrum1)):
        if spectrum1[i] <= 0 or spectrum2[i] <= 0:
            continue
        
        diff = spectrum1[i] - spectrum2[i]
        var = sigma1[i]**2 + sigma2[i]**2
        if var > 0:
            chi2 += (diff**2) / var
            n += 1
    
    if n == 0:
        return 0.0, 0.0
    
    chi2_per_dof = chi2 / n if n > 0 else 0.0
    return chi2, chi2_per_dof

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--on', required=True)
    parser.add_argument('--off', required=True)
    parser.add_argument('--spe', required=True)
    parser.add_argument('--json')
    
    args = parser.parse_args()
    
    # Чтение данных
    header_on, data_on = read_csv(args.on)
    header_off, data_off = read_csv(args.off)
    spe_data = read_spe(args.spe)
    
    # Проверка npsm_enabled
    if header_on['npsm_enabled'] != '1':
        print("Ошибка: --on прогон должен иметь npsm_enabled=1", file=sys.stderr)
        sys.exit(3)
    if header_off['npsm_enabled'] != '0':
        print("Ошибка: --off прогон должен иметь npsm_enabled=0", file=sys.stderr)
        sys.exit(3)
    
    # Извлечение параметров
    energy_keV = float(header_on['energy_keV'])
    n_events_processed = int(header_on['n_events_processed'])
    npsm_bad_s_count = int(header_on['npsm_bad_s_count'])
    npsm_out_of_range_count = int(header_on['npsm_out_of_range_count'])
    
    # Преобразование данных в массивы
    bin_keV = [d['bin_keV'] for d in data_on]
    count_light_on = [d['count_light'] for d in data_on]
    count_light_off = [d['count_light'] for d in data_off]
    
    # Калибровка света в энергию
    # Шкала берётся из поля bin_keV самой таблицы: она РАЗРЕЖЕНА, индекс
    # массива энергией не является.
    bins_on = [d['bin_keV'] for d in data_on]
    bins_off = [d['bin_keV'] for d in data_off]
    calibrated_centers_on, k_on = calibrate_light_to_energy(bins_on, count_light_on, energy_keV)
    calibrated_centers_off, k_off = calibrate_light_to_energy(bins_off, count_light_off, energy_keV)
    
    # Размытие спектров
    sigma_func = calculate_sigma
    
    def sigma_func_adjusted(E):
        return sigma_func(E) * 0.798

    # Общая сетка
    grid = list(range(50, 1001))  # от 50 до 1000 кэВ

    # ⚠ ПОРЯДОК: сначала на общую сетку, ПОТОМ свёртка. Свёртка на исходной
    # шкале расчёта теряет правое крыло пика: при выключенной модели свет
    # тождествен энергии, шкала обрывается ровно на 661.7 кэВ, и размытие
    # некуда положить — спектр за пиком не спадал, ПШПВ не определялась.
    # За пределами расчётной шкалы ставится ноль, а не последнее значение.
    measured_spectrum = interpolate_spectrum(spe_data.counts, spe_data.energy_cal, grid)
    raw_off = np.interp(grid, calibrated_centers_off, count_light_off, left=0.0, right=0.0)
    raw_on = np.interp(grid, calibrated_centers_on, count_light_on, left=0.0, right=0.0)

    convolved_off_interp = np.asarray(convolve_with_gaussian(raw_off, grid, sigma_func))
    convolved_on_interp = np.asarray(convolve_with_gaussian(raw_on, grid, sigma_func))
    convolved_on_adjusted_interp = np.asarray(
        convolve_with_gaussian(raw_on, grid, sigma_func_adjusted))
    
    # Нормировка на площадь пика 662
    peak_energy = 661.658
    
    # Определяем окно ±3σ вокруг 661.658
    sigma_peak = calculate_sigma(peak_energy)
    window_size = int(3 * sigma_peak)
    
    start_idx = max(0, int(peak_energy - window_size) - 50)
    end_idx = min(len(grid), int(peak_energy + window_size) - 50 + 1)
    
    # Вычисляем площадь пика
    def calculate_area(spectrum, start_idx, end_idx):
        # Подложка: среднее из пяти бинов у каждого края
        left_baseline = np.mean(spectrum[start_idx:start_idx+5])
        right_baseline = np.mean(spectrum[end_idx-5:end_idx])
        
        baseline = np.linspace(left_baseline, right_baseline, end_idx - start_idx)
        
        area = 0.0
        for i in range(start_idx, end_idx):
            area += (spectrum[i] - baseline[i - start_idx]) * 1.0  # шаг 1 кэВ
        
        return area
    
    # Нормируем все спектры
    area_off = calculate_area(convolved_off_interp, start_idx, end_idx)
    area_on = calculate_area(convolved_on_interp, start_idx, end_idx)
    area_on_adjusted = calculate_area(convolved_on_adjusted_interp, start_idx, end_idx)
    area_measured = calculate_area(measured_spectrum, start_idx, end_idx)
    
    # Нормируем
    norm_off = convolved_off_interp / area_off if area_off > 0 else convolved_off_interp
    norm_on = convolved_on_interp / area_on if area_on > 0 else convolved_on_interp
    norm_on_adjusted = convolved_on_adjusted_interp / area_on_adjusted if area_on_adjusted > 0 else convolved_on_adjusted_interp
    norm_measured = measured_spectrum / area_measured if area_measured > 0 else measured_spectrum
    
    # Определяем ПШПВ для каждого варианта
    fwhm_off = calculate_fwhm(norm_off, peak_energy, grid)
    fwhm_on = calculate_fwhm(norm_on, peak_energy, grid)
    fwhm_on_adjusted = calculate_fwhm(norm_on_adjusted, peak_energy, grid)
    
    # Определяем положение пика
    def find_peak_position_in_spectrum(spectrum, energy):
        start_idx = max(0, int(energy - 60) - 50)
        end_idx = min(len(grid), int(energy + 60) - 50 + 1)
        
        peak_val = -1
        peak_pos = -1
        
        for i in range(start_idx, end_idx):
            if spectrum[i] > peak_val:
                peak_val = spectrum[i]
                peak_pos = i
        
        return grid[peak_pos] if peak_pos != -1 else None
    
    pos_off = find_peak_position_in_spectrum(norm_off, peak_energy)
    pos_on = find_peak_position_in_spectrum(norm_on, peak_energy)
    pos_on_adjusted = find_peak_position_in_spectrum(norm_on_adjusted, peak_energy)
    
    # ПШПВ измеренного пика берётся ИЗ ТАБЛИЦЫ ПИКОВ САМОГО ФАЙЛА (#CAL-0):
    # прибор её уже посчитал, свой пересчёт допустим только как проверка и
    # только с предъявлением обеих величин рядом (урок W-063: самодельный
    # пересчёт дал +12 %).
    tbl = (spe_data.extras or {}).get("lsrm_peaks_table") or []
    ref_peak = next((p for p in tbl if abs(float(p['energy_keV']) - peak_energy) < 1.0), None)
    if ref_peak is None:
        print("Ошибка: в таблице пиков файла нет линии 662 кэВ", file=sys.stderr)
        sys.exit(3)
    fwhm_measured = float(ref_peak['fwhm_keV'])
    fwhm_measured_own = calculate_fwhm(norm_measured, peak_energy, grid)
    print(f"ПШПВ измеренного пика: {fwhm_measured:.3f} кэВ (таблица пиков файла); "
          f"свой пересчёт для проверки: "
          f"{('%.3f кэВ' % fwhm_measured_own) if fwhm_measured_own else 'не определён'}")
    
    # Вычисляем χ² по континууму 200…450 кэВ
    cont_start_idx = 200 - 50
    cont_end_idx = 450 - 50 + 1
    
    # ⚠ Пуассоновской является СЫРАЯ статистика, а не нормированная. Если
    # norm = raw / A, то σ_norm = √raw / A = √(norm / A). Брать √(norm), как
    # было, значит завысить погрешность в разы и получить бессмысленно малый χ².
    def sigma_norm(arr, area):
        return np.sqrt(np.clip(np.asarray(arr, dtype=float), 0, None) / area) if area > 0 else np.ones_like(arr)

    sl = slice(cont_start_idx, cont_end_idx)
    sig_meas = sigma_norm(norm_measured[sl], area_measured)

    chi2_a, chi2_per_dof_a = calculate_chi2(
        norm_off[sl], norm_measured[sl], sigma_norm(norm_off[sl], area_off), sig_meas)

    chi2_b, chi2_per_dof_b = calculate_chi2(
        norm_on[sl], norm_measured[sl], sigma_norm(norm_on[sl], area_on), sig_meas)

    chi2_v, chi2_per_dof_v = calculate_chi2(
        norm_on_adjusted[sl], norm_measured[sl],
        sigma_norm(norm_on_adjusted[sl], area_on_adjusted), sig_meas)
    
    # Критерии
    criteria = {}
    
    # Г-1: положение пика 662 во всех трёх вариантах
    delta_off = abs(pos_off - peak_energy) if pos_off is not None else float('inf')
    delta_on = abs(pos_on - peak_energy) if pos_on is not None else float('inf')
    delta_v = abs(pos_on_adjusted - peak_energy) if pos_on_adjusted is not None else float('inf')
    
    criteria['G1'] = {
        'threshold': 0.25 * fwhm_measured,
        'value': max(delta_off, delta_on, delta_v),
        'pass': max(delta_off, delta_on, delta_v) <= 0.25 * fwhm_measured
    }
    
    # Г-2: ПШПВ варианта (в)
    criteria['G2'] = {
        'threshold': 41.85 * 0.05,
        'value': abs(fwhm_on_adjusted - 41.85),
        'pass': abs(fwhm_on_adjusted - 41.85) <= 41.85 * 0.05
    }
    
    # Г-3: ПШПВ варианта (б)
    criteria['G3'] = {
        'threshold': None,
        'value': fwhm_on,
        'pass': fwhm_on > fwhm_measured
    }
    
    # Г-4: невязка по континууму 200…450 кэВ.
    # Улучшением считается ТОЛЬКО значимое уменьшение: разность χ² должна
    # превышать статистическую погрешность самой величины χ², которая для n
    # степеней свободы равна sqrt(2n). Простое «стало меньше» критерием не
    # является — на таком пороге любая пара прогонов даёт «улучшение».
    n_dof_g4 = (chi2_a / chi2_per_dof_a) if chi2_per_dof_a else 0.0
    g4_sigma = math.sqrt(2.0 * n_dof_g4) if n_dof_g4 > 0 else float('inf')
    criteria['G4'] = {
        'threshold': g4_sigma,
        'value': chi2_a - chi2_v,          # положительное = (в) лучше
        'pass': (chi2_a - chi2_v) > g4_sigma
    }
    
    # Г-5: счётчики модели в обоих прогонах
    criteria['G5'] = {
        'threshold': "bad_s=0 и out_of_range=0",
        'value': f"{npsm_bad_s_count} и {npsm_out_of_range_count}",
        'pass': npsm_bad_s_count == 0 and npsm_out_of_range_count == 0
    }
    
    # Вывод таблицы
    print("Вариант\tПозиция пика (кэВ)\tПШПВ (кэВ)\tχ²/n")
    print(f"(а)\t{pos_off:.2f}\t{fwhm_off:.3f}\t{chi2_per_dof_a:.4f}")
    print(f"(б)\t{pos_on:.2f}\t{fwhm_on:.3f}\t{chi2_per_dof_b:.4f}")
    print(f"(в)\t{pos_on_adjusted:.2f}\t{fwhm_on_adjusted:.3f}\t{chi2_per_dof_v:.4f}")
    
    # Критерии
    print("\nКРИТЕРИИ")
    for key, crit in criteria.items():
        if key == 'G1':
            print(f"{key}: {crit['value']:.3f} кэВ ≤ {crit['threshold']:.3f} кэВ → {'PASS' if crit['pass'] else 'FAIL'}")
        elif key == 'G2':
            print(f"{key}: ПШПВ(в) {fwhm_on_adjusted:.3f} против измеренных {fwhm_measured:.3f}; "
                  f"|разность| = {crit['value']:.3f} кэВ ≤ {crit['threshold']:.3f} кэВ → {'PASS' if crit['pass'] else 'FAIL'}")
        elif key == 'G3':
            print(f"{key}: {crit['value']:.3f} кэВ > {fwhm_measured:.3f} кэВ → {'PASS' if crit['pass'] else 'FAIL'}")
        elif key == 'G4':
            print(f"{key}: χ²(а) {chi2_a:.1f} − χ²(в) {chi2_v:.1f} = {crit['value']:.1f} "
                  f"при пороге значимости {crit['threshold']:.1f} → "
                  f"{'PASS (улучшение значимо)' if crit['pass'] else 'FAIL (значимого улучшения нет)'}")
        elif key == 'G5':
            print(f"{key}: {crit['value']} → {'PASS' if crit['pass'] else 'FAIL'}")
    
    # Требует толкования
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ")
    # Счёт проверяется по СЫРЫМ отсчётам расчёта: нормированная величина —
    # доля площади пика, и порог «500 отсчётов» к ней неприменим.
    cont_sum = float(np.sum(raw_on[cont_start_idx:cont_end_idx]))
    if cont_sum < 500:
        print(f"Меньше 500 отсчётов расчёта в полосе континуума (200…450 кэВ): {cont_sum:.0f}")
    
    # Погрешность ширины пика ≈ ПШПВ/√(2N) по числу СЫРЫХ событий в пике,
    # а не √ПШПВ: корень от ширины к статистике отношения не имеет и давал
    # ложное срабатывание на каждом прогоне.
    n_peak = float(np.sum(raw_on[start_idx:end_idx]))
    d_fwhm = (fwhm_on_adjusted / math.sqrt(2.0 * n_peak)) if n_peak > 0 else float('inf')
    if abs(fwhm_off - fwhm_on_adjusted) < 2.0 * d_fwhm:
        print(f"ПШПВ (а) {fwhm_off:.3f} и (в) {fwhm_on_adjusted:.3f} различаются на "
              f"{abs(fwhm_off - fwhm_on_adjusted):.3f} кэВ при погрешности {d_fwhm:.3f} — "
              f"различие незначимо")
    
    if pos_off is not None and pos_on is not None:
        delta_pos = abs(pos_off - pos_on)
        if delta_pos > 2.0:
            print(f"Расхождение положения пика (а) и (б) больше 2 кэВ: {delta_pos:.3f} кэВ")
    
    # Итог
    improvement = "NO"
    if criteria['G4']['pass']:
        improvement = "YES"
    elif criteria['G1']['pass'] and criteria['G2']['pass'] and criteria['G3']['pass'] and criteria['G5']['pass']:
        improvement = "INDISTINGUISHABLE"
    
    print(f"\nSTAGE5_ACCEPTANCE={'PASS' if all(crit['pass'] for crit in criteria.values()) else 'FAIL'}")
    print(f"STAGE5_IMPROVEMENT={improvement}")
    
    # Сохранение в JSON
    if args.json:
        result = {
            "variants": [
                {"variant": "a", "peak_position": pos_off, "fwhm": fwhm_off, "chi2_per_dof": chi2_per_dof_a},
                {"variant": "b", "peak_position": pos_on, "fwhm": fwhm_on, "chi2_per_dof": chi2_per_dof_b},
                {"variant": "v", "peak_position": pos_on_adjusted, "fwhm": fwhm_on_adjusted, "chi2_per_dof": chi2_per_dof_v}
            ],
            "criteria": criteria,
            "improvement": improvement
        }
        
        # numpy-типы (bool_, float64) стандартным сериализатором не берутся.
        def _plain(o):
            if isinstance(o, (bool, np.bool_)):
                return bool(o)
            return float(o)

        with open(args.json, 'w', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2, default=_plain)

if __name__ == "__main__":
    main()
