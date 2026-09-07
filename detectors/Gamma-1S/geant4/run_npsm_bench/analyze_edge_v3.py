import sys
import os
import csv
import json
import argparse
import numpy as np
from collections import OrderedDict
from scipy.interpolate import interp1d

# Для сериализации numpy типов в JSON
def numpy_default(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

# Чтение CSV файла с ключами и таблицами
def read_calculation_csv(filename):
    header = {}
    histogram = []
    spectrum = []
    
    with open(filename, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Парсим заголовки. ⚠ Условие «пока в строке есть запятая» негодно с обеих
    # сторон: первая строка файла — комментарий `# npsm_bench …` без запятой,
    # и цикл обрывался на ней, оставляя шапку ПУСТОЙ; а строки обеих таблиц
    # запятые содержат и попали бы в шапку. Признак конца шапки — заголовок
    # таблицы, а полем шапки считается только двухколоночная строка.
    i = 0
    while i < len(lines):
        s = lines[i]
        if s.startswith('n_compt') or s.startswith('bin_keV'):
            break
        if not s.startswith('#'):
            parts = s.split(',')
            if len(parts) == 2:
                header[parts[0]] = parts[1]
        i += 1
    
    # Пропускаем пустые строки
    while i < len(lines) and not lines[i]:
        i += 1
    
    # Парсим гистограмму
    if i < len(lines) and lines[i].startswith('n_compt'):
        i += 1
        while i < len(lines) and ',' in lines[i]:
            # Заголовок СЛЕДУЮЩЕЙ таблицы обрывает эту: обе трёхколоночные,
            # и без явной проверки 'bin_keV' уходило в int().
            if lines[i].startswith('bin_keV'):
                break
            parts = lines[i].split(',')
            if len(parts) == 3:
                histogram.append([int(parts[0]), int(parts[1]), int(parts[2])])
            else:
                break
            i += 1
    
    # Пропускаем пустые строки
    while i < len(lines) and not lines[i]:
        i += 1
    
    # Парсим спектр
    if i < len(lines) and lines[i].startswith('bin_keV'):
        i += 1
        while i < len(lines) and ',' in lines[i]:
            parts = lines[i].split(',')
            if len(parts) == 3:
                spectrum.append([float(parts[0]), int(parts[1]), int(parts[2])])
            else:
                break
            i += 1
    
    return header, histogram, spectrum

# Чтение измеренного спектра
def read_measured_spe(filename):
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
        print("Ошибка: не удалось импортировать read_lsrm_spe из gamma.io.lsrm_spe", file=sys.stderr)
        sys.exit(3)

# Функция для вычисления FWHM_W(E) и sigma_W(E)
def fwhm_w(e):
    return -6.2914323508675 + 1.7442782325498 * np.sqrt(e) + 0.0049483216398 * e

def sigma_w(e):
    return fwhm_w(e) / (2 * np.sqrt(2 * np.log(2)))

# Функция для свертки с Гауссом
def convolve_with_gaussian(data, sigma):
    # Создаем ядро Гаусса
    kernel_size = max(1, int(4 * sigma + 0.5))
    x = np.arange(-kernel_size, kernel_size + 1)
    kernel = np.exp(-0.5 * (x / sigma)**2)
    kernel /= np.sum(kernel)
    
    # Свертка
    result = np.convolve(data, kernel, mode='same')
    return result

# Функция для интерполяции на общую сетку
def interpolate_on_grid(x, y, grid):
    f = interp1d(x, y, kind='linear', fill_value=0.0, bounds_error=False)
    return f(grid)

# Функция для поиска пика полного поглощения
def find_absorption_peak(counts, bin_keV, e0):
    # Ограничиваем область поиска
    lo = 0.3 * e0
    hi = 1.05 * e0
    
    # Фильтруем бины
    valid_indices = np.where((bin_keV >= lo) & (bin_keV <= hi) & (bin_keV >= 1.0))[0]
    
    if len(valid_indices) == 0:
        return None, None
    
    # Сортируем по энергии
    sorted_indices = valid_indices[np.argsort(bin_keV[valid_indices])]
    
    # Находим связную группу
    max_count = np.max(counts[sorted_indices])
    threshold = max_count / 2.0
    
    # Ищем максимальную значимую группу
    groups = []
    current_group = []
    
    for i in sorted_indices:
        if counts[i] >= threshold:
            current_group.append(i)
        else:
            if len(current_group) > 0:
                groups.append(current_group)
                current_group = []
    
    if len(current_group) > 0:
        groups.append(current_group)
    
    # Выбираем значимую группу
    if not groups:
        return None, None
    
    group_sums = [np.sum(counts[group]) for group in groups]
    max_sum_idx = np.argmax(group_sums)
    
    if group_sums[max_sum_idx] < 0.2 * np.max(group_sums):
        return None, None
    
    # Выбираем группу с наибольшим bin_keV
    selected_group = groups[max_sum_idx]
    
    # Центр тяжести
    total_weight = np.sum(counts[selected_group])
    weighted_sum = np.sum(counts[selected_group] * bin_keV[selected_group])
    
    if total_weight > 0:
        center_of_mass = weighted_sum / total_weight
        return center_of_mass, max_count
    else:
        return None, None

# Функция для вычисления производной с помощью полиномиальной регрессии
def local_poly_derivative(data, window_size):
    n = len(data)
    derivative = np.zeros(n)
    
    for i in range(n):
        # Определяем окно
        start = max(0, i - window_size)
        end = min(n, i + window_size + 1)
        
        x_window = np.arange(start, end)
        y_window = data[start:end]
        
        if len(y_window) >= 3:
            # Подгонка параболы
            coeffs = np.polyfit(x_window, y_window, 2)
            derivative[i] = 2 * coeffs[0] * i + coeffs[1]
        else:
            derivative[i] = 0
    
    return derivative

# Функция для поиска минимума производной
def find_min_derivative(derivative, window_size):
    n = len(derivative)
    min_indices = []
    
    for i in range(window_size, n - window_size):
        # Ищем минимум в окне
        window_deriv = derivative[i - window_size:i + window_size + 1]
        min_idx_in_window = np.argmin(window_deriv)
        actual_idx = i - window_size + min_idx_in_window
        
        if actual_idx == i:
            min_indices.append(i)
    
    # Убираем дубликаты и сортируем
    unique_min_indices = list(set(min_indices))
    unique_min_indices.sort()
    
    if not unique_min_indices:
        return None, False
    
    # Проверяем, есть ли несколько минимумов
    min_values = [derivative[i] for i in unique_min_indices]
    sorted_min_values = sorted(min_values)
    
    # ⚠ Производная в минимуме ОТРИЦАТЕЛЬНА, поэтому деление на неё переворачивает
    # знак, и прежнее условие срабатывало всегда — край не находился никогда.
    # Неоднозначность: второй по глубине минимум отличается от первого менее чем
    # на 20 % ПО МОДУЛЮ. И возвращать надо индекс САМОГО ГЛУБОКОГО минимума, а не
    # первого по порядку.
    deepest = int(min(unique_min_indices, key=lambda i: derivative[i]))
    ambiguous = False
    if len(sorted_min_values) > 1 and abs(sorted_min_values[0]) > 0:
        ambiguous = abs(sorted_min_values[1]) > 0.8 * abs(sorted_min_values[0])

    return deepest, ambiguous

# Функция для уточнения минимума параболой
def refine_minimum(derivative, i_min, step):
    if i_min == 0 or i_min == len(derivative) - 1:
        return i_min
    
    d_i_minus_1 = derivative[i_min - 1]
    d_i = derivative[i_min]
    d_i_plus_1 = derivative[i_min + 1]
    
    denominator = d_i_minus_1 - 2 * d_i + d_i_plus_1
    
    if abs(denominator) < 1e-10:
        return i_min
    
    x_star = i_min + 0.5 * step * (d_i_minus_1 - d_i_plus_1) / denominator
    return x_star

# Функция для бутстрепа
def bootstrap_shift(counts_a, counts_b, w, i0, i1, grid0, n_bootstrap=200, seed=12345):
    """Погрешность сдвига розыгрышем счёта по Пуассону.

    ⚠ Прежде здесь вызывалась заглушка compute_shift, возвращавшая 0.0 «для
    примера», и погрешность выходила нулевой на любых данных. Нулевая
    погрешность опаснее отсутствующей: по ней можно объявить значимым любое
    различие. Теперь в каждом повторе полностью повторяется процедура поиска
    края — производная, минимум в окне, уточнение параболой.
    """
    rng = np.random.default_rng(seed)
    a = np.clip(np.asarray(counts_a, dtype=float), 0, None)
    b = np.clip(np.asarray(counts_b, dtype=float), 0, None)

    shifts = []
    for _ in range(n_bootstrap):
        s = compute_shift(rng.poisson(a), rng.poisson(b), w, i0, i1, grid0)
        if s is not None:
            shifts.append(s)

    return float(np.std(shifts)) if len(shifts) > 1 else float('nan')

def compute_shift(counts_a, counts_b, w, i0, i1, grid0):
    """Сдвиг между двумя кривыми на общей сетке (шаг 1 кэВ), метод перегиба."""
    da = local_poly_derivative(counts_a, w)
    db = local_poly_derivative(counts_b, w)
    ma, _ = find_min_derivative(da[i0:i1 + 1], w)
    mb, _ = find_min_derivative(db[i0:i1 + 1], w)
    if ma is None or mb is None:
        return None
    pa = grid0 + refine_minimum(da, ma + i0, 1.0)
    pb = grid0 + refine_minimum(db, mb + i0, 1.0)
    return pa - pb

# Основная функция
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['internal', 'measured'], required=True)
    parser.add_argument('--run', type=str, help='Input CSV file for internal mode')
    parser.add_argument('--calc', type=str, help='Calculation CSV file for measured mode')
    parser.add_argument('--spe', type=str, help='Measured SPE file for measured mode')
    parser.add_argument('--sigma-intr', type=float, default=1.62, help='Internal sigma in percent')
    parser.add_argument('--json', type=str, help='Output JSON file')
    parser.add_argument('--boot-seed', type=int, default=12345, help='Bootstrap seed')
    
    args = parser.parse_args()
    
    # Установка кодировки stdout
    sys.stdout.reconfigure(encoding="utf-8")
    
    # Проверка входных файлов
    if args.mode == 'internal':
        if not args.run:
            print("Ошибка: для режима internal требуется --run", file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(args.run):
            print(f"Ошибка: файл {args.run} не существует", file=sys.stderr)
            sys.exit(1)
    elif args.mode == 'measured':
        if not args.calc:
            print("Ошибка: для режима measured требуется --calc", file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(args.calc):
            print(f"Ошибка: файл {args.calc} не существует", file=sys.stderr)
            sys.exit(1)
        if not args.spe:
            print("Ошибка: для режима measured требуется --spe", file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(args.spe):
            print(f"Ошибка: файл {args.spe} не существует", file=sys.stderr)
            sys.exit(1)
    
    # Чтение данных
    if args.mode == 'internal':
        header, histogram, spectrum = read_calculation_csv(args.run)
        
        # Прогон с ВЫКЛЮЧЕННОЙ моделью подаётся сюда намеренно — это мутация
        # К-4: свет тождествен энергии, сдвиг обязан выйти нулевым. Отвергать
        # такой файл нельзя, иначе мутационную проверку негде выполнить.
        npsm_on = header.get('npsm_enabled', '0') == '1'
        if not npsm_on:
            print("[К-4] прогон с выключенной моделью: сдвиг обязан быть ~0")
        
        # Извлечение данных
        bin_keV = np.array([row[0] for row in spectrum])
        count_edep = np.array([row[1] for row in spectrum])
        count_light = np.array([row[2] for row in spectrum])
        
        # Проверка, что есть данные
        if len(bin_keV) == 0:
            print("Ошибка: пустая таблица спектра", file=sys.stderr)
            sys.exit(1)
        
        # Извлечение энергии
        e0 = float(header.get('energy_keV', '0'))
        if e0 <= 0:
            print("Ошибка: неверная энергия источника", file=sys.stderr)
            sys.exit(1)
            
        # Вычисление теоретического положения края
        e_edge_theor = e0 * (2 * e0 / 511) / (1 + 2 * e0 / 511)
        
        if e_edge_theor < 60:
            print("Ошибка: формула разрешения вне своей области применимости", file=sys.stderr)
            sys.exit(3)
            
        # Калибровка света
        peak_pos, max_count = find_absorption_peak(count_light, bin_keV, e0)
        
        if peak_pos is None:
            print("Ошибка: не удалось найти пик полного поглощения", file=sys.stderr)
            sys.exit(1)
            
        k = e0 / peak_pos
        centers_light = bin_keV * k
        
        # Вычисление sigma_intr
        s_i = args.sigma_intr * e_edge_theor / 100.0
        
        # Определение целевой ширины
        sigma_target = max(sigma_w(e_edge_theor), s_i)
        
        # Добавки к ширине для уравнивания: обе кривые доводятся до
        # sigma_target. Свет уже несёт s_i, энергия не несёт ничего.
        add_light = np.sqrt(max(sigma_target**2 - s_i**2, 0.0))
        add_edep = sigma_target
        equalized = "ok" if s_i < sigma_w(e_edge_theor) else "fallback"

        # Общая сетка
        lo = max(50.0, e_edge_theor - 6 * sigma_w(e_edge_theor))
        hi = min(1.35 * e_edge_theor, 1.05 * e0)
        
        if hi - lo < 6 * sigma_w(e_edge_theor):
            print("Ошибка: недостаточно данных для окна поиска", file=sys.stderr)
            sys.exit(3)
            
        grid = np.arange(lo, hi + 1.0, 1.0)
        
        # ⚠ ПОРЯДОК: сначала на общую сетку, ПОТОМ свёртка. Свёртка на исходной
        # шкале теряет правое крыло: при выключенной модели свет тождествен
        # энергии, шкала обрывается на энергии источника, размытию некуда лечь.
        raw_edep = interpolate_on_grid(bin_keV, count_edep, grid)
        raw_light = interpolate_on_grid(centers_light, count_light, grid)

        interp_edep = convolve_with_gaussian(raw_edep, add_edep)
        interp_light = convolve_with_gaussian(raw_light, add_light)
        
        # Вычисление производных
        w = max(3, int(np.round(sigma_w(e_edge_theor))))
        deriv_edep = local_poly_derivative(interp_edep, w)
        deriv_light = local_poly_derivative(interp_light, w)
        
        # Поиск минимума производной
        window_lo = e_edge_theor - 3 * sigma_w(e_edge_theor)
        window_hi = e_edge_theor + 3 * sigma_w(e_edge_theor)
        
        # Определяем индексы окна
        grid_indices = np.where((grid >= window_lo) & (grid <= window_hi))[0]
        
        if len(grid_indices) == 0:
            print("Ошибка: окно поиска не попадает в данные", file=sys.stderr)
            sys.exit(3)
            
        # Вычисляем положения края
        pos_a = None
        pos_b = None
        
        # ⚠ Минимум производной ищется ТОЛЬКО в окне края. Окно вычислялось, но
        # в функцию не передавалось — минимум брался по всему спектру, где он
        # приходится на спад пика полного поглощения, а не на край.
        # refine_minimum возвращает ДРОБНЫЙ ИНДЕКС сетки, не энергию:
        # переводим как E = grid[0] + индекс (шаг сетки ровно 1 кэВ).
        i0, i1 = int(grid_indices[0]), int(grid_indices[-1])

        min_idx_light, ambiguous = find_min_derivative(deriv_light[i0:i1 + 1], w)
        if min_idx_light is not None:
            pos_a = grid[0] + refine_minimum(deriv_light, min_idx_light + i0, 1.0)

        min_idx_edep, _ = find_min_derivative(deriv_edep[i0:i1 + 1], w)
        if min_idx_edep is not None:
            pos_b = grid[0] + refine_minimum(deriv_edep, min_idx_edep + i0, 1.0)

        if pos_a is None or pos_b is None:
            print("Ошибка: край не найден в окне поиска", file=sys.stderr)
            sys.exit(3)
            
        # Вычисление сдвига
        shift = pos_a - pos_b

        # К-6: тот же сдвиг на НЕУРАВНЕННЫХ кривых (свет несёт s_i, энергия
        # ничего). Разность shift и shift_raw показывает вклад неравенства
        # разрешений — это прямая проверка гипотезы Д2 плана закрытия.
        d_lr = local_poly_derivative(raw_light, w)
        d_er = local_poly_derivative(raw_edep, w)
        m_lr, _ = find_min_derivative(d_lr[i0:i1 + 1], w)
        m_er, _ = find_min_derivative(d_er[i0:i1 + 1], w)
        if m_lr is not None and m_er is not None:
            shift_raw = (refine_minimum(d_lr, m_lr + i0, 1.0)
                         - refine_minimum(d_er, m_er + i0, 1.0))
        else:
            shift_raw = float('nan')

        # Бутстреп
        shift_err = bootstrap_shift(interp_light, interp_edep, w, i0, i1, grid[0], seed=args.boot_seed)
        
        # Вывод CSV строки
        print(f"{args.mode},{e0:.3f},{e_edge_theor:.3f},{pos_a:.3f},{pos_b:.3f},"
              f"{shift:.3f},{shift_err:.3f},{shift_raw:.3f},{k:.5f},{equalized}")
        
        # Критерии
        print("КРИТЕРИИ")
        # К-1 проверяет МЕТОД: найден ли край там, где его помещает теория.
        # Прежняя запись сравнивала с порогом сам СДВИГ — это другая величина,
        # и критерий проходил всегда, ничего не проверяя.
        k1_dev = abs(pos_b - e_edge_theor)
        print(f"К-1: край по энергии {pos_b:.2f} против теории {e_edge_theor:.2f}, "
              f"расхождение {k1_dev:.2f} кэВ = {100*k1_dev/e_edge_theor:.1f} % "
              f"-> {'PASS' if k1_dev <= 0.08 * e_edge_theor else 'FAIL'} (порог 8 %)")
        print(f"К-2: знак сдвига {'PASS' if shift < 0 else 'FAIL'} (shift = {shift:.3f})")
        if not npsm_on:
            print(f"К-4: мутация, |{shift:.3f}| <= 3 кэВ -> "
                  f"{'PASS' if abs(shift) <= 3 else 'FAIL'}")
        else:
            print("К-4: не проверяется (модель включена; нужен прогон npsm=off)")
        # ⚠ shift_raw здесь ЗАТИРАЛСЯ значением pos_a - pos_b, то есть самим
        # shift, и К-6 всегда печатал ноль. Величина считается выше, на
        # неуравненных кривых, и трогать её тут нельзя.
        print(f"К-6: с уравниванием {shift:.3f}, без него {shift_raw:.3f}, "
              f"разница {abs(shift - shift_raw):.3f} кэВ")
        
        # Требует толкования
        print("ТРЕБУЕТ ТОЛКОВАНИЯ")
        if ambiguous:
            print("край размыт, перегиб неоднозначен")
        if e_edge_theor < 60:
            print("формула разрешения вне своей области применимости")
        if s_i >= sigma_w(e_edge_theor):
            print("целевое разрешение задано sigma_intr, а не прибором")
        
        # Итог
        acceptance = "PASS" if (abs(shift) <= 0.08 * e_edge_theor and shift < 0 and abs(shift) <= 3) else "FAIL"
        print(f"EDGE_V3_SHIFT={shift:.3f} keV")
        print(f"EDGE_V3_ACCEPTANCE={acceptance}")
        
    elif args.mode == 'measured':
        header, histogram, spectrum = read_calculation_csv(args.calc)
        
        # Проверка npsm_enabled
        if header.get('npsm_enabled', '0') != '0':
            print("Ошибка: для measured режима требуется npsm_enabled = 0", file=sys.stderr)
            sys.exit(3)
            
        # Извлечение данных
        bin_keV = np.array([row[0] for row in spectrum])
        count_edep = np.array([row[1] for row in spectrum])
        
        # Проверка, что есть данные
        if len(bin_keV) == 0:
            print("Ошибка: пустая таблица спектра", file=sys.stderr)
            sys.exit(1)
            
        # Извлечение энергии
        e0 = float(header.get('energy_keV', '0'))
        if e0 <= 0:
            print("Ошибка: неверная энергия источника", file=sys.stderr)
            sys.exit(1)
            
        # Вычисление теоретического положения края
        e_edge_theor = e0 * (2 * e0 / 511) / (1 + 2 * e0 / 511)
        
        if e_edge_theor < 60:
            print("Ошибка: формула разрешения вне своей области применимости", file=sys.stderr)
            sys.exit(3)
            
        # Чтение измеренного спектра
        measured_data = read_measured_spe(args.spe)
        
        counts = np.asarray(measured_data.counts, dtype=float)
        energy_cal = measured_data.energy_cal
        peaks_table = (measured_data.extras or {}).get('lsrm_peaks_table', [])
        
        # Создание шкалы энергии
        n_channels = len(counts)
        channel_indices = np.arange(n_channels)
        # ⚠ np.polyval ждёт коэффициенты по УБЫВАНИЮ степени, а energy_cal
        # задан по возрастанию — без разворота шкала выходит бессмысленной.
        energy = np.polyval(np.asarray(energy_cal)[::-1], channel_indices)
        
        # Ширина канала
        channel_width = np.gradient(energy)
        
        # Интерполяция на общую шкалу
        grid = np.arange(max(50.0, e_edge_theor - 6 * sigma_w(e_edge_theor)), 
                         min(1.35 * e_edge_theor, 1.05 * e0) + 1.0, 1.0)
        
        # Интерполяция плотности
        density = counts / channel_width
        
        # ПОРЯДОК: сначала на сетку, потом свёртка (иначе теряется крыло).
        interp_edep = convolve_with_gaussian(
            interpolate_on_grid(bin_keV, count_edep, grid), sigma_w(e_edge_theor))
        
        # Вычисление производной
        w = max(3, int(np.round(sigma_w(e_edge_theor))))
        deriv_edep = local_poly_derivative(interp_edep, w)
        
        # Поиск минимума производной
        window_lo = e_edge_theor - 3 * sigma_w(e_edge_theor)
        window_hi = e_edge_theor + 3 * sigma_w(e_edge_theor)
        
        grid_indices = np.where((grid >= window_lo) & (grid <= window_hi))[0]
        
        if len(grid_indices) == 0:
            print("Ошибка: окно поиска не попадает в данные", file=sys.stderr)
            sys.exit(3)
            
        # Окно края — как в режиме internal: минимум ищется ТОЛЬКО в нём,
        # индекс среза смещается обратно, результат переводится в энергию.
        i0, i1 = int(grid_indices[0]), int(grid_indices[-1])

        min_idx_edep, ambiguous = find_min_derivative(deriv_edep[i0:i1 + 1], w)
        pos_b = (grid[0] + refine_minimum(deriv_edep, min_idx_edep + i0, 1.0)
                 if min_idx_edep is not None else None)

        measured_interp = interpolate_on_grid(energy, density, grid)
        deriv_measured = local_poly_derivative(measured_interp, w)

        min_idx_measured, _ = find_min_derivative(deriv_measured[i0:i1 + 1], w)
        pos_a = (grid[0] + refine_minimum(deriv_measured, min_idx_measured + i0, 1.0)
                 if min_idx_measured is not None else None)

        if pos_a is None or pos_b is None:
            print("Ошибка: край не найден в окне поиска", file=sys.stderr)
            sys.exit(3)
            
        # Вычисление сдвига
        shift = pos_a - pos_b

        # shift_raw считается ЗДЕСЬ, до печати: ниже он использовался раньше,
        # чем присваивался, и режим measured падал на UnboundLocalError.
        deriv_edep_raw = local_poly_derivative(
            interpolate_on_grid(bin_keV, count_edep, grid), w)
        m_er, _ = find_min_derivative(deriv_edep_raw[i0:i1 + 1], w)
        shift_raw = (pos_a - (grid[0] + refine_minimum(deriv_edep_raw, m_er + i0, 1.0))
                     if m_er is not None else float('nan'))

        # Бутстреп
        shift_err = bootstrap_shift(measured_interp, interp_edep, w, i0, i1, grid[0], seed=args.boot_seed)

        # Вывод CSV строки
        print(f"{args.mode},{e0:.3f},{e_edge_theor:.3f},{pos_a:.3f},{pos_b:.3f},"
              f"{shift:.3f},{shift_err:.3f},{shift_raw:.3f},1.00000,measured")
        
        # Критерии
        print("КРИТЕРИИ")
        k1_dev = abs(pos_b - e_edge_theor)
        print(f"К-1: край расчёта {pos_b:.2f} против теории {e_edge_theor:.2f}, "
              f"расхождение {k1_dev:.2f} кэВ = {100*k1_dev/e_edge_theor:.1f} % "
              f"-> {'PASS' if k1_dev <= 0.08 * e_edge_theor else 'FAIL'} (порог 8 %)")
        print(f"К-2: знак сдвига {'PASS' if shift < 0 else 'FAIL'} (shift = {shift:.3f})")
        print("К-4: в режиме measured не проверяется (мутация — режим internal)")
        
        print(f"К-6: с уравниванием {shift:.3f}, без него {shift_raw:.3f}, "
              f"разница {abs(shift - shift_raw):.3f} кэВ")
        
        # Требует толкования
        print("ТРЕБУЕТ ТОЛКОВАНИЯ")
        if ambiguous:
            print("край размыт, перегиб неоднозначен")
        if e_edge_theor < 60:
            print("формула разрешения вне своей области применимости")
        
        # Итог
        acceptance = "PASS" if (abs(shift) <= 0.08 * e_edge_theor and shift < 0 and abs(shift) <= 3) else "FAIL"
        print(f"EDGE_V3_SHIFT={shift:.3f} keV")
        print(f"EDGE_V3_ACCEPTANCE={acceptance}")

if __name__ == "__main__":
    main()
