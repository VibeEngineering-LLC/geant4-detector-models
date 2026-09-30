# -*- coding: utf-8 -*-
r"""#GS-28: остаточный сдвиг шкалы GS2020 по ФОРМЕ групп линий (табличные энергии и доли IAEA/ENSDF, амплитуда на
компоненту, общий сдвиг на группу). Численная схема — донор Gamma-1S/analysis/deconv.py:241-290 (веса 1/σ, lsq_linear
trf с границами, ковариация через SVD), обобщена на N компонент; импорт донора невозможен (требует каталог сборки
Гамма-1С). Запуск: python gs2020_calib_mplet.py"""
import os, sys, csv, json, math
import numpy as np
from scipy.optimize import lsq_linear
from scipy.special import erfc
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs2020_calib as C

LINES_DIR = r"C:\g4work\gs2020\lines"
# GS28_SCALE=file: подгонка на ШКАЛЕ ФАЙЛА (без поправки по реперам) — источник узлов совместной калибровки (gs2020_calib)
FILE_SCALE = os.environ.get("GS28_SCALE", "own") == "file"
SCAN = 15.0 if FILE_SCALE else 8.0
OUT_JSON = os.path.join(os.path.dirname(C.OUT_JSON), "cal_mplet_file.json" if FILE_SCALE else "cal_mplet.json")

# Константы ширины линии
K_FWHM = 0.6061
P_FWHM = 0.6690
CONST_2SQRTPI = math.sqrt(2 * math.pi)

def fwhm(E):
    return K_FWHM * (E ** P_FWHM)

def sigma(E):
    return fwhm(E) / 2.354820045

def load(nuc, kind):
    """Читает CSV с линиями, возвращает список (E, I)."""
    path = os.path.join(LINES_DIR, f"{nuc}_{kind}.csv")
    lines = []
    if not os.path.exists(path):
        return lines
    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                e_str = row.get('energy', '').strip()
                i_str = row.get('intensity', '').strip()
                if not e_str or not i_str:
                    continue
                E = float(e_str)
                I = float(i_str)
                lines.append((E, I))
            except ValueError:
                continue
    return lines

# Определение групп
GROUPS_DEF = [
    {"name": "X-K", "lo": 55, "hi": 110, "comp": [("Bi K", "214pb", "x"), ("Pb K", "208tl", "x"), ("Po K", "214bi", "x"), ("Th K", "228ac", "x")], "cont": "quad"},
    {"name": "238", "lo": 212, "hi": 268, "comp": [("Pb-212", "212pb", "g"), ("Ra-224", "224ra", "g"), ("Pb-214", "214pb", "g")], "cont": "lin+step"},
    {"name": "352", "lo": 272, "hi": 378, "comp": [("Pb-214", "214pb", "g"), ("Ac-228", "228ac", "g"), ("Pb-212", "212pb", "g"), ("Tl-208", "208tl", "g")], "cont": "lin+step"},   # +Pb-212 300,1 и Tl-208 277,4 — в окне
    {"name": "609", "lo": 555, "hi": 645, "comp": [("Tl-208", "208tl", "g"), ("Bi-214", "214bi", "g")], "cont": "lin+step"},
    {"name": "911", "lo": 875, "hi": 1000, "comp": [("Ac-228", "228ac", "g"), ("Bi-214", "214bi", "g")], "cont": "lin+step"},
    # #DBG-1 (синтетика): окно <2 ПШПВ даёт вырождение «сдвиг ↔ наклон континуума» — окна расширены до ≈±2 ПШПВ
    {"name": "1461", "lo": 1300, "hi": 1640, "comp": [("K-40", "40k", "g"), ("Ac-228", "228ac", "g"), ("Bi-214", "214bi", "g")], "cont": "lin+step"},
    {"name": "1764", "lo": 1580, "hi": 1950, "comp": [("Bi-214", "214bi", "g")], "cont": "lin+step"},
    {"name": "2614", "lo": 2380, "hi": 2850, "comp": [("Tl-208", "208tl", "g"), ("Bi-214", "214bi", "g")], "cont": "lin+step"},
]

def get_lines_in_window(comp_name, nuc, kind, lo, hi):
    """Фильтрует линии компоненты по окну и интенсивности."""
    all_lines = load(nuc, kind)
    if not all_lines:
        return []
    
    # Предварительный отбор по энергии для оценки максимума интенсивности
    # Окно расширяется на 3 сигмы краев для учета размытия
    lo_ext = lo - 3 * sigma(lo)
    hi_ext = hi + 3 * sigma(hi)
    
    candidates = []
    for E, I in all_lines:
        if lo_ext <= E <= hi_ext and I > 0:
            candidates.append((E, I))
            
    if not candidates:
        return []
        
    max_I = max(I for _, I in candidates)
    threshold = 0.02 * max_I
    
    # Финальный отбор по строгим границам окна и порогу интенсивности
    final_lines = [(E, I) for E, I in candidates if lo <= E <= hi and I >= threshold]
    return final_lines

def build_matrix(x, delta, group_lines, cont_type):
    """
    Строит матрицу плана M.
    x: энергии каналов.
    delta: сдвиг шкалы.
    group_lines: список списков линий [(E, I), ...] для каждой компоненты.
    cont_type: 'quad' или 'lin+step'.
    Возвращает (M, comp_names, n_comp).
    """
    M_cols = []
    comp_names = []
    
    # Колонки компонент
    for lines in group_lines:
        if not lines:
            continue
        col = np.zeros_like(x)
        for E_tab, I_rel in lines:
            E_peak = E_tab + delta
            sig = sigma(E_tab) # Сигма зависит от табличной энергии (или пиковой? Обычно от пиковой, но разность мала. По ТЗ: sigma(E_k))
            # Формула из ТЗ: I_k * dE / (sigma * sqrt(2pi)) * exp(...)
            # Но в матрице плана для линейной модели y = M p, где p - амплитуда.
            # Обычно модель: Sum A_i * G(x; E_i+delta, sigma_i).
            # ТЗ говорит: "амплитуда — своя у каждой компоненты".
            # Колонка компоненты: сумма гауссов с весами I_k.
            # Нормировка: интеграл от суммы должен быть связан с амплитудой?
            # В ТЗ: "Колонка компоненты: Σ_k I_k·dE/(σ_k√(2π))·exp..."
            # Это выглядит как вклад в интенсивность на канал.
            # dE = np.gradient(x). Но gradient возвращает массив той же длины.
            # Для точности лучше использовать среднее расстояние или просто 1, если шкала линейна.
            # В ТЗ явно указано dE.
            
            # Вычисляем гауссиану для этой линии
            gauss = np.exp(-0.5 * ((x - E_peak) / sig)**2)
            col += I_rel * gauss / (sig * CONST_2SQRTPI)
        
        M_cols.append(col)
        comp_names.append(f"Comp_{len(comp_names)}") # Имена будут подставлены позже из контекста

    # Колонки континуума
    x0 = (x[0] + x[-1]) / 2
    span = x[-1] - x[0]
    if span == 0: span = 1e-9
    t = (x - x0) / span
    
    if cont_type == "quad":
        M_cols.append(np.ones_like(x))
        M_cols.append(t)
        M_cols.append(t**2)
    elif cont_type == "lin+step":
        M_cols.append(np.ones_like(x))
        M_cols.append(t)
        # Ступенька: E_s - средняя энергия линий группы + delta
        # Нужно вычислить среднюю энергию всех линий в группе
        all_E = []
        for lines in group_lines:
            for E, _ in lines:
                all_E.append(E)
        if all_E:
            E_mean_group = np.mean(all_E)
            # sigma_s = max sigma_k
            sigs = [sigma(E) for E in all_E]
            sig_s = max(sigs) if sigs else 1.0
            step_col = 0.5 * erfc((x - (E_mean_group + delta)) / (sig_s * math.sqrt(2)))
            M_cols.append(step_col)
            
    M = np.column_stack(M_cols)
    return M, comp_names

def solve_for_delta(x, y, group_lines, cont_type, delta_val):
    """Решает задачу наименьших квадратов для фиксированного delta."""
    dE = np.gradient(x)
    
    # Строим матрицу
    M_raw, _ = build_matrix(x, delta_val, group_lines, cont_type)
    
    # Умножаем колонки компонент на dE? 
    # В ТЗ: "Колонка компоненты: Σ_k I_k·dE/(σ_k√(2π))·exp..."
    # Это значит, что сама колонка уже включает dE.
    # Но build_matrix выше не умножал на dE явно в цикле, а добавлял гауссианы.
    # Нужно поправить: каждая колонка компоненты должна быть умножена на dE?
    # Или dE входит в нормировку? 
    # "Колонка компоненты: ... I_k * dE / ..." -> Да, dE часть формулы колонки.
    
    # Пересчитаем M с учетом dE для колонок компонент
    n_comp = len(group_lines)
    M = np.zeros_like(M_raw)
    
    # Копируем и масштабируем колонки компонент
    for i in range(n_comp):
        M[:, i] = M_raw[:, i] * dE
        
    # Континуум (последние колонки) не масштабируются dE в формуле ТЗ?
    # В ТЗ про континуум сказано отдельно: "колонки 1, t, t^2". Без dE.
    n_cont = M.shape[1] - n_comp
    if n_cont > 0:
        M[:, n_comp:] = M_raw[:, n_comp:]
        
    # Веса
    sy = np.sqrt(np.maximum(y, 1))
    Mw = M / sy[:, None]
    yw = y / sy
    
    # Границы
    lb = -np.inf * np.ones(M.shape[1])   # #DBG-1: континуум свободен (было 0 у всех — наклон фона изображался сдвигом)
    ub = np.inf * np.ones(M.shape[1])
    
    # Компоненты >= 0
    for i in range(n_comp):
        lb[i] = 0
        
    # Ступенька (если есть) >= 0. Она последняя в lin+step.
    if cont_type == "lin+step":
        lb[-1] = 0
        
    # Решение
    try:
        res = lsq_linear(Mw, yw, bounds=(lb, ub), method="trf", tol=1e-10, max_iter=400)
        p = res.x
        cost = 2.0 * res.cost   # res.cost у lsq_linear = ½·Σ(yw − Mw·p)²
        
        # Chi2 nu
        n_data = len(y)
        n_params = M.shape[1]
        nu = n_data - n_params
        if nu <= 0: return None, None, None, None
        
        chi2_nu = cost / nu
        
        # Ковариация через SVD
        U, s, Vt = np.linalg.svd(Mw, full_matrices=False)
        # Фильтр малых сингулярных значений для стабильности
        s_inv = 1.0 / s
        cov = (Vt.T * s_inv**2) @ Vt * chi2_nu
        
        # Значимости компонент
        sigmas_p = np.sqrt(np.diag(cov))
        significances = []
        for i in range(n_comp):
            if sigmas_p[i] > 0:
                significances.append(p[i] / sigmas_p[i])
            else:
                significances.append(0.0)
                
        return p, chi2_nu, significances, cov
        
    except Exception as e:
        return None, None, None, None

def main():
    tags = ["bg", "bgw", "kcl", "sample"]
    results = {} # Для JSON и сводки
    
    # Структура для хранения результатов по группам и тегам
    all_results = {tag: {} for tag in tags}
    
    print("Начало калибровки GS2020 по форме групп...")
    
    for tag in tags:
        d = C.read(tag)
        counts = d["counts"]
        x_all = C.s.channel_to_energy(np.arange(len(counts)), d["coeffs"]) if FILE_SCALE else C.energy_axis(tag, len(counts))
        
        for grp in GROUPS_DEF:
            name = grp["name"]
            lo, hi = grp["lo"], grp["hi"]
            cont_type = grp["cont"]
            
            # 1. Сбор линий для каждой компоненты в окне
            group_lines_list = [] # Список списков (E, I)
            comp_names_real = []
            valid_group = True
            
            for cname, nuc, kind in grp["comp"]:
                lines = get_lines_in_window(cname, nuc, kind, lo, hi)
                if not lines:
                    print(f"  {tag} {name}: компонента {cname}: нет линий в окне")
                    # Если у компоненты нет линий, она не вносит вклад. 
                    # Но ТЗ говорит "Компонента без линий в окне — пропускается".
                    # Значит, мы просто не добавляем её в модель?
                    # Да, иначе колонка будет нулевой или бессмысленной.
                    continue
                group_lines_list.append(lines)
                comp_names_real.append(cname)
            
            if not group_lines_list:
                print(f"  {tag} {name}: нет ни одной компоненты с линиями")
                all_results[tag][name] = {"weak": True, "edge": False, "delta_keV": 0, "sigma_keV": 0, "chi2_nu": 0, "E_mean": (lo+hi)/2, "components": {}}
                continue
                
            # 2. Определение окна данных
            mask = (x_all >= lo) & (x_all <= hi)
            x_win = x_all[mask]
            y_win = counts[mask]
            
            if len(x_win) < 10:
                 print(f"  {tag} {name}: слишком мало точек в окне")
                 all_results[tag][name] = {"weak": True, "edge": False, "delta_keV": 0, "sigma_keV": 0, "chi2_nu": 0, "E_mean": (lo+hi)/2, "components": {}}
                 continue

            # 3. Скан delta
            deltas = np.arange(-SCAN, SCAN + 0.05, 0.05)
            chi2_vals = []
            best_res = None
            
            for delta in deltas:
                p, chi2_nu, sigs, cov = solve_for_delta(x_win, y_win, group_lines_list, cont_type, delta)
                if p is not None:
                    chi2_vals.append(chi2_nu)
                    # Сохраняем результат для лучшего delta (пока просто последний валидный, потом найдем min)
                    # Чтобы не хранить все, будем искать минимум на лету или после
                else:
                    chi2_vals.append(np.inf)
            
            if not any(c < np.inf for c in chi2_vals):
                 print(f"  {tag} {name}: решение не найдено")
                 all_results[tag][name] = {"weak": True, "edge": False, "delta_keV": 0, "sigma_keV": 0, "chi2_nu": 0, "E_mean": (lo+hi)/2, "components": {}}
                 continue

            chi2_arr = np.array(chi2_vals)
            min_idx = np.argmin(chi2_arr)
            delta_best = deltas[min_idx]
            chi2_min = chi2_arr[min_idx]
            
            # Повторное решение для лучшего delta, чтобы получить параметры и ковариацию
            p_best, chi2_nu_best, sigs_best, cov_best = solve_for_delta(x_win, y_win, group_lines_list, cont_type, delta_best)
            
            if p_best is None:
                 print(f"  {tag} {name}: ошибка при финальном расчете")
                 all_results[tag][name] = {"weak": True, "edge": False, "delta_keV": 0, "sigma_keV": 0, "chi2_nu": 0, "E_mean": (lo+hi)/2, "components": {}}
                 continue

            # 4. Оценка погрешности delta
            # Точки, где chi2 <= chi2_min + chi2_min/nu
            # nu нужно взять из последнего решения
            n_data = len(y_win)
            n_params = len(p_best)
            nu = n_data - n_params
            
            threshold = chi2_min + chi2_min / nu if nu > 0 else chi2_min + 1e-6
            
            indices_in_range = np.where(chi2_arr <= threshold)[0]
            
            if len(indices_in_range) == 0:
                # Если минимум острый и не попадает в диапазон (редко), берем шаг
                sigma_delta = 0.05 / 2
            else:
                delta_range = deltas[indices_in_range]
                width = delta_range[-1] - delta_range[0]
                sigma_delta = max(width / 2, 0.05) # Не меньше шага
            
            # Проверка на КРАЙ
            is_edge = abs(delta_best) >= SCAN - 0.1
            
            # Проверка на СЛАБО
            # Значимость компоненты p_i / sqrt(cov_ii). 
            # sigs_best уже содержит эти значения для компонент.
            # "Если ни одна компонента не значима (все < 5)"
            # значимость СУММЫ компонент: совпадающие линии (238,6/241,0/242,0) по отдельности незначимы, вместе — сильный пик
            nc = len(comp_names_real); vt = float(cov_best[:nc, :nc].sum()) if cov_best is not None else 0.0
            # годна, если значима ЛЮБАЯ компонента или сумма (дисперсию суммы раздувает прижатая к нулю слабая колонка)
            is_weak = not (max(sigs_best) >= 5 or (vt > 0 and float(p_best[:nc].sum()) / math.sqrt(vt) >= 5))
            
            # Формирование строки вывода
            comp_str_parts = []
            comp_json_data = {}
            for i, cname in enumerate(comp_names_real):
                sig_val = sigs_best[i] if i < len(sigs_best) else 0
                amp_val = p_best[i] if i < len(p_best) else 0
                # Погрешность амплитуды из ковариации
                err_amp = math.sqrt(cov_best[i, i]) if cov_best is not None and i < cov_best.shape[0] else 0
                
                comp_str_parts.append(f"{cname} {sig_val:.1f}")
                comp_json_data[cname] = {"amp": float(amp_val), "sig": float(sig_val)} # sig здесь как значимость по ТЗ вывода? Или sigma? В JSON просили amp, sig. Обычно sig=sigma. Но в тексте "значимость". Давайте запишем значимость как 'sig' для соответствия запросу "p_i/sqrt(cov_ii)" в описании значимости, но в JSON поле называется sig. Уточню: в JSON примере {"amp","sig"}. В тексте вывода "a.a sigma". Вероятно, имеется в виду погрешность или значимость. 
                # Перечитаю ТЗ: "компоненты: <имя> a.a σ; ... (для компоненты — значимость p_i/σ_i с одним знаком)".
                # В JSON: {"amp","sig"}. Обычно sig=sigma. Но контекст вывода подразумевает значимость. 
                # Однако, в научном контексте sig часто sigma. 
                # Посмотрю на пример вывода: "a.a σ". Это может быть "значение погрешность". 
                # Но дальше сказано "значимость ... с одним знаком".
                # Давайте запишем в JSON 'sig' как значимость (S/N), так как это ключевой метрический показатель качества подгонки компоненты.
                
            comp_str = ", ".join(comp_str_parts)
            
            flags = ""
            if is_edge: flags += " [КРАЙ]"
            if is_weak: flags += " [СЛАБО]"
            
            print(f"{tag} {name} окно {lo}–{hi}: δ = {delta_best:+.2f} ± {sigma_delta:.2f} кэВ, χ²/ν {chi2_nu_best:.2f}, компоненты: {comp_str}{flags}")
            
            # Сохранение в структуру
            all_results[tag][name] = {
                "lo": lo,
                "hi": hi,
                "delta_keV": float(delta_best),
                "sigma_keV": float(sigma_delta),
                "chi2_nu": float(chi2_nu_best),
                "edge": bool(is_edge),
                "weak": bool(is_weak),
                "E_mean": float(np.mean([np.mean([l[0] for l in lines]) for lines in group_lines_list])), # Средняя энергия линий группы
                "components": comp_json_data,
                # центр узла: табличная энергия, взвешенная по вкладу линий в модель (амплитуда × интенсивность)
                "E_node": float(sum(max(p_best[i], 0) * I * E for i, ls in enumerate(group_lines_list) for E, I in ls)
                                / max(sum(max(p_best[i], 0) * I for i, ls in enumerate(group_lines_list) for E, I in ls), 1e-30))
            }

    # Сводка по группам
    print("\n--- Сводка по сдвигам ---")
    interpretation_issues = []
    
    for grp in GROUPS_DEF:
        name = grp["name"]
        deltas_grp = {}
        sigmas_grp = {}
        
        for tag in tags:
            res = all_results[tag].get(name)
            if res and not res.get("weak", True):
                deltas_grp[tag] = res["delta_keV"]
                sigmas_grp[tag] = res["sigma_keV"]
                
                # Проверка на проблемы для интерпретации
                if res["chi2_nu"] > 3:
                    interpretation_issues.append(f"{tag} {name}: χ²/ν = {res['chi2_nu']:.2f}")
                if abs(res["delta_keV"]) > 0.25 * fwhm((grp["lo"]+grp["hi"])/2):
                    interpretation_issues.append(f"{tag} {name}: |δ| велико")
                if res["edge"]:
                    interpretation_issues.append(f"{tag} {name}: КРАЙ")
        
        # Вывод строки сводки
        parts = []
        for tag in tags:
            if tag in deltas_grp:
                parts.append(f"{tag} δ={deltas_grp[tag]:+.2f}")
            else:
                parts.append(f"{tag} -")
        
        print(f"{name}: {'; '.join(parts)}")
        
        # Разности kcl-bgw, sample-bgw
        if "bgw" in deltas_grp:
            if "kcl" in deltas_grp:
                diff_k = deltas_grp["kcl"] - deltas_grp["bgw"]
                err_k = math.sqrt(sigmas_grp["kcl"]**2 + sigmas_grp["bgw"]**2)
                print(f"  kcl − bgw = {diff_k:+.2f} ± {err_k:.2f}")
            if "sample" in deltas_grp:
                diff_s = deltas_grp["sample"] - deltas_grp["bgw"]
                err_s = math.sqrt(sigmas_grp["sample"]**2 + sigmas_grp["bgw"]**2)
                print(f"  sample − bgw = {diff_s:+.2f} ± {err_s:.2f}")

    # Блок ТРЕБУЕТ ТОЛКОВАНИЯ
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    if interpretation_issues:
        for issue in set(interpretation_issues): # Уникальные проблемы
            print(f"  - {issue}")
    else:
        print("  нет")

    # Запись JSON
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, ensure_ascii=False, indent=1)
    print(f"\nРезультаты записаны в {OUT_JSON}")

if __name__ == "__main__":
    main()
