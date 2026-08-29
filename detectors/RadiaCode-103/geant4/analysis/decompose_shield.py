# -*- coding: utf-8 -*-
"""
Понуклидное разложение спектра в свинцовом домике при ИСХОДНЫХ активностях,
взятых из открытого фона, плюс энергетическая зависимость ослабления по каждому
нуклиду. Ничего не подгоняется: активности переносятся как есть, а ослабление
берётся из отношения МК-шаблонов «с домиком / без домика», посчитанных одной и
той же программой в одной и той же геометрии.

Выход: текстовые таблицы в stdout и PNG с двумя панелями.

Коэффициент пропускания T(E) — не подгоночная кривая и не аналитическая
формула ослабления, а прямое отношение двух расчётов, поэтому в нём уже сидят
рассеяние в свинце, вклад через открытый верх и геометрия полости.
"""

import os, sys, io, contextlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc
import fit_physical_chains as fpc
import predict_shield as ps

matplotlib.rcParams["font.size"] = 9
matplotlib.rcParams["figure.dpi"] = 130
matplotlib.rcParams["axes.grid"] = True
matplotlib.rcParams["grid.alpha"] = 0.3
OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PNG = os.path.join(OUT_DIR, "RC103_shield_decomposition.png")

# Палитра — проверенные слоты категориальной темы скилла dataviz, взяты в
# фиксированном порядке fpc.NUCS. Прогон validate_palette.js --mode light:
# все жёсткие проверки пройдены (худшая пара при дальтонизме ΔE 9,1 при пороге 8;
# при нормальном зрении 19,6 при пороге 15). Прежняя подобранная на глаз палитра
# проверку провалила по четырём пунктам из пяти.
# Космика намеренно НЕ берёт категориальный слот: это не нуклид, а другая
# физика, и нейтральный серый отделяет её от рядов.
COLORS = {
    "K40":   "#2a78d6",
    "Ra226": "#eb6834",
    "Pb214": "#1baf7a",
    "Bi214": "#eda100",
    "Pb212": "#e87ba4",
    "Ac228": "#008300",
    "Bi212": "#4a3aa7",
    "Tl208": "#e34948",
    "mu":    "#8c8c8c",
}
RU = {"K40": "K-40", "Ra226": "Ra-226", "Pb214": "Pb-214", "Bi214": "Bi-214",
      "Pb212": "Pb-212", "Ac228": "Ac-228", "Bi212": "Bi-212", "Tl208": "Tl-208",
      "mu": "космика"}

def fit_open():
    """Подгонка открытого фона для получения исходных активностей."""
    cnt_o, e_o, live_o, _ = ps.read_meas(ftc.MEAS_NAME, ftc.CAL_ROOM)
    cps_o, var_o, hm_o = ps.load_cps(ps.FMT_OPEN, ftc.MUON_CSV)
    
    if not hm_o:
        sys.exit("нет мюонного шаблона для открытого фона")
        
    names, A_o, V_o = ps.columns_on_grid(cps_o, var_o, e_o, hm_o)
    sel = (e_o >= ftc.E_LO) & (e_o < ftc.E_HI)
    
    # Подгонка в отсчетах: матрица * live_time, измерение - сырые счеты
    # ftc.fit возвращает кортеж (amp, sd, pred, chi2ndf, shape)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        # live_o — СКАЛЯР (живое время измерения), индексировать его нельзя.
        # Сигнатура позиционная: fit(A_counts, meas_counts, weights, names,
        # title, note, var_counts=...).
        result = ftc.fit(
            A_o[sel] * live_o,
            cnt_o[sel],
            1.0 / np.sqrt(np.maximum(cnt_o[sel], 1.0)),
            names, "", "",
            var_counts=V_o[sel] * live_o * live_o,
        )

    amp, sd, pred, chi2ndf, shape = result
    
    return (dict(zip(names, amp)), names)

def nuclide_activities(amp_map):
    """Разворачивает цепочечные амплитуды в понуклидные активности."""
    required_keys = ["K40", "Ra226_chain", "Th232_chain", "mu"]
    for key in required_keys:
        if key not in amp_map:
            print(f"Отсутствует ключ '{key}'. Имеющиеся ключи: {list(amp_map.keys())}")
            sys.exit(1)

    A_K40   = amp_map["K40"]
    A_Ra226 = amp_map["Ra226_chain"]
    # Вклад Pb-214 и Bi-214 зависит от эманации радона (F_RN)
    A_Pb214 = A_Bi214 = A_Ra226 * (1 - ps.F_RN)
    
    A_Ac228 = amp_map["Th232_chain"]
    # Вклад Pb-212 и Bi-212 зависит от эманации торона (F_TN) и отношения цепочек
    A_Pb212 = A_Bi212 = A_Ac228 * ps.R_TH * (1 - ps.F_TN)
    
    # Tl-208 — ветвление распада Bi-212
    A_Tl208 = fpc.TL208_BRANCH * A_Bi212
    
    A_mu    = amp_map["mu"]
    
    acts = {
        "K40": A_K40,
        "Ra226": A_Ra226,
        "Pb214": A_Pb214,
        "Bi214": A_Bi214,
        "Ac228": A_Ac228,
        "Pb212": A_Pb212,
        "Bi212": A_Bi212,
        "Tl208": A_Tl208,
        "mu": A_mu
    }
    
    # Возвращаем в порядке fpc.NUCS + ["mu"]
    ordered_acts = {}
    for nuc in fpc.NUCS:
        ordered_acts[nuc] = acts[nuc]
    ordered_acts["mu"] = acts["mu"]
    
    return ordered_acts

def per_nuclide_spectra(fmt, muon_csv, e_grid, acts):
    """Собирает вклад каждого нуклида отдельно в имп/с на сетке измерения."""
    spectra = {}
    
    # Обработка нуклидов
    for nuc in fpc.NUCS:
        template_path = os.path.join(ftc.TEMPLATE_DIR, fmt % nuc)
        if not os.path.exists(template_path):
            print(f"Файл шаблона не найден: {template_path}")
            sys.exit(1)
            
        _meta, arr, _cnt = ftc.read_template(template_path)  # ридер отдаёт КОРТЕЖ
        # Сворачивание спектра в геометрию детектора "103"
        col = ftc.rcspec.fold(arr, "103", tail_T=ps.TAIL_T)
        # Перевод на сетку измерения и масштабирование активностью
        spec = ps.to_meas_grid(col, e_grid) * acts[nuc]
        spectra[nuc] = spec
        
    # Обработка мюонов
    if not muon_csv or not os.path.exists(muon_csv):
        print(f"Файл мюонного шаблона не найден или пуст: {muon_csv}")
        sys.exit(1)
        
    _meta_mu, arr_mu, _cnt_mu = ftc.read_template(muon_csv)
    col_mu = ftc.rcspec.fold(arr_mu, "103", tail_T=ps.TAIL_T)
    spec_mu = ps.to_meas_grid(col_mu, e_grid) * acts["mu"]
    spectra["mu"] = spec_mu
    
    return spectra

def transmission(e_grid):
    """Коэффициент пропускания по энергии: отношение шаблонов с домиком и без."""
    T_dict = {}
    
    for nuc in fpc.NUCS:
        # Шаблон без домика (открытый фон)
        path_open = os.path.join(ftc.TEMPLATE_DIR, ps.FMT_OPEN % nuc)
        # Шаблон с домиком
        path_shield = os.path.join(ftc.TEMPLATE_DIR, ps.FMT_SHIELD % nuc)
        
        if not os.path.exists(path_open) or not os.path.exists(path_shield):
            print(f"Отсутствуют шаблоны для {nuc}: {path_open}, {path_shield}")
            sys.exit(1)
            
        _m1, arr_open, _c1 = ftc.read_template(path_open)
        _m2, arr_shield, _c2 = ftc.read_template(path_shield)
        
        # Сворачивание в геометрию детектора
        col_open = ftc.rcspec.fold(arr_open, "103", tail_T=ps.TAIL_T)
        col_shield = ftc.rcspec.fold(arr_shield, "103", tail_T=ps.TAIL_T)
        
        # Перевод на сетку измерения
        spec_open = ps.to_meas_grid(col_open, e_grid)
        spec_shield = ps.to_meas_grid(col_shield, e_grid)
        
        # Сглаживание скользящим средним (25 каналов) для подавления шума перед делением
        kernel = np.ones(25) / 25.0
        spec_open_smooth = np.convolve(spec_open, kernel, mode="same")
        spec_shield_smooth = np.convolve(spec_shield, kernel, mode="same")
        
        # Отношение только там, где знаменатель значим
        threshold = 1e-3 * np.max(spec_open_smooth)
        T = np.full_like(e_grid, np.nan)
        mask = spec_open_smooth > threshold
        T[mask] = spec_shield_smooth[mask] / spec_open_smooth[mask]
        
        T_dict[nuc] = T
        
    return T_dict

def report(acts, spec_shield, meas_shield, live_s, e_s):
    """Печатает три таблицы: активности, вклад нуклидов, вклад по полосам."""
    
    # 1. Исходные активности (модельная оценка)
    print("\n" + "="*60)
    print("МОДЕЛЬНАЯ ОЦЕНКА ИСХОДНЫХ АКТИВНОСТЕЙ (из открытого фона)")
    print("="*60)
    print(f"{'Нуклид':<10} {'Активность':>12} {'Единица'}")
    print("-"*40)
    for nuc in fpc.NUCS:
        unit = "Бк/кг"
        print(f"{RU[nuc]:<10} {acts[nuc]:>12.4e} {unit}")
    print(f"{RU['mu']:<10} {acts['mu']:>12.4e} с^-1")
    
    # 2. Вклад нуклидов в спектр в домике (полный интеграл)
    print("\n" + "="*60)
    print("ВКЛАД НУКЛИДОВ В СПЕКТР В ДОМИКЕ (интегрально)")
    print("="*60)
    
    # Массивы приходят УЖЕ отобранными по E_LO..E_HI (отбор сделан в main),
    # повторно применять маску нельзя — длины не совпадут.
    # Величины заданы как имп/с НА КАНАЛ, а не как плотность на кэВ, поэтому
    # полный счёт есть простая сумма по каналам; домножать на ширину канала
    # было бы ошибкой — счёт вырос бы в среднюю ширину (~2,7) раз.
    delta_e = 1.0
        
    contributions = {}
    total_model = 0.0
    
    for nuc in fpc.NUCS + ["mu"]:
        spec_nuc = spec_shield[nuc]
        # Интеграл по спектру (сумма cps * ширина канала)
        integral = np.sum(spec_nuc) * delta_e
        contributions[nuc] = integral
        total_model += integral
        
    # Сортировка по убыванию вклада
    sorted_nucs = sorted(contributions.keys(), key=lambda k: contributions[k], reverse=True)
    
    print(f"{'Нуклид':<10} {'Счет, имп/с':>12} {'Доля, %':>10}")
    print("-"*40)
    for nuc in sorted_nucs:
        frac = contributions[nuc] / total_model * 100 if total_model > 0 else 0
        print(f"{RU[nuc]:<10} {contributions[nuc]:>12.5f} {frac:>10.2f}")
        
    # Сумма модели и измерение
    meas_integral = np.sum(meas_shield)
    print("-"*40)
    print(f"{'Сумма':<10} {total_model:>12.5f} {'100.00':>10}")
    print(f"{'Измерено':<10} {meas_integral:>12.5f} {'-':>10}")

    # 3. Вклад по полосам
    print("\n" + "="*60)
    print("ВКЛАД ПО ЭНЕРГЕТИЧЕСКИМ ПОЛОСАМ")
    print("="*60)
    
    # Заголовки
    header = f"{'Полоса':<15}"
    for nuc in fpc.NUCS + ["mu"]:
        header += f" {RU[nuc]:>6}"
    header += f" {'Изм,имп/с':>10} {'Мод,имп/с':>10}"
    print(header)
    print("-"*len(header))
    
    # ftc.BANDS — СПИСОК пар (lo, hi), не словарь; имя полосы строим сами.
    for (e_lo, e_hi) in ftc.BANDS:
        band_name = "%d-%d" % (e_lo, e_hi)
        sel_band = (e_s >= e_lo) & (e_s < e_hi)
        if not np.any(sel_band):
            continue
            
        # Интегралы в полосе
        band_integrals = {}
        total_band_model = 0.0
        
        for nuc in fpc.NUCS + ["mu"]:
            spec_nuc = spec_shield[nuc][sel_band]
            integral = np.sum(spec_nuc)
            band_integrals[nuc] = integral
            total_band_model += integral
            
        meas_band = np.sum(meas_shield[sel_band])
        
        row = f"{band_name:<15}"
        for nuc in fpc.NUCS + ["mu"]:
            if total_band_model > 1e-9:
                frac = band_integrals[nuc] / total_band_model * 100
                row += f" {frac:>6.1f}"
            else:
                row += f" {'-':>6}"
                
        row += f" {meas_band:>10.5f} {total_band_model:>10.5f}"
        print(row)

def draw(acts, spec_shield, meas_shield, e_s, T_dict):
    """Рисует две панели: разложение спектра и пропускание."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 10), 
                                   gridspec_kw={"height_ratios": [3, 2]})
    
    sel = (e_s >= ftc.E_LO) & (e_s < ftc.E_HI)
    e_plot = e_s[sel]
    meas_plot = meas_shield[sel]
    
    # --- Верхняя панель: Разложение спектра ---
    
    # Измерение
    ax1.step(e_plot, meas_plot, where="mid", color="#111111", lw=1.0, label="измерено, домик")
    
    # Стопка вкладов нуклидов
    # Сортируем нуклиды по полному вкладу (возрастание), чтобы мелкие были снизу
    contributions = {}
    for nuc in fpc.NUCS + ["mu"]:
        spec_nuc = spec_shield[nuc]
        contributions[nuc] = np.sum(spec_nuc)
        
    sorted_nucs = sorted(contributions.keys(), key=lambda k: contributions[k])
    
    # Строим стопку вручную через fill_between для контроля порядка и цветов
    y_stack = np.zeros_like(e_plot)
    for nuc in sorted_nucs:
        spec_nuc = spec_shield[nuc]
        ax1.fill_between(e_plot, y_stack, y_stack + spec_nuc, 
                         color=COLORS[nuc], linewidth=0, label=RU[nuc])
        y_stack += spec_nuc
        
    # Полная модель сверху
    ax1.plot(e_plot, y_stack, color="#d62728", lw=1.4, label="модель, сумма")
    
    # Настройки осей
    min_meas = np.min(meas_plot[meas_plot > 0])
    max_meas = np.max(meas_plot)
    ax1.set_ylim(0.5 * min_meas, 2 * max_meas)
    ax1.set_yscale("log")
    ax1.set_xlim(ftc.E_LO, ftc.E_HI)
    ax1.set_ylabel("скорость счёта, имп/с на канал")
    ax1.set_title("RadiaCode-103: разложение спектра в свинцовом домике при исходных активностях")
    ax1.set_xlabel("") # Пустая подпись X для верхней панели
    
    # Легенда вне поля
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=5, frameon=False)
    
    # --- Нижняя панель: Пропускание T(E) ---
    
    for nuc in fpc.NUCS + ["mu"]:
        if nuc in T_dict:
            T_vals = T_dict[nuc][sel]
            ax2.plot(e_plot, T_vals, color=COLORS[nuc], lw=1.2, label=RU[nuc])
            
    # Опорная линия 1.0
    ax2.axhline(1.0, color="#666666", lw=0.8, ls="--")
    
    # Настройки осей
    ax2.set_yscale("log")
    ax2.set_xlim(ftc.E_LO, ftc.E_HI)
    ax2.set_ylabel("пропускание домика T(E) = с домиком / без домика")
    ax2.set_xlabel("энергия, кэВ")
    
    # Форматирование Y-оси обычными числами
    ax2.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    
    # Легенда справа вне поля
    ax2.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), frameon=False)
    
    fig.tight_layout()
    fig.savefig(OUT_PNG)
    print(f"\nГрафик сохранен: {OUT_PNG}")

def main():
    # 1. Подгонка открытого фона
    amp_map, names = fit_open()
    
    # 2. Расчет понуклидных активностей
    acts = nuclide_activities(amp_map)
    
    # 3. Чтение данных для домика
    cnt_s, e_s, live_s, _ = ps.read_meas(ps.MEAS_SHIELD, None)
    sel_s = (e_s >= ftc.E_LO) & (e_s < ftc.E_HI)
    
    # Проверка мюонного шаблона для домика
    if not ps.MUON_SHIELD_CSV or not os.path.exists(ps.MUON_SHIELD_CSV):
        sys.exit("Без космической компоненты разложение в домике бессмысленно: файл мюонов не найден.")
        
    # 4. Расчет спектров по нуклидам в домике
    spec_shield = per_nuclide_spectra(ps.FMT_SHIELD, ps.MUON_SHIELD_CSV, e_s[sel_s], acts)
    
    # Измеренный спектр в имп/с, ОТОБРАННЫЙ по тому же интервалу, что и модель:
    # дальше и report, и draw работают на одной сетке e_s[sel_s].
    meas_shield = (cnt_s / live_s)[sel_s]
    
    # 5. Расчет пропускания (нужна полная сетка или та же sel_s? 
    # Функция transmission принимает e_grid, лучше передать отфильтрованную для соответствия графикам)
    T_dict = transmission(e_s[sel_s])
    
    # 6. Отчет и график
    report(acts, spec_shield, meas_shield, live_s, e_s[sel_s])
    draw(acts, spec_shield, meas_shield, e_s[sel_s], T_dict)

if __name__ == "__main__":
    main()
