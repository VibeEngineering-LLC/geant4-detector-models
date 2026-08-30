# -*- coding: utf-8 -*-
"""
predict_shield.py

Проверка модели ослабления гамма-фона свинцовым домиком.
Активности нуклидов берутся из подгонки фона в открытой комнате и используются без изменений для предсказания спектра в домике.
Свободных параметров в предсказании нет — это проверка модели, а не подгонка.
Дополнительно выполняется свободная подгонка спектра в домике для диагностики расхождений между известными активностями и теми, что "видит" детектор за защитой.
"""

import os, sys, io, contextlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc
import fit_physical_chains as fpc

# Константы модуля
TAIL_T = float(os.environ.get("G4M_TAIL_T", "0.75"))
# F_RN, R_TH, F_TN — ДЕЙСТВУЮЩАЯ постановка контура, зафиксированная по открытому фону.
# Подбирать их здесь запрещено.
F_RN = float(os.environ.get("G4M_FIX_FRN", "0.10"))
R_TH = float(os.environ.get("G4M_FIX_RTH", "1.00"))
F_TN = float(os.environ.get("G4M_FIX_FTN", "0.00"))

FMT_OPEN = "rc103_field_room_%s.csv"
FMT_SHIELD = os.environ.get("G4M_SHIELD_FMT", "rc103_field_room_shield_%s.csv")
# Имя записано escape-последовательностями намеренно: файл лежит у оператора и его имя по-русски,
# а исходник обязан оставаться ASCII в этой строке.
MEAS_SHIELD = os.environ.get("G4M_MEAS_SHIELD", "\u0424\u043e\u043d \u0434\u043e\u043c\u0438\u043a 23 \u0434\u043d\u044f.xml")
MUON_SHIELD_CSV = os.environ.get("G4M_MUON_SHIELD_CSV", "")

BANDS = ftc.BANDS


def _check_posture(meta, tag, posture):
    """Посадка (P-005) — только у шаблонов С ДОМИКОМ. Смешать в одном
    предсказании файлы с РАЗНОЙ посадкой — молчаливый W-035: числа сойдутся
    по форме файла и разъедутся по физике. FATAL, не предупреждение."""
    if meta.get("shield", 0.0) != 1.0:
        return posture
    p = (meta.get("stand_mm"), meta.get("screen_up"))
    if posture is not None and p != posture:
        print(f"FATAL: посадка '{tag}' {p} не совпадает с прежней {posture} "
              f"— шаблоны для разных постановок смешивать нельзя (P-005/W-035).")
        sys.exit(1)
    return p if posture is None else posture


def load_cps(fmt, muon_csv):
    """
    Загрузка шаблонов для заданного формата имени файла.

    Аргументы:
        fmt: строка формата имени шаблона (с %s для нуклида).
        muon_csv: путь к мюонному шаблону (пустая строка = мюонов нет).

    Возвращает:
        (cps, var, has_muon, posture)

    posture — (stand_mm, screen_up) из шапки CSV, если шаблон построен с
    домиком (shield=1), иначе None. P-005 (30.08.2026): посадка прибора в
    полости — параметр, а не константа кода; смешать в одном предсказании
    шаблоны, посчитанные для РАЗНЫХ посадок, — тот же класс дефекта W-035
    («сравниваемые постановки обязаны различаться ровно одной вещью»), только
    молчаливый — числа сойдутся по форме файла и разъедутся по физике.
    FATAL здесь, а не предупреждение: посадка входит в состав каждого числа
    итоговой таблицы, ошибку нельзя откатить постфактум.
    """
    cps = {}
    var = {}
    posture = None

    # Загрузка нуклидов
    for nuc in fpc.NUCS:
        path = os.path.join(ftc.TEMPLATE_DIR, fmt % nuc)
        if not os.path.exists(path):
            print(f"Файл шаблона не найден: {path}")
            sys.exit(1)
            
        meta, arr, cnt_mc = ftc.read_template(path)
        # Сворачиваем спектр с учетом хвоста
        cps[nuc] = ftc.rcspec.fold(arr, "103", tail_T=TAIL_T)
        # Вычисляем дисперсию (ftc.template_variance делает свой fold без хвоста, это принятое приближение)
        var[nuc] = ftc.template_variance(cnt_mc, float(meta.get("t_run_s", 0.0)))
        posture = _check_posture(meta, nuc, posture)

    has_muon = False
    
    # Загрузка мюонного фона, если указан и файл существует
    if muon_csv and os.path.exists(muon_csv):
        meta_mu, arr_mu, cnt_mu = ftc.read_template(muon_csv)
        cps["mu"] = ftc.rcspec.fold(arr_mu, "103", tail_T=TAIL_T)
        var["mu"] = ftc.template_variance(cnt_mu, float(meta_mu.get("n_events", 0.0)))
        has_muon = True
        # Справочная сверка с потоком PDG обязана считаться для ТОГО диска,
        # которым посчитан шаблон: в ftc зашито значение для R=300 мм, а здесь
        # радиус другой. Само число лежит в метаданных шаблона.
        pdg = float(meta_mu.get("pdg_expected_per_s", 0.0))
        if pdg > 0:
            ftc.MUON_PDG_PER_S = pdg
        posture = _check_posture(meta_mu, "mu", posture)

    return cps, var, has_muon, posture


def to_meas_grid(col, e_meas):
    """Перевод модели на сетку измерения."""
    return ftc.fl.rebin_model_to_meas(np.arange(len(col)) + 0.5, col, e_meas)


def read_meas(name, cal):
    """
    Чтение файла измерения.
    
    Аргументы:
        name: имя файла в MEAS_DIR.
        cal: коэффициенты калибровки (None для приборной калибровки из файла).
        
    Возвращает:
        (cnt, e, live, cal_tuple)
    """
    smp = ftc.read_rcxml.read(os.path.join(ftc.MEAS_DIR, name))[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    
    if cal is not None:
        e = sum(c * ch**i for i, c in enumerate(cal))
        cal_tuple = tuple(cal)
    else:
        # Приборная калибровка из файла
        cal = smp.coef
        e = sum(c * ch**i for i, c in enumerate(cal))
        cal_tuple = tuple(cal)
        
    return cnt, np.asarray(e, dtype=float), float(smp.live), cal_tuple


def columns_on_grid(cps, var, e_meas, has_muon):
    """
    Построение матриц столбцов и дисперсий на сетке измерения.
    
    Аргументы:
        cps: словарь спектров нуклидов.
        var: словарь дисперсий нуклидов.
        e_meas: энергетическая сетка измерения.
        has_muon: флаг наличия мюонного компонента.
        
    Возвращает:
        (names, A, V)
    """
    # Устанавливаем флаг модуля fpc для включения/исключения мюонного столбца
    fpc.NO_MUON = not has_muon
    
    names, cols, vars_ = fpc.build_columns(cps, var, F_RN, R_TH, F_TN)
    
    # Переводим каждый столбец и дисперсию на сетку измерения
    A_cols = []
    V_cols = []
    for i, name in enumerate(names):
        col_rebinned = to_meas_grid(cols[i], e_meas)
        var_rebinned = to_meas_grid(vars_[i], e_meas)
        A_cols.append(col_rebinned)
        V_cols.append(var_rebinned)
        
    # Собираем матрицы: каналы x столбцы
    A = np.column_stack(A_cols)
    V = np.column_stack(V_cols)
    
    return names, A, V


def bands_table(title, meas_cps, pred_cps, e_meas):
    """
    Печать таблицы сравнения по энергетическим полосам.
    
    Аргументы:
        title: заголовок таблицы.
        meas_cps: массив измеренных cps (имп/с) на сетке e_meas.
        pred_cps: массив предсказанных cps (имп/с) на сетке e_meas.
        e_meas: энергетическая сетка.
        
    Возвращает:
        список кортежей (lo, hi, m, p, ratio)
    """
    print(f"\n{title}")
    print("-" * 60)
    print("%-13s | %12s | %12s | %7s" % ("полоса, кэВ", "измерено", "модель", "м/и"))
    print("-" * 60)
    
    results = []
    
    # По полосам BANDS
    for lo, hi in BANDS:
        mask = (e_meas >= lo) & (e_meas < hi)
        m_sum = np.sum(meas_cps[mask])
        p_sum = np.sum(pred_cps[mask])
        
        if m_sum > 0:
            ratio = p_sum / m_sum
        else:
            ratio = float('nan')
            
        print("%4d-%4d кэВ | %12.6f | %12.6f | %7.3f" % (lo, hi, m_sum, p_sum, ratio))
        results.append((lo, hi, m_sum, p_sum, ratio))
        
    # Итог по всему интервалу ftc.E_LO..ftc.E_HI
    mask_total = (e_meas >= ftc.E_LO) & (e_meas < ftc.E_HI)
    m_total = np.sum(meas_cps[mask_total])
    p_total = np.sum(pred_cps[mask_total])
    
    if m_total > 0:
        ratio_total = p_total / m_total
    else:
        ratio_total = float('nan')
        
    print("%4d-%4d кэВ | %12.6f | %12.6f | %7.3f" % (ftc.E_LO, ftc.E_HI, m_total, p_total, ratio_total))
    results.append((ftc.E_LO, ftc.E_HI, m_total, p_total, ratio_total))
    
    print("-" * 60)
    return results


def main():
    # 1. Открытый фон
    cnt_o, e_o, live_o, cal_o = read_meas(ftc.MEAS_NAME, ftc.CAL_ROOM)
    print(f"\nОткрытый фон: {ftc.MEAS_NAME}")
    print(f"Живое время: {live_o/3600:.2f} ч ({live_o/86400:.2f} сут)")
    print(f"Полный счёт: {np.sum(cnt_o):.0f}, имп/с: {np.sum(cnt_o)/live_o:.4f}")
    
    # 2. Шаблоны без домика
    cps_o, var_o, hm_o, _posture_o = load_cps(FMT_OPEN, ftc.MUON_CSV)
    if not hm_o:
        print("ВНИМАНИЕ: Мюонный шаблон не найден или не указан. Работа продолжается без космической компоненты.")
    
    # 3. Подгонка по открытому фону
    names, A_o, V_o = columns_on_grid(cps_o, var_o, e_o, hm_o)
    sel = (e_o >= ftc.E_LO) & (e_o < ftc.E_HI)
    y_o = cnt_o[sel]
    A_o_sel = A_o[sel, :]
    V_o_sel = V_o[sel, :]
    
    # Веса пуассоновские
    w = 1.0 / np.sqrt(np.maximum(y_o, 1.0))
    
    title_open = "Открытый фон — источник известных активностей"
    note_open = f"F_RN={F_RN}, R_TH={R_TH}, F_TN={F_TN}, TAIL_T={TAIL_T}"
    
    # Вызов fit печатает таблицу амплитуд сама
    amp, sd, pred_o, chi2_o, shape_o = ftc.fit(
        A_o_sel * live_o, 
        y_o, 
        w, 
        names, 
        title_open, 
        note_open, 
        var_counts=V_o_sel * live_o * live_o
    )
    
    # 4. Активности в физических единицах
    # ⚠ Это МОДЕЛЬНАЯ ОЦЕНКА, а не измерение: перевод в Бк/кг идёт через
    # расчёт поля в помещении, геометрия которого принята, а не измерена.
    print("\nИзвестные активности (модельная оценка по открытому фону), Бк/кг:")
    A_K, A_Ra, A_Th = amp[0], amp[1], amp[2]
    A_Pb214 = A_Bi214 = A_Ra * (1 - F_RN)
    A_Pb212 = A_Bi212 = A_Th * R_TH * (1 - F_TN)
    A_Tl208 = fpc.TL208_BRANCH * A_Bi212
    print(f"  K-40                    {A_K:9.3f}")
    print(f"  Ra-226                  {A_Ra:9.3f}")
    print(f"  Pb-214 = Bi-214         {A_Pb214:9.3f}")
    print(f"  Ac-228                  {A_Th:9.3f}")
    print(f"  Pb-212 = Bi-212         {A_Pb212:9.3f}")
    print(f"  Tl-208                  {A_Tl208:9.3f}")
    if hm_o:
        print(f"  мюоны, с^-1             {amp[3]:9.3f}")
    
    # 5. Измерение в домике
    cnt_s, e_s, live_s, cal_s = read_meas(MEAS_SHIELD, None)
    print(f"\nИзмерение в домике: {MEAS_SHIELD}")
    print(f"Живое время: {live_s/3600:.2f} ч ({live_s/86400:.2f} сут)")
    print(f"Полный счёт: {np.sum(cnt_s):.0f}, имп/с: {np.sum(cnt_s)/live_s:.4f}")
    print(f"Коэффициенты калибровки: {cal_s[0]:.6f}, {cal_s[1]:.6f}, {cal_s[2]:.6f}")
    
    # 6. Шаблоны с домиком
    cps_s, var_s, hm_s, posture_s = load_cps(FMT_SHIELD, MUON_SHIELD_CSV)
    if not hm_s:
        print("\nВНИМАНИЕ: Мюонный шаблон с домиком не найден или не указан.")
        print("Предсказание строится БЕЗ космической компоненты. Ожидается дефицит в жёсткой части спектра.")
    # Постановка расчёта (P-005) — обязательная оговорка при публикации:
    # None означает шаблоны БЕЗ posture-полей (собраны ДО 30.08.2026 —
    # прежнее допущение "центр габарита в (0,0,0), экраном вниз").
    if posture_s is not None:
        stand, screen = posture_s
        print(f"\nПостановка расчёта (P-005): опора {stand} мм, "
              f"экран {'ВВЕРХ' if screen == 1.0 else 'вниз'}.")
    else:
        print("\nПостановка расчёта: посадка НЕ указана в шаблонах — это "
              "старые файлы (до 30.08.2026, допущение P-005, экраном вниз, "
              "кристалл в центре полости). Якорь на них НЕ актуален.")

    # 7. Предсказание
    names_s, A_s, V_s = columns_on_grid(cps_s, var_s, e_s, hm_s)
    sel_s = (e_s >= ftc.E_LO) & (e_s < ftc.E_HI)
    y_s = cnt_s[sel_s]
    A_s_sel = A_s[sel_s, :]
    V_s_sel = V_s[sel_s, :]
    
    # Перенос амплитуд по именам
    amp_map = dict(zip(names, amp))
    amp_s_list = []
    for n in names_s:
        if n not in amp_map:
            print(f"Ошибка: имя '{n}' из шаблонов с домиком отсутствует в результатах подгонки открытого фона.")
            sys.exit(1)
        amp_s_list.append(amp_map[n])
        
    amp_s = np.array(amp_s_list)
    
    # Предсказанный счёт
    pred_s_counts = (A_s_sel @ amp_s) * live_s
    
    # 8. Сравнение
    meas_cps_s = y_s / live_s
    pred_cps_s = pred_s_counts / live_s
    
    bands_table("Сравнение: Измерение в домике vs Предсказание", meas_cps_s, pred_cps_s, e_s[sel_s])
    
    # Chi2/ndf и невязка формы. ftc.metrics возвращает ДВА значения:
    # (chi2/ndf, невязка формы). Число параметров ноль — в предсказании не
    # подбиралось ни одного, поэтому ndf = число каналов.
    chi2ndf_s, shape_s = ftc.metrics(pred_s_counts, y_s, 0)
    print(f"\nchi2/ndf = {chi2ndf_s:.3f} (свободных параметров НЕТ, ndf = {len(y_s)})")
    print(f"невязка формы = {shape_s:.4f}")
    
    # 9. Ослабление
    print("\nОслабление фона домиком (во сколько раз счёт упал):")
    print("%-13s | %9s | %9s | %7s" % ("полоса, кэВ", "измерено", "модель", "и/м"))
    print("-" * 50)

    # ⚠ pred_o и pred_s_counts живут на ОТОБРАННЫХ сетках (sel и sel_s), а
    # cnt_o/cnt_s — на полных. Маски обязаны строиться каждая на своей сетке,
    # иначе индексация молча съедет.
    e_o_sel = e_o[sel]
    e_s_sel = e_s[sel_s]

    # Для каждой полосы считаем ослабление
    for lo, hi in BANDS:
        # Измеренное ослабление: открытый фон / домик (оба в имп/с)
        mask_o = (e_o >= lo) & (e_o < hi)
        mask_s = (e_s >= lo) & (e_s < hi)
        
        meas_open_sum = np.sum(cnt_o[mask_o]) / live_o
        meas_shield_sum = np.sum(cnt_s[mask_s]) / live_s
        
        if meas_shield_sum > 0:
            meas_attenuation = meas_open_sum / meas_shield_sum
        else:
            meas_attenuation = float('inf')
            
        # Модельное ослабление: предсказание открытого фона / предсказание домика
        pred_open_sum = np.sum(pred_o[(e_o_sel >= lo) & (e_o_sel < hi)]) / live_o
        pred_shield_sum = np.sum(pred_s_counts[(e_s_sel >= lo) & (e_s_sel < hi)]) / live_s

        if pred_shield_sum > 0:
            model_attenuation = pred_open_sum / pred_shield_sum
        else:
            model_attenuation = float('inf')
            
        # Отношение ослаблений
        if model_attenuation > 0 and meas_attenuation != float('inf'):
            att_ratio = meas_attenuation / model_attenuation
        else:
            att_ratio = float('nan')
            
        print("%4d-%4d кэВ | %8.2f | %8.2f | %7.3f" % (lo, hi, meas_attenuation, model_attenuation, att_ratio))
        
    # Итог по всему интервалу
    mask_o_total = (e_o >= ftc.E_LO) & (e_o < ftc.E_HI)
    mask_s_total = (e_s >= ftc.E_LO) & (e_s < ftc.E_HI)
    
    meas_open_total = np.sum(cnt_o[mask_o_total]) / live_o
    meas_shield_total = np.sum(cnt_s[mask_s_total]) / live_s
    
    if meas_shield_total > 0:
        meas_att_total = meas_open_total / meas_shield_total
    else:
        meas_att_total = float('inf')
        
    pred_open_total = np.sum(pred_o) / live_o
    pred_shield_total = np.sum(pred_s_counts) / live_s
    
    if pred_shield_total > 0:
        model_att_total = pred_open_total / pred_shield_total
    else:
        model_att_total = float('inf')
        
    if model_att_total > 0 and meas_att_total != float('inf'):
        att_ratio_total = meas_att_total / model_att_total
    else:
        att_ratio_total = float('nan')
        
    print("%4d-%4d кэВ | %8.2f | %8.2f | %7.3f" % (ftc.E_LO, ftc.E_HI, meas_att_total, model_att_total, att_ratio_total))
    
    # 10. Свободная подгонка в домике (ДИАГНОСТИКА)
    print("\n" + "="*60)
    print("ДИАГНОСТИКА: Свободная подгонка спектра в домике")
    print("="*60)
    
    w_s = 1.0 / np.sqrt(np.maximum(y_s, 1.0))
    title_shield_fit = "Свободная подгонка в домике (ДИАГНОСТИКА)"
    note_shield_fit = f"F_RN={F_RN}, R_TH={R_TH}, F_TN={F_TN}, TAIL_T={TAIL_T}"
    
    amp_free, sd_free, pred_free, chi2_free, shape_free = ftc.fit(
        A_s_sel * live_s, 
        y_s, 
        w_s, 
        names_s, 
        title_shield_fit, 
        note_shield_fit, 
        var_counts=V_s_sel * live_s * live_s
    )
    
    # Отношение свободных амплитуд к известным
    print("\nОтношение свободных амплитуд (домик) к известным (открытый фон):")
    for i, name in enumerate(names_s):
        known_amp = amp_map[name]
        free_amp = amp_free[i]
        if known_amp != 0:
            ratio = free_amp / known_amp
        else:
            ratio = float('nan')
        print(f"{name}: {ratio:.4f}")


if __name__ == "__main__":
    main()
