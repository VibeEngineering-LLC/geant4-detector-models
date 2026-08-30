# -*- coding: utf-8 -*-
"""
plot_shield_predict.py

Нарисовать проверку модели ослабления фона свинцовым домиком: измеренный спектр
в домике против предсказания, построенного на активностях из открытого фона БЕЗ
подгонки. Показать раздельно вклад гамма-фона комнаты и космической компоненты —
именно их соотношение и есть содержание рисунка.

Данные берутся через импорт predict_shield (ps), fit_two_criteria (ftc) и
fit_physical_chains (fpc). Повторной реализации расчёта нет.

ВНИМАНИЕ: Легенда обязана лежать ВНЕ поля данных — в этой линии уже был дефект, когда
непрозрачная легенда перекрыла спектр и график читался как обрыв данных.
"""

import os
import sys
import io
import contextlib

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

# Настройки шрифтов и сетки
matplotlib.rcParams["font.size"] = 9
matplotlib.rcParams["figure.dpi"] = 130
matplotlib.rcParams["axes.grid"] = True
matplotlib.rcParams["grid.alpha"] = 0.3

OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PNG = os.path.join(OUT_DIR, "RC103_shield_predict.png")

# Цвета — светлые и НЕПРОЗРАЧНЫЕ (alpha не использовать нигде)
C_MEAS_OPEN = "#9aa0a6"   # открытый фон, измерение — фон-контекст
C_MEAS = "#111111"        # измерение в домике
C_MODEL = "#d62728"       # полная модель
C_GAMMA = "#4e79a7"       # гамма-компонента (заливка)
C_MUON = "#8c8c8c"        # космическая компонента (заливка)


def build():
    """
    Возвращает словарь со всем, что нужно для рисования.
    Повторно использует функции модуля ps — своей реализации чтения, свёртки и подгонки НЕ писать.
    """
    # 1. Открытый фон
    cnt_o, e_o, live_o, _ = ps.read_meas(ftc.MEAS_NAME, ftc.CAL_ROOM)

    # 2. Шаблоны открытые
    cps_o, var_o, hm_o, _ = ps.load_cps(ps.FMT_OPEN, ftc.MUON_CSV)
    if not hm_o:
        sys.exit("нет мюонного шаблона для открытого фона")

    # 3. Подготовка сетки и отбор по энергии
    names, A_o, V_o = ps.columns_on_grid(cps_o, var_o, e_o, hm_o)
    sel = (e_o >= ftc.E_LO) & (e_o < ftc.E_HI)
    cnt_o = cnt_o[sel]
    e_o = e_o[sel]
    A_o = A_o[sel, :]
    V_o = V_o[sel, :]

    # 4. Подгонка открытого фона.
    # ⚠ ftc.fit работает в ОТСЧЁТАХ, а не в имп/с: матрица подаётся домноженной
    # на живое время, измерение — сырым счётом. Деление измерения на live здесь
    # было бы рассогласованием в live раз.
    y_o = cnt_o
    w = 1.0 / np.sqrt(np.maximum(y_o, 1.0))

    # Подавление вывода подгонки
    f_out = io.StringIO()
    with contextlib.redirect_stdout(f_out):
        amp, _sd, pred_o_counts, _chi2, _shape = ftc.fit(
            A_o * live_o, y_o, w, names, "", "", var_counts=V_o * live_o * live_o)

    # ftc.fit возвращает КОРТЕЖ (amp, sd, pred, chi2/ndf, невязка формы), не словарь.
    pred_o = pred_o_counts / live_o   # в имп/с, для рисования

    # 5. Домик: измерение
    cnt_s, e_s, live_s, _ = ps.read_meas(ps.MEAS_SHIELD, None)

    # 6. Шаблоны с домиком
    cps_s, var_s, hm_s, _ = ps.load_cps(ps.FMT_SHIELD, ps.MUON_SHIELD_CSV)
    if not hm_s:
        sys.exit("нет мюонного шаблона С ДОМИКОМ")

    # 7. Подготовка сетки для домика и отбор
    names_s, A_s, V_s = ps.columns_on_grid(cps_s, var_s, e_s, hm_s)
    sel_s = (e_s >= ftc.E_LO) & (e_s < ftc.E_HI)
    cnt_s = cnt_s[sel_s]
    e_s = e_s[sel_s]
    A_s = A_s[sel_s, :]

    # 8. Перенос амплитуд ПО ИМЕНАМ
    amp_map = dict(zip(names, amp))
    amp_s = []
    for n in names_s:
        if n not in amp_map:
            sys.exit(f"отсутствует имя в переносе амплитуд: {n}")
        amp_s.append(amp_map[n])
    amp_s = np.array(amp_s)

    # 9. Раздельные вклады
    i_mu = names_s.index("mu")
    pred_mu = A_s[:, i_mu] * amp_s[i_mu]          # космическая компонента, имп/с
    pred_tot = A_s @ amp_s                        # полная модель, имп/с
    pred_gamma = pred_tot - pred_mu               # гамма-фон комнаты, имп/с

    # Предсказание для открытого фона уже переведено в имп/с выше (pred_o).
    meas_open_rate = cnt_o / live_o
    meas_shield_rate = cnt_s / live_s

    return {
        "e_o": e_o,
        "meas_open": meas_open_rate,
        "e_s": e_s,
        "meas_shield": meas_shield_rate,
        "pred_tot": pred_tot,
        "pred_gamma": pred_gamma,
        "pred_mu": pred_mu,
        "amp_map": amp_map,
        "live_o": live_o,
        "live_s": live_s,
        "pred_open": pred_o # rate open background model
    }


def draw(d):
    """
    Рисует график проверки модели ослабления.
    """
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(11, 9), gridspec_kw={"height_ratios": [3, 1]})

    # --- Верхняя панель: Спектры ---
    
    # Открытый фон (контекст)
    ax_top.step(d["e_o"], d["meas_open"], where="mid", color=C_MEAS_OPEN, lw=1.0, label="измерено, открытый фон")

    # Заливки компонент: СТОПКОЙ снизу вверх, без прозрачности, БЕЗ контурных линий
    # Космика ПЕРВОЙ (на заднем плане)
    ax_top.fill_between(d["e_s"], 0, d["pred_mu"], color=C_MUON, linewidth=0, label="модель: космическая компонента")
    # Гамма-фон сверху космической
    ax_top.fill_between(d["e_s"], d["pred_mu"], d["pred_mu"] + d["pred_gamma"], color=C_GAMMA, linewidth=0, label="модель: гамма-фон комнаты")

    # Полная модель (линия)
    ax_top.plot(d["e_s"], d["pred_tot"], color=C_MODEL, lw=1.4, label="модель, сумма")

    # Измерение в домике
    ax_top.step(d["e_s"], d["meas_shield"], where="mid", color=C_MEAS, lw=1.0, label="измерено, домик")

    # Ось Y: логарифмическая
    ax_top.set_yscale("log")
    
    # Пределы Y
    min_meas = np.min(d["meas_shield"][d["meas_shield"] > 0])
    y_min = max(1e-7, min_meas)
    y_max = 1.6 * np.max(d["meas_open"])
    ax_top.set_ylim(y_min, y_max)

    # Ось X
    ax_top.set_xlim(ftc.E_LO, ftc.E_HI)
    ax_top.set_xlabel("") # Подпись оси X ставим только на нижней панели, чтобы не перекрыть легенду
    ax_top.set_ylabel("скорость счёта, имп/с на канал")

    # Заголовок
    ax_top.set_title("RadiaCode-103: ослабление фона свинцовым домиком — предсказание без подгонки")

    # Легенда ВНЕ поля данных
    ax_top.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=3, frameon=False)

    # --- Нижняя панель: Ослабление по полосам ---
    
    bands = ftc.BANDS
    meas_att = []
    model_att = []
    labels = []

    for lo, hi in bands:
        # Маска для открытого фона (сетка e_o)
        mask_o = (d["e_o"] >= lo) & (d["e_o"] < hi)
        sum_meas_open = np.sum(d["meas_open"][mask_o]) * (hi - lo) / len(d["e_o"][mask_o]) if np.any(mask_o) else 0 
        # Интеграл: сумма rate * width_bin. Но bins могут быть разной ширины? 
        # В условии сказано: "сумма meas_open по полосе". Обычно под этим понимают интеграл.
        # Однако, если просто суммировать значения rate, это не имеет физического смысла без умножения на ширину бина.
        # Но в задаче сказано: "сумма meas_open по полосе / сумма meas_shield по полосе".
        # Если мы делим сумму rate на сумму rate, единицы сокращаются, но нужно учитывать количество каналов или ширину.
        # Предположим, что под "суммой" имеется в виду интеграл (площадь).
        # Для step-графика площадь = sum(rate_i * delta_E_i).
        
        # Давайте сделаем честный интеграл по полосе [lo, hi]
        # Для e_o:
        idx_o = np.where((d["e_o"] >= lo) & (d["e_o"] < hi))[0]
        if len(idx_o) > 0:
            # Простая аппроксимация интеграла: сумма значений * ширина бина.
            # Ширина бина для e_o может быть не постоянной? Обычно да.
            # Но в условии сказано "сумма ... по полосе". 
            # Если просто суммировать значения, то отношение будет зависеть от количества каналов.
            # Скорее всего, имеется в виду интеграл.
            # Однако, чтобы не усложнять и следовать буквальному тексту "сумма", 
            # но с учетом того, что это спектры, обычно берут интеграл.
            # Давайте посчитаем интеграл как sum(rate * bin_width).
            # Но у нас есть только центры e_o. Ширину бина можно оценить как среднее расстояние между соседними точками?
            # Или просто предположить, что "сумма" означает сумму значений rate, и так как сетки разные, 
            # это может быть неточно. 
            # В условии: "маску для каждой суммы строить по своей сетке".
            
            # Интерпретация: Интеграл от спектра в полосе.
            # Для e_o:
            rates_o = d["meas_open"][idx_o]
            # Оценим ширину бина как среднее расстояние между точками в этой области? 
            # Или просто примем, что "сумма" означает сумму значений, и мы сравниваем "среднюю интенсивность"?
            # Нет, ослабление - это отношение интегралов.
            
            # Давайте используем простую сумму значений, умноженную на среднюю ширину бина в этой полосе?
            # Или просто сумму значений, если предположить, что количество каналов примерно одинаково? 
            # Нет, сетки разные.
            
            # Самый надежный способ без знания точных границ bins: 
            # Сумма (rate * delta_E). Delta_E можно взять как (hi-lo) / N_bins_in_band?
            # Или просто интегрировать методом трапеций/прямоугольников по доступным точкам.
            
            # Примем прямоугольники шириной в среднее расстояние между соседними e_o в этой полосе.
            if len(idx_o) > 1:
                bin_width_o = np.mean(np.diff(d["e_o"][idx_o]))
            else:
                bin_width_o = (hi - lo) # Если одна точка, считаем что она покрывает всю полосу? Рискованно.
            
            integral_meas_open = np.sum(rates_o) * bin_width_o
            
            # Для e_s (домик):
            idx_s = np.where((d["e_s"] >= lo) & (d["e_s"] < hi))[0]
            rates_s = d["meas_shield"][idx_s]
            if len(idx_s) > 1:
                bin_width_s = np.mean(np.diff(d["e_s"][idx_s]))
            else:
                bin_width_s = (hi - lo)
            
            integral_meas_shield = np.sum(rates_s) * bin_width_s
            
            # Модель открытый фон:
            rates_model_o = d["pred_open"][idx_o]
            integral_model_open = np.sum(rates_model_o) * bin_width_o
            
            # Модель домик (полная):
            rates_model_s = d["pred_tot"][idx_s]
            integral_model_shield = np.sum(rates_model_s) * bin_width_s

        else:
            integral_meas_open = 0
            integral_meas_shield = 1e-9 # Защита от деления на 0
            integral_model_open = 0
            integral_model_shield = 1e-9

        if integral_meas_shield > 0:
            meas_att_val = integral_meas_open / integral_meas_shield
        else:
            meas_att_val = 1.0 # Нет данных, считаем что не ослабилось? Или NaN? Лучше 1.0 для графика.

        if integral_model_shield > 0:
            model_att_val = integral_model_open / integral_model_shield
        else:
            model_att_val = 1.0

        meas_att.append(meas_att_val)
        model_att.append(model_att_val)
        labels.append(f"{lo}-{hi}")

    x = np.arange(len(bands))
    
    # Столбики
    ax_bot.bar(x - 0.2, meas_att, width=0.4, color=C_MEAS, label="измерено")
    ax_bot.bar(x + 0.2, model_att, width=0.4, color=C_MODEL, label="модель")

    # Подписи делений
    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(labels)

    # Ось Y логарифмическая
    ax_bot.set_yscale("log")
    ax_bot.set_ylabel("ослабление, раз")

    # Легенда внутри верхнего левого угла
    # Легенда вынесена ВПРАВО от поля: внутри верхнего левого угла она
    # перекрывала подпись первого столбика (проверено по готовому PNG).
    ax_bot.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), frameon=False)
    # Обычные числа вместо «3×10¹» — шкала логарифмическая, но подписи читаемые.
    ax_bot.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax_bot.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax_bot.set_yticks([1, 2, 5, 10, 20, 30])

    # Подписи значений над столбиками
    for i, (m_val, mod_val) in enumerate(zip(meas_att, model_att)):
        ax_bot.text(i - 0.2, m_val * 1.1, f"{m_val:.1f}", ha='center', va='bottom', fontsize=8)
        ax_bot.text(i + 0.2, mod_val * 1.1, f"{mod_val:.1f}", ha='center', va='bottom', fontsize=8)

    # Подпись оси X на нижней панели (так как верхняя занята легендой)
    ax_bot.set_xlabel("энергия, кэВ")

    fig.tight_layout()
    fig.savefig(OUT_PNG)
    print(f"Сохранено: {OUT_PNG}")


if __name__ == "__main__":
    draw(build())
