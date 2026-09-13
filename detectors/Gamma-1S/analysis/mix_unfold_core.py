"""Ядро разбора смеси AmTiCsEu на Гамма-1С (11.09.2026).

Одна реализация подгонки нетто-спектра суммой МК-шаблонов: ею пользуются и
командная строка mix_unfold_g1s.py, и выгрузка веб-страницы
(web-th232/export_amticseu_data.py). Две реализации одного расчёта — дефект,
поэтому расчёт живёт только здесь. Сгенерировано qwen3-coder:30b по спеке
_spec_mix_unfold_core.md, принято с правкой amplitude_errors.
"""
import os
import sys
import numpy as np
# scipy.optimize.nnls больше не вызывается здесь: подгонку ведёт crit_bench (D-020). Оставленный
# мёртвый импорт создавал бы впечатление, что расчёт живёт в этом файле.

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mix_unfold_g1s as g1s   # read_lsrm_spe, read_template, make_fwhm, broaden, recalibrate_energy, TAIL_T
# Критерии подгонки — из стенда #CRIT-1, импортом, а не копией (§33): правка донора обязана
# расходиться на всех потребителей. Рабочий критерий A2 и вторая мера E1 выбраны решением D-020.
import crit_bench as cb        # fit_A1 (дисперсия), fit_A2 (рабочий критерий), fit_E1 (вторая мера)

def channel_edges(e):
    """Возвращает границы каналов по центрам."""
    return np.concatenate(([e[0] - 0.5*(e[1]-e[0])], 0.5*(e[:-1]+e[1:]), [e[-1] + 0.5*(e[-1]-e[-2])]))

def rebin_counts(counts, src_edges, dst_edges):
    """Перенос гистограммы с одной шкалы на другую с сохранением счета."""
    if len(src_edges) != len(counts)+1:
        raise ValueError("Несовпадение длин массивов: длина границ должна быть на 1 больше длины отсчетов")
    if not np.all(np.diff(src_edges) > 0):
        raise ValueError("Шкала источника не возрастает")
    if not np.all(np.diff(dst_edges) > 0):
        raise ValueError("Шкала приемника не возрастает")

    cum = np.concatenate(([0.0], np.cumsum(np.asarray(counts, dtype=float))))
    v = np.interp(dst_edges, src_edges, cum, left=0.0, right=cum[-1])
    return np.diff(v)

def unfold(spe, bg, templates, fwhm_points, lo=40.0, hi=1500.0,
          recalibrate=False, tail=None, bg_energy_of_ch=None, verbose=True, extra_refs=None,
          conv="energy", blur=1.0, ch_offset=0.0, light_scale=None):
    """conv="channel" — свёртка шаблонов в КАНАЛАХ (g1s.broaden_ch), blur — множитель ширины
    (для шаблонов npsm=on: собственный разброс уже в свете, этап 5 — 0,798).
    light_scale — пара (a, b) прямой «канал = a + b·свет» (analysis/light_scale_g1s.py). Задана —
    гистограммы шаблонов npsm=on читаются как СВЕТ и переводятся в энергию шкалы пробы
    (analysis/light_mode.py); не задана — прежнее поведение, свет трактуется как энергия."""
    """Основная функция разбора смеси AmTiCsEu на Гамма-1С."""
    if isinstance(spe, str):
        spec = g1s.read_lsrm_spe(spe)
    else:
        spec = spe

    if isinstance(bg, str):
        bgs = g1s.read_lsrm_spe(bg)
    else:
        bgs = bg

    if tail is not None:
        g1s.TAIL_T = float(tail)

    # Шаг c: границы каналов пробы
    if recalibrate:
        e_of_ch, n_refs = g1s.recalibrate_energy(spec, verbose=verbose, extra_refs=extra_refs,
                                                 ch_offset=ch_offset)
        e = np.array([e_of_ch(i) for i in range(spec.n_channels)])
    else:
        n_refs = 0
        e = np.array([spec.channel_to_energy(i) for i in range(spec.n_channels)])
    ch_edges = channel_edges(e)

    # Шаг b: масштабирование фона
    # k_bg — тот самый множитель, которым фон приведён к живому времени пробы. Дисперсия
    # масштабированного фона равна k²·фон = k_bg·bg_scaled, поэтому критерию нужен именно он,
    # и берётся он отсюда, а не вычисляется вторично (два источника одной величины разошлись бы).
    k_bg = spec.live_time / bgs.live_time
    if bg_energy_of_ch is None:
        bg_scaled = np.array(bgs.counts) * k_bg
    else:
        eb = np.asarray(bg_energy_of_ch, dtype=float)
        if len(eb) != len(bgs.counts):
            raise ValueError("Длина bg_energy_of_ch не совпадает с числом каналов фона")
        bg_on = rebin_counts(np.asarray(bgs.counts, dtype=float), channel_edges(eb), ch_edges)
        bg_scaled = bg_on * k_bg

    # Шаг d: чистые отсчёты
    y = np.maximum(np.array(spec.counts) - bg_scaled, 0)

    # Шаг e: свёртка шаблонов
    names, cols, n_events = [], [], []
    fwhm_func = g1s.make_fwhm(fwhm_points)
    for name, path in templates:
        hist, n_events_i, npsm = g1s.read_template(path)
        if light_scale is not None:
            # Шаблон npsm=on даёт СВЕТ в своих единицах. Переводим: канал = a + b·свет, затем
            # энергия шкалы пробы. Без этого свет трактуется как энергия и шаблон уезжает вниз.
            import light_mode
            if npsm != 1:
                raise SystemExit("ОТКАЗ: задана шкала света, но шаблон %s посчитан без "
                                 "непропорциональности (npsm_enabled=0)" % name)
            e_of = lambda c: float(np.interp(c, np.arange(len(e), dtype=float), e))
            hist, st = light_mode.light_to_energy(hist, light_scale[0], light_scale[1], e_of, len(e))
            if st["sum_used"] <= 0:
                raise SystemExit("ОТКАЗ: шаблон %s — в шкалу пробы не попал ни один бин света" % name)
            if verbose:
                print("  %s: свет → энергия, бинов %d из %d, счёта %.1f%%"
                      % (name, st["n_used"], st["n_in"], 100.0 * st["sum_used"] / st["sum_in"]))
        if conv not in ("energy", "channel"):   # B-5: опечатка молча уходила в энергию
            raise ValueError("conv должен быть 'energy' или 'channel', получено %r" % (conv,))
        if conv == "channel":
            col = g1s.broaden_ch(hist, n_events_i, e, fwhm_func,
                                 lambda c: float(spec.channel_to_energy(c + ch_offset)), blur)
        else:
            col = g1s.broaden(hist, n_events_i, ch_edges, lambda E: blur * fwhm_func(E))
        names.append(name)
        cols.append(col)
        n_events.append(n_events_i)

    # Шаг f: подгонка. Критерий — A2, вторая мера — E1 (решение D-020). Прежний критерий
    # (веса 1/√max(y,1), нетто обрезано нулём) убран БЕЗ переключателя «как было»; он остаётся
    # доступен только в стенде crit_bench.fit_A0 — там же, где меряются все восемь критериев.
    sel = (e >= lo) & (e <= hi)
    sigma = np.sqrt(np.maximum(y, 1))
    sigma[sigma == 0] = 1
    A = np.array(cols).T[sel] / sigma[sel][:, None]
    b = y[sel] / sigma[sel]
    cols_arr = np.array(cols)
    counts_arr = np.array(spec.counts, dtype=float)
    net = counts_arr - bg_scaled          # нетто БЕЗ обрезки нулём — его требует критерий
    # Вес канала 1/√(p + k·b + Σ aₖ²·colₖ/n_eventsₖ): дисперсия нетто плюс дисперсия шаблонов.
    # Веса зависят от искомых амплитуд, поэтому внутри шесть итераций.
    coef, sd, _ = cb.fit_A2(cols_arr, counts_arr, bg_scaled, k_bg, n_events, sel)
    activities = coef / spec.live_time
    model = np.sum(cols_arr * coef[:, None], axis=0)
    # χ² — по дисперсии ТОГО ЖЕ критерия, которым шла подгонка: иначе χ²/n.d.f. относится к
    # одному критерию, а амплитуды к другому. Дисперсия измерения берётся у стенда (fit_A1),
    # слагаемое шаблонов — той же формулой, что внутри A2.
    var_a1 = cb.fit_A1(cols_arr, counts_arr, bg_scaled, k_bg, n_events, sel)[2]["var"]
    var_tpl = (cols_arr[:, sel].T / np.asarray(n_events, dtype=float)[None, :]) @ (coef ** 2)
    chi2 = float(np.sum((model[sel] - net[sel]) ** 2 / (var_a1 + var_tpl)))
    # Единая метрика сравнения критериев (дисперсия A1) — ею меряются все восемь в #CRIT-1.
    chi2_ref = float(np.sum((model[sel] - net[sel]) ** 2 / var_a1))
    # Дисперсия критерия НА ВСЕХ каналах — нужна долям полос. Считать вклад полосы по одним
    # весам, а делить на χ² по другим нельзя: 12.09.2026 это дало долю 120 % у полосы 40–90
    # (поймано прогоном командной строки сразу после переноса критерия).
    var = np.maximum(counts_arr + k_bg * bg_scaled, 1.0) \
        + (cols_arr.T / np.asarray(n_events, dtype=float)[None, :]) @ (coef ** 2)
    ndof = int(np.sum(sel)) - len(coef)   # честное ν: параметры ЭТОГО критерия (D-020, п. 4)

    # Вторая мера — E1: минимум полной вариации между нормированными формами при условии
    # Σ модель = Σ нетто, точное ЛП; погрешности — параметрический бутстрап. Считается ВСЕГДА
    # и публикуется рядом с A2: расхождение критериев по Am-241 (14,6 %) — факт о модели.
    coef_e1, sd_e1, ex_e1 = cb.fit_E1(cols_arr, counts_arr, bg_scaled, k_bg, n_events, sel)
    model_e1 = np.sum(cols_arr * coef_e1[:, None], axis=0)
    e1 = {"coef": coef_e1, "activities": coef_e1 / spec.live_time, "sd": sd_e1,
          "model": model_e1, "tv": float(ex_e1["tv"]),
          "chi2_ref": float(np.sum((model_e1[sel] - net[sel]) ** 2 / var_a1))}

    return {
        "spec": spec,
        "bg": bgs,
        "e": e,
        "ch_edges": ch_edges,
        "counts": np.array(spec.counts, dtype=float),
        "bg_scaled": bg_scaled,
        "y": y,
        "sigma": sigma,
        "sel": sel,
        "names": names,
        "cols": np.array(cols),
        "n_events": n_events,
        "net": net,
        "var": var,
        "coef": coef,
        "sd": sd,
        "activities": activities,
        "model": model,
        "crit": "A2",
        "chi2": chi2,
        "chi2_ref": chi2_ref,
        "ndof": ndof,
        "e1": e1,
        "A": A,
        "b": b,
        "n_refs": n_refs,
        "fwhm": fwhm_func,
        "live_s": float(spec.live_time)
    }

def amplitude_errors(A, coef):
    """Статистическая погрешность амплитуд NNLS в единицах coef: ковариация
    (A_act.T @ A_act)⁻¹ по активным столбцам (coef > 0). У прижатых к нулю
    столбцов погрешность не определена — nan. Вырожденная матрица — громкий
    отказ, а не тихий nan (прежний вариант генерации глотал LinAlgError)."""
    coef = np.asarray(coef, dtype=float)
    act = coef > 0
    result = np.full(len(coef), np.nan)
    if not act.any():
        return result
    try:
        cov = np.linalg.inv(A[:, act].T @ A[:, act])
    except np.linalg.LinAlgError as ex:
        raise ValueError("ОТКАЗ: матрица плана вырождена, погрешности не определены: %s" % ex)
    result[act] = np.sqrt(np.diag(cov))
    return result

def band_shares(r, bands=((40, 90), (90, 200), (200, 500), (500, 1000), (1000, 1500))):
    """Доли хи-квадрат по полосам."""
    result = []
    for lo_b, hi_b in bands:
        mask = (r["e"] >= lo_b) & (r["e"] <= hi_b)
        if mask.sum() == 0:
            continue
        # Та же дисперсия и то же нетто, что у χ² критерия: иначе доля полосы считается в
        # одной метрике, а нормируется на другую, и сумма долей расходится с сотней.
        chi2_band = float(np.sum((r["model"] - r["net"])[mask] ** 2 / r["var"][mask]))
        share_pct = chi2_band / r["chi2"] * 100
        result.append({"lo": lo_b, "hi": hi_b, "chi2": chi2_band, "share_pct": share_pct})
    return result

def zone(r, lo=50.0, hi=70.0):
    """Информация о конкретной полосе."""
    mask = (r["e"] >= lo) & (r["e"] <= hi)
    if mask.sum() == 0:
        raise ValueError("В окне не попало ни одного канала")
    sum_meas = float(np.sum(r["y"][mask]))
    sum_model = float(np.sum(r["model"][mask]))
    if sum_meas == 0:
        raise ValueError("Сумма измерений в окне равна нулю")
    ratio = sum_model / sum_meas

    model_counts = np.array([float(np.sum(r["cols"][i][mask]) * r["coef"][i]) for i in range(len(r["cols"]))])
    model_shares_pct = model_counts / np.sum(model_counts) * 100
    template_shares_pct = np.array([float(np.sum(r["cols"][i][mask])) for i in range(len(r["cols"]))])
    template_shares_pct = template_shares_pct / np.sum(template_shares_pct) * 100

    peak_meas_keV = r["e"][mask][np.argmax(r["y"][mask])]
    peak_model_keV = r["e"][mask][np.argmax(r["model"][mask])]

    return {
        "sum_meas": sum_meas,
        "sum_model": sum_model,
        "ratio": ratio,
        "model_counts": model_counts,
        "model_shares_pct": model_shares_pct,
        "template_shares_pct": template_shares_pct,
        "peak_meas_keV": peak_meas_keV,
        "peak_model_keV": peak_model_keV,
        "names": r["names"]
    }

def net_area(e, col, E0, F, win=1.25, side=(1.6, 3.0)):
    """Чистая площадь пика."""
    d = np.abs(e - E0)
    inwin = d <= win * F
    sb = (d >= side[0] * F) & (d <= side[1] * F)
    if sb.sum() < 3 or inwin.sum() == 0:
        raise ValueError("Недостаточно точек для подложки или окно пустое")
    k, c = np.polyfit(e[sb], col[sb], 1)
    return float(np.sum(col[inwin] - (k * e[inwin] + c)))

def peak_row(table, E0, tol):
    """Поиск строки с ближайшим пиком."""
    if not table:
        return None
    min_diff = min(abs(row["energy_keV"] - E0) for row in table)
    if min_diff <= tol:
        return next(row for row in table if abs(row["energy_keV"] - E0) == min_diff)
    return None

def selftest():
    """Самопроверка."""
    # T1
    counts = np.array([1, 2, 3, 4])
    src_edges = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    dst_edges = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    result = rebin_counts(counts, src_edges, dst_edges)
    if not np.allclose(result, counts):
        print("T1 FAIL: тождественная шкала не совпала")
        return 1

    # T2
    try:
        dst_edges = np.array([4.0, 3.0, 2.0, 1.0, 0.0])
        rebin_counts(counts, src_edges, dst_edges)
        print("T2 FAIL: убывающая шкала не бросила ошибку")
        return 1
    except ValueError:
        pass

    # T3
    e = np.arange(0, 400, 1.0) + 0.5
    sigma = 10 / 2.3548
    col = 1000 * np.exp(-0.5 * ((e - 200) / sigma)**2) / (sigma * np.sqrt(2 * np.pi)) + (5 + 0.02 * e)
    try:
        area = net_area(e, col, 200, 10)
        if abs(area - 1000) > 20:
            print("T3 FAIL: площадь не совпала с точностью 2%")
            return 1
    except Exception as ex:
        print(f"T3 FAIL: ошибка в net_area: {ex}")
        return 1

    # T4
    A = np.eye(3)
    coef = [1, 0, 2]
    err = amplitude_errors(A, coef)
    expected = [1.0, np.nan, 1.0]
    if not np.allclose(err, expected, equal_nan=True):
        print("T4 FAIL: погрешности не совпали")
        return 1

    # T5
    try:
        r = {
            "e": np.array([100., 200.]),
            "y": np.array([0., 0.]),
            "model": np.array([0., 0.]),
            "sigma": np.array([1., 1.]),
            "cols": [np.array([0., 0.]), np.array([0., 0.])],
            "coef": np.array([0., 0.])
        }
        zone(r, lo=50.0, hi=70.0)
        print("T5 FAIL: не бросила ValueError при пустом окне")
        return 1
    except ValueError:
        pass

    print("SELFTEST OK")
    return 0

if __name__ == "__main__":
    sys.exit(selftest())
