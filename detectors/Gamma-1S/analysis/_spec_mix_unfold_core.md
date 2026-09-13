Напиши ОДИН модуль Python 3 целиком. Верни ТОЛЬКО код, без объяснений и без
markdown-заборов. Комментарии и docstring — на РУССКОМ языке. Никаких голых
`except`, никаких `except: pass`: любая ошибка входа — громкий отказ
(`raise SystemExit("ОТКАЗ: ...")` или `raise ValueError(...)` с русским текстом).

ИМЯ ФАЙЛА: mix_unfold_core.py (лежит в той же папке, что mix_unfold_g1s.py).

ШАПКА МОДУЛЯ (docstring): ядро разбора смеси AmTiCsEu на Гамма-1С. Одна
реализация подгонки нетто-спектра суммой МК-шаблонов, которой пользуются и
командная строка mix_unfold_g1s.py, и выгрузка веб-страницы
(web-th232/export_amticseu_data.py). Две реализации одного расчёта — дефект,
поэтому расчёт живёт только здесь.

ИМПОРТЫ (ровно так):
    import os
    import sys
    import numpy as np
    from scipy.optimize import nnls
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import mix_unfold_g1s as g1s   # read_lsrm_spe, read_template, make_fwhm, broaden, recalibrate_energy, TAIL_T

ФУНКЦИИ:

1. channel_edges(e)
   e — массив энергий центров каналов (возрастает). Вернуть границы каналов,
   длина len(e)+1, СТРОГО этой формулой:
       np.concatenate(([e[0] - 0.5*(e[1]-e[0])], 0.5*(e[:-1]+e[1:]), [e[-1] + 0.5*(e[-1]-e[-2])]))

2. rebin_counts(counts, src_edges, dst_edges)
   Перенос гистограммы с одной энергетической шкалы на другую С СОХРАНЕНИЕМ
   СЧЁТА. counts — отсчёты по каналам источника (длина N), src_edges — их
   границы в кэВ (длина N+1, строго возрастают), dst_edges — границы каналов
   приёмника (длина M+1). Внутри канала источника счёт считается равномерным по
   энергии. Алгоритм:
       cum = np.concatenate(([0.0], np.cumsum(np.asarray(counts, dtype=float))))
       v = np.interp(dst_edges, src_edges, cum, left=0.0, right=cum[-1])
       return np.diff(v)
   Перед этим проверить: len(src_edges) == len(counts)+1, иначе ValueError;
   np.all(np.diff(src_edges) > 0) и np.all(np.diff(dst_edges) > 0), иначе
   ValueError («шкала не возрастает»).

3. unfold(spe, bg, templates, fwhm_points, lo=40.0, hi=1500.0,
          recalibrate=False, tail=None, bg_energy_of_ch=None, verbose=True)
   spe, bg — путь к .spe (str) либо уже прочитанный объект спектра; строку
   читать через g1s.read_lsrm_spe(path).
   templates — список пар (имя, путь_к_csv).
   fwhm_points — путь к CSV точек ПШПВ (передаётся в g1s.make_fwhm).
   tail — если не None: присвоить g1s.TAIL_T = float(tail) ДО свёртки шаблонов.
   bg_energy_of_ch — None либо массив энергий центров каналов ФОНА по его
   собственной проверенной шкале.

   Шаги СТРОГО в этом порядке и этими выражениями (результат обязан совпасть
   бит в бит с прежней командной строкой при bg_energy_of_ch=None):
     a) spec, bgs — объекты спектров.
     b) если bg_energy_of_ch is None:
            bg_scaled = np.array(bgs.counts) * (spec.live_time / bgs.live_time)
        иначе:
            eb = np.asarray(bg_energy_of_ch, dtype=float)
            если len(eb) != len(bgs.counts) — ValueError;
            bg_on = rebin_counts(np.asarray(bgs.counts, dtype=float), channel_edges(eb), ch_edges_пробы)
            bg_scaled = bg_on * (spec.live_time / bgs.live_time)
        (для ветки «иначе» границы каналов пробы нужны заранее — поэтому шаг c
        выполнить ДО шага b; значения от этого не меняются.)
     c) энергии каналов пробы:
            если recalibrate: e_of_ch, n_refs = g1s.recalibrate_energy(spec, verbose=verbose)
                              e = np.array([e_of_ch(i) for i in range(spec.n_channels)])
            иначе: n_refs = 0
                   e = np.array([spec.channel_to_energy(i) for i in range(spec.n_channels)])
            ch_edges = channel_edges(e)
     d) y = np.maximum(np.array(spec.counts) - bg_scaled, 0)
     e) fwhm_func = g1s.make_fwhm(fwhm_points); для каждого (имя, путь):
            hist, n_events, npsm = g1s.read_template(путь)
            col = g1s.broaden(hist, n_events, ch_edges, fwhm_func)
        собрать names (список имён), cols (список массивов), n_events (список int).
     f) sel = (e >= lo) & (e <= hi)
        sigma = np.sqrt(np.maximum(y, 1)); sigma[sigma == 0] = 1
        A = np.array(cols).T[sel] / sigma[sel][:, None]
        b = y[sel] / sigma[sel]
        coef, _ = nnls(A, b)
        activities = coef / spec.live_time
        model = np.sum(np.array(cols) * coef[:, None], axis=0)
        chi2 = float(np.sum(((model[sel] - y[sel]) / sigma[sel]) ** 2))
        ndof = int(np.sum(sel)) - len(templates)
   Вернуть dict с ключами: "spec", "bg", "e", "ch_edges", "counts"
   (np.array(spec.counts, dtype=float)), "bg_scaled", "y", "sigma", "sel",
   "names", "cols" (np.array(cols), форма (число шаблонов, число каналов)),
   "n_events", "coef", "activities", "model", "chi2", "ndof", "A", "b",
   "n_refs", "fwhm" (fwhm_func), "live_s" (float(spec.live_time)).

4. amplitude_errors(A, coef)
   Статистическая погрешность амплитуд NNLS в единицах coef. Ковариация —
   обратная к A_act.T @ A_act, где A_act — столбцы A с coef > 0. Для
   прижатых к нулю столбцов (coef == 0) погрешность не определена — nan.
   Вернуть np.array длины len(coef).

5. band_shares(r, bands=((40, 90), (90, 200), (200, 500), (500, 1000), (1000, 1500)))
   r — результат unfold. Для каждой полосы (lo_b, hi_b):
       mask = (r["e"] >= lo_b) & (r["e"] <= hi_b)
       если mask.sum() == 0 — полосу пропустить;
       chi2_band = float(np.sum(((r["model"] - r["y"])[mask] / r["sigma"][mask]) ** 2))
       share_pct = chi2_band / r["chi2"] * 100
   Вернуть список dict {"lo": lo_b, "hi": hi_b, "chi2": chi2_band, "share_pct": share_pct}.

6. zone(r, lo=50.0, hi=70.0)
   mask = (r["e"] >= lo) & (r["e"] <= hi); если mask.sum() == 0 — ValueError.
   sum_meas = float(np.sum(r["y"][mask])); sum_model = float(np.sum(r["model"][mask]))
   ratio = sum_model / sum_meas если sum_meas != 0, иначе ValueError.
   model_counts: для каждого i — float(np.sum(r["cols"][i][mask]) * r["coef"][i])
       (ВКЛАД шаблона в модель — с его амплитудой);
   model_shares_pct = model_counts / sum(model_counts) * 100.
   template_shares_pct: доли НЕМАСШТАБИРОВАННЫХ шаблонов np.sum(r["cols"][i][mask])
       (так печатала прежняя командная строка; величина — вероятность на
       распад, как «вклад» она смысла не имеет, оставлена только для
       сопоставимости со старым выводом).
   Также: peak_meas_keV и peak_model_keV — энергии максимумов y и model внутри
   mask (r["e"][mask][np.argmax(...)]).
   Вернуть dict со всеми этими полями и "names".

7. net_area(e, col, E0, F, win=1.25, side=(1.6, 3.0))
   Чистая площадь пика на спектре col (отсчёты по каналам с энергиями e):
       d = np.abs(e - E0)
       inwin = d <= win * F
       sb = (d >= side[0] * F) & (d <= side[1] * F)
       если sb.sum() < 3 или inwin.sum() == 0 — ValueError;
       k, c = np.polyfit(e[sb], col[sb], 1)
       return float(np.sum(col[inwin] - (k * e[inwin] + c)))
   В docstring: подложка линейная, МНК по двум боковым полосам
   [E0−3,0F, E0−1,6F] и [E0+1,6F, E0+3,0F]; окно E0 ± 1,25F.

8. peak_row(table, E0, tol)
   table — список dict строк таблицы пиков прибора (ключ "energy_keV").
   Вернуть строку с наименьшим |energy_keV − E0|, если это расстояние <= tol,
   иначе None. Пустая таблица — None.

9. selftest()
   Самопроверки, печать "SELFTEST OK" и return 0, либо печать причины и return 1:
   T1 rebin_counts на тождественной шкале возвращает исходные counts
      (np.allclose); сумма сохраняется при переносе на шкалу, сдвинутую на
      0,37 канала внутри диапазона (np.isclose сумм, допуск 1e-9 отн., когда
      приёмник целиком накрывает источник).
   T2 rebin_counts на убывающей шкале бросает ValueError.
   T3 net_area: e = np.arange(0, 400, 1.0) + 0.5; col = 1000·гаусс(центр 200,
      sigma 10/2.3548, нормирован на единичную площадь) + (5 + 0.02·e);
      net_area(e, col, 200, 10) должен дать 1000 с точностью 2 %
      (окно ±1,25 ПШПВ содержит ~99,9 % гаусса).
   T4 amplitude_errors: A = единичная матрица 3×3, coef = [1, 0, 2] →
      [1, nan, 1].
   T5 zone бросает ValueError, если в окно не попал ни один канал
      (подать r с e = np.array([100., 200.]) и остальными полями нужной длины).

10. if __name__ == "__main__": sys.exit(selftest())
