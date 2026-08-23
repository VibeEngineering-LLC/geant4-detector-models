# -*- coding: utf-8 -*-
"""Подгонка K-40/Ra-226/Th-232/мюон ПО ЛИНИЯМ — замена вырожденной подгонки
континуума в fit_room_field.py (см. план, разделы «Одним множителем не
обойтись» и «Подгонка 4 амплитуд»: корреляция K/Ra/Th по 7 широким полосам
1,000 — комптоновский континуум неразличим по родительскому нуклиду там, где
сосредоточена вся статистика измерения, поэтому подгонка континуума в любом
варианте (равный вес/1/sqrt(cps)/1/sqrt(counts)) даёт три разных ответа).

МЕТОД. Площадь пика за вычетом ЛИНЕЙНОЙ подложки по двум боковым окнам — тот
же геометрический приём, что уже в peaks.py для Cs-137/K-40 пробы, — считается
ОТДЕЛЬНО на измеренном спектре (реальная смесь) и на модельной ЕДИНИЧНОЙ
(1 Бк/кг) кривой серии. Каждая кривая wallfield с gSeries=S по построению
несёт ТОЛЬКО линии этой серии (см. wallfield.cc: gSeries>=0 -> амплитуда
чужих линий = 0) — модельная площадь на "чужой" линии физический нуль, на
"своей" — чистая, без примеси. Амплитуда серии = площадь(измерено, cps) /
площадь(модель, cps на 1 Бк/кг), для Ra/Th — несколько линий, усреднение
весом 1/sigma^2, разброс между линиями печатается явно (диагностика метода,
не то, что можно спрятать усреднением).

Мюон — континуум, не линия: окно ВЫШЕ самой жёсткой линии Th-232 (2614,5),
где K/Ra/Th по построению (GPS Arb — один квант на событие, без каскадного
суммирования) не дают вклада; из измеренного счёта в этом окне вычитается
ПРЕДСКАЗАНИЕ по уже подобранным a_K/a_Ra/a_Th (должно быть ~0, это проверка,
не допущение), остаток делится на модельный отклик мюонной кривой в том же
окне.

Зависит от уже посчитанного: results/m200/background/bg_cyl_field_{K,Ra,Th}.csv
(rcspec.rdir) и wf_{K,Ra,Th}.csv + cosmicmu.csv (paths.build) — все читаются
существующими функциями fit_room_field.py, здесь не дублируются.

Запуск: python fit_lines.py
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for _d in ("analysis", "drivers"):
    _p = os.path.join(HERE, "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths  # noqa: E402
import rcspec  # noqa: E402
import read_rcxml  # noqa: E402
import fit_room_field as frf  # noqa: E402  — переиспользуем read_wallfield/read_cosmicmu/пути

# Диагностические линии — ТЕ ЖЕ энергии и выходы, что в LINES из wallfield.cc
# (сверено построчно 12.08.2026, не переписано на глаз): K-40 одна линия,
# Ra-226 и Th-232 — по несколько, для явной проверки согласия между ними.
DIAG_LINES = {
    "K":  [1460.8],
    "Ra": [609.3, 1764.5],
    "Th": [583.2, 911.2, 2614.5],
}
MODEL = "103"


def line_net_area(e, y, E0, nsig=2.5, gap=1.6, model=MODEL):
    """Площадь пика в (e,y) за вычетом линейной подложки по боковым окнам.

    Геометрия окон — БУКВАЛЬНО та же, что peaks.peak_area (тот же nsig=2.5,
    gap=1.6, тот же коэффициент покрытия 0.9876) — не независимо изобретённая,
    чтобы результат был сравним с уже принятой в проекте конвенцией.
    var = gross + |cont| — тот же приближённый (не строго propagated из
    отдельных дисперсий подложки) учёт, что в peaks.py; годится для порядка
    величины ошибки амплитуды, не для метрологического протокола.
    """
    # Ширина ДЛЯ ОКНА — из ИЗМЕРЕНИЯ (rcspec.fwhm_win), не из модельной кривой:
    # та выше 1 МэВ завышает, окно 2614 растягивалось за край шкалы и давало
    # ОТРИЦАТЕЛЬНОЕ нетто (21.08, указание оператора «ширину для окон из
    # измерения»). Свёртка модели по-прежнему на rcspec.fwhm — другая задача.
    sigma = rcspec.fwhm_win(E0) / 2.35482
    mp = (e >= E0 - nsig * sigma) & (e <= E0 + nsig * sigma)
    ml = (e >= E0 - (gap + 2.2) * sigma) & (e <= E0 - gap * sigma * 1.6)
    mh = (e >= E0 + gap * sigma * 1.6) & (e <= E0 + (gap + 2.2) * sigma)
    if mp.sum() < 1 or ml.sum() < 2 or mh.sum() < 2:
        return None
    yl, yh = y[ml].mean(), y[mh].mean()
    xl, xh = e[ml].mean(), e[mh].mean()
    base = yl + (yh - yl) * (e[mp] - xl) / (xh - xl)
    gross = y[mp].sum()
    cont = base.sum()
    area = (gross - cont) / 0.9876
    var = gross + abs(cont)
    return dict(area=area, sd=math.sqrt(max(var, 0.0)) / 0.9876, sigma=sigma,
               lo=e[mp][0], hi=e[mp][-1], cont=cont, gross=gross)


def model_cps_curve(series):
    """-> (e_model, cps) для серии K/Ra/Th: та же нормировка rate/t_run, что
    fit_room_field.fit() (тождество Ф=4N/S), свёрнутая с разрешением."""
    c = frf.CYL_M200
    r, hz = c["r"] / 10, 0.5 * (c["z1"] - c["z0"]) / 10
    area = 2 * math.pi * r * (r + 2 * hz)
    e_flu, flu = frf.read_wallfield(frf.wf_csv_path(series))
    fluence_total = flu.sum()
    rate = fluence_total * area / 4.0
    meta, hist_counts = rcspec.read_spec(frf.bg_csv_path(series))
    n = int(meta["N_primaries"])
    t_run = n / rate
    cps = rcspec.fold(hist_counts, MODEL) / t_run
    counts_folded = rcspec.fold(hist_counts, MODEL)   # для Пуассон-sigma МК, в счётах
    e_model = np.arange(len(cps)) + 0.5
    return e_model, cps, e_model, counts_folded, t_run


def rebin_model_to_meas(e_model, cps_model, e_meas):
    """Модель (1 кэВ/канал, ПЛОТНОСТЬ cps/кэВ) -> real-канальная сетка
    измерения (2,4..3,2 кэВ/канал у RC-103, 1024 канала, КВАДРАТИЧНАЯ
    калибровка — не равномерная!).

    НАЙДЕНО 12.08.2026: точечная интерполяция np.interp(e_meas, e_model,
    cps_model) (использовалась и в fit_room_field.py, и в первой версии
    "проверки" ниже в этом файле) берёт модельную плотность ТОЛЬКО в одной
    точке на канал измерения, а не интегрирует по всей ширине реального
    канала — систематически НЕДОСЧИТЫВАЕТ модель во столько раз, во сколько
    реальный канал шире модельного (здесь ~2,4-3,2x). Объясняет и провал
    "проверки" (модель/изм ~0.28 вместо ~1), и, ретроспективно, разъезд
    подгонки по континууму в fit_room_field.py (NNLS раздувал K до 1065,
    компенсируя систематический недосчёт модели, а не находя физический
    ответ). Правильно — сохраняющий поток ремешок (cumulative-sum), тот же
    приём, что peaks.rebin_to(), обобщённый на модельную сетку.
    """
    cum = np.concatenate(([0.0], np.cumsum(cps_model)))
    edges_model = np.arange(len(cps_model) + 1, dtype=float)   # 0,1,2,...  (кэВ, 1 кэВ/бин)

    def edges_of(e):
        return np.concatenate(([e[0] - (e[1] - e[0]) / 2],
                               (e[:-1] + e[1:]) / 2,
                               [e[-1] + (e[-1] - e[-2]) / 2]))

    edges_meas = edges_of(e_meas)
    cum_at_meas = np.interp(edges_meas, edges_model, cum, left=0.0, right=cum[-1])
    return np.diff(cum_at_meas)


def model_mu_curve():
    if not os.path.exists(frf.COSMICMU_CSV):
        return None
    meta_mu, hist_mu = frf.read_cosmicmu(frf.COSMICMU_CSV)
    n_mu = int(meta_mu["N_primaries"])
    cps_per_primary = rcspec.fold(hist_mu, MODEL) / n_mu
    e_model = np.arange(len(cps_per_primary)) + 0.5
    return e_model, cps_per_primary


def fit_series_amplitude(series, e_meas, cps_meas_area_fn):
    """Амплитуда a_S [Бк/кг] по своим линиям серии; возвращает (a, sd, per_line)."""
    e_model, cps_model, e_model_c, counts_model_c, t_run = model_cps_curve(series)
    per_line = []
    for E0 in DIAG_LINES[series]:
        rm = line_net_area(e_model, cps_model, E0)
        rmc = line_net_area(e_model_c, counts_model_c, E0)   # для МК-sigma в счётах
        rd = cps_meas_area_fn(E0)
        if rm is None or rd is None or rm["area"] <= 0:
            per_line.append((E0, None))
            continue
        sd_model_cps = (rmc["sd"] / t_run) if rmc is not None else 0.0
        a = rd["area"] / rm["area"]
        rel2 = (rd["sd"] / rd["area"]) ** 2 if rd["area"] else 0.0
        rel2 += (sd_model_cps / rm["area"]) ** 2 if rm["area"] else 0.0
        sd = abs(a) * math.sqrt(rel2)
        per_line.append((E0, dict(a=a, sd=sd, area_meas=rd["area"], area_model=rm["area"])))
    ok = [(E0, r) for E0, r in per_line if r is not None and r["sd"] > 0]
    if not ok:
        return None, None, per_line
    w = np.array([1.0 / r["sd"] ** 2 for _, r in ok])
    a_arr = np.array([r["a"] for _, r in ok])
    a = float((w * a_arr).sum() / w.sum())
    sd = float(1.0 / math.sqrt(w.sum()))
    return a, sd, per_line


# Правило контура: «образец и фон калибруются отдельно» (оператор, 18.08). Штатные
# коэффициенты XML этого файла дают невязку 5,38 кэВ rms на якорях Pb Ka1+K-40
# (см. #SHIELD-21, метод донора SpectraVibe — SNIP + якоря в каналах). Шкала здесь
# перестроена ТЕМИ ЖЕ коэффициентами, что применены в plot_bg_compare.py для того же
# измерения (не тот же ФАЙЛ — это «фон комнаты», в plot_bg_compare — «фон домика»,
# но обе перестроены одним методом на своих якорях).
# P-016: ЛИНЕЙНАЯ калибровка по ДВУМ якорям (Pb Ka1 + K-40) промахивалась на
# Tl-208 2614,5 на -116 кэВ — шкала RC-103 квадратична, экстраполяция за пределы
# якорей негодна. Третий якорь найден SNIP-поиском: ch=952,02, выступ 14,8 при
# подложке 5,9 (6 сигма); у донора та же линия на ch=953 — сходится.
# Три якоря: Pb Ka1 74,97 (ch 32,05, вес 0,5) + K-40 1460,82 (ch 558,20) +
# Tl-208 2614,51 (ch 952,02). rms штатной шкалы 6,10 кэВ -> новой 0,00.
CAL_ROOM = [-3.711311, 2.444318, 0.000321]   # E(ch) фона комнаты, кэВ (#SHIELD-21)


def fit():
    smp = read_rcxml.read(frf.MEASURED_BG)[0]
    _ch = np.arange(len(smp.counts))
    e_meas = sum(c * _ch ** i for i, c in enumerate(CAL_ROOM))
    counts_meas = smp.counts   # СЫРЫЕ отсчёты — для честной Пуассон-sigma

    def meas_area_cps(E0):
        r = line_net_area(e_meas, counts_meas, E0)
        if r is None:
            return None
        return dict(area=r["area"] / smp.live, sd=r["sd"] / smp.live,
                    lo=r["lo"], hi=r["hi"], cont=r["cont"] / smp.live)

    print("=== Подгонка K-40/Ra-226/Th-232 по диагностическим линиям ===")
    print("измеренный спектр: %s, живое %.0f с\n" % (frf.MEASURED_BG, smp.live))

    amps = {}
    for series in ("K", "Ra", "Th"):
        a, sd, per_line = fit_series_amplitude(series, e_meas, meas_area_cps)
        amps[series] = (a, sd)
        print("--- %s ---" % series)
        for E0, r in per_line:
            if r is None:
                print("  %8.1f кэВ: окно вне спектра или площадь <= 0" % E0)
                continue
            print("  %8.1f кэВ: измерено %.5f±%.5f cps, модель(1 Бк/кг) %.4e cps/(Бк/кг)"
                  " -> a=%.2f±%.2f Бк/кг"
                  % (E0, r["area_meas"], 0.0, r["area_model"], r["a"], r["sd"]))
        if a is None:
            print("  ИТОГО %s: ни одна линия не дала площадь > 0 — подгонка невозможна" % series)
        else:
            print("  ИТОГО %s (взвешенное среднее по своим линиям): a=%.2f ± %.2f Бк/кг"
                  % (series, a, sd))
        print()

    # --- мюонный континуум: окно выше самой жёсткой линии Th-232 -----------
    a_K, sd_K = amps["K"]
    a_Ra, sd_Ra = amps["Ra"]
    a_Th, sd_Th = amps["Th"]

    mu_curve = model_mu_curve()
    if mu_curve is None:
        print("[mu] cosmicmu.csv не найден — мюонная амплитуда не считается")
        a_mu = sd_mu = None
    else:
        e_mu, cps_mu_per = mu_curve
        E_hi_th = max(DIAG_LINES["Th"])
        sigma_hi = rcspec.fwhm(E_hi_th, MODEL) / 2.35482
        # 2 сигма, а не 4: после перекалибровки разрешения (FWHM(662)=9,83 % вместо
        # 8,4 %) сигма на 2614,5 выросла до ~86 кэВ, и отступ 4 сигма уводил нижнюю
        # границу окна (2958 кэВ) ВЫШЕ конца спектра (2833) — окно вырождалось в
        # пустое, мюонная амплитуда не определялась вовсе. При 2 сигма хвост
        # Tl-208 даёт <2,3 % площади линии, это допустимая примесь; проверка
        # «предсказание K+Ra+Th в окне ~0» ниже остаётся в силе и её видно.
        lo_win = E_hi_th + 2.0 * sigma_hi
        hi_win = min(e_meas.max(), e_mu.max()) - 5.0
        if hi_win <= lo_win:
            print("[mu] окно [%.0f, %.0f] пусто — спектр измерения короче, "
                  "чем нужно для мюонного окна" % (lo_win, hi_win))
            a_mu = sd_mu = None
        else:
            mwin = (e_meas >= lo_win) & (e_meas < hi_win)
            meas_win = counts_meas[mwin].sum() / smp.live
            meas_win_sd = math.sqrt(counts_meas[mwin].sum()) / smp.live

            def model_in_win(series, a):
                if a is None:
                    return 0.0
                e_model, cps_model, *_ = model_cps_curve(series)
                pred = rebin_model_to_meas(e_model, cps_model, e_meas[mwin])
                return a * pred.sum()

            pred_krath = model_in_win("K", a_K) + model_in_win("Ra", a_Ra) + model_in_win("Th", a_Th)
            net_mu_win = meas_win - pred_krath

            mu_model_win = rebin_model_to_meas(e_mu, cps_mu_per, e_meas[mwin]).sum()
            print("[mu] окно %.0f..%.0f кэВ: измерено %.5f±%.5f cps, "
                  "предсказание K+Ra+Th в этом окне %.5f cps (проверка ~0)"
                  % (lo_win, hi_win, meas_win, meas_win_sd, pred_krath))
            if mu_model_win > 0 and net_mu_win > 0:
                a_mu = net_mu_win / mu_model_win
                sd_mu = meas_win_sd / mu_model_win
                print("  чистый мюонный остаток %.5f cps / модельный отклик %.4e "
                      "-> a_mu=%.2f ± %.2f (условные ед. потока через диск)"
                      % (net_mu_win, mu_model_win, a_mu, sd_mu))
                r_disk_cm2 = math.pi * 7.0 ** 2
                expect_mu = 1.0 / 60.0 * r_disk_cm2
                print("  сверка порядка величины: поток PDG ~1/см²/мин * площадь диска "
                      "%.0f см² = %.2f 1/с (отношение к подобранному %.2f)"
                      % (r_disk_cm2, expect_mu, a_mu / expect_mu if expect_mu else float("nan")))
            else:
                print("  чистый мюонный остаток <= 0 (%.5f cps) — амплитуда мюона не определена "
                      "этим окном" % net_mu_win)
                a_mu = sd_mu = None

    # --- итоговая проверка по широким полосам (та же разбивка, что раньше) --
    print("\n=== Проверка: предсказание по подобранным линейным амплитудам, "
          "широкие полосы ===")
    BANDS = [(20, 100), (100, 300), (300, 600), (600, 750),
             (750, 1400), (1400, 1550), (1550, 2700), (2700, 3200)]
    e_model_K, cps_K, *_ = model_cps_curve("K")
    e_model_Ra, cps_Ra, *_ = model_cps_curve("Ra")
    e_model_Th, cps_Th, *_ = model_cps_curve("Th")
    pred_total = np.zeros_like(e_meas)
    for a, e_m, cps_m in ((a_K, e_model_K, cps_K), (a_Ra, e_model_Ra, cps_Ra),
                          (a_Th, e_model_Th, cps_Th)):
        if a is None:
            continue
        pred_total += a * rebin_model_to_meas(e_m, cps_m, e_meas)
    if a_mu is not None:
        e_mu, cps_mu_per = mu_curve
        pred_total += a_mu * rebin_model_to_meas(e_mu, cps_mu_per, e_meas)
    cps_meas_all = counts_meas / smp.live
    print("%-12s %10s %10s %10s" % ("полоса,кэВ", "измерено", "модель", "модель/изм"))
    for lo, hi in BANDS:
        m = (e_meas >= lo) & (e_meas < hi)
        ym, pm = cps_meas_all[m].sum(), pred_total[m].sum()
        print("%5d-%-6d %10.5f %10.5f %10.3f" % (lo, hi, ym, pm, pm / ym if ym > 0 else float("nan")))

    print("\nИТОГОВЫЕ АМПЛИТУДЫ:")
    print("  K   %8.2f ± %5.2f Бк/кг" % (a_K, sd_K) if a_K is not None else "  K   не определено")
    print("  Ra  %8.2f ± %5.2f Бк/кг" % (a_Ra, sd_Ra) if a_Ra is not None else "  Ra  не определено")
    print("  Th  %8.2f ± %5.2f Бк/кг" % (a_Th, sd_Th) if a_Th is not None else "  Th  не определено")
    if a_mu is not None:
        print("  mu  %8.2f ± %5.2f (условные ед.)" % (a_mu, sd_mu))
    return dict(K=(a_K, sd_K), Ra=(a_Ra, sd_Ra), Th=(a_Th, sd_Th),
               mu=(a_mu, sd_mu) if a_mu is not None else None)


if __name__ == "__main__":
    fit()
