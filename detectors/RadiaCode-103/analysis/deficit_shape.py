# -*- coding: utf-8 -*-
"""ФОРМА недостающего фона: разностный спектр «измерение минус модель».

bg_budget.py даёт отношение модель/измерение по широким полосам — видно, ЧТО
не сходится, но не видно, КАКОГО спектра не хватает. Здесь тот же бюджет
считается в узких окнах (50 кэВ) и с разбивкой модели по слагаемым
K/Ra/Th/мюоны/Pb-210, а печатается РАЗНОСТЬ (измерение минус модель).
Форма остатка — портрет недостающего источника: плоская плотность дефицита
означает жёсткий спектр, спад — мягкий; постоянное отношение дефицита к
какой-то одной компоненте означает, что не хватает просто её масштаба, а не
нового источника.

Ничего не подгоняется: ни одного коэффициента, подобранного по измерению.
Только арифметика над готовыми массивами модели и измеренным спектром."""
import math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths, rcspec, read_rcxml
import fit_lines as fl
import run_bg_shield as bg

# --- измерение --------------------------------------------------------------
meas = str(paths.measured("RadiaCode-103") / "Фон домик 23 дня.xml")
smp = read_rcxml.read(meas)[0]
# P-007: последний канал — переполнение (всё выше шкалы), сравнивать его как
# энергетический интервал нельзя.
# P-018: НЕ smp.energy — штатная шкала прибора врёт (rms 6,1 кэВ на якорях).
# Та же калибровка, что в bg_budget.CAL_BG: реперы Pb-210 46,539 и суммарный
# XRF-комплекс свинца (эффективный центроид 76,68) + K-40 1460,82 + Tl-208
# 2614,51, квадратичная по 4 якорям, rms 0,68 кэВ.
CAL_BG = [0.888094236, 2.484743967, 0.000221908]
_ch = np.arange(len(smp.counts))
e_meas = sum(c * _ch ** i for i, c in enumerate(CAL_BG))[:-1]
cps_meas = (smp.counts / smp.live)[:-1]
# Статистика измерения на канал: пуассоновская, sqrt(N)/T. Пустой канал
# считаем как 1 отсчёт, иначе получилась бы нулевая погрешность там, где её
# быть не может.
d_cps = (np.sqrt(np.maximum(smp.counts, 1.0)) / smp.live)[:-1]

# --- модель: ЕРН + мюоны РАЗДЕЛЬНО из b_room_prior.csv ----------------------
# Формат: E_keV,cps_total,cps_K,cps_Ra,cps_Th,cps_mu; мюоны уже в total
# (#SHIELD-13). Здесь total не используется вовсе — берутся поля 2..5 по
# отдельности, чтобы видеть вклад каждого слагаемого в форму остатка.
g = os.path.join(bg.out_dir(50.0), "b_room_prior.csv")
if not os.path.exists(g):
    raise SystemExit("нет %s — сначала run_bg_shield.py" % g)
full = {k: np.zeros(rcspec.NBINS) for k in ("K", "Ra", "Th", "mu")}
for line in open(g, encoding="utf-8"):
    if line.startswith("#") or not line[:1].isdigit():
        continue
    p = line.split(",")
    i = int(float(p[0]))
    if 0 <= i < rcspec.NBINS:
        full["K"][i] = float(p[2])
        full["Ra"][i] = float(p[3])
        full["Th"][i] = float(p[4])
        full["mu"][i] = float(p[5])

# --- Pb-210 в свинце домика -------------------------------------------------
# Прогоны pbself: 1 распад Pb-210 (цепочка до Bi-210, P-010) на историю,
# равномерно по объёму ячейки. При активности a Бк/кг ячейки массой m идёт
# a*m распадов/с; N историй = T = N/(a*m) секунд, cps = fold(hist)*a*m/N.
# Массы ячеек — из геометрии слоя, сумма 210,3 кг сверена против
# massPb_kg=210.259 в шапке прогона (0,02 %).
PB_CELLS = {"sh0_Pb_bot": 35.5, "sh0_Pb_xhi": 54.6, "sh0_Pb_xlo": 54.6,
            "sh0_Pb_yhi": 32.8, "sh0_Pb_ylo": 32.8}
pb_per_bqkg = np.zeros(rcspec.NBINS)   # cps на 1 Бк/кг (вся защита)
for cell, mass in PB_CELLS.items():
    p = os.path.join(str(paths.build("RadiaCode-103")), "pb210_%s.csv" % cell)
    meta, h = rcspec.read_spec(p)
    pb_per_bqkg += rcspec.fold(h, "103") * mass / float(meta["N_primaries"])

# Литературная активность Pb-210 в обычном коммерческом свинце (GeMSE),
# НЕ подгонка: активность свинца ЭТОГО домика не измерена.
PB210_BQKG = 91.0
full["Pb"] = pb_per_bqkg * PB210_BQKG

# --- на канальную сетку прибора ---------------------------------------------
COMPS = ("K", "Ra", "Th", "mu", "Pb")
e_mod = np.arange(rcspec.NBINS) + 0.5
comp = {k: fl.rebin_model_to_meas(e_mod, full[k], e_meas) for k in COMPS}

# Окна по 50 кэВ. Ниже 100 кэВ картину определяет XRF-комплекс свинца, выше
# 1200 кэВ статистика измерения в окне шириной 50 кэВ уже мала.
WIN = [(lo, lo + 50) for lo in range(100, 1200, 50)]


def _div(a, b):
    """Отношение с nan вместо деления на ноль — таблица не должна падать."""
    return a / b if b else float("nan")


def main():
    print("измерение: %s" % os.path.basename(meas))
    print("модель: %s" % g)
    print("Pb-210: %.0f Бк/кг — литературная (GeMSE), НЕ подгонка" % PB210_BQKG)
    print("живое время измерения: %.2f сут (канал переполнения исключён, P-007)"
          % (smp.live / 86400))
    print("НИЧЕГО не подгоняется: ни одной нормировки на измерение.\n")

    rows = []
    for lo, hi in WIN:
        m = (e_meas >= lo) & (e_meas < hi)
        a = float(cps_meas[m].sum())
        parts = {k: float(comp[k][m].sum()) for k in COMPS}
        mod = sum(parts.values())
        dfc = a - mod
        err = math.sqrt(float((d_cps[m] ** 2).sum()))
        rows.append((lo, hi, a, mod, parts, dfc, err, dfc / (hi - lo)))

    # --- Таблица A: бюджет и дефицит по окнам -------------------------------
    print("ТАБЛИЦА A. Бюджет по окнам 50 кэВ (cps в окне)")
    print("%-13s %9s %9s %9s %9s %9s %9s %9s %9s %9s %7s %9s"
          % ("окно, кэВ", "измерено", "модель", "K", "Ra", "Th", "мю",
             "Pb210", "ДЕФИЦИТ", "погр.деф", "м/и", "cps/кэВ"))
    for lo, hi, a, mod, parts, dfc, err, dens in rows:
        print("%-13s %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e"
              " %7.3f %9.3e"
              % ("%d-%d" % (lo, hi), a, mod, parts["K"], parts["Ra"],
                 parts["Th"], parts["mu"], parts["Pb"], dfc, err,
                 _div(mod, a), dens))

    # --- Таблица B: форма дефицита ------------------------------------------
    print("\nТАБЛИЦА B. Плотность дефицита, нормированная на первое окно")
    print("Постоянство столбца xN означает ЖЁСТКИЙ (плоский) источник,"
          " спад — МЯГКИЙ.")
    print("%-13s %9s %9s" % ("окно, кэВ", "cps/кэВ", "xN"))
    dens0 = rows[0][7]
    for lo, hi, a, mod, parts, dfc, err, dens in rows:
        print("%-13s %9.3e %9.3f"
              % ("%d-%d" % (lo, hi), dens, _div(dens, dens0)))

    # --- Таблица C: на какую компоненту похож дефицит ------------------------
    print("\nТАБЛИЦА C. Отношение дефицита к каждой компоненте модели")
    print("Если столбец ПОСТОЯНЕН по энергии — дефицит имеет ФОРМУ этой"
          " компоненты,")
    print("то есть не хватает просто её масштаба, а не нового источника.")
    print("%-13s %9s %9s %9s %9s %9s"
          % ("окно, кэВ", "деф/K", "деф/Ra", "деф/Th", "деф/mu", "деф/Pb"))
    for lo, hi, a, mod, parts, dfc, err, dens in rows:
        print("%-13s %9.3f %9.3f %9.3f %9.3f %9.3f"
              % ("%d-%d" % (lo, hi),
                 _div(dfc, parts["K"]), _div(dfc, parts["Ra"]),
                 _div(dfc, parts["Th"]), _div(dfc, parts["mu"]),
                 _div(dfc, parts["Pb"])))

    tot_d = sum(r[5] for r in rows)
    tot_e = math.sqrt(sum(r[6] ** 2 for r in rows))
    print("\nсуммарный дефицит 100-1200 кэВ: %.4e +- %.1e cps" % (tot_d, tot_e))
    print("живое время измерения: %.2f сут" % (smp.live / 86400))

    out = os.path.join(rcspec.RESULTS, "deficit_shape.csv")
    with open(out, "w", encoding="utf-8") as f:
        f.write("E_lo,E_hi,cps_meas,cps_model,cps_deficit,d_stat,"
                "dens_per_keV\n")
        for lo, hi, a, mod, parts, dfc, err, dens in rows:
            f.write("%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e\n"
                    % (lo, hi, a, mod, dfc, err, dens))
    print("CSV: %s" % out)


if __name__ == "__main__":
    main()
