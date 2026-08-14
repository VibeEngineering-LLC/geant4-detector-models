# -*- coding: utf-8 -*-
"""МДА по ISO 11929 для проекта свинцовой защиты — задача №5.

Формулы ПЕРЕНЕСЕНЫ в этот репозиторий (в geant4-detector-models их не было),
с юнит-тестом, сверяющим результат с независимой реализацией SpectraVibe
(`scripts/gamma/math/iso_11929_thresholds.py`, SPECTRAVIBE_ROOT) на общих
входах — чтобы не разъезжались две копии одной формулы (правило проекта,
см. план раздел 5).

КОНВЕНЦИЯ gross/bg (сверено с докстрингом SpectraVibe дословно, ISO
11929-1:2019 §5.4.3): u(n_net|y=0) = sqrt(N_gross + N_bg) — ДВЕ независимые
пуассоновские выборки одного и того же истинного фона (проба И отдельное
контрольное измерение фона, ОБЕ приведены к одному времени набора — правило
оператора от 12.08.2026, см. план раздел 0 "Общее методическое правило").
Для теоретического МДА (нет реального измерения, есть МОДЕЛЬНЫЙ прогноз
уровня фона) N_gross = N_bg = предсказанный счёт в окне за t_набора — то
есть при y=0 сигнала нет, и "проба" статистически неотличима от "фона".

    y* = k_(1-α) · sqrt(N_gross + N_bg) / (ε · I · m · t)
    y# = (k_(1-α) + k_(1-β)) · sqrt(N_gross + N_bg) / (ε · I · m · t)

При α=β=0.05 (KTA, симметрично): y# = 2·y*, тот же частный случай, что у
SpectraVibe (docstring, "y# ≈ 2·y*").

Запуск:
    python mda_shield.py test    — юнит-тест против SpectraVibe
    python mda_shield.py grid    — МДА(pb, t_набора) по сетке run_shield_grid.py
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

K_95 = 1.6449   # ISO 11929 Table A.1, α=β=0.05 — ТА ЖЕ константа, что SpectraVibe


def decision_threshold(gross_counts, bg_counts, efficiency, branching_ratio,
                       mass_kg, live_time_s, k=K_95):
    """y* [Бк/кг], ISO 11929-1:2019 §5.4.3. None при некорректных входах."""
    if efficiency <= 0 or branching_ratio <= 0 or mass_kg <= 0 or live_time_s <= 0:
        return None
    sigma_0 = math.sqrt(max(0.0, gross_counts) + max(0.0, bg_counts)) / (
        efficiency * branching_ratio * mass_kg * live_time_s)
    return k * sigma_0


def detection_limit(gross_counts, bg_counts, efficiency, branching_ratio,
                    mass_kg, live_time_s, k_alpha=K_95, k_beta=K_95):
    """y# [Бк/кг], ISO 11929-1:2019 §5.4.4, низкостатистическое приближение
    ũ(y#)≈ũ(0) (тот же частный случай, что у SpectraVibe — обоснование
    приближения там же: остаток фона доминирует над искомым сигналом на
    уровне порога, что для нашей задачи — предсказание МДА, а не проверка
    конкретного найденного пика — выполняется по построению)."""
    if efficiency <= 0 or branching_ratio <= 0 or mass_kg <= 0 or live_time_s <= 0:
        return None
    sigma_0 = math.sqrt(max(0.0, gross_counts) + max(0.0, bg_counts)) / (
        efficiency * branching_ratio * mass_kg * live_time_s)
    return (k_alpha + k_beta) * sigma_0


def time_to_target_mda(target_bqkg, bg_cps, efficiency, branching_ratio,
                       mass_kg, k_alpha=K_95, k_beta=K_95):
    """Обратная задача: t [с] такое, что y#(t) = target_bqkg.

    y#(t) = (k_a+k_b)*sqrt(2*bg_cps*t)/(eps*I*m*t) = (k_a+k_b)*sqrt(2*bg_cps/t)/(eps*I*m)
    -> t = [(k_a+k_b)*sqrt(2*bg_cps) / (target_bqkg*eps*I*m)]^2
    ЧИСЛЕННАЯ инверсия (не приближение «фон >> сигнала», как требует план,
    раздел 5) — но в НАШЕЙ формуле (gross=bg=bg_cps*t) сигнал изначально не
    входит в дисперсию вовсе (низкостатистическое приближение УЖЕ встроено
    в саму y#, см. docstring detection_limit), поэтому обратная задача здесь
    решается АНАЛИТИЧЕСКИ точно относительно ЭТОЙ формулы, без дополнительного
    приближения поверх неё.
    """
    if bg_cps <= 0 or efficiency <= 0 or branching_ratio <= 0 or mass_kg <= 0 or target_bqkg <= 0:
        return None
    denom = target_bqkg * efficiency * branching_ratio * mass_kg
    t = ((k_alpha + k_beta) * math.sqrt(2.0 * bg_cps) / denom) ** 2
    return t


# ─────────────────────────────────────────────────────────────────────────
def _spectravibe_module():
    # Личный путь НЕ хардкодить (правило репозитория, common/py/paths.py) —
    # require_spectravibe() уже даёт понятную ошибку вместо хардкода.
    root = str(paths.require_spectravibe("юнит-тест mda_shield.py против SpectraVibe"))
    p = os.path.join(root, "scripts", "gamma", "math")
    if p not in sys.path:
        sys.path.insert(0, p)
    import iso_11929_thresholds as sv  # noqa
    return sv


def unit_test():
    """Сверка с независимой реализацией SpectraVibe на общих входах — не
    просто «числа похожи», а np.isclose с жёстким rtol, на НЕСКОЛЬКИХ точках
    (низкая/высокая статистика, разные ε/I/m/t), как требует правило проекта
    (interpretation-is-hypothesis: формула обязана сходиться на нескольких
    точках, не на одной)."""
    sv = _spectravibe_module()

    cases = [
        # gross, bg,   eff,      I,     m_kg,   t_s
        (83.41, 83.41, 6.635e-4, 0.851, 0.1001, 3600.0),
        (10.0,  10.0,  6.635e-4, 0.851, 0.1001, 3600.0),
        (5000.0, 4800.0, 1.3397e-5, 0.1055, 0.1001, 612250.0),
        (0.0,   0.0,   6.635e-4, 0.851, 0.1001, 3600.0),
        (1.0,   500.0, 3.6491e-4, 0.4549, 0.1001, 86400.0),
    ]
    n_ok = 0
    for gross, bg, eff, I, m, t in cases:
        ys_mine = decision_threshold(gross, bg, eff, I, m, t)
        yd_mine = detection_limit(gross, bg, eff, I, m, t)
        ys_sv = sv.decision_threshold(gross, bg, eff, I, m, t, regime="KTA")
        yd_sv = sv.detection_limit(gross, bg, eff, I, m, t, regime="KTA")
        ok_s = (ys_mine is None and ys_sv is None) or (
            ys_mine is not None and ys_sv is not None
            and np.isclose(ys_mine, ys_sv, rtol=1e-9))
        ok_d = (yd_mine is None and yd_sv is None) or (
            yd_mine is not None and yd_sv is not None
            and np.isclose(yd_mine, yd_sv, rtol=1e-9))
        status = "OK" if (ok_s and ok_d) else "РАСХОДИТСЯ"
        print("gross=%-7.2f bg=%-7.2f eff=%-10.4e I=%-6.4f m=%-6.4f t=%-9.1f "
              "y*mine=%s y*sv=%s y#mine=%s y#sv=%s  [%s]"
              % (gross, bg, eff, I, m, t, ys_mine, ys_sv, yd_mine, yd_sv, status))
        if ok_s and ok_d:
            n_ok += 1
        else:
            raise SystemExit("РАСХОЖДЕНИЕ с SpectraVibe на случае gross=%s bg=%s" % (gross, bg))
    print("\n%d/%d случаев сошлись с SpectraVibe (rtol=1e-9)." % (n_ok, len(cases)))

    # y# = 2*y* при alpha=beta=0.05 (KTA) — явная сверка частного случая,
    # тот же, что задокументирован в SpectraVibe docstring.
    ys = decision_threshold(100.0, 100.0, 6.635e-4, 0.851, 0.1001, 3600.0)
    yd = detection_limit(100.0, 100.0, 6.635e-4, 0.851, 0.1001, 3600.0)
    assert np.isclose(yd, 2 * ys, rtol=1e-9), "y# != 2*y* при alpha=beta=0.05"
    print("y# = 2*y* при alpha=beta=0.05 — подтверждено.")


# ─────────────────────────────────────────────────────────────────────────
EPS_P_662 = 6.635e-4     # m200/organic 0.50, план "Опорные числа" (±2.7%)
I_CS137 = 0.851
MASS_KG = 0.1001
K_BASELINE_BQKG = 250.0

# N_MIN — эвристический порог, НЕ формальная цитата ISO 11929 (там нет явно
# заданного минимума). Формула y*/y# использует гауссово приближение
# Пуассона (k·sqrt(N) как квантиль) — приближение тем хуже, чем меньше N;
# инженерное правило (Currie-школа счётной статистики, общепринятое, но не
# численно фиксированное стандартом) — не доверять точечному значению ниже
# ~20 отсчётов. Найдено 12.08.2026: вопрос оператора "реальны ли МДА<50 при
# таком cps" — ячейки с N<N_MIN больше НЕ печатаются как число, только как
# явный прочерк с самим N, чтобы не выдавать оптимистичную точку за факт.
N_MIN = 20.0

BUILD = str(paths.build("RadiaCode-103"))
GRID_DIR = os.path.join(BUILD, "shield_grid")


def read_cps_window(path, lo=632.2, hi=691.2):
    e, c = [], []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line[:1].isdigit():
            continue
        a, b = line.split(",")
        e.append(float(a))
        c.append(float(b))
    e, c = np.array(e), np.array(c)
    m = (e >= lo) & (e < hi)
    return c[m].sum()


def grid():
    """МДА(pb, t_набора) по уже посчитанной сетке run_shield_grid.py."""
    pb_list = []
    for fn in os.listdir(GRID_DIR):
        if fn.startswith("b_room_gamma_pb") and fn.endswith(".csv"):
            pb_list.append(float(fn[len("b_room_gamma_pb"):-4]))
    pb_list.sort()
    if not pb_list:
        raise SystemExit("нет b_room_gamma_pb*.csv в %s — сначала run_shield_grid.py" % GRID_DIR)

    TIMES = [("1 час", 3600.0), ("4 часа", 14400.0), ("1 сутки", 86400.0),
             ("7 суток", 604800.0)]

    print("=== МДА(pb, t) по Cs-137 (662 кэВ), ISO 11929 KTA, k=%.4f ===" % K_95)
    print("(gross=bg=предсказанный счёт в окне — теоретическое МДА, не измерение)")
    print("(ячейка с N<%.0f отсчётов — гауссово приближение Пуассона ненадёжно," % N_MIN)
    print(" точка НЕ печатается, см. N_MIN в коде)\n")
    header = "%6s" % "pb,мм"
    for name, _ in TIMES:
        header += " %14s" % name
    print(header)
    t_min_row = {}
    for pb in pb_list:
        b_gamma = read_cps_window(os.path.join(GRID_DIR, "b_room_gamma_pb%.0f.csv" % pb))
        mu_path = os.path.join(GRID_DIR, "b_room_mu_pb%.0f.csv" % pb)
        b_mu = read_cps_window(mu_path) if os.path.exists(mu_path) else 0.0
        s_k40_path = os.path.join(GRID_DIR, "s_K40_per_bqkg_pb%.0f.csv" % pb)
        s_k40 = (K_BASELINE_BQKG * read_cps_window(s_k40_path)
                if os.path.exists(s_k40_path) else 0.0)
        bg_cps_total = b_gamma + b_mu + s_k40   # ВЕСЬ фон под окном при Cs137=0
        t_min_row[pb] = N_MIN / bg_cps_total if bg_cps_total > 0 else float("inf")

        row = "%6.0f" % pb
        for _name, t in TIMES:
            n = bg_cps_total * t
            if n < N_MIN:
                row += " %14s" % ("N=%.1f<%.0f" % (n, N_MIN))
                continue
            yd = detection_limit(n, n, EPS_P_662, I_CS137, MASS_KG, t)
            row += " %14.2f" % yd if yd is not None else " %14s" % "-"
        print(row)

    print("\n=== минимальное время набора для достоверной МДА (N>=%.0f) ===" % N_MIN)
    for pb in pb_list:
        tm = t_min_row[pb]
        print("  pb=%3.0f мм: t_min = %.0f с (%.1f мин)" % (pb, tm, tm / 60.0))

    print("\n⚠️ Использует B_room_ПОЛНЫЙ (гамма+мюон, задачи №4/№6, replay")
    print("   biased по объёму кристалла — задача №13, закрыта 12.08.2026) —")
    print("   БЕЗ поправки на дефицит континуума №12 (~25-35%, недооценка")
    print("   фона, МДА здесь ОПТИМИСТИЧНЕЕ реальности). a_mu не проверен")
    print("   независимо (задача №3). Сетка — 9 точек (5..100мм, решение")
    print("   оператора 12.08), без 2D Cd/Cu (задача №4 продолжается).")
    print("   σ — ТОЛЬКО пуассоновская счётная статистика: дрейф усиления/")
    print("   разрешения прибора на многосуточном наборе и нестационарность")
    print("   фона помещения за неделю НЕ смоделированы (вопрос оператора")
    print("   12.08) — реальная МДА может быть хуже показанной.")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "test"
    if cmd == "test":
        unit_test()
    elif cmd == "grid":
        grid()
    else:
        raise SystemExit("команды: test | grid")
