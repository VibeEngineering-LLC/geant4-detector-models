# -*- coding: utf-8 -*-
"""Плотность мюонного отклика в узком окне — без поканального счёта.

ЗАЧЕМ. Мюонный отклик в рабочем окне 632-691 кэВ набирает единицы отсчётов
(5-13 при 500 000 первичных, найдено 13.08.2026 — см. план, «Находка 2»),
и поканальные числа там статистический шум, а не сигнал: скачок 3,48e-4
(pb=20) -> 8,65e-4 (pb=30) объяснялся 5 против 12 отсчётов. Прямое лечение
статистикой не проходит: при насыщенном диске (находка 3) даже 10^7
первичных дают в окне порядка 19 отсчётов.

ЧТО ДЕЛАЕМ. Мюонный спектр депозита в диапазоне 0-3200 кэВ — ГЛАДКИЙ
континуум без линий (распределение длин хорд через кристалл, свёрнутое с
dE/dx; пик полного пролёта лежит выше, в области 4-10 МэВ). Поэтому
плотность в узком окне оценивается по ВСЕЙ накопленной статистике спектра,
а не по попавшим именно в окно каналам: бины укрупняются адаптивно до
заданного минимума отсчётов, дальше плотность интерполируется.

ЧЕГО НЕ ДЕЛАЕМ. Никакой параметрической формы не навязывается — это не
подгонка модели, а перегруппировка тех же отсчётов. Интеграл сохраняется
тождественно (проверяется юнит-тестом self_test()).

ГРАНИЦА ПРИМЕНИМОСТИ. Приём законен ровно потому, что спектр гладкий на
масштабе окна. К спектрам С ЛИНИЯМИ (гамма-компоненты B_room, S_K40,
S_Cs137) НЕ применять — там укрупнение размажет пик по континууму.

Запуск: python mu_smooth.py <spectrum.csv> [<E_lo> <E_hi>]
        python mu_smooth.py test
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

MIN_COUNTS_PER_BIN = 25.0     # эвристика, не норматив: при 25 отсчётах
                              # пуассоновская погрешность бина 20 %


def read_spectrum(path):
    """-> (e_kev, counts, n_primaries). Читает CSV shieldrun (E_keV,counts)."""
    e, c, n = [], [], None
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                n = int(line.split("=")[1].strip())
            continue
        if not line[:1].isdigit():
            continue
        a, b = line.split(",")
        e.append(float(a))
        c.append(float(b))
    return np.array(e), np.array(c), n


def adaptive_edges(e, c, min_counts=MIN_COUNTS_PER_BIN):
    """Границы укрупнённых бинов: набираем каналы, пока не наберётся
    min_counts. Последний бин может недобрать — он сливается с предыдущим."""
    if len(e) == 0:
        return np.array([]), []
    step = e[1] - e[0] if len(e) > 1 else 1.0
    lo = e[0] - step / 2.0
    edges, acc, groups, cur = [lo], 0.0, [], []
    for i in range(len(e)):
        cur.append(i)
        acc += c[i]
        if acc >= min_counts:
            edges.append(e[i] + step / 2.0)
            groups.append(cur)
            acc, cur = 0.0, []
    if cur:                      # хвост недобрал — приклеиваем к последнему
        if groups:
            groups[-1].extend(cur)
            edges[-1] = e[cur[-1]] + step / 2.0
        else:
            groups.append(cur)
            edges.append(e[cur[-1]] + step / 2.0)
    return np.array(edges), groups


def density(e, c, min_counts=MIN_COUNTS_PER_BIN):
    """-> (e_centres, dens, dens_sd) — плотность отсчётов на кэВ и её
    пуассоновская погрешность, по адаптивно укрупнённым бинам."""
    edges, groups = adaptive_edges(e, c, min_counts)
    if len(groups) == 0:
        return np.array([]), np.array([]), np.array([])
    cen, den, sd = [], [], []
    for k, g in enumerate(groups):
        n = float(sum(c[i] for i in g))
        w = edges[k + 1] - edges[k]
        cen.append(0.5 * (edges[k] + edges[k + 1]))
        den.append(n / w)
        sd.append(np.sqrt(max(n, 1.0)) / w)
    return np.array(cen), np.array(den), np.array(sd)


_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy 2.x переименовал


def window_counts(e, c, e_lo, e_hi, min_counts=MIN_COUNTS_PER_BIN,
                  _checked=False):
    """Оценка числа отсчётов в окне [e_lo, e_hi) по СГЛАЖЕННОЙ плотности.

    ⚠️ `c` — СЫРЫЕ ОТСЧЁТЫ, не cps и не плотность. Функция НЕЛИНЕЙНА по
    масштабу входа: укрупнение набирает бины до `min_counts` отсчётов, и у
    входа, домноженного на константу, группировка получится другой. Порядок
    работы — сначала окно на отсчётах, потом нормировка:
        val, sd = window_counts(e, folded_counts, lo, hi)
        cps = val * a_mu / n_primaries          # а НЕ window_counts(e, cps...)
    Поймано 13.08.2026: подача cps дала 1,33e-3 вместо 4,97e-4 (расхождение
    2,7×) и погрешность 460 % — при том, что калибровка на тех же данных,
    но на сырых отсчётах, сходилась с измерением.

    -> (value, sd). Плотность интерполируется линейно по центрам укрупнённых
    бинов и интегрируется по окну.

    Погрешность НЕ интегрируется вместе с плотностью (складывать σ линейно
    означало бы считать соседние бины полностью коррелированными и завышать
    результат). Берётся относительная: оценка опирается на N_used — сырые
    отсчёты тех укрупнённых бинов, которые перекрывают окно, — поэтому
    σ(N_win) = N_win / sqrt(N_used). Это и есть выигрыш приёма: N_used на
    порядок-два больше числа отсчётов, попавших в само окно.
    """
    total = float(np.sum(c))
    if not _checked and 0 < total < min_counts:
        raise ValueError(
            "window_counts: суммарный вход %.4g меньше порога укрупнения %.0f — "
            "похоже, поданы cps/плотность вместо СЫРЫХ ОТСЧЁТОВ (функция "
            "нелинейна по масштабу, см. докстринг). Считайте окно на отсчётах "
            "и нормируйте результат; если вход действительно счётный и просто "
            "очень бедный, передайте _checked=True." % (total, min_counts))
    edges, groups = adaptive_edges(e, c, min_counts)
    if len(groups) == 0:
        return 0.0, 0.0
    cen, den, _sd = density(e, c, min_counts)
    grid = np.linspace(e_lo, e_hi, 201)
    d = np.interp(grid, cen, den, left=den[0], right=den[-1])
    val = float(_trapz(d, grid))

    n_used = 0.0
    for k, g in enumerate(groups):
        if edges[k + 1] > e_lo and edges[k] < e_hi:      # бин перекрывает окно
            n_used += float(sum(c[i] for i in g))
    sd = val / np.sqrt(n_used) if n_used > 0 else float("inf")
    return val, float(sd)


def raw_window_counts(e, c, e_lo, e_hi):
    """Поканальный счёт в окне — для сравнения со сглаженным."""
    m = (e >= e_lo) & (e < e_hi)
    n = float(c[m].sum())
    return n, float(np.sqrt(max(n, 1.0)))


def self_test():
    """Интеграл обязан сохраняться: сумма по укрупнённым бинам = сумма по
    исходным каналам. Плюс проверка на заведомо гладком спектре."""
    ok = True
    rng = np.random.default_rng(20260813)

    e = np.arange(0.5, 3200.5, 1.0)
    true_dens = 3.0 * np.exp(-e / 1500.0)         # гладкая, без линий
    c = rng.poisson(true_dens).astype(float)

    edges, groups = adaptive_edges(e, c)
    n_grouped = sum(sum(c[i] for i in g) for g in groups)
    if abs(n_grouped - c.sum()) > 1e-9:
        print("ПРОВАЛ: интеграл не сохранён: %.6f против %.6f"
              % (n_grouped, c.sum())); ok = False
    else:
        print("ok: интеграл сохранён (%.0f отсчётов)" % c.sum())

    idx = np.concatenate([np.array(g) for g in groups])
    if len(idx) != len(e) or len(set(idx.tolist())) != len(e):
        print("ПРОВАЛ: разбиение не покрывает все каналы ровно один раз")
        ok = False
    else:
        print("ok: разбиение покрывает все %d каналов ровно один раз" % len(e))

    lo, hi = 632.2, 691.2
    exact = float(_trapz(3.0 * np.exp(-np.linspace(lo, hi, 201) / 1500.0),
                         np.linspace(lo, hi, 201)))
    sm, _ = window_counts(e, c, lo, hi)
    raw, raw_sd = raw_window_counts(e, c, lo, hi)
    print("окно %.1f-%.1f: истина %.2f | сглажено %.2f | поканально %.2f ± %.2f"
          % (lo, hi, exact, sm, raw, raw_sd))
    if abs(sm - exact) > 3.0 * max(raw_sd, 1.0):
        print("ПРОВАЛ: сглаженная оценка ушла дальше 3 сигма от истины")
        ok = False
    else:
        print("ok: сглаженная оценка в пределах 3 сигма от истины")

    # разреженный случай — ради чего всё и затевалось
    c_sparse = rng.poisson(true_dens * 0.02).astype(float)
    sm_s, _ = window_counts(e, c_sparse, lo, hi)
    raw_s, raw_s_sd = raw_window_counts(e, c_sparse, lo, hi)
    exact_s = exact * 0.02
    print("разреженный (x0,02): истина %.3f | сглажено %.3f | поканально %.3f ± %.3f"
          % (exact_s, sm_s, raw_s, raw_s_sd))
    print("  всего отсчётов в спектре: %.0f" % c_sparse.sum())

    print("ИТОГ:", "ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ" if ok else "ЕСТЬ ПРОВАЛЫ")
    return ok


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        sys.exit(0 if self_test() else 1)
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    path = sys.argv[1]
    lo = float(sys.argv[2]) if len(sys.argv) > 3 else 632.2
    hi = float(sys.argv[3]) if len(sys.argv) > 3 else 691.2
    e, c, n = read_spectrum(path)
    sm, sm_sd = window_counts(e, c, lo, hi)
    raw, raw_sd = raw_window_counts(e, c, lo, hi)
    print("%s" % os.path.basename(path))
    print("  N_primaries = %s, всего отсчётов %.0f" % (n, c.sum()))
    print("  окно %.1f-%.1f кэВ: поканально %.1f ± %.1f | сглажено %.2f ± %.2f"
          % (lo, hi, raw, raw_sd, sm, sm_sd))
