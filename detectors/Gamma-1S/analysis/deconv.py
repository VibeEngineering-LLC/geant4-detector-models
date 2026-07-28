# -*- coding: utf-8 -*-
"""Связанная деконволюция группы линий: активность одна, площади привязаны.

ЗАЧЕМ. Отбор по разделённости оставил у Th-232 ОДНУ годную линию — 2614,5 кэВ.
Запаса нет: если с ней что-то не так, торий проверить нечем. Линии 238,6, 583,2 и
911,2 отброшены не потому, что плохи, а потому что оконный съём площади не умеет
их разделить. Умеет подгонка.

МЕТОД — по ЛСРМ, «Алгоритмические основы SpectraLine», формула 5.2-7. В χ² к
обычному члену по каналам добавляется СВЯЗЫВАЮЩИЙ:

    χ² = Σ_i w_i (y_i − континуум_i − Σ_k S_k·ψ(i, E_k, σ_k))²
       + Σ_k w_k (A·I_k·ε_k·t − S_k)²          ← площади не свободны

Здесь связь взята в предельно жёсткой форме: площади ПОДСТАВЛЕНЫ, а не
оштрафованы, то есть S_k = A·I_k·ε_k·t тождественно. Тогда у группы из трёх линий
остаётся ОДИН свободный параметр активности вместо трёх площадей, и задача из
плохо обусловленной становится хорошо обусловленной.

ЧТО ЗАФИКСИРОВАНО И ПОЧЕМУ:
  положения    E_k из библиотеки плюс ОДИН общий сдвиг Δ на всю группу —
               калибровочный дефект общий, а не свой у каждой линии;
  ПШПВ         по закону σ(E) = FWHM₆₆₂·√(E/661,657)/2,355, не свободна.
               У ЛСРМ есть прямое правило: при dS/S > 0,1 полуширину снимают с
               подгонки, а в блендес перекрытием больше половины ПШПВ свободная
               полуширина уводит решение всегда;
  интенсивности I_k — из спектра испускания ТОГО ЖЕ прогона Geant4 (*_emit.csv),
               а не из справочника. Это правило проекта: числовые данные берутся
               из той же базы, что и транспорт;
  континуум    линейный по энергии в пределах области. Комптоновская ступенька не
               вводится: у ЛСРМ её снимают при dS/S > 0,05, а здесь области
               широкие и ступенька не отделяется от наклона.

МОДЕЛЬ ЛИНЕЙНА по (A, c₀, c₁) при фиксированном Δ — значит решается одним МНК без
итераций, а Δ перебирается по сетке. Ковариация даёт погрешность активности.

ЧЕГО ЭТОТ МЕТОД НЕ ДАЁТ, читать обязательно. Связь через A·I_k·ε_k использует
РАСЧЁТНУЮ эффективность, поэтому внутригрупповые отношения площадей заданы
моделью, а не измерены. Проверяется абсолютный масштаб — активность против
паспорта, — а не форма кривой внутри группы. Для отношений формы годятся только
одиночные линии. Так же устроено и у ЛСРМ: их связывающий член тоже опирается на
их собственную ε.

СОСТОЯНИЕ: ПРОТОТИП, ЧИСЛА НЕ ПУБЛИКОВАТЬ. Первый прогон на записи Th-232 в
маринелли (паспорт 3104 Бк):

    238,6 кэВ   A/пасп 1,034   4 линии в группе
    583,2 кэВ   A/пасп 0,798   3 линии
    911,2 кэВ   A/пасп 0,802   5 линий
   2614,5 кэВ   A/пасп 0,631   1 линия

Обнадёживает то, что 583 и 911 — обе ранее отброшенные как нераздельные — дали
0,798 и 0,802, то есть сошлись между собой на 0,5 % и с маринелльным множителем
0,78, полученным по чистым линиям совсем другим способом. Ради этого метод и
затевался.

Но ДВА признака говорят, что он ещё не готов:
  1. Одиночная линия 2614,5 обязана дать то же, что оконный съём (0,796), а даёт
     0,631. У одиночной линии деконволюции нечего делать, значит расхождение
     20 % — это дефект нормировки, а не физика.
  2. χ²/dof выходит 0,00–0,01. Так не бывает: веса взяты как 1/(|y|+1) вместо
     обратной дисперсии, поэтому и χ², и ковариация, и погрешность dA
     недостоверны.

Пока эти два пункта не закрыты, отбор линий по чистоте (kit_recalc) остаётся в
силе, а числа из этого модуля никуда не идут.

    python detectors/Gamma-1S/analysis/deconv.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
FWHM662 = 49.9
SIG = 2.0 * math.sqrt(2.0 * math.log(2.0))      # ПШПВ = 2,355·σ
EMIT_HALF = 3.0            # кэВ, ширина «самой линии» в спектре испускания
SHIFT_GRID = np.arange(-12.0, 12.01, 0.5)       # кэВ, перебор общего сдвига


def sigma(E):
    return FWHM662 * math.sqrt(E / 661.657) / SIG


def load_hist(path):
    h, N = {}, None
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"):
            if "N_primaries" in ln:
                N = int(ln.split("=")[1])
            continue
        if ln[:1].isdigit():
            e, c = ln.split(",")
            h[float(e)] = int(c)
    return h, N


def group_lines(base, E0, half, min_frac=0.02):
    """Линии внутри области из спектра испускания: [(E, выход на распад)].

    Порог 2 % от полного выхода области: слабее этого линия в подгонке ничего не
    меняет, а параметров не добавляет — их и так подставляем, но список короче
    читается.
    """
    p = os.path.join(BUILD, base + "_emit.csv")
    if not os.path.exists(p):
        return None
    emit, N = load_hist(p)
    if not N:
        return None
    tot = sum(c for e, c in emit.items() if abs(e - E0) <= half)
    if tot <= 0:
        return None
    # группируем соседние каналы спектра испускания в линии
    peaks = []
    for e in sorted(e for e in emit if abs(e - E0) <= half):
        c = emit[e]
        if c <= 0:
            continue
        if peaks and e - peaks[-1][1] <= EMIT_HALF:
            w = peaks[-1][0] + c
            peaks[-1] = (w, (peaks[-1][1] * peaks[-1][0] + e * c) / w)
        else:
            peaks.append((c, e))
    out = [(e, c / N) for c, e in peaks if c > min_frac * tot]
    return sorted(out, key=lambda x: -x[1])


def eps_at(gtag, E):
    """Эффективность ППП по сетке моноэнергий: линейно по log-log между узлами."""
    import glob
    import re
    pts = []
    for p in glob.glob(os.path.join(BUILD, "grid", gtag + "_E*.csv")):
        m = re.search(r"_E(\d+\.\d)\.csv$", p)
        if not m:
            continue
        h, N = load_hist(p)
        if not N:
            continue
        Eg = float(m.group(1))
        gross = sum(c for e, c in h.items() if abs(e - Eg) <= 6.0)
        side = sum(c for e, c in h.items() if Eg - 30 <= e <= Eg - 10)
        n = math.floor(Eg + 6 - 0.5) - math.ceil(Eg - 6 - 0.5) + 1
        ns = math.floor(Eg - 10 - 0.5) - math.ceil(Eg - 30 - 0.5) + 1
        net = gross - side / ns * n
        if net > 0:
            pts.append((Eg, net / N))
    if len(pts) < 4:
        return None
    pts.sort()
    xs = np.log([p[0] for p in pts])
    ys = np.log([p[1] for p in pts])
    return float(np.exp(np.interp(math.log(E), xs, ys)))


def fit_group(sp, bg, lines, eps, live, lo, hi, shift):
    """МНК по каналам области [lo, hi] при заданном сдвиге. -> (A, dA, chi2/dof)."""
    ch = np.arange(len(sp.n), dtype=float)
    en = sp.energy(ch)
    m = (en >= lo) & (en <= hi)
    if m.sum() < 8:
        return None
    x, y = en[m], sp.n[m].astype(float)
    if bg is not None:
        bch = np.arange(len(bg.n), dtype=float)
        ybg = np.interp(x, bg.energy(bch), bg.n.astype(float))
        y = y - ybg * (live / bg.live)
    # площадь линии на канал: гауссиана, нормированная на ширину канала
    dE = np.gradient(x)
    shape = np.zeros_like(x)
    for E, I in lines:
        e = eps.get(round(E, 1))
        if not e:
            continue
        s = sigma(E)
        shape += (I * e * live * dE
                  / (s * math.sqrt(2 * math.pi))
                  * np.exp(-0.5 * ((x - (E + shift)) / s) ** 2))
    if shape.max() <= 0:
        return None
    # линейная модель: y = A*shape + c0 + c1*(x - x0)
    x0 = 0.5 * (lo + hi)
    M = np.vstack([shape, np.ones_like(x), (x - x0) / max(hi - lo, 1.0)]).T
    w = 1.0 / np.maximum(np.abs(y) + 1.0, 1.0)
    Mw = M * w[:, None]
    yw = y * w
    sol, *_ = np.linalg.lstsq(Mw, yw, rcond=None)
    resid = yw - Mw @ sol
    dof = max(1, len(x) - M.shape[1])
    chi2 = float((resid ** 2).sum()) / dof
    try:
        cov = np.linalg.inv(Mw.T @ Mw) * chi2
        dA = float(math.sqrt(max(cov[0, 0], 0.0)))
    except np.linalg.LinAlgError:
        return None
    return float(sol[0]), dA, chi2


def deconvolve(sp, bg, base, gtag, E0, span=1.6):
    """Активность по группе вокруг E0. -> dict с результатом или None."""
    half = span * FWHM662 * math.sqrt(E0 / 661.657)
    lines = group_lines(base, E0, half)
    if not lines:
        return None
    eps = {}
    for E, _I in lines:
        e = eps_at(gtag, E)
        if e:
            eps[round(E, 1)] = e
    if not eps:
        return None
    best = None
    for d in SHIFT_GRID:
        r = fit_group(sp, bg, lines, eps, sp.live, E0 - half, E0 + half, d)
        if r and r[0] > 0 and (best is None or r[2] < best[3]):
            best = (r[0], r[1], float(d), r[2])
    if best is None:
        return None
    A, dA, shift, chi2 = best
    return dict(A=A, dA=dA, shift=shift, chi2=chi2, lines=lines,
                half=half, n_lines=len(lines))
