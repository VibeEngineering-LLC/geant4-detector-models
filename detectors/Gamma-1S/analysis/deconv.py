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

Связь взята в предельно жёсткой форме: площади ПОДСТАВЛЕНЫ, а не оштрафованы,
S_k = A·I_k·ε_k·t тождественно. У группы из трёх линий остаётся ОДИН свободный
параметр вместо трёх площадей, и задача из плохо обусловленной становится хорошо
обусловленной.

НОРМИРОВКА — ЧЕРЕЗ ВТОРУЮ ТАКУЮ ЖЕ ПОДГОНКУ, и это главное в модуле.

Произведение I_k·ε_k нигде не выписывается. Вместо него та же самая подгонка
делается второй раз — по МОДЕЛЬНОМУ спектру полного распада, уширенному до
разрешения прибора: тот же участок, тот же набор линий, та же форма, тот же
континуум, те же веса. Тогда

    A·t = N_распадов · S_измер / S_модель

и из отношения выпадает всё, что одинаково с двух сторон: выход линии, ППП,
каскадное суммирование, вклад соседей в бленде, континуум под пиком, доля пика,
срезанная краями участка. Остаётся только то, ради чего расчёт и делался.

ПОЧЕМУ НЕ ЧЕРЕЗ СЕТКУ МОНОЭНЕРГИЙ, как было в первой редакции. Там ε бралась с
сетки моноэнергий и умножалась на выход линии из спектра испускания. Измерено на
Tl-208 2614,5 кэВ в маринелли: сетка не знает о совпадениях и даёт ε на 12 %
больше, чем прогон цепочки (2,526e-3 против 2,231e-3 на распад), а до оконной
величины 1,917e-3 не хватает ещё столько же — на разницу в ширине окна. Отсюда и
бралось расхождение 20 % на ОДИНОЧНОЙ линии, где деконволюции делать нечего.
Это была не «нормировка», а недостающая физика плюс несогласованность окон.
Правило самосогласованности («площади с обеих сторон снимаются одним и тем же
алгоритмом») записано в шапке kit_recalc; здесь оно проведено до конца.

ЧТО ЗАФИКСИРОВАНО И ПОЧЕМУ:
  положения    E_k из библиотеки плюс ОДИН общий сдвиг Δ на всю группу —
               калибровочный дефект общий, а не свой у каждой линии. У модели
               сдвига нет и быть не может: у неё линии стоят на истинных
               энергиях, поэтому модельная подгонка идёт при Δ = 0;
  ПШПВ         НЕ по закону корня, а по калибровке прибора ПШПВ² = a + b·E,
               снятой с трёх сильных ОДИНОЧНЫХ линий маринелли (662, 1461,
               2614,5). Закон √E из одной точки на этом приборе врёт: на 583
               завышает (46,8 против измеренных 40,4), на 2614 занижает (99
               против 105). Окну это почти безразлично, подгонке — нет: на
               одиночной 2614,5 закон корня стоил 9 % площади. Независимый
               контроль калибровки — линия 238,6, в подгонку не входившая:
               измерено 23,8, модель 23,6 кэВ. Свободной полуширина не
               делается: у ЛСРМ прямое правило снимать её с подгонки при
               dS/S > 0,1, а в бленде с перекрытием больше половины ПШПВ
               свободная ширина уводит решение всегда (gamma.peaks.
               fit_stability);
  интенсивности I_k — из спектра испускания ТОГО ЖЕ прогона (*_emit.csv), а не из
               справочника: правило проекта — числа из той же базы, что транспорт;
  континуум    линейный плюс ступенька erfc под группой. Ступенька обязательна
               там, где выше группы излучения почти нет (2614,5 кэВ — верх
               спектра): прямая через такой перепад садится на склон и забирает
               часть пика.

ЧИСЛЕННАЯ ЧАСТЬ ВЗЯТА У SPECTRAVIBE (gamma/peaks/coupled_multiplet.py), потому
что там она сделана правильно, а у меня в первой редакции — нет:
  веса         1/σ при σ = √(отсчёты), то есть обратная СИГМА. Стояло 1/(|y|+1) —
               обратная ДИСПЕРСИЯ. Из-за этого каналы пика получали ничтожный вес,
               а χ²/dof выходил 0,00–0,01, чего не бывает;
  границы      амплитуда и ступенька ≥ 0 через scipy.optimize.lsq_linear (trf);
  ковариация   через SVD матрицы плана, а не через inv(XᵀX): у нормальной матрицы
               число обусловленности в квадрате, и на тесных мультиплетах
               погрешности искажаются.

ЧЕГО ЭТОТ МЕТОД НЕ ДАЁТ, читать обязательно. Отношения площадей внутри группы
заданы моделью, а не измерены: связь через одну активность именно это и означает.
Проверяется АБСОЛЮТНЫЙ масштаб — активность против паспорта, — а не форма кривой
внутри группы. Для отношений формы годятся только одиночные чистые линии. Так же
устроено и у ЛСРМ.

СОСТОЯНИЕ ПОСЛЕ ПРОВЕРКИ 28.07.2026.

Нормировка ПОДТВЕРЖДЕНА: на всех одиночных линиях, где деконволюции нечего
делать, она сходится с оконным съёмом — 662: 0,773 против 0,775; 1461: 0,712
против 0,720; 2614,5 у Петри 1,196 против 1,195, у «Денты» 1,079 против 1,090
(у маринелли 0,769 против 0,796 — 3 %, χ²/dof 3,6: при её статистике вылезает
несовершенство гауссовой формы). Оба дефекта первой редакции закрыты, χ²/dof
стал осмысленным (0,8–2 там, где форма описывается).

НЕ ЗАКРЫТО: группы 583 и 911 систематически ВЫШЕ одиночной 2614,5 того же
нуклида — на +14/+12 % в маринелли, +24/+20 % в «Денте», +26/+2 % в Петри.
Это не веса и не нормировка (проверено), не 511 кэВ из комнаты (маскирование
511±1,5σ с обеих сторон ничего не меняет) — это несовпадение формы континуума
модели и измерения в середине спектра. Пока оно не разобрано, для АКТИВНОСТИ
остаётся в силе отбор линий по чистоте (kit_recalc); числа деконволюции — для
отработки метода.

    python detectors/Gamma-1S/analysis/deconv.py          # весь комплект
    python detectors/Gamma-1S/analysis/deconv.py --one Th-232
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402
import kit_recalc as kr  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
FWHM662 = kr.FWHM662
SIG = 2.0 * math.sqrt(2.0 * math.log(2.0))      # ПШПВ = 2,355·σ
EMIT_HALF = 3.0            # кэВ, ширина «самой линии» в спектре испускания
SHIFT_GRID = np.arange(-12.0, 12.01, 0.5)       # кэВ, перебор общего сдвига
SPAN = 1.6                 # полуширина участка подгонки в ПШПВ
MIN_FRAC = 0.02            # порог интенсивности линии внутри участка

# Опорные ОДИНОЧНЫЕ линии для калибровки ПШПВ² = a + b·E. Только маринелли:
# у мелких кювет статистики на полувысоту не хватает (у «Денты» на 2614
# автомат намерил 66 кэВ — это шум, а не разрешение).
_FWHM_PROBES = (("*M_cs_*", 661.657), ("*M_k_*", 1460.822),
                ("*Th232*", 2614.511))


# Объявление наблюдаемой — что именно за число лежит в таблице. Без него
# таблицу нельзя сравнивать ни с какой другой: за один вечер 30.07.2026
# подмена определения стоила вывода четыре раза (method-rules §5).
OBS = {
    "quantity":
        "активность по группе линий; восстановленная связанной деконволюцией окна",
    "area":
        "площади линий группы из совместной подгонки; не по отдельности",
    "window":
        "окно группы охватывает все линии бленда",
    "shelf":
        "подложка подгоняется вместе с линиями группы",
    "blurred":
        "измерение как есть; модель размыта приборной ПШПВ",
}


def _stamp(inputs=None):
    return stamp.lines("detectors/Gamma-1S/analysis/deconv.py", OBS,
                       inputs=inputs,
                       geometry_dir=str(paths.geometry("Gamma-1S")),
                       names=stamp.SRC_LISTS["Gamma-1S"],
                       repo_dir=str(paths.REPO))


def _fwhm_calibrate():
    kd = paths.kit_dir("Marinelli_1L")
    pts = []
    if kd:
        for mask, E in _FWHM_PROBES:
            fs = sorted(str(p) for p in kd.rglob(mask))
            if not fs:
                continue
            sp, _bg = bm.read(fs[0])
            f = bm.fwhm_at(sp, E)
            if f:
                pts.append((E, f))
    if len(pts) < 3:
        return None
    E = np.array([p[0] for p in pts])
    F = np.array([p[1] for p in pts])
    (a, b), *_ = np.linalg.lstsq(np.vstack([np.ones_like(E), E]).T,
                                 F ** 2, rcond=None)
    return float(a), float(b)


_FWHM_AB = _fwhm_calibrate()


def fwhm(E):
    """ПШПВ прибора, кэВ: калибровка ПШПВ² = a + b·E, при её отсутствии —
    закон корня от 662 (тогда числа хуже, см. шапку)."""
    if _FWHM_AB:
        a, b = _FWHM_AB
        v = a + b * E
        if v > 1.0:
            return math.sqrt(v)
    return FWHM662 * math.sqrt(E / 661.657)


def sigma(E):
    return fwhm(E) / SIG


def group_lines(base, E0, half, min_frac=MIN_FRAC):
    """Линии, видимые внутри участка, из спектра испускания: [(E, выход)].

    Берутся не только линии ВНУТРИ участка [E0−half, E0+half], но и соседи
    СНАРУЖИ в пределах 3σ от края: их склон достаёт до участка, и без них
    подгонка садится на чужой хвост. Измерено на Ra-226: линия 295,2 стоит в
    шести кэВ за краем участка 351,9 — без неё χ²/dof 14–39 и завышение до
    24 %, с ней участок описывается. Модельной подгонке эти же линии нужны по
    той же причине — обе стороны остаются согласованными.

    Порог 2 % от полного выхода: слабее этого линия в подгонке ничего не
    меняет. Соседние каналы спектра испускания сливаются в одну линию, если
    отстоят меньше чем на EMIT_HALF.
    """
    p = os.path.join(BUILD, base + "_emit.csv")
    if not os.path.exists(p):
        return None
    emit, N = kr.load_hist(p)
    if not N:
        return None
    reach = half + 3.0 * sigma(E0)
    tot = sum(c for e, c in emit.items() if abs(e - E0) <= reach)
    if tot <= 0:
        return None
    peaks = []
    for e in sorted(e for e in emit if abs(e - E0) <= reach):
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


def _erfc(z):
    from scipy.special import erfc as _e
    return _e(z)


def _design(x, dE, lines, shift, E_step, sig_step):
    """Матрица плана: [форма группы, 1, наклон, ступенька]."""
    shape = np.zeros_like(x)
    for E, I in lines:
        s = sigma(E)
        shape += (I * dE / (s * math.sqrt(2 * math.pi))
                  * np.exp(-0.5 * ((x - (E + shift)) / s) ** 2))
    x0 = 0.5 * (x[0] + x[-1])
    span = max(x[-1] - x[0], 1.0)
    step = 0.5 * _erfc((x - E_step) / (sig_step * math.sqrt(2.0)))
    return np.vstack([shape, np.ones_like(x), (x - x0) / span, step]).T


def _solve(M, y, sy):
    """Взвешенный МНК с границами: амплитуда и ступенька ≥ 0.

    Веса, границы и ковариация — как в gamma/peaks/coupled_multiplet.py:
    w = 1/σ при σ = √отсчётов, решатель lsq_linear(method='trf'),
    ковариация через SVD матрицы плана.
    """
    w = 1.0 / sy
    Mw = M * w[:, None]
    yw = y * w
    n = M.shape[1]
    lb = np.full(n, -np.inf)
    ub = np.full(n, np.inf)
    lb[0] = 0.0        # амплитуда группы
    lb[3] = 0.0        # ступенька континуума
    try:
        from scipy.optimize import lsq_linear
        res = lsq_linear(Mw, yw, bounds=(lb, ub), method="trf",
                         tol=1e-10, max_iter=400)
        p, ok = res.x, bool(res.success)
    except ImportError:
        p, *_ = np.linalg.lstsq(Mw, yw, rcond=None)
        p = np.where((np.arange(n) == 0) | (np.arange(n) == 3),
                     np.maximum(p, 0.0), p)
        ok = True
    resid = yw - Mw @ p
    dof = max(1, len(y) - n)
    chi2 = float((resid ** 2).sum()) / dof
    try:
        _U, s, Vt = np.linalg.svd(Mw, full_matrices=False)
        if s.size == 0 or s[-1] <= 0:
            return None
        cov = (Vt.T * (1.0 / (s * s))) @ Vt * chi2
        dS = float(math.sqrt(max(cov[0, 0], 0.0)))
    except np.linalg.LinAlgError:
        return None
    return float(p[0]), dS, chi2, ok, np.asarray(p, dtype=float)


def _fit_measured(sp, bg, lines, lo, hi, shift):
    """Подгонка ИЗМЕРЕННОГО спектра за вычетом фона. -> (S, dS, chi2/dof)."""
    ch = np.arange(len(sp.n), dtype=float)
    en = sp.energy(ch)
    m = (en >= lo) & (en <= hi)
    if m.sum() < 8:
        return None
    x = en[m]
    gross = sp.n[m].astype(float)
    ybg = np.zeros_like(gross)
    if bg is not None:
        bch = np.arange(len(bg.n), dtype=float)
        ybg = np.interp(x, bg.energy(bch),
                        bg.n.astype(float)) * (sp.live / bg.live)
    y = gross - ybg
    # Дисперсия чистой площади — сумма дисперсий пробы и приведённого фона,
    # а не «отсчёты после вычитания»: после вычета y может быть около нуля,
    # тогда как неопределённость там наибольшая.
    sy = np.sqrt(np.maximum(gross, 1.0) + np.maximum(ybg, 0.0))
    E_step = 0.5 * (lines[0][0] + lines[-1][0]) + shift
    sig_step = max(sigma(E) for E, _ in lines)
    M = _design(x, np.gradient(x), lines, shift, E_step, sig_step)
    r = _solve(M, y, sy)
    if r is None:
        return None
    # Кривые для рисования — из ТЕХ ЖЕ колонок, что и решение, чтобы картинка
    # не могла разойтись с числом (см. analysis/spectra_figs.py).
    return r + (dict(x=x, y=y, sy=sy, dE=np.gradient(x), shift=shift,
                     model=M @ r[4], cont=M[:, 1:] @ r[4][1:]),)


def _fit_model(arr, lines, lo, hi, bin_keV=1.0):
    """Та же подгонка по УШИРЕННОМУ модельному спектру. Сдвига у модели нет."""
    i0, i1 = int(round(lo / bin_keV)), int(round(hi / bin_keV))
    i0, i1 = max(0, i0), min(len(arr), i1)
    if i1 - i0 < 8:
        return None
    x = np.arange(i0, i1, dtype=float) * bin_keV
    y = arr[i0:i1].astype(float)
    sy = np.sqrt(np.maximum(y, 1.0))
    E_step = 0.5 * (lines[0][0] + lines[-1][0])
    sig_step = max(sigma(E) for E, _ in lines)
    M = _design(x, np.full_like(x, bin_keV), lines, 0.0, E_step, sig_step)
    return _solve(M, y, sy)


_BROAD = {}


def _broadened(base):
    # Уширение — ТЕМ ЖЕ калиброванным законом, каким подгоняется измерение:
    # на этом держится сокращение нормировки.
    if base not in _BROAD:
        hist, N = kr.load_hist(os.path.join(BUILD, base + ".csv"))
        _BROAD[base] = (bm.broaden(hist, fwhm_of=fwhm), N)
    return _BROAD[base]


def deconvolve(sp, bg, base, E0, geom=None, rho_src=None, span=SPAN):
    """Активность по группе вокруг E0, Бк. -> dict или None.

    A·t = N_распадов · S_измер / S_модель, где обе площади сняты ОДНОЙ и той же
    подгонкой. Если плотность источника отличается от плотности прогона, вводится
    отношение поправок самопоглощения f(mu·ρ·d) — единственный множитель, который
    из отношения не выпадает.
    """
    half = span * fwhm(E0)
    lines = group_lines(base, E0, half)
    if not lines:
        return None
    lines = sorted(lines)
    lo, hi = E0 - half, E0 + half
    arr, N = _broadened(base)
    if not N:
        return None
    fm = _fit_model(arr, lines, lo, hi)
    if not fm or fm[0] <= 0:
        return None
    best = None
    for d in SHIFT_GRID:
        r = _fit_measured(sp, bg, lines, lo, hi, float(d))
        if r and r[0] > 0 and (best is None or r[2] < best[2]):
            best = (r[0], r[1], r[2], float(d), r[5])
    if best is None:
        return None
    S, dS, chi2, shift, fit = best
    corr = 1.0
    if geom and rho_src:
        rho_run, mat_run = kr.RUNRHO[geom]
        dmm = kr.GRIDS[geom][2]
        key = min(kr.MU_O, key=lambda k: abs(k - E0))
        mu_run = kr.MU_O[key] if mat_run == "OISN16" else kr.MU_W[key]
        mu_src = kr.MU_O[key] if rho_src > 1.3 else kr.MU_W[key]
        corr = (kr.fx(mu_src * rho_src * dmm / 10)
                / kr.fx(mu_run * rho_run * dmm / 10))
    A = N * S / fm[0] / sp.live / corr
    dA = A * math.hypot(dS / S, fm[1] / fm[0])
    # Вклад каждой линии по отдельности — для врезок деконволюции на странице.
    x, dEx = fit["x"], fit["dE"]
    fit["parts"] = []
    for E, I in lines:
        s = sigma(E)
        fit["parts"].append((E, I, S * I * dEx / (s * math.sqrt(2 * math.pi))
                             * np.exp(-0.5 * ((x - (E + shift)) / s) ** 2)))
    return dict(A=A, dA=dA, shift=shift, chi2=chi2, chi2_model=fm[2],
                n_lines=len(lines), lines=lines, half=half,
                S=S, S_model=fm[0], N=N, fit=fit)


# ---------------------------------------------------------------------------

def _run():
    only = None
    if "--one" in sys.argv:
        only = sys.argv[sys.argv.index("--one") + 1]
    print("Связанная деконволюция комплекта: активность из ОДНОЙ подгонки на\n"
          "группу, нормировка второй такой же подгонкой по модельному спектру.\n"
          "Сверка — с оконным съёмом (kit_recalc) на тех же линиях.\n")
    print("%-13s %-8s %8s %5s %9s %9s %8s %8s %7s"
          % ("геометрия", "нуклид", "E, кэВ", "линий", "A деконв.", "A оконн.",
             "деконв.", "оконн.", "χ²/dof"))
    rows = []
    for geom, mask, nuc, aspec, dpct, d0, mass, vol in kr.VOLUME_RECORDS:
        if only and nuc != only:
            continue
        kd = paths.kit_dir(geom)
        files = sorted(str(p) for p in kd.rglob(mask)) if kd else []
        if not files:
            continue
        sp, bg, _cal = bm.read_checked(files[0])
        import re
        txt = open(files[0], encoding="utf-8", errors="replace").read()
        md = re.search(r"<StartTime>(\d{4}-\d{2}-\d{2})", txt)
        md = md.group(1) if md else None
        A0 = aspec * mass / 1000.0 * kr.decay_factor(nuc, d0, md)
        rho = mass / vol
        R = float(sp.n.sum()) / sp.live
        pile = math.exp(2 * kr.TAU_SHAPE * R)
        lines, ckey = kr.VLINES[nuc]
        base = kr.RUNBASE.get((geom, ckey))
        if not base:
            continue
        for E in lines:
            r = deconvolve(sp, bg, base, E, geom=geom, rho_src=rho)
            if not r:
                continue
            A = r["A"] * pile
            dA = r["dA"] * pile
            # Оконный столбец обязан ВОСПРОИЗВОДИТЬ kit_recalc до последней
            # цифры, иначе сверка ничего не проверяет. Отсюда два условия:
            # ПШПВ по закону корня (у kit_recalc он) и матрица источника по
            # его же правилу — у лёгких засыпок (ρ ≤ 1,3) состав в файлах не
            # записан и берётся ВОДА, а не ОИСН-16 с её 71 % железа. Оба
            # места я сначала взял по-своему, и таблица разошлась: ПШПВ дала
            # до 12 % на Ra-226 351,9, матрица — ровные 0,8 % на всех лёгких.
            fw = FWHM662 * math.sqrt(E / 661.657)
            frac, _dirt = kr.purity(base, E, fw)
            nr = bm.net_rate(sp, bg, E, fw, roi=1.0, side=1.0)
            key = min(kr.MU_O, key=lambda k: abs(k - E))
            mu_src = kr.MU_O[key] if rho > 1.3 else kr.MU_W[key]
            eps = kr.eps_per_decay(geom, ckey, E, fw, rho, mu_src)
            win = (nr[0] * pile / eps / A0) if (nr and eps and nr[0] > 0) else None
            print("%-13s %-8s %8.1f %5d %9.1f %9s %8.3f %8s %7.2f"
                  % (geom, nuc, E, r["n_lines"], A,
                     "%.1f" % (win * A0) if win else "—",
                     A / A0, "%.3f" % win if win else "—", r["chi2"]))
            rows.append((geom, nuc, E, r["n_lines"], A, dA, A0, A / A0,
                         win if win else 0.0, frac if frac else 0.0,
                         r["chi2"], r["chi2_model"], r["shift"]))
        # Активность нуклида по правилу ЛСРМ — по ВСЕМ линиям: в связанной
        # подгонке бленд разобран, и отбрасывать линию за неразделённость
        # незачем. В этом смысл затеи. НО: пока не разобрано несовпадение
        # формы континуума (см. шапку), эти средние — отработка метода,
        # опубликованная активность считается по kit_recalc.
        sel = [r for r in rows if r[0] == geom and r[1] == nuc]
        av = kr.lsrm_average([(r[4], r[5]) for r in sel])
        if av:
            m, dm, kind, n = av
            print("   %-10s по %d линиям: A = %.4g ± %.2g Бк, "
                  "A/пасп = %.3f ± %.3f (%s)"
                  % (nuc, n, m, dm, m / sel[0][6], dm / sel[0][6], kind))
    if rows:
        out = os.path.abspath(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..", "results", "deconv_lines.csv"))
        csvio.write(
            out,
            ["geometry", "nuclide", "E_keV", "n_lines", "A_Bq", "dA_Bq",
             "A_pass_Bq", "ratio", "ratio_window", "purity", "chi2_meas",
             "chi2_model", "shift_keV"],
            [(r[0], r[1], "%.3f" % r[2], "%d" % r[3], "%.2f" % r[4],
              "%.2f" % r[5], "%.2f" % r[6], "%.4f" % r[7], "%.4f" % r[8],
              "%.3f" % r[9], "%.3f" % r[10], "%.3f" % r[11],
              "%+.2f" % r[12]) for r in rows],
            comments=[
                "Связанная деконволюция: активность по группе линий.",
                "A = N_распадов*S_измер/S_модель/t — обе площади сняты",
                "  одной и той же подгонкой, см. шапку deconv.py.",
                "ratio_window — то же отношение оконным съёмом; у чистых",
                "  линий (purity>=0.95) обязано совпадать.",
            ],
        stamp=_stamp())
        print("\nтаблица: %s (%d строк)" % (out, len(rows)))


if __name__ == "__main__":
    _run()
