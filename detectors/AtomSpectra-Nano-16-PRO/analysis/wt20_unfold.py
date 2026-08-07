# -*- coding: utf-8 -*-
"""Разложение измеренного спектра по шаблонам нуклидов ряда тория.

Метод шаблонный: измеренный спектр представляется суммой откликов ОТДЕЛЬНЫХ
нуклидов, посчитанных Монте-Карло в той же геометрии, плюс измеренный фон.
Шаблон нормирован на ОДИН распад, поэтому коэффициенты разложения — сразу
активности в беккерелях, без промежуточной «эффективности по линии».

    измеренное(E) = Σ aᵢ · t · Tᵢ(E) + b · фон(E)

где aᵢ — активность нуклида i, Бк; t — время набора, с; Tᵢ — шаблон, отсчёты на
распад в канале E; b — множитель фона (ожидается около единицы).

Побочные пики отдельными компонентами НЕ вводятся: вылет аннигиляционных
квантов, обратное рассеяние, характеристический рентген вольфрама и суммирование
каскада возникают в переносе сами и лежат внутри шаблона своего нуклида. Тем
методика отличается от разложения на аналитические компоненты, где каждый такой
пик приходится перечислять руками.

Приборное разрешение навешивается ЗДЕСЬ: ПШПВ² = f0 + f1·E подбирается вместе с
активностями (перебором по сетке), потому что паспортного разрешения у прибора
нет, а измеренные ширины отдельных линий на мультиплетах ненадёжны.

    python analysis/wt20_unfold.py <спектр.xml> <каталог шаблонов> [каталог вывода]
"""
import csv
import io
import math
import os
import re
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.environ.get("SPECTRAVIBE_ROOT")
if not _ROOT:
    raise SystemExit("не задана переменная окружения SPECTRAVIBE_ROOT")
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

from gamma.io.atomspectra_xml import read_atomspectra_xml        # noqa: E402

E_MAX = 3700.0          # верх энергетической сетки, кэВ
E_STEP = 2.0            # шаг сетки разложения, кэВ
E_FIT = (60.0, 3000.0)  # окно подгонки и показа (верх — по директиве оператора)

# Ступенчатая подгонка с жёсткого конца (2614 -> 911 -> мягкий участок) введена
# 07.08.2026 и в тот же день ОТМЕНЕНА как умолчание решением оператора. Код
# оставлен за флагом WT20_STEPWISE=1, метод описан в `docs/wt20-method.md`
# §4.1а. Вернуться к нему предполагается после того, как будет закрыт избыток
# модельного континуума: именно он мешал и выбору ширины, и оценке масштаба.
E_HARD = 2614.51
E_MID = 911.20
W_SOFT = (25.0, 150.0)
STEPWISE = os.environ.get("WT20_STEPWISE") == "1"
GLOBAL_FIT = not STEPWISE

# Границы окна переопределяются переменными окружения — так разброс по окну
# считается прогоном одного и того же кода (analysis/wt20_window_scan.py).
if os.environ.get("WT20_FIT_LO"):
    E_FIT = (float(os.environ["WT20_FIT_LO"]),
             float(os.environ.get("WT20_FIT_HI", E_FIT[1])))

# Удельная активность Th-232, Бк/г. ВЫЧИСЛЕНА: A = ln2·N_A/(T½·M) при
# T½ = 4,41797·10¹⁷ с и M = 232,038054 а.е.м. (API МАГАТЭ, ENSDF).
SPEC_ACT_TH232 = 0.6931472 * 6.02214076e23 / (4.41796963644288e17 * 232.038054)

# Порядок и подписи компонент — по месту в ряду. XI и XW — не члены цепочки,
# это флуоресценция иода в кристалле и K-серия вольфрама, вносимая самим
# источником. Их шаблоны — линейчатые, из каталога `reference/roi/`
# (`analysis/build_xray_templates.py`); в отличие от нуклидных компонент, их
# коэффициент подгонки — темп эмиссии рентгена, эмиссий/с, а не активность.
ORDER = [
    ("Th232", "Th-232", "#7a5c3a"),
    ("Ra228", "Ra-228", "#a08a5c"),
    ("Ac228", "Ac-228", "#d81b8c"),
    ("Th228", "Th-228", "#8a6d3b"),
    ("Ra224", "Ra-224", "#c98b1e"),
    ("Rn220", "Rn-220", "#6b8f3a"),
    ("Po216", "Po-216", "#9bb06a"),
    ("Pb212", "Pb-212", "#b07d2a"),
    ("Bi212", "Bi-212", "#2f6b34"),
    ("Tl208", "Tl-208", "#c8cf7a"),
    ("Po212", "Po-212", "#8fa0a8"),
    ("XI",    "X-I (флуоресценция иода в CsI)",   "#2E7D32"),
    ("XW",    "X-W (K-серия вольфрама электрода)", "#C03535"),
    ("XD1",   "X-K дочерних (ветвь A1)",           "#7b3fb3"),
    ("XD2",   "X-K дочерних (ветвь A2)",           "#3b6ec0"),
]

ORDER_MAP = [(k, lab) for k, lab, _ in ORDER]

# Ветвление Bi-212: доля распадов, идущих через Tl-208. Число берётся из
# библиотеки МАГАТЭ при запуске (поле decay_% в файле линий), а не пишется сюда.
LIB = os.path.normpath(os.path.join(_HERE, "..", "reference", "nuclide-lines"))


def branch_to_tl208():
    """Доля распадов Bi-212 по альфа-ветви (на Tl-208), из файла линий Tl-208.

    В библиотеке МАГАТЭ поле `decay_%` строк Tl-208 — это доля РОДИТЕЛЬСКОГО
    распада, приводящая к этому нуклиду. Берётся оттуда, а не пишется числом:
    35,94 % — величина оценённая и может уточняться.
    """
    p = os.path.join(LIB, "212bi_gammas.csv")
    for r in csv.DictReader(io.open(p, encoding="utf-8")):
        d, pc = r.get("decay"), r.get("decay_%")
        if d and d.strip().upper().startswith("A") and pc:
            try:
                return float(pc) / 100.0
            except ValueError:
                pass
    raise SystemExit("в %s нет альфа-ветви Bi-212" % p)


def read_template(path):
    """Спектр Geant4: шапка + E_keV,counts. -> (dict шапки, E[], counts[])."""
    head, e, c = {}, [], []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("# ").split("=", 1)
                head[k.strip()] = v.strip()
            continue
        if not ln or ln.startswith("E_keV"):
            continue
        a, b = ln.split(",")
        e.append(float(a))
        c.append(float(b))
    return head, np.array(e), np.array(c)


def poly(coefs, x):
    out = np.zeros_like(np.asarray(x, dtype=float))
    for k, c in enumerate(coefs):
        out = out + c * np.asarray(x, dtype=float) ** k
    return out


def read_correction(path):
    """Поправки калибровки из wt20_calibration.py. -> dict имя->коэффициенты."""
    out = {}
    if not os.path.exists(path):
        return out
    for row in csv.reader(io.open(path, encoding="utf-8")):
        if not row or row[0].startswith("#") or row[0] == "спектр":
            continue
        out[row[0]] = [float(v) for v in row[1:] if v not in ("", None)]
    return out


def rebin_to_grid(counts, cal, corr, grid_edges):
    """Отсчёты по каналам -> отсчёты по энергетической сетке.

    Канал считается равномерно заполненным по своему энергетическому интервалу
    [E(ch−½), E(ch+½)]; его отсчёты делятся между ячейками сетки по перекрытию.
    Так сохраняется ПОЛНОЕ число отсчётов — проверяется в конце.
    """
    n = len(counts)
    ch = np.arange(n)
    lo = poly(cal, ch - 0.5)
    hi = poly(cal, ch + 0.5)
    if corr:
        lo = poly(corr, lo)
        hi = poly(corr, hi)
    out = np.zeros(len(grid_edges) - 1)
    g0, g1 = grid_edges[0], grid_edges[-1]
    step = grid_edges[1] - grid_edges[0]
    for i in range(n):
        c = counts[i]
        if c <= 0:
            continue
        a, b = lo[i], hi[i]
        if b <= g0 or a >= g1 or b <= a:
            continue
        a, b = max(a, g0), min(b, g1)
        ka = int((a - g0) / step)
        kb = int((b - g0) / step)
        if ka == kb:
            out[ka] += c
            continue
        w = b - a
        for k in range(ka, min(kb + 1, len(out))):
            left = max(a, g0 + k * step)
            right = min(b, g0 + (k + 1) * step)
            if right > left:
                out[k] += c * (right - left) / w
    return out


# Форма аппаратной линии. TAIL_L и TAIL_R — точки перехода от гауссианы к
# экспоненте, в единицах сигмы; None означает чистую гауссиану.
#
# Значения по умолчанию взяты из самого файла замера: там записана калибровка
# ПШПВ, сделанная в программе (узел SimpleSqrtFwhmCalibration), с моделью пика
# ExpGaussExp и параметрами ExpGaussExpLeftTail = 1,10, ExpGaussExpRightTail =
# 1,70 при хи²/ndf = 2,89 по восьми опорным пикам. Это НЕ паспорт производителя,
# а подгонка в программе; кем и когда сделана, из файла не видно. Берётся как
# эмпирическое описание формы линии ЭТОГО прибора, с той же силой, что любая
# другая подгонка.
TAIL_L = 1.10
TAIL_R = 1.70


def line_shape(x, e, s, tail_l=TAIL_L, tail_r=TAIL_R):
    """ExpGaussExp: гауссиана, переходящая в экспоненты за |x−e| > tail·сигма.

    Сшивка непрерывна и по значению, и по производной — это и есть смысл
    параметра перехода: exp(t²/2 ∓ t·u) при u = (x−e)/сигма совпадает с
    exp(−u²/2) в точке u = ∓t вместе с наклоном.
    """
    u = (x - e) / s
    g = np.exp(-0.5 * u ** 2)
    if tail_l:
        m = u < -tail_l
        if m.any():
            g[m] = np.exp(0.5 * tail_l ** 2 + tail_l * u[m])
    if tail_r:
        m = u > tail_r
        if m.any():
            g[m] = np.exp(0.5 * tail_r ** 2 - tail_r * u[m])
    return g


def fwhm_from_file(path):
    """Кривая ПШПВ(E) из самого файла замера — узел SimpleSqrtFwhmCalibration.

    Восемь опорных точек калибровки превращаются в интерполятор: ПШПВ в
    каналах умножается на локальную производную dE/dch заводской шкалы, между
    точками — линейная интерполяция по log E, вне — степенная экстраполяция по
    краевой паре.

    Двучлен ПШПВ² = f0 + f1·E эти точки не описывает: подгонка по всем восьми
    даёт f0 = −204, то есть ширина обращается в нуль около 62 кэВ, а невязка
    на мягком крае достигает +7,8 кэВ. Поэтому берётся сама кривая.

    Кривая задаёт ХОД ширины по энергии, но не сами значения: она снята по
    другому спектру и в этом замере не воспроизводится. Прямой промер пиков
    измеренного спектра даёт ПШПВ 74 кэВ на 2614,5 против 92,4 по кривой и
    44 против 41,7 на 911,2 — то есть ход и величина расходятся. Поэтому
    кривая берётся как опорная форма, а поверх неё подгоняются два параметра
    (см. fit_width): множитель k и показатель наклона p,

        ПШПВ(E) = k · ПШПВ_файла(E) · (E / 662)^p.

    Отпускать ширину произвольным двучленом ПШПВ² = f0 + f1·E нельзя: при
    свободной паре подгонка уходила на ПШПВ(662) = 26,7 кэВ, вдвое ниже
    собственного разрешения прибора.
    """
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    i = s.find("<SimpleSqrtFwhmCalibration>")
    if i < 0:
        return None
    blk = s[i:s.find("</SimpleSqrtFwhmCalibration>", i)]
    pk = re.findall(r"<Channel>([\d.]+)</Channel>\s*"
                    r"<Energy>([\d.eE+-]+)</Energy>\s*"
                    r"<FWHM>([\d.eE+-]+)</FWHM>", blk)
    if len(pk) < 3:
        return None
    ch = np.array([float(a) for a, _, _ in pk])
    en = np.array([float(b) for _, b, _ in pk])
    fw_ch = np.array([float(c) for _, _, c in pk])
    dEdch = np.polyval(np.polyder(np.polyfit(ch, en, 3)), ch)
    fw = fw_ch * dEdch
    order = np.argsort(en)
    en, fw = en[order], fw[order]
    lg_e, lg_f = np.log(en), np.log(fw)
    k_lo = (lg_f[1] - lg_f[0]) / (lg_e[1] - lg_e[0])
    k_hi = (lg_f[-1] - lg_f[-2]) / (lg_e[-1] - lg_e[-2])

    def curve(e):
        e = np.asarray(e, float)
        le = np.log(np.maximum(e, 1e-6))
        out = np.interp(le, lg_e, lg_f)
        out = np.where(le < lg_e[0], lg_f[0] + k_lo * (le - lg_e[0]), out)
        out = np.where(le > lg_e[-1], lg_f[-1] + k_hi * (le - lg_e[-1]), out)
        return np.exp(out)

    return curve


def broaden(raw_e, raw_c, grid_centres, f0, f1, tail_l=TAIL_L, tail_r=TAIL_R,
            fwhm_fn=None):
    """Свёртка линейного спектра с аппаратной формой.

    Ширина берётся из fwhm_fn(E), если она задана (кривая прибора из файла
    замера), иначе из двучлена ПШПВ² = f0 + f1·E. Ядро нормируется на единицу
    площади, поэтому полное число отсчётов шаблона не зависит от того, есть
    хвосты или нет: хвосты только перекладывают отсчёты из пика в подножие, а
    это ровно тот эффект, который проверяется.
    """
    out = np.zeros(len(grid_centres))
    step = grid_centres[1] - grid_centres[0]
    lo0 = grid_centres[0] - 0.5 * step
    # с хвостами ядро тянется дальше гауссова: экспонента с показателем tail
    # спадает в e раз на сигму, восьми сигм хватает на четыре порядка
    reach = 8.0 if (tail_l or tail_r) else 4.0
    for e, c in zip(raw_e, raw_c):
        if c <= 0:
            continue
        if fwhm_fn is not None:
            fw = float(fwhm_fn(e))
        else:
            # Вырожденная пара (0, 0) означает, что вызывающий забыл передать
            # кривую ширины: молча свернуть с ядром в один килоэлектронвольт —
            # значит выдать линейчатый спектр за уширенный. Так уже случилось
            # с заливкой компонент на рисунке.
            if f0 == 0.0 and f1 == 0.0:
                raise ValueError("broaden: не задана ширина линии — ни fwhm_fn, "
                                 "ни коэффициенты f0, f1")
            fw = math.sqrt(max(f0 + f1 * e, 1.0))
        fw = max(fw, 1.0)
        s = fw / 2.3548
        k0 = int((e - reach * s - lo0) / step)
        k1 = int((e + reach * s - lo0) / step) + 1
        k0 = max(k0, 0)
        k1 = min(k1, len(out))
        if k1 <= k0:
            continue
        x = grid_centres[k0:k1]
        g = line_shape(x, e, s, tail_l, tail_r)
        ssum = g.sum()
        if ssum <= 0:
            continue
        out[k0:k1] += c * g / ssum
    return out


def nnls_fit(A, y, w):
    """Взвешенный МНК с неотрицательными коэффициентами."""
    from scipy.optimize import nnls
    Aw = A * w[:, None]
    yw = y * w
    x, rnorm = nnls(Aw, yw)
    return x, rnorm


def broaden_var(raw_e, raw_c_norm, n_prim, grid_centres, f0, f1,
                tail_l=TAIL_L, tail_r=TAIL_R):
    """Дисперсия уширенного шаблона от статистики Монте-Карло, на распад².

    Шаблон — линейная комбинация пуассоновских счётчиков: T_i = Σ_j c_j·K_ij/N,
    где K — нормированное ядро формы линии. Отсюда Var[T_i] = Σ_j c_j·K_ij²/N².
    При хранящемся нормированном c_j/N это Σ_j (c_j/N)·K_ij²/N — то же ядро,
    но ВОЗВЕДЁННОЕ В КВАДРАТ, и деление на N один раз.
    """
    out = np.zeros(len(grid_centres))
    step = grid_centres[1] - grid_centres[0]
    lo0 = grid_centres[0] - 0.5 * step
    reach = 8.0 if (tail_l or tail_r) else 4.0
    for e, c in zip(raw_e, raw_c_norm):
        if c <= 0:
            continue
        fw = math.sqrt(max(f0 + f1 * e, 1.0))
        s = fw / 2.3548
        k0 = max(int((e - reach * s - lo0) / step), 0)
        k1 = min(int((e + reach * s - lo0) / step) + 1, len(out))
        if k1 <= k0:
            continue
        x = grid_centres[k0:k1]
        g = line_shape(x, e, s, tail_l, tail_r)
        ssum = g.sum()
        if ssum <= 0:
            continue
        out[k0:k1] += c * (g / ssum) ** 2
    return out / max(n_prim, 1.0)


def lsrm_solve(A, VarA, y_raw, bg_scaled, n_iter=3):
    """Решение разложения по ЛСРМ «Алгоритмические основы», гл. 12.

    Минимизируется хи²-функционал (формула 12.1) с весами, включающими ПОЛНУЮ
    дисперсию каждой точки:

        W_i = 1 / { (S_i + f_i·t/t_f) + Σ_k (Ã_k · σr_ik)² },

    где первая скобка — счётная дисперсия образца и вычтенного фона, вторая —
    дисперсия матрицы чувствительности (у нас — статистика Монте-Карло
    шаблонов) при текущей оценке активностей Ã_k. Оценки входят в веса,
    поэтому решение итерационное: на первом проходе Ã_k = 0 (§12.1).

    Условие минимума даёт нормальную систему D·A = b с
        D_lk = Σ_i r_il·r_ik·W_i,   b_l = Σ_i S_i^f·r_il·W_i,
    случайная неопределённость — из диагонали обратной матрицы:
        σ(A_k) = √(D⁻¹)_kk                                  (формула 12.1-2).

    Система решается КАК ЛИНЕЙНАЯ, без ограничения знака. Отрицательный
    коэффициент по гл. 14 означает «нуклид не идентифицирован» — такая
    компонента исключается из набора, и система решается заново.

    A      — матрица столбцов-шаблонов (отсчёты за время набора на 1 Бк);
    VarA   — их дисперсии той же размерности (отсчёты² на 1 Бк²);
    y_raw  — измеренный спектр (до вычета фона);
    bg_scaled — фон, приведённый к времени образца.

    Возвращает (x, sigma, chi2_per_dof, active_mask).
    """
    m_pts, n_cmp = A.shape
    resid = y_raw - bg_scaled
    # Компонента, чей шаблон в окне подгонки пуст (все линии за окном),
    # делает нормальную матрицу вырожденной — исключается сразу, с тем же
    # смыслом «в этом окне не определяется».
    active = A.any(axis=0)
    x = np.zeros(n_cmp)
    sigma = np.zeros(n_cmp)
    chi2_dof = float("nan")
    for _ in range(6):                       # внешние: исключение отрицательных
        x_try = np.zeros(n_cmp)
        for _ in range(n_iter):              # внутренние: веса от оценок Ã
            var = y_raw + bg_scaled
            if x_try.any():
                var = var + (VarA * (x_try ** 2)[None, :]).sum(axis=1)
            W = 1.0 / np.maximum(var, 1.0)
            Aa = A[:, active]
            D = (Aa * W[:, None]).T @ Aa
            b = (Aa * W[:, None]).T @ resid
            try:
                sol = np.linalg.solve(D, b)
            except np.linalg.LinAlgError:
                sol = np.linalg.lstsq(D, b, rcond=None)[0]
            x_try = np.zeros(n_cmp)
            x_try[active] = sol
        neg = active & (x_try < 0)
        if not neg.any():
            x = x_try
            try:
                Dinv = np.linalg.inv(D)
            except np.linalg.LinAlgError:
                Dinv = np.linalg.pinv(D)
            sigma = np.zeros(n_cmp)
            sigma[active] = np.sqrt(np.maximum(np.diag(Dinv), 0.0))
            r = (A[:, active] @ x[active]) - resid
            chi2 = float((r ** 2 * W).sum())
            dof = max(1, m_pts - int(active.sum()))
            chi2_dof = chi2 / dof
            # §12.1: если хи² превышает число степеней свободы, случайная
            # неопределённость раздувается в корень из хи²/dof раз
            if chi2_dof > 1.0:
                sigma = sigma * math.sqrt(chi2_dof)
            break
        active[neg] = False
        if not active.any():
            break
    return x, sigma, chi2_dof, active


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, tdir = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(src)
    os.makedirs(outdir, exist_ok=True)

    spec = read_atomspectra_xml(src)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    t_smp = float(spec.real_time)

    corr = read_correction(os.path.join(outdir, "calibration_fitted.csv"))
    if not corr:
        print("ВНИМАНИЕ: поправок калибровки нет — работаем по заводской шкале")

    edges = np.arange(0.0, E_MAX + E_STEP, E_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])

    y = rebin_to_grid(np.asarray(spec.counts, float), list(spec.energy_cal),
                      corr.get("sample"), edges)
    print("образец: %d отсчётов в сетке из %d (в файле %d)"
          % (int(y.sum()), len(centres), int(np.sum(spec.counts))))

    if bg is not None:
        ybg = rebin_to_grid(np.asarray(bg.counts, float),
                            list(bg.energy_cal), corr.get("background"), edges)
        ybg = ybg * (t_smp / float(bg.real_time))     # к времени образца
        print("фон: %.0f отсчётов, приведён к %.0f с (было %.0f с)"
              % (ybg.sum(), t_smp, bg.real_time))
    else:
        ybg = np.zeros_like(y)

    # --- шаблоны --------------------------------------------------------------
    templates, names, colours, stamps = [], [], [], set()
    head0 = {}
    for key, label, colour in ORDER:
        p = os.path.join(tdir, "%s.csv" % key)
        if not os.path.exists(p):
            continue
        head, e, c = read_template(p)
        head0 = head0 or head
        n = float(head["N_primaries"])
        if c.sum() == 0:
            print("  %s: ни одного события — компонента пропущена" % label)
            continue
        stamps.add(head.get("src_sha1", "?"))
        templates.append((key, label, colour, e, c / n, n))
        names.append(label)
        colours.append(colour)
    if len(stamps) > 1:
        raise SystemExit("шаблоны разных ревизий: %s" % ", ".join(stamps))
    print("шаблонов %d, штамп %s" % (len(templates), stamps.pop() if stamps
                                     else "?"))

    m = (centres >= E_FIT[0]) & (centres <= E_FIT[1])
    yy = y[m]
    wgt = 1.0 / np.sqrt(np.maximum(yy, 1.0))          # пуассоновские веса

    # --- сборка компонент подгонки -------------------------------------------
    # Свободная подгонка ПОНУКЛИДНО не годится: шаблоны Th-232, Th-228, Ra-224
    # набраны единицами событий (2, 10 и 100 при 200 тыс. розыгрышей — у этих
    # звеньев почти нет гамма-выхода), и подгонка тянет их множители в сотни
    # килобеккерелей, подгоняя ими шум. Поэтому ряд собирается в ДВЕ подцепочки,
    # внутри которых равновесие обеспечено самими периодами полураспада:
    #
    #   A1 — верхняя часть: Th-232 -> Ra-228 -> Ac-228 (сигнал даёт Ac-228);
    #   A2 — нижняя часть: Th-228 -> Ra-224 -> Rn-220 -> Po-216 -> Pb-212 ->
    #        Bi-212 -> (Tl-208 | Po-212). Все периоды от 0,15 с до 3,6 суток,
    #        равновесие достигается за недели.
    #
    # Ветвление Bi-212 -> Tl-208 берётся из библиотеки МАГАТЭ, а не пишется
    # числом. Отношение A1/A2 остаётся СВОБОДНЫМ: именно оно показывает,
    # нарушено ли равновесие между Ra-228 и Th-228.
    br_tl = branch_to_tl208()
    print("ветвление Bi-212 -> Tl-208: %.2f %% (МАГАТЭ)" % (100.0 * br_tl))
    # Основные подцепочки — как раньше. XI и XW — самостоятельные компоненты
    # (флуоресценция иода в кристалле и K-серия вольфрама, вносимая источником);
    # их коэффициенты подгонки — темп эмиссии рентгена, эмиссий/с. В подгонке
    # они снимают вклад, который иначе уходил бы в подложку и тянул A1 вверх.
    GROUP = {
        "A1 (Ra-228 -> Ac-228)": {"Th232": 1.0, "Ra228": 1.0, "Ac228": 1.0},
        "A2 (Th-228 -> Tl-208)": {"Th228": 1.0, "Ra224": 1.0, "Rn220": 1.0,
                                  "Po216": 1.0, "Pb212": 1.0, "Bi212": 1.0,
                                  "Tl208": br_tl, "Po212": 1.0 - br_tl},
    }
    tmap = {k: (e, c) for k, _, _, e, c, _ in templates}
    nmap = {k: n for k, _, _, _, _, n in templates}
    XGROUPS = [("XI (флуоресценция иода в CsI)", "XI"),
               ("XW (K-серия вольфрама электрода)", "XW"),
               ("XD1 (K-серия дочерних, ветвь A1)", "XD1"),
               ("XD2 (K-серия дочерних, ветвь A2)", "XD2")]
    for gname, k in XGROUPS:
        if k in tmap:
            GROUP[gname] = {k: 1.0}

    # --- ширина линии: двучлен ПШПВ² = f₀ + f₁·E, перебор по сетке ----------
    # Оба параметра перебираются ВМЕСТЕ с активностями, критерий — хи² модели и
    # измерения (директива оператора, 07.08.2026). Две другие параметризации
    # пробовались и отвергнуты:
    #   — кривая из файла замера как есть: снята по другому спектру, на
    #     2614,51 кэВ даёт 92 кэВ против 64-80 по промеру этого спектра;
    #   — та же кривая с подгоняемыми множителем и наклоном: критерий по
    #     невязке уводил ширину на 94 кэВ, компенсируя ею избыток модельного
    #     континуума вместо того, чтобы измерять разрешение.
    # Кривая из файла оставлена как независимая справка: печатается рядом с
    # подогнанной, расхождение видно числом.
    fwc0 = fwhm_from_file(src)
    if fwc0 is not None:
        print("справка, калибровка ПШПВ из файла замера: %.1f кэВ на 662, "
              "%.1f на 2614" % (float(fwc0(661.657)), float(fwc0(2614.51))))

    # Опорные линии ширины — только ОДИНОЧНЫЕ. Линия 911,20 сюда не годится:
    # рядом 964,77 и 968,97, при разрешении прибора они сливаются в один
    # бугор, и промер по полувысоте даёт 60 кэВ вместо собственной ширины
    # линии — подгонка по такому якорю растягивала всю кривую.
    W_ANCHORS = [238.63, 583.19, 2614.51]

    def measure_fwhm(v, e0, guess):
        """ПШПВ пика ПОДГОНКОЙ ПРОФИЛЯ в окне, кэВ.

        Промер по полувысоте над медианой краёв окна здесь непригоден: на
        мягких линиях пик сидит на крутом спаде континуума, и «подложка»,
        взятая одним числом, срезает пик несимметрично. На 238,63 кэВ такой
        промер давал 22 кэВ (9,2 % от энергии) при ожидаемых для CsI(Tl)
        13-15 %, и подогнанная по нему кривая ширины выходила заниженной.

        Здесь в окне подгоняется модель «линия + ЛИНЕЙНАЯ подложка»:

            v(E) = A · L((E − E₀)/σ) + a + b·(E − E₀),

        где L — та же форма ExpGaussExp, что в свёртке шаблонов. Линейные
        параметры A, a, b при каждом (E₀, σ) находятся точно, методом
        наименьших квадратов; перебираются только два нелинейных. ПШПВ
        возвращается как 2,3548·σ.
        """
        half = max(2.5 * guess, 40.0)
        sel = (centres > e0 - half) & (centres < e0 + half)
        x, yv = centres[sel], v[sel].astype(float)
        if len(x) < 9:
            return None
        w = 1.0 / np.sqrt(np.maximum(np.abs(yv), 1.0))
        best = None
        for s in np.linspace(0.25 * guess, 1.2 * guess, 40) / 2.3548:
            for de in np.linspace(-0.4 * guess, 0.4 * guess, 9):
                L = line_shape(x, e0 + de, s)
                M = np.column_stack([L, np.ones_like(x), x - e0])
                Mw = M * w[:, None]
                try:
                    c, *_ = np.linalg.lstsq(Mw, yv * w, rcond=None)
                except np.linalg.LinAlgError:
                    continue
                if c[0] <= 0:                     # амплитуда линии должна быть
                    continue                      # положительной, иначе это не пик
                r = (M @ c - yv) * w
                q = float((r ** 2).sum())
                if best is None or q < best[0]:
                    best = (q, s)
        if best is None:
            return None
        return 2.3548 * best[1]

    # Ширина — канонический двучлен ПШПВ² = f₀ + f₁·E. Ниже по коду параметры
    # называются (k, p) только потому, что так подписаны аргументы кэша; их
    # смысл — f₀ и f₁.
    def fwc_kp(e, f0, f1):
        return np.sqrt(np.maximum(f0 + f1 * np.asarray(e, float), 1.0))

    # Свёртка считается ОДИН раз на полной сетке для каждой пары (k, p), дальше
    # берутся срезы. Считать broaden на обрезанной подсетке нельзя: ядро
    # нормируется на свою сумму, и у края подсетки оно теряет хвосты — пик
    # выходит выше. Из-за этого ступень 1 садилась на завышенную модель, и пик
    # 2614 при сверке на полной сетке давал 1,29 вместо единицы.
    _gcache = {}

    def gfull(gname, f0, f1):
        key = (gname, round(f0, 6), round(f1, 6))
        if key not in _gcache:
            acc = np.zeros(len(centres))
            for nk, wgt_k in GROUP[gname].items():
                if nk in tmap:
                    e, c = tmap[nk]
                    acc += wgt_k * broaden(e, c, centres, f0, f1)
            _gcache[key] = acc * t_smp
        return _gcache[key]

    # Дисперсия столбцов от статистики Монте-Карло — для весов ЛСРМ §12.1.
    # Синтетические линейчатые шаблоны (XI/XW/XD, N_primaries = 1) статистики
    # не несут: их дисперсия нулевая, а не 1/1.
    _vcache = {}

    def gvar(gname, f0, f1):
        key = (gname, round(f0, 6), round(f1, 6))
        if key not in _vcache:
            acc = np.zeros(len(centres))
            for nk, wgt_k in GROUP[gname].items():
                if nk in tmap and nmap.get(nk, 1.0) > 1.0:
                    e, c = tmap[nk]
                    acc += wgt_k ** 2 * broaden_var(e, c, nmap[nk],
                                                    centres, f0, f1)
            _vcache[key] = acc * t_smp ** 2
        return _vcache[key]

    def build(k, p):
        return np.column_stack([gfull(g, k, p)[m] for g in GROUP])

    def build_var(k, p):
        return np.column_stack([gvar(g, k, p)[m] for g in GROUP])

    gnames = list(GROUP)
    resid_full = y - ybg
    gA1 = next(g for g in gnames if g.startswith("A1"))
    gA2 = next(g for g in gnames if g.startswith("A2"))

    if GLOBAL_FIT:
        # --- единая подгонка по всему окну (умолчание) ----------------------
        # Решение — по ЛСРМ «Алгоритмические основы», гл. 12: нормальная
        # система D·A = b с итерационными весами, включающими дисперсию
        # шаблонов, и σ(A) из диагонали D⁻¹ (см. lsrm_solve). ПШПВ² = f₀ + f₁·E
        # перебирается по сетке ВМЕСТЕ с активностями, критерий — хи²/dof.
        # ФОН НЕ ПОДГОНЯЕТСЯ: он измерен тем же прибором и приведён к времени
        # образца. Свободный множитель фона в пробной подгонке уходил на 4,8 —
        # затыкал фоном нехватку континуума, то есть лечил симптом.
        best = None
        for kk in np.arange(-600.0, 601.0, 100.0):        # f₀
            for pp in np.arange(1.0, 5.01, 0.2):          # f₁
                if kk + pp * E_FIT[0] <= 4.0:
                    continue
                A = build(kk, pp)
                VarA = build_var(kk, pp)
                x, sig, chi2, act = lsrm_solve(A, VarA, yy, ybg[m])
                if not act.any():
                    continue
                if best is None or chi2 < best[0]:
                    best = (chi2, kk, pp, x, A, sig, act)
        chi2, k_fw, p_fw, x, A, x_sigma, x_active = best
        width_rows = []
        print("\n[единая подгонка по окну %.0f-%.0f кэВ, решение нормальной "
              "системы по ЛСРМ гл. 12]" % E_FIT)
    else:
        # --- ступенчатая подгонка с жёсткого конца (умолчание) --------------
        ms = (centres >= W_SOFT[0]) & (centres <= W_SOFT[1])
        gs = centres[ms]
        rs = resid_full[ms]
        ws = 1.0 / np.sqrt(np.maximum(y[ms], 1.0))

        def net_in(v, grid, e0, fw):
            """Нетто в окне ±1 ПШПВ с подложкой по внешним подокнам.

            Считается ОДИНАКОВО у измеренного и у модели — иначе коэффициент
            вбирает разницу подложек, а не пиков.
            """
            sel = (grid >= e0 - fw) & (grid <= e0 + fw)
            bl = (grid >= e0 - 2.0 * fw) & (grid < e0 - fw)
            br = (grid > e0 + fw) & (grid <= e0 + 2.0 * fw)
            if sel.sum() < 3 or bl.sum() < 2 or br.sum() < 2:
                return None
            b = 0.5 * (v[bl].mean() + v[br].mean())
            return float(v[sel].sum() - b * sel.sum())

        def scale_by_peak(T, r, grid, e0, fw):
            """Коэффициент, при котором нетто модели равно нетто измеренного.

            Подгонка по сырому сигналу здесь не годится: в окне жёсткой линии
            подложка (наложения, космика, суммирование) сравнима с самим пиком
            и в модели описана иначе, чем в измерении. Взвешенный МНК по всему
            окну прижимал коэффициент к подложке, и пик 2614 выходил впятеро
            ниже измеренного. Приравнивание НЕТТО ставит пик на место по
            построению — а расхождение подложек остаётся видимым в невязке.
            """
            nm = net_in(T, grid, e0, fw)
            ny = net_in(r, grid, e0, fw)
            if not nm or nm <= 0 or ny is None:
                return 0.0
            return max(0.0, ny / nm)

        def band(e0, fw):
            """Маска ±2,5 ПШПВ вокруг линии — вмещает пик и оба подокна."""
            return (centres >= e0 - 2.5 * fw) & (centres <= e0 + 2.5 * fw)

        def steps12(k, p):
            """Ступени 1 и 2 при данной ширине: (a1, a2)."""
            fwh, fwm = float(fwc_kp(E_HARD, k, p)), float(fwc_kp(E_MID, k, p))
            mh, mm_ = band(E_HARD, fwh), band(E_MID, fwm)
            gh, gm = centres[mh], centres[mm_]
            rh, rm = resid_full[mh], resid_full[mm_]
            Th2, Th1 = gfull(gA2, k, p)[mh], gfull(gA1, k, p)[mh]
            Tm2, Tm1 = gfull(gA2, k, p)[mm_], gfull(gA1, k, p)[mm_]
            # Ступени зацеплены: в окно жёсткой линии попадает вклад верхней
            # ветви (суммирование каскада Ac-228), а в окно 911 — вклад нижней.
            # Без итерации пик 2614 выходил на 1,29 измеренного вместо единицы.
            # Пяти проходов хватает: поправка на пятом ниже 0,1 %.
            a1 = a2 = 0.0
            for _ in range(5):
                a2 = scale_by_peak(Th2, rh - a1 * Th1, gh, E_HARD, fwh)
                a1 = scale_by_peak(Tm1, rm - a2 * Tm2, gm, E_MID, fwm)
            return a1, a2

        # --- ширина: по совпадению формы измеренного и ПОЛНОЙ модели ---------
        # Сравнивать надо не с одиночной гауссианой, а с моделью: в окно любой
        # линии попадают соседи, и подгонка изолированного профиля приписывает
        # их вклад ширине самой линии. Промер профиля на 238,63 кэВ давал
        # 19,0 кэВ (относительное разрешение 8,0 %), тогда как из 583,19 и
        # 2614,51 следует около 9,5 %: разницу создают линии 240,99 (Ra-224),
        # 270,85 (Ac-228) и 277,37 кэВ (Tl-208), попадающие в то же окно.
        # Модель их содержит, поэтому мерой служит невязка «измеренное минус
        # модель». Масштаб при этом задан нетто опорных пиков и от ширины почти
        # не зависит, так что остаётся именно ширина.
        wsel = np.zeros(len(centres), bool)
        for e0 in W_ANCHORS + [E_MID]:
            g0 = float(fwc0(e0)) if fwc0 is not None else math.sqrt(2.0 * e0)
            wsel |= (centres >= e0 - 2.0 * g0) & (centres <= e0 + 2.0 * g0)
        r_obs = resid_full[wsel]
        w_obs = 1.0 / np.sqrt(np.maximum(y[wsel], 1.0))
        best = None
        for kk in np.arange(-600.0, 601.0, 100.0):        # f₀
            for pp in np.arange(1.0, 5.01, 0.2):          # f₁
                if kk + pp * W_SOFT[0] <= 4.0:
                    continue
                a1c, a2c = steps12(kk, pp)
                if a1c <= 0 or a2c <= 0:
                    continue
                mod = (a1c * gfull(gA1, kk, pp) + a2c * gfull(gA2, kk, pp))[wsel]
                q = float((((mod - r_obs) * w_obs) ** 2).sum()) / max(
                    1, int(wsel.sum()) - 4)
                if best is None or q < best[0]:
                    best = (q, kk, pp, a1c, a2c)
        if best is None:
            raise SystemExit("ступенчатая подгонка не сошлась: нетто опорного "
                             "пика в модели или в измерении неположительно")
        chi2_hm, k_fw, p_fw, a1, a2 = best

        # Промер ширины опорных линий — уже ПРОВЕРКА, а не вход подгонки:
        # сравнивается ширина линии в измеренном спектре и в модели, снятая
        # одинаково.
        width_rows = []
        mdl_w = a1 * gfull(gA1, k_fw, p_fw) + a2 * gfull(gA2, k_fw, p_fw)
        for e0 in W_ANCHORS:
            g0 = float(fwc_kp(e0, k_fw, p_fw))
            fo = measure_fwhm(resid_full, e0, g0)
            fm = measure_fwhm(mdl_w, e0, g0)
            if fo and fm:
                width_rows.append((e0, fo, fm))
        print("\nширина линии: ПШПВ² = %.0f + %.2f·E" % (k_fw, p_fw))
        print("  %8s %11s %11s %8s" % ("E, кэВ", "изм., кэВ", "модель", "Δ, %"))
        for e0, fo, fm in width_rows:
            print("  %8.2f %11.1f %11.1f %8.1f"
                  % (e0, fo, fm, 100.0 * (fm - fo) / fo))

        # ступень 3: рентгеновские компоненты на мягком участке при
        # зафиксированных A1 и A2
        xg = [g for g in gnames if g.startswith(("XI", "XW", "XD"))]
        xvals = {g: 0.0 for g in xg}
        if xg:
            base_s = (a1 * gfull(gA1, k_fw, p_fw)
                      + a2 * gfull(gA2, k_fw, p_fw))[ms]
            As = np.column_stack([gfull(g, k_fw, p_fw)[ms] for g in xg])
            xs, _ = nnls_fit(As, rs - base_s, ws)
            xvals = dict(zip(xg, xs))

        x = np.array([{**{gA1: a1, gA2: a2}, **xvals}.get(g, 0.0)
                      for g in gnames])
        A = build(k_fw, p_fw)
        # у ступенчатой ветки строгих σ нет (задача №35) — нули как признак
        x_sigma = np.zeros(len(x))
        x_active = x > 0
        print("\n[ступенчатая подгонка: %.0f -> %.0f -> %.0f-%.0f кэВ]"
              % (E_HARD, E_MID, W_SOFT[0], W_SOFT[1]))
        print("хи²/n по окнам опорных линий = %.1f" % chi2_hm)
        r_all = (A @ x + ybg[m] - yy) * wgt
        chi2 = float((r_all ** 2).sum()) / max(1, len(yy) - len(x))

    def fw_at(e):
        return float(fwc_kp(e, k_fw, p_fw))

    fw662 = fw_at(661.657)
    print("\nПШПВ² = %.0f + %.2f·E  (ПШПВ 662 кэВ = %.1f кэВ, %.2f %%), "
          "хи²/n по полосе %.0f-%.0f = %.1f"
          % (k_fw, p_fw, fw662, 100.0 * fw662 / 661.657,
             E_FIT[0], E_FIT[1], chi2))

    model = A @ x + ybg[m]
    # Идентификация — по гл. 14 ЛСРМ: параметр достоверности δa = σ(A)/A.
    # При δa < 1 нуклид считается идентифицированным (при доверительной
    # вероятности весов это допускает вероятностную трактовку, §14.1);
    # исключённая из решения компонента (a < 0) — не идентифицирована.
    print("\n--- коэффициенты разложения (A ± σ, δa = σ/A) ---")
    rows = []
    a_by_gname = {}
    for i, gname in enumerate(GROUP):
        unit = "эмиссий/с" if gname.startswith(("XI", "XW")) else "Бк"
        if not x_active[i]:
            verdict = "не идентифицирован (a < 0)"
            print("  %-38s %s" % (gname, verdict))
        else:
            da = x_sigma[i] / x[i] if x[i] > 0 else float("inf")
            verdict = ("идентифицирован" if da < 1.0 else
                       "не идентифицирован (δa ≥ 1)")
            if x_sigma[i] > 0:
                print("  %-38s %10.0f ± %6.0f %s  δa=%.3f  %s"
                      % (gname, x[i], x_sigma[i], unit, da, verdict))
            else:
                print("  %-38s %10.0f %s" % (gname, x[i], unit))
        rows.append((gname, x[i], x_sigma[i]))
        a_by_gname[gname] = x[i]
    a1 = next((v for k, v in a_by_gname.items() if k.startswith("A1")), 0.0)
    a2 = next((v for k, v in a_by_gname.items() if k.startswith("A2")), 0.0)
    if a2 > 0:
        s1 = next((s for (k, v, s) in rows if k.startswith("A1")), 0.0)
        s2 = next((s for (k, v, s) in rows if k.startswith("A2")), 0.0)
        ratio_v = a1 / a2
        # неопределённость отношения — по независимым σ; ковариация A1, A2
        # из D⁻¹ здесь не учитывается, что помечается явно
        dr = (ratio_v * math.hypot(s1 / a1 if a1 else 0,
                                   s2 / a2 if a2 else 0)
              if (s1 or s2) else 0.0)
        print("  отношение A1/A2 = %.3f%s (равновесие ряда -> 1,000)"
              % (ratio_v, " ± %.3f (без ковариации)" % dr if dr else ""))

    # --- удельная активность и сверка с номиналом этикетки ------------------
    # Масса пачки берётся ИЗ ШАПКИ ШАБЛОНА, а не пересчитывается здесь: считает
    # её геометрия по построенным телам, и второй счёт в другом месте — это
    # ровно тот случай, когда числа расходятся молча.
    mass_g = float(head0.get("wt20_mass_g", "0").split()[0])
    # Полусумму A1 и A2 здесь считать НЕЛЬЗЯ. Прежняя редакция брала
    # 0,5·(A1+A2) с оговоркой «если равновесен» и делила на номинал Th-232,
    # хотя строкой выше сама печатала A1/A2 = 0,60. Средним двух неравных
    # активностей подменялась величина, которой в этом случае нет: при
    # нарушенном равновесии активность ряда одним числом не описывается.
    # Ниже печатаются ОБЕ ветви порознь, каждая со своей долей номинала, а
    # сводное число даётся только при сошедшемся равновесии.
    ratio = a1 / a2 if a2 > 0 else float("nan")
    equilibrium = abs(ratio - 1.0) <= 0.10       # заведомо мягкий порог
    print("\n--- удельная активность ---")
    if mass_g > 0:
        # Номинал этикетки: 2 % масс. ThO2, доля тория в ThO2 0,878809,
        # удельная активность Th-232 4072 Бк/г (вычислена из T1/2 МАГАТЭ).
        th_g = mass_g * 0.02 * 0.878809
        a_nom = th_g * SPEC_ACT_TH232
        print("  масса пачки %.1f г (из шапки шаблона)" % mass_g)
        print("  ПО ЭТИКЕТКЕ (2 %% ThO2): тория %.3f г -> %.0f Бк на звено"
              % (th_g, a_nom))
        for gname, a, _s in rows:
            if gname.startswith(("XI", "XW")):
                # эмиссии/с — не активность цепочки, к массе не приводится
                print("  %-38s %8.0f эмиссий/с" % (gname, a))
                continue
            print("  %-38s %8.0f Бк | %6.0f Бк/кг | %.3f номинала"
                  % (gname, a, 1000.0 * a / mass_g,
                     a / a_nom if a_nom else 0))
        if equilibrium:
            a_meas = 0.5 * (a1 + a2)
            print("  равновесие сошлось (A1/A2 = %.3f), ряд в целом:" % ratio)
            print("    %.0f Бк, %.0f Бк/кг, %.3f номинала, "
                  "эквивалент %.2f %% масс. ThO2"
                  % (a_meas, 1000.0 * a_meas / mass_g, a_meas / a_nom,
                     2.0 * a_meas / a_nom))
        else:
            print("  РАВНОВЕСИЕ НЕ СОШЛОСЬ: A1/A2 = %.3f." % ratio)
            print("  Сводной «активности ряда» и эквивалентного содержания")
            print("  ThO2 в этом случае нет — две ветви называются порознь.")
            print("  Избыток Th-228 над Ra-228 в замкнутой системе после")
            print("  химической очистки тория невозможен: там Ra-228 идёт")
            print("  впереди. Значит либо радий потерян на переделе, либо")
            print("  занижена A1, которая держится на линиях Ac-228.")

    # --- сверка по опорным пикам: попала ли модель в измеренное -------------
    # Смысл ступенчатой подгонки в том, что жёсткий пик обязан вписаться. Это
    # проверяется числом, а не глазом по рисунку: берётся площадь в окне
    # ±1 ПШПВ за вычетом подложки по внешним подокнам — одинаково у измеренного
    # и у модели.
    def peak_net(v, e0):
        fw = fw_at(e0)
        sel = (centres >= e0 - fw) & (centres <= e0 + fw)
        bl = (centres >= e0 - 2.0 * fw) & (centres < e0 - fw)
        br = (centres > e0 + fw) & (centres <= e0 + 2.0 * fw)
        if sel.sum() < 3 or bl.sum() < 2 or br.sum() < 2:
            return None
        b = 0.5 * (v[bl].mean() + v[br].mean())
        return float(v[sel].sum() - b * sel.sum())

    yfull = y - ybg
    mfull = np.zeros(len(centres))
    for gi, gname in enumerate(GROUP):
        mfull += x[gi] * gfull(gname, k_fw, p_fw)
    print("\n--- сверка по опорным пикам ---")
    print("  %8s %12s %12s %8s %8s %8s"
          % ("E, кэВ", "измерено", "модель", "мод/изм", "ПШПВ изм", "ПШПВ мод"))
    peak_rows = []
    for e0 in (238.63, 583.19, 911.20, 968.97, 1620.50, 2614.51):
        a_meas = peak_net(yfull, e0)
        a_mod = peak_net(mfull, e0)
        if a_meas is None or a_mod is None:
            continue
        rel = a_mod / a_meas if a_meas else float("nan")
        g = fw_at(e0)
        fy = measure_fwhm(yfull, e0, g) or float("nan")
        fm = measure_fwhm(mfull, e0, g) or float("nan")
        print("  %8.2f %12.0f %12.0f %8.3f %8.1f %8.1f"
              % (e0, a_meas, a_mod, rel, fy, fm))
        peak_rows.append((e0, a_meas, a_mod, rel, fy, fm))

    with io.open(os.path.join(outdir, "unfold_activities.csv"), "w",
                 encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["компонента", "активность_Бк", "сигма_Бк"])
        for gname, a, s in rows:
            w.writerow([gname, "%.5g" % a, "%.4g" % s if s else ""])
        w.writerow(["ПШПВ_662_кэВ", "%.3g" % fw662, ""])
        w.writerow(["ПШПВ2_f0", "%.4g" % k_fw, ""])
        w.writerow(["ПШПВ2_f1", "%.4g" % p_fw, ""])
        w.writerow(["хи2_на_канал", "%.4g" % chi2, ""])
        for e0, am, mo, rel, fy, fm in peak_rows:
            w.writerow(["пик_%.0f_модель_на_измерено" % e0, "%.3f" % rel, ""])

    with io.open(os.path.join(outdir, "peak_check.csv"), "w",
                 encoding="utf-8", newline="") as f:
        f.write("# сверка модели и измеренного по опорным пикам,\n")
        f.write("# нетто в окне +-1 ПШПВ с вычетом подложки по внешним подокнам;\n")
        f.write("# ПШПВ снята по полувысоте над подложкой одинаково у обоих\n")
        f.write("E_кэВ;измерено;модель;модель_на_измерено;"
                "ПШПВ_изм_кэВ;ПШПВ_мод_кэВ\n")
        for e0, am, mo, rel, fy, fm in peak_rows:
            f.write("%.2f;%.0f;%.0f;%.4f;%.4g;%.4g\n"
                    % (e0, am, mo, rel, fy, fm))

    with io.open(os.path.join(outdir, "width_fit.csv"), "w",
                 encoding="utf-8", newline="") as f:
        f.write("# подгонка ширины линии по одиночным опорным пикам\n")
        f.write("# ПШПВ(E) = k * ПШПВ_файла(E) * (E/662)^p\n")
        f.write("# k;%.5g\n" % k_fw)
        f.write("# p;%.5g\n" % p_fw)
        f.write("E_кэВ;ПШПВ_измерена_кэВ;ПШПВ_подгонка_кэВ;отклонение_%\n")
        for e0, fw, pred in width_rows:
            f.write("%.2f;%.4g;%.4g;%.3g\n"
                    % (e0, fw, pred, 100.0 * (pred - fw) / fw))

    # Вклад ОТДЕЛЬНЫХ нуклидов при найденных активностях — для рисунка.
    names = []
    parts = []
    for gi, (gname, members) in enumerate(GROUP.items()):
        for k, wk in members.items():
            if k not in tmap:
                continue
            e, c = tmap[k]
            # ШИРИНА ТА ЖЕ, что в подгонке. Прежняя редакция звала broaden без
            # fwhm_fn, и при подогнанной кривой f0 = f1 = 0 давали ядро шириной
            # 1 кэВ: на рисунке заливка компонент выходила иглами при том, что
            # сама модель была уширена правильно.
            v = wk * broaden(e, c, centres, 0.0, 0.0,
                             fwhm_fn=lambda q: fwc_kp(q, k_fw, p_fw)
                             )[m] * t_smp * x[gi]
            if v.sum() <= 0:
                continue
            names.append(dict(ORDER_MAP).get(k, k))
            parts.append(v)
    names.append("фон")
    parts.append(ybg[m])

    # Пишется через io.open с явным UTF-8: np.savetxt кодирует шапку системной
    # кодировкой, и кириллические имена компонент выходили нечитаемыми.
    tab = np.column_stack([centres[m], yy, model, *parts])
    with io.open(os.path.join(outdir, "unfold_spectrum.csv"), "w",
                 encoding="utf-8", newline="") as f:
        f.write("E_keV,измерено,модель," + ",".join(names) + "\n")
        for row in tab:
            f.write(",".join("%.6g" % v for v in row) + "\n")
    print("\nзаписано: %s" % os.path.join(outdir, "unfold_spectrum.csv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
