"""Трек 1c: мост методик — размытие МК-спектра и съём площадей фитом.

Две наблюдаемые, которые весь проект сравнивал напрямую, — разные:
наша eps — сумма депозитов в окне ±1 ПШПВ (после уширения) либо ±6 кэВ
(сетки), аттестованная — площадь гауссова фита со ступенькой на реальном
спектре. Здесь обе стороны приводятся к ОДНОЙ процедуре:

1. Модельный спектр распада (chain_Th232.csv, бины 1 кэВ) размывается
   гауссианой с приборной шириной. Ширина — НЕ модельный закон
   49,9*sqrt(E) (он врёт до +53%/-7%, задача 114), а in-situ полином по
   9 пикам самого спектра Th232_420-7-17_Маринелли_0cm (fwhm_check.py).
2. Измеренный спектр читается штатным читателем SpectraVibe.
3. Обе стороны фитуются ОДНОЙ функцией: сумма гауссиан + ступенька
   (высота на амплитуду — из пика-образа .cpt прибора) + линейный фон,
   зонами вокруг групп линий; ширина и положения плавают, как у
   СпектраЛайн (Minimize=FWHM,Position,Linear).
4. Сравниваются: (а) фит измерения против площадей сеанса СпектраЛайн
   29.07.2026 — валидация нашего фита; (б) фит модели против оконного
   съёма модели — цена конвенции площади для eps; (в) отношение
   изм/модель в ЕДИНОЙ конвенции — очищенный от моста уровень.

Высоты ступеньки из .cpt (доля амплитуды): 239 кэВ 0,00183; 662 0,00035;
2613 0,00139 — вшиты константами (файл прибора, в репо не тянем).
"""
import math
import os
import sys

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fwhm_check import MAR_POLY, poly_sqrtE  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
SPE = os.path.join("detectors", "Gamma-1S", "raw_lsrm", "Work", "BG",
                   "Gamma-1S", "Spe - поверки", "Поверка 2024", "Маринелли",
                   "Th232_420-7-17_Маринелли_0cm.spe")
TLIVE_HINT = 11359.164

# Зоны фита: (границы кэВ, [линии], h_step из .cpt ближайшего образа)
ZONES = [
    ((520.0, 650.0), [583.187, 609.312], 0.0006),
    ((840.0, 1050.0), [911.204, 964.766, 968.971], 0.0009),
    ((2380.0, 2900.0), [2614.511], 0.0014),
]
# 609,3 — Bi-214 фона в измерении; в модельном спектре чистого Th-232 его
# нет, но лишняя свободная линия там просто занулится.


def fwhm_of(E):
    return poly_sqrtE(MAR_POLY, E)


def blur(edges_counts, emax=3200):
    """Размытие спектра 1-кэВ бинов in-situ гауссианой."""
    src = np.zeros(emax)
    for e, c in edges_counts:
        i = int(e)
        if 0 <= i < emax:
            src[i] += c
    out = np.zeros(emax)
    xs = np.arange(emax, dtype=float)
    for i in np.nonzero(src)[0]:
        s = fwhm_of(max(xs[i], 20.0)) / 2.3548
        lo, hi = max(0, int(xs[i] - 6 * s)), min(emax, int(xs[i] + 6 * s))
        g = np.exp(-0.5 * ((xs[lo:hi] - xs[i]) / s) ** 2)
        out[lo:hi] += src[i] * g / (s * math.sqrt(2 * math.pi))
    return xs, out


def load_mc(path):
    ec = []
    N = None
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"):
            if "N_primaries" in ln:
                N = int(ln.split("=")[1])
            continue
        if ln and ln[0].isdigit():
            e, c = ln.split(",")
            ec.append((float(e), float(c)))
    return ec, N


def zone_model(lines, hstep):
    """Модель зоны: гауссианы (площадь, сдвиг общий, масштаб ПШПВ общий)
    + ступеньки (высота = hstep*амплитуда) + линейный фон."""
    def f(x, *p):
        # p = [A1..An, shift, wscale, b0, b1]
        n = len(lines)
        shift, wscale, b0, b1 = p[n], p[n + 1], p[n + 2], p[n + 3]
        y = b0 + b1 * (x - x.mean())
        for A, E0 in zip(p[:n], lines):
            s = wscale * fwhm_of(E0) / 2.3548
            mu = E0 + shift
            amp = A / (s * math.sqrt(2 * math.pi))
            y = y + amp * np.exp(-0.5 * ((x - mu) / s) ** 2)
            y = y + amp * hstep * 0.5 * erfc((x - mu) / (s * math.sqrt(2)))
        return y
    return f


def fit_zone(xs, ys, zone, lines, hstep, sig=None):
    """Фит зоны. xs/ys — ПЛОТНОСТЬ (отсчёты на кэВ) на любой сетке:
    для измерения это канальная сетка (пересыпать 3-кэВ каналы в 1-кэВ
    бины нельзя — две трети бинов пустые, и нули душат фит), для
    модельного спектра — его родные 1-кэВ бины."""
    lo, hi = zone
    m = (xs >= lo) & (xs <= hi)
    x, y = xs[m], ys[m]
    s = sig[m] if sig is not None else np.sqrt(np.maximum(y, 1.0))
    p0 = [max(y.sum() * (x[1] - x[0]) / max(len(lines), 1), 1.0)] \
        * len(lines) + [0.0, 1.0, float(np.median(y)), 0.0]
    bounds = ([0.0] * len(lines) + [-30.0, 0.5, -np.inf, -np.inf],
              [np.inf] * len(lines) + [30.0, 2.0, np.inf, np.inf])
    p, _ = curve_fit(zone_model(lines, hstep), x, y, p0=p0, sigma=s,
                     bounds=bounds, maxfev=20000)
    return p


def spectrum_density(counts, cal):
    """Канальный спектр -> (энергии каналов, плотность имп/кэВ, сигма)."""
    ch = np.arange(len(counts), dtype=float)
    en = sum(c * ch ** k for k, c in enumerate(cal))
    w = np.gradient(en)                       # кэВ на канал
    y = np.asarray(counts, dtype=float)
    dens = y / w
    sig = np.sqrt(np.maximum(y, 1.0)) / w
    return en, dens, sig


if __name__ == "__main__":
    root = paths.require_spectravibe("мост методик: фит измеренного спектра")
    sys.path.insert(0, os.path.join(str(root), "scripts"))
    from gamma.io.lsrm_spe import read_lsrm_spe

    sp = read_lsrm_spe(os.path.join(str(root), SPE))
    xs_m, y_meas, sig_m = spectrum_density(sp.counts, list(sp.energy_cal))
    tlive = float(sp.live_time)

    ec, N = load_mc(os.path.join(BUILD, "chain_Th232.csv"))
    xs_c, y_mc = blur(ec)

    # Оконный съём модели (старая конвенция ±1 ПШПВ с полками) — для цены
    # конвенции берём то же, что kit_recalc: area_sim по уширенному спектру.
    import kit_recalc as kr
    hist = {float(e): c for e, c in ec}

    print("Мост методик: фит гаусс+ступенька+фон, ширина in-situ, зоны как"
          " у прибора.\n")
    print("%9s %12s %12s %10s | %12s %12s %8s" %
          ("линия", "изм.фит", "SL-сеанс", "фит/SL",
           "МК.фит", "МК.окно", "фит/окно"))
    rows = []
    SL = {583.187: 233661.0, 911.204: 107813.0, 964.766: 77700.0,
          968.971: 4080.0, 2614.511: 54025.0, 609.312: None}
    for zone, lines, hstep in ZONES:
        pm = fit_zone(xs_m, y_meas, zone, lines, hstep, sig=sig_m)
        pc = fit_zone(xs_c, y_mc, zone, lines, hstep)
        for i, E in enumerate(lines):
            a_meas = pm[i]
            a_mc = pc[i]
            fw = kr.FWHM662 * math.sqrt(E / 661.657)
            a_win = kr.area_sim(hist, E, fwhm=fw, key="bridge")
            sl = SL.get(E)
            print("%9.1f %12.0f %12s %10s | %12.1f %12.1f %8.3f"
                  % (E, a_meas,
                     "%.0f" % sl if sl else "—",
                     "%.3f" % (a_meas / sl) if sl else "—",
                     a_mc, a_win,
                     a_mc / a_win if a_win > 0 else float("nan")))
            rows.append((E, a_meas, sl, a_mc, a_win))

    # Очищенный от моста уровень на 2614,5: изм/модель В ОДНОЙ конвенции
    E = 2614.511
    r = [r for r in rows if r[0] == E][0]
    eps_fit = r[3] / N                       # eps на распад, конвенция фита
    rate_fit = r[1] / tlive
    A_fit = rate_fit / eps_fit
    print("\n2614,5 в единой конвенции фита: A = %.0f Бк, A/паспорт = %.3f"
          % (A_fit, A_fit / 3104.0))
    print("(было в оконной конвенции: 0,796; вклад моста = %.3f)"
          % (A_fit / 3104.0 / 0.796))

    out = os.path.join(RESULTS := os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results")),
        "blur_fit.csv")
    csvio.write(
        out,
        ["E_keV", "area_meas_fit", "area_SL_session", "area_mc_fit",
         "area_mc_window"],
        [("%.3f" % e, "%.1f" % am, "" if sl is None else "%.1f" % sl,
          "%.1f" % ac, "%.1f" % aw) for e, am, sl, ac, aw in rows],
        comments=[
            "Мост методик (Трек 1c): измеренный и размытый модельный"
            " спектры фитуются одной функцией (гаусс+ступенька из"
            " пика-образа+линейный фон), ширина in-situ.",
            "area_meas_fit против area_SL_session — валидация фита;"
            " area_mc_fit против area_mc_window — цена конвенции площади.",
        ])
    print("таблица: %s" % out)
