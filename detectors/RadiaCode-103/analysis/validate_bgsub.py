# -*- coding: utf-8 -*-
"""Проверка поправки k(E) на искажение фона пробой — по измерениям.

Поправка k(E) до сих пор нигде не сверялась с данными: считалась, но не
проверялась. А проверить её можно, потому что у оператора есть ПАРА наборов БЕЗ
свинцового домика, снятых в одном месте:

    «Фон 7 дней без домика»                  — пустой сосуд, 612 250 с
    «Черника авторская маринелли без домика» — та же проба, 186 739 с

Приложение само связало их (BackgroundSpectrumFile), то есть пара штатная.

ИДЕЯ. Собственное излучение пробы берётся из НАБОРА В ДОМИКЕ, где фон помещения
подавлен в 19 раз, — а не из расчёта. Тогда

    k(E) = [ (проба без домика) - (своё излучение) ] / (пустой сосуд без домика)

измеряется вообще без модельных чисел: три набора и один шаг арифметики. Модель
входит только как то, что проверяется.

Почему это корректно: множитель 0,837, которым иначе пришлось бы нормировать
расчётный вклад пробы, сам получен из набора в домике, так что расчётный путь
численно совпал бы с измеренным по построению — и проверка стала бы кольцевой.
Расчётный вклад всё равно печатается рядом, для контроля.

ЧЕГО ЭТА ПРОВЕРКА НЕ МОЖЕТ. Наборы разнесены по времени (7 суток против 2,2), а
радон в помещении меняется. Поэтому линии радонового ряда (352, 609, 1120, 1764)
могут разойтись сами по себе; линия Tl-208 2614 кэВ от радона не зависит и
служит контролем.

Запуск:  python validate_bgsub.py m500
"""
import sys
import os

import numpy as np

# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec
import read_rcxml

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'common', 'py'))
import paths  # noqa: E402


BASE = str(paths.measured("RadiaCode-103"))  # измеренные спектры прибора
SAMPLE = os.path.join(BASE, "Черника_авторская_маринелли_без_домика.xml")
BG = os.path.join(BASE, "Фон 7 дней без домика.xml")

# Состав пробы. Cs-137 — аттестованная привязка (приложение сверено с АТОМТЕХ),
# K-40 и Sr-90 — из независимых данных оператора по этой партии ягод.
MASS_KG = 0.246
ACT = {"Cs137": 3340.0 * MASS_KG, "K40": 240.0 * MASS_KG, "Sr90": 50.0 * MASS_KG}
K_TOT = 0.837        # нормировка полного отклика, см. normalization.py
NUCDIR = "organic_0.49"
BG_SAMPLE = "bg_cyl_organic_0.49.csv"
BG_EMPTY = "bg_cyl_air_0.00.csv"

BANDS = [(20, 100), (100, 300), (300, 600), (600, 750), (750, 1400),
         (1400, 1550), (1550, 2700)]


def edges(spec):
    e = spec.energy
    return np.concatenate(([e[0] - (e[1] - e[0]) / 2],
                           (e[:-1] + e[1:]) / 2,
                           [e[-1] + (e[-1] - e[-2]) / 2]))


def to_grid(hist1kev, spec):
    """Гистограмму 1 кэВ/канал — на каналы измеренного спектра."""
    src = np.arange(len(hist1kev) + 1, dtype=float)      # границы, кэВ
    cum = np.concatenate(([0.0], np.cumsum(hist1kev)))
    return np.diff(np.interp(edges(spec), src, cum))


def rebin_meas(spec, other):
    """Спектр other — на каналы spec (у наборов разная калибровка)."""
    cum = np.concatenate(([0.0], np.cumsum(other.counts)))
    return np.diff(np.interp(edges(spec), edges(other), cum))


def pick(path, needle):
    for s in read_rcxml.read(path):
        if needle in s.name:
            return s
    raise SystemExit("не нашёл набор «%s» в %s" % (needle, path))


def own_measured(target):
    """Собственное излучение пробы, имп/с на канал target — из набора в домике.

    В свинцовом домике фон помещения подавлен в 19 раз, поэтому разность
    «проба в домике минус холостой набор в домике» — это почти чистое излучение
    самой пробы. Остаточное искажение: свинец подсвечивает мягкую область своими
    K-линиями 72–88 кэВ, но они есть в обоих наборах и вычитаются.
    """
    from fit_peak import SAMPLE as S_SHIELD, BG as B_SHIELD
    if not (os.path.exists(S_SHIELD) and os.path.exists(B_SHIELD)):
        return None
    s = read_rcxml.read(S_SHIELD)[0]
    b = read_rcxml.read(B_SHIELD)[0]
    net = s.counts - rebin_meas(s, b) * (s.live / b.live)
    # приводим к каналам целевого набора через накопленную сумму
    cum = np.concatenate(([0.0], np.cumsum(net)))
    return np.diff(np.interp(edges(target), edges(s), cum)) / s.live


def fit662(spec, bgscaled, kbg=1.0):
    """Площадь пика 662 подгонкой «гауссиана + линейная подложка»."""
    from fit_peak import peak_model
    from scipy.optimize import curve_fit
    net = spec.counts - bgscaled
    e, dE = spec.energy, np.gradient(spec.energy)
    m = (e > 560) & (e < 790)
    x, y = e[m], net[m] / dE[m]
    # Дисперсия ЧИСТОГО счёта: N_проба + k²·N_фон_сырой, где k —
    # множитель приведения фона по времени. Через уже
    # приведённый фон это N_проба + k·bg_приведённый. Стояло без k,
    # то есть при k < 1 погрешность занижалась, и сверка «модель
    # против аттестации» показывала согласие лучше действительного.
    # Эталон вывода — detectors/Gamma-1S/analysis/export_curves.py.
    err = np.sqrt(np.maximum(spec.counts[m] + kbg * bgscaled[m], 1)) / dE[m]
    p, cov = curve_fit(peak_model, x, y,
                       p0=[max(net[m].sum(), 1.0), 661.657,
                           0.084 * 662 / 2.355, y.min(), 0.0],
                       sigma=err, absolute_sigma=True)
    return p[0], np.sqrt(cov[0, 0])


def peak_cross_check(smp_open, bg_open):
    """Один и тот же образец, снятый БЕЗ домика и В домике, обязан дать одну
    активность. Проверка ловит смещение от вычитания фона: без домика фон под
    пиком в 19 раз выше, поэтому если вычитание кривит площадь — увидим здесь."""
    from fit_peak import SAMPLE as S_SHIELD, BG as B_SHIELD
    if not (os.path.exists(S_SHIELD) and os.path.exists(B_SHIELD)):
        print("\nнет набора в домике — сверка по пику пропущена")
        return
    smp_sh = read_rcxml.read(S_SHIELD)[0]
    bg_sh = read_rcxml.read(B_SHIELD)[0]

    k_o = smp_open.live / bg_open.live
    a_o, d_o = fit662(smp_open, rebin_meas(smp_open, bg_open) * k_o, k_o)
    k_s = smp_sh.live / bg_sh.live
    a_s, d_s = fit662(smp_sh, rebin_meas(smp_sh, bg_sh) * k_s, k_s)
    r_o, dr_o = a_o / smp_open.live, d_o / smp_open.live
    r_s, dr_s = a_s / smp_sh.live, d_s / smp_sh.live

    print("\nсверка по площади пика 662: одна проба, два набора")
    print("%-16s %10s %12s %14s" % ("набор", "живое, с", "площадь", "имп/с"))
    print("%-16s %10d %6.0f±%-5.0f %.5f±%.5f"
          % ("без домика", smp_open.live, a_o, d_o, r_o, dr_o))
    print("%-16s %10d %6.0f±%-5.0f %.5f±%.5f"
          % ("в домике", smp_sh.live, a_s, d_s, r_s, dr_s))
    d = np.hypot(dr_o, dr_s)
    print("отношение без домика / в домике: %.3f ± %.3f (обязано быть 1)"
          % (r_o / r_s, d / r_s))
    print("расхождение %.1f сигма" % (abs(r_o - r_s) / d))


def main():
    v = rcspec.vessel()
    smp = pick(SAMPLE, "без домика")
    bg = pick(BG, "Фон")
    print("проба:  %s" % smp)
    print("фон:    %s" % bg)
    print("привязка фона в файле пробы: %s" % smp.bgname)

    # Собственное излучение пробы — ИЗМЕРЕННОЕ, из набора в домике.
    own_g = own_measured(smp)
    if own_g is None:
        print("\nнет набора в домике — проверку сделать нечем")
        return

    # то же по расчёту, для контроля (в полном счёте совпадёт по построению)
    own_mc = np.zeros(rcspec.NBINS)
    print("\nрасчётный вклад пробы (прогоны полного распада x %.3f), для сверки:"
          % K_TOT)
    for nuc, a in ACT.items():
        p = rcspec.rdir("nuclides", NUCDIR, "nuc_%s.csv" % nuc, v=v)
        if not os.path.exists(p):
            print("  %-6s нет прогона — пропущен" % nuc)
            continue
        meta, h = rcspec.read_spec(p)
        n = float(meta["N_primaries"])
        own_mc += h / n * a * K_TOT
        print("  %-6s %6.1f Бк -> %.4f имп/с" % (nuc, a, h[20:].sum() / n * a * K_TOT))
    own_mc_g = to_grid(rcspec.fold(own_mc), smp)
    thr0 = smp.energy > 20
    print("  итого расчёт %.3f имп/с, измерено в домике %.3f имп/с"
          % (own_mc_g[thr0].sum(), own_g[thr0].sum()))

    # измеренные скорости счёта
    s_rate = smp.counts / smp.live
    b_rate = rebin_meas(smp, bg) / bg.live

    thr = smp.energy > 20
    print("\nполные скорости счёта выше 20 кэВ:")
    print("  пустой сосуд       %.3f имп/с" % b_rate[thr].sum())
    print("  проба в сосуде     %.3f имп/с" % s_rate[thr].sum())
    print("  превышение         %.3f имп/с" % (s_rate[thr] - b_rate[thr]).sum())
    print("  своё излучение     %.3f имп/с (измерено в домике)" % own_g[thr].sum())
    print("  => на фон остаётся %+.3f имп/с, то есть k_полное = %.3f"
          % ((s_rate[thr] - b_rate[thr] - own_g[thr]).sum(),
             (s_rate[thr] - own_g[thr]).sum() / b_rate[thr].sum()))

    peak_cross_check(smp, bg)

    # модельный k(E)
    pe = rcspec.rdir("background", BG_EMPTY, v=v)
    ps = rcspec.rdir("background", BG_SAMPLE, v=v)
    if not (os.path.exists(pe) and os.path.exists(ps)):
        print("\nмодельных прогонов фона ещё нет:\n  %s\n  %s" % (pe, ps))
        return
    me, he = rcspec.read_spec(pe)
    ms, hs = rcspec.read_spec(ps)
    he = rcspec.fold(he / float(me["N_primaries"]))
    hs = rcspec.fold(hs / float(ms["N_primaries"]))
    he_g, hs_g = to_grid(he, smp), to_grid(hs, smp)

    print("\n%-14s %10s %10s %9s %9s %8s" %
          ("полоса, кэВ", "фон,имп/с", "проба-своё", "k изм.", "k модель",
           "разн."))
    for lo, hi in BANDS:
        m = (smp.energy >= lo) & (smp.energy < hi)
        b = b_rate[m].sum()
        s = (s_rate[m] - own_g[m]).sum()
        if b <= 0:
            continue
        kmeas = s / b
        kmod = hs_g[m].sum() / he_g[m].sum() if he_g[m].sum() > 0 else np.nan
        # статистика: отсчёты в обоих наборах
        ds = np.sqrt(smp.counts[m].sum()) / smp.live
        db = np.sqrt(rebin_meas(smp, bg)[m].sum()) / bg.live
        dk = kmeas * np.hypot(ds / max(s, 1e-9), db / b)
        print("%-14s %10.4f %10.4f %6.3f±%.3f %9.3f %+7.1f %%"
              % ("%d..%d" % (lo, hi), b, s, kmeas, dk, kmod,
                 100 * (kmeas / kmod - 1) if kmod == kmod else float("nan")))

    print("\nдоля собственного излучения пробы в полосе (для чтения таблицы):")
    for lo, hi in BANDS:
        m = (smp.energy >= lo) & (smp.energy < hi)
        if s_rate[m].sum() > 0:
            print("  %-12s %5.1f %%" % ("%d..%d" % (lo, hi),
                                        100 * own_g[m].sum() / s_rate[m].sum()))


if __name__ == "__main__":
    main()
