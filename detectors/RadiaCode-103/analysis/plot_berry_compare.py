# -*- coding: utf-8 -*-
"""Картинка №2: спектр пробы «черника» — измерение против расчёта.

Сравнивается НЕТТО (проба минус измеренный фон домика) с модельным откликом
пробы. Так из сравнения уходит вопрос о фоне (модель его занижает) и остаётся
ровно то, что проверяется: отклик прибора на известную активность.

Cs-137 берётся с ЭТАЛОНА — 834 Бк по приложению RadiaCode (скриншот оператора).
K-40 нормируется на активность, ИЗМЕРЕННУЮ по его же линии 1460,8 той же
моделью: это метрология (как активности бетона), а не подгонка формы.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths        # noqa: E402
import rcspec       # noqa: E402
import read_rcxml   # noqa: E402
import fit_lines as fl   # noqa: E402

APP_BQ = 834.0        # Бк Cs-137 по приложению (скриншот, +-1,15 %)
MASS_KG = 0.246
E_CS, Y_CS = 661.657, 0.851
E_K, Y_K = 1460.82, 0.1066


def load_model(nuc):
    """-> cps на 1 Бк активности пробы (свёрнуто с разрешением прибора)."""
    p = os.path.join(str(paths.build("RadiaCode-103")), "berry_%s.csv" % nuc)
    meta, h = rcspec.read_spec(p)
    return rcspec.fold(h, "103") / float(meta["N_primaries"])


def net_area(e, y, e0, half_sig=2.5, gap=3.2, side=3.0):
    # gap отсчитывается от ЦЕНТРА пика: при gap=1.6 сигма боковое окно
    # начиналось на 1,6 сигма, то есть ВНУТРИ окна пика (±2,5 сигма) — подложка
    # считалась по хвосту самого пика и завышалась, площадь занижалась на ~18 %.
    # gap обязан быть больше half_sig.
    """Нетто-площадь пика с линейной подложкой по боковым окнам."""
    sig = rcspec.fwhm(e0, "103") / 2.35482
    mp = (e >= e0 - half_sig * sig) & (e <= e0 + half_sig * sig)
    ml = (e >= e0 - (gap + side) * sig) & (e < e0 - gap * sig)
    mh = (e > e0 + gap * sig) & (e <= e0 + (gap + side) * sig)
    yl, yh = y[ml].mean(), y[mh].mean()
    xl, xh = e[ml].mean(), e[mh].mean()
    base = yl + (yh - yl) * (e[mp] - xl) / (xh - xl)
    return y[mp].sum() - base.sum(), mp


# Правило контура (оператор, 18.08): «всегда образец и фон калибруются отдельно».
# Записанные в XML коэффициенты RC-103 реальность НЕ описывают: у файлов «черника» и
# «Фон 7 дней без домика» они ОДИНАКОВЫ, а линия K-40 стоит на 15 кэВ в разных местах —
# усиление дрейфует между измерениями. Шкалы перестроены по якорям В КАНАЛАХ методом
# SpectraVibe (03_ecal_rebuild.py, SNIP + полином): проба — по Cs-137 661,657 и Pb Ka1
# 74,97 (rms штатной 2,59 кэВ), фон домика — по Pb Ka1 + K-40 1460,82 + Tl-208 2614,51
# (rms штатной 3,41 кэВ). Проверка: пик Cs-137 после поправки встаёт на 661,7 кэВ.
# P-017: шкала RC-103 НЕ описывается квадратичным полиномом — между якорями 75 и
# 1461 кэВ она провисает (Bi-214 609,3 промахивалась на -15,1 кэВ), потому что у
# CsI(Tl) отклик непропорционален и сильнее всего ниже 200 кэВ. ПОЛИНОМ 4-й СТЕПЕНИ
# (разрешён оператором 18.08). ФОРМА нелинейности снята на фоне КОМНАТЫ по шести
# якорям: Pb Ka 74,97 · Bi-214 609,32 · дублет Ac-228 911+969 (эфф. 933,14) ·
# дублет Bi-214 1120+1155 (эфф. 1123,73) · K-40 1460,82 · Tl-208 2614,51 — rms 1,94 кэВ
# при одной степени свободы (против 8,4 у квадратичной), шкала монотонна.
# Нелинейность — свойство ПРИБОРА, поэтому на другие файлы переносится ФОРМА, а
# свободны только усиление и смещение (дрейф между измерениями).
# проба: 2 якоря (Pb Ka 74,97 ch 30,27 и Cs-137 661,657 ch 262,68) на 2 свободных
# параметра -> невязка ноль ПО ПОСТРОЕНИЮ, но ФОРМА не подгонялась (взята с комнаты).
CAL_SAMPLE = [-18.3177083934, 3.1958256172, -0.0039870149, 7.4289e-06, -4e-09]
CAL_BG = [-20.5583527552, 3.2593915726, -0.0040663178, 7.5767e-06, -4.1e-09]


def ecal(coef, n):
    ch = np.arange(n)
    return sum(c * ch ** i for i, c in enumerate(coef))


def load_meas():
    b = str(paths.measured("RadiaCode-103"))
    smp = read_rcxml.read(os.path.join(
        b, "RC103 черника маринелли авторская домик 246 гр.xml"))[0]
    bg = read_rcxml.read(os.path.join(b, "Фон домик 23 дня.xml"))[0]
    # последний канал — переполнение (P-007)
    e = ecal(CAL_SAMPLE, len(smp.counts))[:-1]
    e_b = ecal(CAL_BG, len(bg.counts))[:-1]
    cps_s = (smp.counts / smp.live)[:-1]
    cps_b_own = (bg.counts / bg.live)[:-1]
    # Фон на шкалу пробы — интерполяцией ПЛОТНОСТИ (имп/с на кэВ) с обратным
    # умножением на ширину канала: иначе при растяжении шкалы теряется площадь.
    cps_b = np.interp(e, e_b, cps_b_own / np.gradient(e_b),
                      left=0.0, right=0.0) * np.gradient(e)
    return (e, cps_s, cps_b,
            (np.sqrt(np.maximum(smp.counts, 0)) / smp.live)[:-1], smp.live)


def build():
    e, cps_s, cps_b, sd, live = load_meas()
    net = cps_s - cps_b
    e_mod = np.arange(rcspec.NBINS) + 0.5
    cs = fl.rebin_model_to_meas(e_mod, load_model("Cs137"), e)
    k40 = fl.rebin_model_to_meas(e_mod, load_model("K40"), e)
    # активность K-40 ИЗМЕРЯЕМ по его же линии той же моделью (не подгонка формы)
    a_meas, _ = net_area(e, net, E_K)
    a_mod, _ = net_area(e, k40, E_K)
    bq_k = a_meas / a_mod if a_mod > 0 else 0.0
    return e, cps_s, cps_b, net, sd, cs, k40, bq_k, live


LINES = ((661.7, "Cs-137\n661,7"), (1460.8, "K-40\n1460,8"),
         (477.6, "комптон-край\nCs-137 477"), (184.0, "обратное\nрассеяние ~184"),
         (72.8, "ХРИ Pb\n72,8-88"))


def draw_lines(ax, emin, emax):
    for x, nm in LINES:
        if not (emin <= x <= emax):
            continue
        ax.axvline(x, color="#8d99ae", lw=0.6, ls=":", alpha=0.85)
        ax.annotate(nm, xy=(x, ax.get_ylim()[1]),
                    xytext=(x, ax.get_ylim()[1] * 0.30), fontsize=6.6,
                    rotation=90, color="#495057", ha="center", va="top")


def ratio_panel(axr, e, net, model, emin, emax):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(net > 0, model / net, np.nan)
    k = np.ones(9) / 9
    axr.step(e, np.convolve(np.nan_to_num(r), k, mode="same"), where="mid",
             lw=1.0, color="#c1121f")
    axr.axhline(1.0, color="#1b1b1b", lw=0.8)
    axr.set_xlim(emin, emax)
    axr.set_ylim(0, 2.0)
    axr.set_ylabel("модель / изм.")
    axr.set_xlabel("энергия, кэВ")
    axr.grid(alpha=0.22, lw=0.4)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "berry_compare.png"
    emin = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
    emax = float(sys.argv[3]) if len(sys.argv) > 3 else 1700.0
    # raw=1 — БЕЗ вычитания фона (замечание оператора 17.08): фон домика снят
    # в ДРУГОЙ посадке прибора (лёжа, без сосуда), проба — стоя в маринелли,
    # поэтому вычитаемое не тождественно фону измерения пробы. В режиме raw
    # сравниваются полные спектры, а модельный фон берётся ИЗМЕРЕННЫЙ — так
    # из сравнения уходит и вычитание, и заниженная модель фона.
    raw = len(sys.argv) > 4 and sys.argv[4] == "raw"
    e, cps_s, cps_b, net, sd, cs, k40, bq_k, live = build()
    w0 = (e >= 562.0) & (e < 762.0)
    model = ((cps_s - cps_b)[w0].sum() / cs[w0].sum()) * cs + bq_k * k40
    if raw:
        net = cps_s          # измерение — полный спектр пробы
        model = model + cps_b   # модель пробы + измеренный фон домика
    print("Cs-137: %.0f Бк (эталон приложения)" % APP_BQ)
    print("K-40:   %.0f Бк = %.0f Бк/кг (измерено по линии 1460,8 этой моделью)"
          % (bq_k, bq_k / MASS_KG))
    for nm, lo, hi in (("нетто 20-3000", 20, 3000), ("пик Cs-137", 600, 725),
                       ("пик K-40", 1390, 1530), ("20-300", 20, 300)):
        m = (e >= lo) & (e < hi)
        a, b = net[m].sum(), model[m].sum()
        print("  %-14s изм %.4e  модель %.4e  м/и %.3f"
              % (nm, a, b, b / a if a else float("nan")))

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(11.5, 8.6), dpi=160, sharex=True,
        gridspec_kw=dict(height_ratios=[3, 1], hspace=0.06))
    ax.step(e, net, where="mid", lw=1.1, color="#1b1b1b",
            label=("измерено: черника 246 г, БЕЗ вычета фона" if raw
                   else "измерено: черника 246 г, фон домика вычтен"))
    if raw:
        ax.step(e, cps_b, where="mid", lw=0.9, color="#2a9d8f", alpha=0.9,
                label="измеренный фон домика (для масштаба)")
    ax.fill_between(e, np.maximum(net - sd, 1e-9), net + sd, step="mid",
                    color="#1b1b1b", alpha=0.18, lw=0)
    ax.step(e, model, where="mid", lw=1.3, color="#c1121f",
            label=("модель: Cs-137 834 Бк + K-40 + измеренный фон" if raw
                   else "модель ПОЛНАЯ: Cs-137 834 Бк + K-40"))
    # Активность по НАШЕЙ модели — в том же окне 562-762, что у BecqMoni,
    # операция та же (проба минус фон). Это измерение активности нашей
    # эффективностью, а не подгонка формы: форма при этом не трогается.
    w = (e >= 562.0) & (e < 762.0)
    bq_own = (cps_s - cps_b)[w].sum() / cs[w].sum()
    ax.step(e, bq_own * cs, where="mid", lw=1.2, color="#0a6ebd", alpha=0.95,
            label="модель: Cs-137, %.0f Бк (наша эффективность)" % bq_own)
    ax.step(e, APP_BQ * cs, where="mid", lw=1.0, color="#0a6ebd", alpha=0.55,
            ls="--", label="модель: Cs-137, 834 Бк (калибровка приложения)")
    ax.step(e, bq_k * k40, where="mid", lw=1.0, color="#e07b00", alpha=0.9,
            label="модель: K-40, %.0f Бк/кг (по линии 1461)" % (bq_k / MASS_KG))
    ax.set_yscale("log")
    ax.set_ylabel("скорость счёта, имп/с на канал")
    ax.set_xlim(emin, emax)
    vis = (e >= emin) & (e <= emax) & (net > 0)
    ax.set_ylim(max(1e-7, net[vis].min() * 0.5), net[vis].max() * 3)
    ax.grid(alpha=0.22, lw=0.4)
    ax.legend(fontsize=8.4, loc="upper right", framealpha=0.95)
    ax.set_title("Проба «черника» в маринелли: расчёт против измерения\n"
                 "Активность Cs-137 НЕ подбиралась — взята с прибора (834 Бк)",
                 fontsize=11)
    draw_lines(ax, emin, emax)
    ratio_panel(axr, e, net, model, emin, emax)
    fig.text(0.5, 0.012,
             "Cs-137 нормирован на 834 Бк — показание прибора, НЕ подгонка формы. "
             "K-40 нормирован по площади СВОЕЙ линии 1460,8 этой же моделью.\n"
             "Отношение модель/измерение в пике Cs-137 — проверяемая величина: оно "
             "измеряет расхождение нашей эффективности с заводской калибровкой.",
             ha="center", fontsize=7.2, color="#343a40")
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(out)
    print("записано", out)


if __name__ == "__main__":
    main()
