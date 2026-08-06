# -*- coding: utf-8 -*-
"""Эффект фрезеровки: две расчётные кривые одной ревизии плюс библиотечная.

Обе расчётные кривые посчитаны ОДНОЙ сборкой и отличаются единственным
параметром `/asn16/capWindow`. Поэтому их отношение есть мера самой фрезеровки,
а не суммы всех расхождений между двумя моделями.

Библиотечная кривая наложена третьей: она показывает, какую из двух наших
постановок воспроизводит модель автора. Вывода о самой чужой модели отсюда не
делается — у автора два разных торцевых слоя, и каким из них снята его кривая,
не установлено (см. README, раздел «Слепой тест»).

    python analysis/draw_window_effect.py
"""
import io
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter

_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(_HERE, "..", "results"))
REF = os.path.normpath(os.path.join(_HERE, "..", "reference",
                                    "becqmoni-library-curve.csv"))
OUT = os.path.normpath(os.path.join(_HERE, "..", "drawings",
                                    "nano16pro_window_effect.png"))

C_WIN = "#1f5fa8"
C_NOWIN = "#b5651d"
C_LIB = "#2f7a4a"


def read(path, cols):
    rows, head = [], None
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        p = ln.split(",")
        if head is None:
            head = p
            continue
        r = dict(zip(head, p))
        rows.append([float(r[c]) for c in cols])
    return rows


def stamp_of(path):
    for ln in io.open(path, encoding="utf-8"):
        if ln.startswith("# src_sha1"):
            return ln.split("=", 1)[1].strip()
    return "?"


def main():
    win = read(os.path.join(RES, "eff_point_end10cm.csv"),
               ("E_keV", "eps_peak", "d_eps_peak"))
    nowin = read(os.path.join(RES, "eff_point_end10cm_nowin.csv"),
                 ("E_keV", "eps_peak", "d_eps_peak"))
    lib = read(REF, ("E_keV", "eff_peak", "d_eff_pct"))
    stamp = stamp_of(os.path.join(RES, "eff_point_end10cm.csv"))

    fig = plt.figure(figsize=(11.4, 8.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.1, 1.0], hspace=0.09,
                          left=0.10, right=0.975, top=0.885, bottom=0.085)

    ax = fig.add_subplot(gs[0])
    ax.errorbar([r[0] for r in win], [r[1] for r in win],
                yerr=[r[2] for r in win], fmt="o-", ms=4.2, lw=1.4,
                color=C_WIN, ecolor=C_WIN, capsize=2.2,
                label="расчёт: крышка Al 1,50 с фрезеровкой до 0,60 "
                      "(окно 0,4089 г/см²)")
    ax.errorbar([r[0] for r in nowin], [r[1] for r in nowin],
                yerr=[r[2] for r in nowin], fmt="s--", ms=4.0, lw=1.3,
                color=C_NOWIN, ecolor=C_NOWIN, capsize=2.2,
                label="расчёт: крышка СПЛОШНАЯ 1,50 (0,6518 г/см²)")
    ax.plot([r[0] for r in lib], [r[1] for r in lib], "^:", ms=4.6, lw=1.3,
            color=C_LIB, label="библиотечная BecqMoni (торец, каким снята, "
                               "не установлен)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Эффективность по ППП на 4π")
    ax.grid(True, which="both", alpha=0.25, lw=0.5)
    ax.legend(fontsize=8.6, loc="lower left", framealpha=0.92)
    ax.tick_params(labelbottom=False)

    ax2 = fig.add_subplot(gs[1], sharex=ax)
    gain = [(w[0], 100.0 * (w[1] / n[1] - 1.0)) for w, n in zip(win, nowin)]
    ax2.plot([g[0] for g in gain], [g[1] for g in gain], "o-", ms=4.0, lw=1.4,
             color="#7a2020")
    ax2.axhline(0, color="#888888", lw=0.8, ls=":")
    ax2.set_xscale("log")
    ax2.set_xlabel("Энергия, кэВ")
    ax2.set_ylabel("выигрыш от\nфрезеровки, %")
    ax2.grid(True, which="both", alpha=0.25, lw=0.5)
    # Подписи оси энергий — по УЗЛАМ сетки: мягкий край, где всё и происходит,
    # на декадных подписях неразличим.
    es = [r[0] for r in win]
    ticks = [t for t in (20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000) if min(es) * 0.98 <= t <= max(es) * 1.02]
    for a in (ax, ax2):
        a.xaxis.set_major_locator(FixedLocator(ticks))
        a.xaxis.set_major_formatter(ScalarFormatter())
        a.xaxis.set_minor_formatter(NullFormatter())
    for lab in ax2.get_xticklabels():
        lab.set_rotation(45)
        lab.set_ha("right")
        lab.set_fontsize(8)
    top = max(g[1] for g in gain)
    ax2.annotate("%+.0f %% на %.0f кэВ" % (gain[0][1], gain[0][0]),
                 xy=(gain[0][0], gain[0][1]),
                 xytext=(gain[0][0] * 1.5, gain[0][1] * 0.82),
                 fontsize=8.5, color="#7a2020",
                 arrowprops=dict(arrowstyle="->", color="#7a2020", lw=0.8))

    fig.suptitle("AtomSpectra Nano 16 PRO: что даёт фрезеровка крышки\n"
                 "точечный источник на оси, 10 см от наружной плоскости "
                 "корпуса, торец в пучке", fontsize=12, y=0.965)
    fig.text(0.10, 0.905,
             "Обе расчётные кривые — ОДНА сборка, один штамп %s; различие "
             "только в /asn16/capWindow, поэтому их отношение есть мера самой "
             "фрезеровки." % stamp,
             fontsize=8.4, color="#555555", ha="left")
    fig.text(0.5, 0.012,
             "Библиотечная кривая ложится на расчёт БЕЗ фрезеровки (+3,8 % на "
             "40 кэВ), а не на расчёт с ней (+17,6 %). Вывода о чужой модели "
             "отсюда НЕ следует:\n"
             "у автора два разных торцевых слоя, 0,5398 и 0,7718 г/см², и "
             "наблюдаемый знак согласуется только со вторым — каким торцом "
             "снята его кривая, не установлено.\n"
             "Жёсткий край расходится одинаково в обоих случаях, то есть к "
             "входному окну отношения не имеет.",
             fontsize=8.4, ha="center", color="#555555")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=160)
    print("записано: %s" % OUT)
    print("выигрыш от фрезеровки: %+.1f %% (%.0f кэВ) ... %+.1f %% (%.0f кэВ)"
          % (gain[0][1], gain[0][0], gain[-1][1], gain[-1][0]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
