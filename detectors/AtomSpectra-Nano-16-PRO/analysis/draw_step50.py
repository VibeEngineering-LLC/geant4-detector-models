# -*- coding: utf-8 -*-
"""Кривая на сетке 30-3000 кэВ с шагом 50, по 2e5 историй на узел.

Сетка и число историй заданы не нами: 2e5 — столько историй на точку стоит в
описании библиотечной кривой BecqMoni, шаг 50 кэВ — задание оператора. Смысл
рисунка не в самой кривой (при таком числе историй статистика грубая), а в том,
чтобы показать ход и цену истории при аналоговом розыгрыше.

Ось энергий ЛИНЕЙНАЯ и размечена по узлам: мелкие штрихи через 50 кэВ — это
ровно шаг сетки, подписи через 250. На логарифмической оси шаг сетки не виден,
а он здесь и есть предмет.

    python analysis/draw_step50.py
"""
import io
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(_HERE, "..", "results"))
REF = os.path.normpath(os.path.join(_HERE, "..", "reference",
                                    "becqmoni-library-curve.csv"))
SRC = os.path.join(RES, "eff_step50_200k.csv")
OUT = os.path.normpath(os.path.join(_HERE, "..", "drawings",
                                    "nano16pro_step50_200k.png"))

C_PEAK = "#1f5fa8"
C_TOT = "#7a2020"
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


def head_val(path, key):
    for ln in io.open(path, encoding="utf-8"):
        if ln.startswith("# " + key) or ln.startswith("#@ " + key):
            return ln.split("=", 1)[1].strip()
    return "?"


def main():
    if not os.path.exists(SRC):
        raise SystemExit("нет %s — сначала прогон и export_curve.py" % SRC)
    rows = read(SRC, ("E_keV", "eps_peak", "d_eps_peak", "eps_total",
                      "fep_counts", "N_primaries", "solid_angle_frac"))
    lib = read(REF, ("E_keV", "eff_peak", "d_eff_pct"))
    stamp = head_val(SRC, "src.spectra_sha1")

    es = [r[0] for r in rows]
    peak = [r[1] for r in rows]
    dpeak = [r[2] for r in rows]
    tot = [r[3] for r in rows]
    # Погрешность полной эффективности: счёт событий с сигналом восстанавливается
    # из самой величины, отдельной колонки для него нет.
    dtot = []
    for r in rows:
        n_sig = r[3] * r[5] / r[6]
        dtot.append(r[3] / math.sqrt(n_sig) if n_sig > 0 else 0.0)

    fig, ax = plt.subplots(figsize=(12.6, 7.4))
    ax.errorbar(es, tot, yerr=dtot, fmt="s-", ms=3.4, lw=1.2, color=C_TOT,
                ecolor=C_TOT, capsize=1.8, elinewidth=0.8,
                label="эффективность регистрации (любой депозит)")
    ax.errorbar(es, peak, yerr=dpeak, fmt="o-", ms=3.4, lw=1.2, color=C_PEAK,
                ecolor=C_PEAK, capsize=1.8, elinewidth=0.8,
                label="эффективность по ППП, строгое окно ±1,5 кэВ")
    ax.plot([r[0] for r in lib], [r[1] for r in lib], "^:", ms=4.4, lw=1.1,
            color=C_LIB, label="библиотечная BecqMoni по ППП (2·10⁵ историй/точку)")

    ax.set_yscale("log")
    ax.set_xlim(0, 3050)
    ax.set_xlabel("Энергия, кэВ")
    ax.set_ylabel("Эффективность, отнесённая к 4π")
    ax.xaxis.set_major_locator(MultipleLocator(250))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.grid(True, which="major", axis="both", alpha=0.30, lw=0.6)
    ax.grid(True, which="minor", axis="x", alpha=0.13, lw=0.4)
    ax.grid(True, which="minor", axis="y", alpha=0.13, lw=0.4)
    ax.legend(fontsize=9.2, loc="lower left", framealpha=0.93)

    fig.suptitle("AtomSpectra Nano 16 PRO: сетка 30–3000 кэВ с шагом 50, "
                 "2·10⁵ историй на узел\n"
                 "точечный источник на оси, 10 см от наружной плоскости "
                 "передней крышки, крышка с фрезеровкой", fontsize=12, y=0.975)
    fig.text(0.5, 0.012,
             "Мелкие штрихи по оси энергий — узлы сетки (шаг 50 кэВ). Усы — "
             "статистика прогона; при 2·10⁵ историй она грубая, и кривая "
             "публикуемым результатом не является.\n"
             "Штамп исходников %s. Библиотечная кривая наложена при том же "
             "числе историй на точку: сравниваются оценщики, а не приборы."
             % stamp,
             fontsize=8.6, ha="center", color="#555555")
    fig.subplots_adjust(left=0.075, right=0.985, top=0.885, bottom=0.125)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=160)
    print("записано: %s" % OUT)
    print("узлов %d, статистика по ППП: %.1f %% на %.0f кэВ ... %.1f %% на %.0f кэВ"
          % (len(es), 100.0 * dpeak[0] / peak[0], es[0],
             100.0 * dpeak[-1] / peak[-1], es[-1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
