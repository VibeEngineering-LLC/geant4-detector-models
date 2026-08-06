# -*- coding: utf-8 -*-
"""Разложение измеренного спектра по нуклидам — рисунок с накоплением.

Читает `unfold_spectrum.csv`, который пишет `wt20_unfold.py`: измеренный спектр,
модель и вклад каждой компоненты в тех же каналах. Компоненты рисуются с
накоплением (площадями), измеренное — линией поверх.

    python analysis/draw_wt20_unfold.py <каталог с unfold_spectrum.csv>
"""
import io
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(_HERE, "..", "drawings",
                                    "nano16pro_wt20_unfold.png"))

COLOUR = {
    "Th-232": "#7a5c3a", "Ra-228": "#a08a5c", "Ac-228": "#d81b8c",
    "Th-228": "#8a6d3b", "Ra-224": "#c98b1e", "Rn-220": "#6b8f3a",
    "Po-216": "#9bb06a", "Pb-212": "#b07d2a", "Bi-212": "#2f6b34",
    "Tl-208": "#c8cf7a", "Po-212": "#8fa0a8", "фон": "#8899a6",
}


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    d = sys.argv[1]
    path = os.path.join(d, "unfold_spectrum.csv")
    with io.open(path, encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    e, meas, model = data[:, 0], data[:, 1], data[:, 2]
    comps = [(head[i], data[:, i]) for i in range(3, data.shape[1])]
    comps = [(n, v) for n, v in comps if v.sum() > 0]
    comps.sort(key=lambda t: -t[1].sum())
    total = sum(v.sum() for _, v in comps)

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(15.0, 9.4), sharex=True,
                                  gridspec_kw=dict(height_ratios=[3.4, 1.0],
                                                   hspace=0.06))
    base = np.zeros_like(e)
    for name, v in comps:
        ax.fill_between(e, np.maximum(base, 1e-9), np.maximum(base + v, 1e-9),
                        step="mid", facecolor=COLOUR.get(name, "#999999"),
                        edgecolor="none", alpha=0.92,
                        label="%s — %.2f %%" % (name, 100.0 * v.sum() / total))
        base = base + v
    ax.step(e, np.maximum(meas, 1e-9), where="mid", color="#111111", lw=0.9,
            label="измерено")
    ax.set_yscale("log")
    ax.set_ylim(max(1.0, meas[meas > 0].min() * 0.5), meas.max() * 3.0)
    ax.set_ylabel("отсчётов в канале %.0f кэВ" % (e[1] - e[0]))
    ax.grid(True, which="major", alpha=0.25, lw=0.6)
    ax.grid(True, which="minor", axis="x", alpha=0.12, lw=0.4)
    ax.xaxis.set_major_locator(MultipleLocator(250))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.legend(fontsize=8.4, loc="upper right", ncol=2, framealpha=0.95)

    # --- остаток ------------------------------------------------------------
    sig = np.sqrt(np.maximum(meas, 1.0))
    axr.axhline(0, color="#555555", lw=0.8)
    axr.step(e, (meas - model) / sig, where="mid", color="#7a2020", lw=0.8)
    axr.set_ylim(-8, 8)
    axr.set_ylabel("(изм. − модель)/σ")
    axr.set_xlabel("Энергия, кэВ")
    axr.grid(True, alpha=0.25, lw=0.6)
    axr.xaxis.set_major_locator(MultipleLocator(250))
    axr.xaxis.set_minor_locator(MultipleLocator(50))

    fig.suptitle("AtomSpectra Nano 16 PRO на пачке WT-20: разложение спектра "
                 "по нуклидам ряда тория\nшаблоны посчитаны Монте-Карло в той "
                 "же геометрии, нормировка — на один распад",
                 fontsize=12, y=0.985)
    fig.text(0.5, 0.005,
             "Побочные пики отдельными компонентами НЕ вводятся: вылет "
             "аннигиляционных квантов, обратное рассеяние, рентген вольфрама и "
             "суммирование каскада возникают в переносе и лежат внутри шаблона "
             "своего нуклида.",
             fontsize=8.4, ha="center", color="#555555")
    fig.subplots_adjust(left=0.062, right=0.988, top=0.915, bottom=0.075)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("записано: %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
