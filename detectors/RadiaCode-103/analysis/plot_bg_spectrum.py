# -*- coding: utf-8 -*-
"""Модельный фоновый спектр в том же виде, в каком его показывает прибор:
свёрнутый с разрешением, в импульсах в секунду на канал, логарифмическая шкала.
Нужен для сравнения формы с реальным набором."""
import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
BG = rcspec.rdir("background")
FIGS = rcspec.rdir("figures")

CYL_S = 2 * math.pi * 4.5 * (4.5 + 16.5)

LINES = [(238.6, "Pb-212"), (351.9, "Pb-214"), (583.2, "Tl-208"),
         (609.3, "Bi-214"), (911.2, "Ac-228"), (1120.3, "Bi-214"),
         (1460.8, "K-40"), (1764.5, "Bi-214"), (2614.5, "Tl-208")]


def field_fluence():
    for line in open(os.path.join(RESULTS, "field_spectrum.mac"), encoding="utf-8"):
        m = re.match(r"#\s*FLUENCE_TOTAL_CM2_S\s*=\s*([\d.eE+-]+)", line)
        if m:
            return float(m.group(1))
    raise SystemExit("нет нормировки флюенса")


def main():
    # Приводим модель к дозе конкретного помещения: измерено 0.09 мкЗв/ч,
    # расчётное поле соответствует 0.134 мкЗв/ч.
    scale = 0.09 / 0.134
    phi = field_fluence()
    rate = phi * CYL_S / 4.0

    meta, hist = rcspec.read_spec(os.path.join(BG, "bg_cyl_air_0.00.csv"))
    n = float(meta["N_primaries"])
    cps = hist / n * rate * scale
    folded = rcspec.fold(cps, "103")

    tot = folded[20:].sum()
    print("модель, приведённая к 0.09 мкЗв/ч: %.2f имп/с выше 20 кэВ" % tot)
    print("измерено прибором: 6.36 имп/с")
    print("отношение модель/измерение: %.2f" % (tot / 6.36))

    fig, ax = plt.subplots(figsize=(9, 6))
    E = np.arange(len(folded)) + 0.5
    m = (E > 15) & (E < 2900)
    ax.fill_between(E[m], folded[m], color="#c77b30", alpha=0.35)
    ax.plot(E[m], folded[m], color="#a85c10", lw=1.1)
    for e0, nuc in LINES:
        ax.axvline(e0, color="0.7", lw=0.6, ls=":", zorder=0)
        ax.text(e0, folded[int(e0)] * 2.2, nuc, rotation=90, fontsize=7.5,
                ha="center", va="bottom", color="0.35")
    ax.set_yscale("log")
    ax.set_xlim(0, 2900)
    ax.set_xlabel("энергия, кэВ")
    ax.set_ylabel("имп/с на канал 1 кэВ")
    ax.set_title("Модельный фон в пустом сосуде, приведён к 0,09 мкЗв/ч\n"
                 "свёрнут с разрешением прибора (8,4 %% на 662 кэВ)", fontsize=11)
    ax.grid(True, which="both", alpha=0.25, lw=0.6)
    fig.tight_layout()
    out = os.path.join(FIGS, "background_spectrum.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("график:", out)

    # доли счёта по диапазонам — для сравнения с формой измеренного
    print("\nраспределение счёта по диапазонам (модель):")
    for lo, hi in [(20, 100), (100, 300), (300, 700), (700, 1200),
                   (1200, 1600), (1600, 2900)]:
        print("  %4d..%4d кэВ: %5.1f %%" % (lo, hi, 100 * folded[lo:hi].sum() / tot))


if __name__ == "__main__":
    main()
