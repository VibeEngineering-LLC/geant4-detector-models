# -*- coding: utf-8 -*-
"""Карта функции полного отклика и три её сечения.

Матрица строится `analysis/response_matrix.py` из сетки моноэнергетических
прогонов. Здесь она только показывается: карта — весь отклик разом, сечения —
то, что видит спектрометрист как «отклик на линию».

Верхняя панель: строка — энергия падающего кванта, столбец — канал
энерговыделения, цвет — вероятность на один испущенный квант в 4π
(логарифмическая шкала). Диагональ — пик полного поглощения; ниже неё
комптоновский континуум; параллельные диагонали ниже пика на 511 и 1022 кэВ —
вылет одного и обоих аннигиляционных квантов.

Нижняя панель: три горизонтальных сечения карты.

    python analysis/draw_response.py [<матрица.csv>]
"""
import io
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator

_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(_HERE, "..", "results"))
OUT = os.path.normpath(os.path.join(_HERE, "..", "drawings",
                                    "nano16pro_response.png"))
SLICES = (180.0, 1480.0, 3000.0)


def read_matrix(path):
    head, cols, rows, es = {}, None, [], []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("#@ ").split("=", 1)
                head[k.strip()] = v.strip()
            continue
        p = ln.split(",")
        if cols is None:
            cols = [float(x) for x in p[1:]]
            continue
        es.append(float(p[0]))
        rows.append([float(x) for x in p[1:]])
    return head, cols, es, rows


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        RES, "response_matrix_10keV.csv")
    head, cols, es, rows = read_matrix(src)
    stamp = head.get("src.spectra_sha1", "?")

    # Нули логарифмическая шкала не рисует; порог берётся от минимума
    # ненулевых значений, а не назначается на глаз.
    nz = [v for r in rows for v in r if v > 0]
    vmin, vmax = min(nz), max(nz)

    fig = plt.figure(figsize=(12.8, 10.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], hspace=0.22,
                          left=0.085, right=0.985, top=0.90, bottom=0.085)

    ax = fig.add_subplot(gs[0])
    im = ax.pcolormesh(cols, es, rows, norm=LogNorm(vmin=vmin, vmax=vmax),
                       cmap="magma", shading="auto")
    ax.plot([0, 3050], [0, 3050], color="#8fd0ff", lw=0.7, ls="--", alpha=0.8)
    ax.plot([0, 3050 - 511], [511, 3050], color="#7fe0a0", lw=0.7, ls=":",
            alpha=0.8)
    ax.plot([0, 3050 - 1022], [1022, 3050], color="#7fe0a0", lw=0.7, ls=":",
            alpha=0.6)
    ax.set_xlim(0, 3200)
    ax.set_ylim(min(es), max(es))
    ax.set_xlabel("Энерговыделение, кэВ")
    ax.set_ylabel("Энергия падающего кванта, кэВ")
    ax.xaxis.set_major_locator(MultipleLocator(250))
    ax.yaxis.set_major_locator(MultipleLocator(250))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    cb = fig.colorbar(im, ax=ax, pad=0.012)
    cb.set_label("вероятность на квант в 4π, на канал 10 кэВ", fontsize=9)
    ax.text(0.985, 0.06, "штриховая — пик полного поглощения\n"
                         "точечные — вылет 511 и 1022 кэВ",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.4,
            color="#e8e8e8",
            bbox=dict(facecolor="#00000088", edgecolor="none", pad=3.5))

    ax2 = fig.add_subplot(gs[1])
    for e_want, col in zip(SLICES, ("#1f5fa8", "#b5651d", "#7a2020")):
        i = min(range(len(es)), key=lambda k: abs(es[k] - e_want))
        ax2.step(cols, rows[i], where="mid", lw=1.2, color=col,
                 label="падающий квант %.0f кэВ" % es[i])
    ax2.set_yscale("log")
    ax2.set_xlim(0, 3200)
    ax2.set_ylim(max(vmin, 1e-9), vmax * 1.6)
    ax2.set_xlabel("Энерговыделение, кэВ")
    ax2.set_ylabel("вероятность на квант в 4π")
    ax2.xaxis.set_major_locator(MultipleLocator(250))
    ax2.xaxis.set_minor_locator(MultipleLocator(50))
    ax2.grid(True, which="major", alpha=0.28, lw=0.6)
    ax2.grid(True, which="minor", axis="x", alpha=0.12, lw=0.4)
    ax2.legend(fontsize=9, loc="upper right", framealpha=0.93)

    fig.suptitle("AtomSpectra Nano 16 PRO: функция полного отклика, "
                 "61 узел 30–3000 кэВ с шагом 50\n"
                 "точечный источник на оси, 10 см от наружной плоскости "
                 "передней крышки; свёрнуто с ПШПВ(E) = 41,60·√(E/661,657) кэВ",
                 fontsize=12, y=0.972)
    fig.text(0.5, 0.012,
             "Матрица собрана из спектров энерговыделения тех же прогонов, "
             "отдельного расчёта отклика не делалось. Штамп исходников %s."
             % stamp, fontsize=8.6, ha="center", color="#555555")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("записано: %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
