# -*- coding: utf-8 -*-
"""Передняя крышка Nano 16 PRO: алюминий 1,50 мм с фрезеровкой до 0,60 мм.

Заведён отдельным рисунком по той же причине, что и `draw_stack.py`: входное
окно — место, где кривая на мягком краю решается целиком, и оно уже дважды
меняло определение за один день. На общем разрезе прибора выборка глубиной
0,90 мм неразличима.

Размеры читаются из `geometry/ASN16Detector.hh`; своих констант файл не держит.

    python analysis/draw_cap.py
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

_HERE = os.path.dirname(os.path.abspath(__file__))
HH = os.path.normpath(os.path.join(_HERE, "..", "geometry", "ASN16Detector.hh"))
OUT = os.path.normpath(os.path.join(_HERE, "..", "drawings",
                                    "nano16pro_cap_window.png"))

RHO = {"Al": 2.699, "ПТФЭ": 2.200}
C_AL = "#9aa5ad"
C_AL_THIN = "#c6ced4"
C_CSI = "#d8ad3c"
C_PTFE = "#f4f4f0"
C_FOIL = "#c9cdd1"
C_AIR = "#eef3f7"
EC = "#3a3f44"


def geom_from_header(path):
    src = open(path, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^\s*double\s+(\w+)\s*=\s*([0-9.]+)\s*;", src,
                         re.MULTILINE):
        out[m.group(1)] = float(m.group(2))
    need = ("cryX", "cryY", "cryZ", "ptfe", "alFoil", "wCap", "wCapWin",
            "capWinPad", "bodyX", "bodyY", "wSide", "wBot", "wFront")
    miss = [k for k in need if k not in out]
    if miss:
        raise SystemExit("в %s не найдены поля: %s" % (path, ", ".join(miss)))
    return out


def ru(x, nd=2):
    return ("%.*f" % (nd, x)).replace(".", ",")


def main():
    g = geom_from_header(HH)
    winX = g["cryX"] + 2 * g["capWinPad"]
    winY = g["cryY"] + 2 * g["capWinPad"]
    xCav = 0.5 * g["bodyX"] - g["wSide"]
    yFoilT = 0.5 * g["cryY"] + g["ptfe"] + g["alFoil"]
    yCavB = (yFoilT + g["wFront"]) - g["bodyY"] + g["wBot"]
    depth = g["wCap"] - g["wCapWin"]

    sd_win = (g["ptfe"] * RHO["ПТФЭ"] + g["alFoil"] * RHO["Al"]
              + g["wCapWin"] * RHO["Al"]) / 10.0
    sd_rim = (g["ptfe"] * RHO["ПТФЭ"] + g["alFoil"] * RHO["Al"]
              + g["wCap"] * RHO["Al"]) / 10.0

    fig = plt.figure(figsize=(13.0, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.13,
                          left=0.05, right=0.975, top=0.84, bottom=0.10)

    # ---------------- вид со стороны источника ------------------------------
    ax = fig.add_subplot(gs[0])
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Вид со стороны источника (плоскость X–Y)", fontsize=10)
    ax.add_patch(Rectangle((-xCav, yCavB), 2 * xCav, yFoilT - yCavB,
                           facecolor=C_AL, edgecolor=EC, lw=1.0))
    ax.add_patch(Rectangle((-0.5 * winX, -0.5 * winY), winX, winY,
                           facecolor=C_AL_THIN, edgecolor="#7a2020", lw=1.3,
                           ls="--"))
    ax.add_patch(Rectangle((-0.5 * g["cryX"], -0.5 * g["cryY"]),
                           g["cryX"], g["cryY"], facecolor="none",
                           edgecolor="#8a6d3b", lw=1.0, ls=":"))
    ax.text(0, yFoilT + 1.4, "окно фрезеровки %s × %s мм (штриховая рамка)"
            % (ru(winX), ru(winY)), fontsize=8.5, ha="center",
            color="#7a2020")
    ax.text(0, 0, "кристалл\n%s × %s" % (ru(g["cryX"]), ru(g["cryY"])),
            fontsize=8, ha="center", va="center", color="#4a3405")
    ax.annotate("", xy=(0.5 * g["cryX"], -0.5 * winY - 1.2),
                xytext=(0.5 * winX, -0.5 * winY - 1.2),
                arrowprops=dict(arrowstyle="<->", color="#7a2020", lw=0.9))
    ax.text(0.5 * (0.5 * g["cryX"] + 0.5 * winX), -0.5 * winY - 3.2,
            "припуск %s" % ru(g["capWinPad"]), fontsize=7.5, ha="center",
            color="#7a2020")
    ax.text(0, yCavB - 7.4, "пластина крышки по сечению полости: %s × %s мм"
            % (ru(2 * xCav), ru(yFoilT - yCavB)), fontsize=8, ha="center",
            color="#3a3f44")
    ax.set_xlim(-xCav - 4, xCav + 4)
    ax.set_ylim(yCavB - 10.5, yFoilT + 5)

    # ---------------- разрез по оси -----------------------------------------
    ax2 = fig.add_subplot(gs[1])
    ax2.set_aspect("equal")
    ax2.axis("off")
    ax2.set_title("Разрез по оси кристалла (плоскость Y–Z): что стоит в пучке",
                  fontsize=10)
    z0 = 0.0                       # внутренняя грань крышки
    zOut = z0 + g["wCap"]          # наружная плоскость корпуса
    yT = 0.5 * winY + 3.0
    # тело крышки
    ax2.add_patch(Rectangle((z0, -yT), g["wCap"], 2 * yT,
                            facecolor=C_AL, edgecolor=EC, lw=1.0))
    # выборка: воздух от внутренней грани на глубину depth
    ax2.add_patch(Rectangle((z0, -0.5 * winY), depth, winY,
                            facecolor=C_AIR, edgecolor="#7a2020", lw=1.1))
    # обёртка и кристалл слева от крышки
    ax2.add_patch(Rectangle((z0 - g["alFoil"], -0.5 * winY - 1.0),
                            g["alFoil"], winY + 2.0, facecolor=C_FOIL,
                            edgecolor=EC, lw=0.6))
    ax2.add_patch(Rectangle((z0 - g["alFoil"] - g["ptfe"],
                             -0.5 * winY - 1.0), g["ptfe"], winY + 2.0,
                            facecolor=C_PTFE, edgecolor=EC, lw=0.6))
    ax2.add_patch(Rectangle((z0 - g["alFoil"] - g["ptfe"] - 3.0,
                             -0.5 * g["cryY"]), 3.0, g["cryY"],
                            facecolor=C_CSI, edgecolor=EC, lw=0.8))
    ax2.text(z0 - g["alFoil"] - g["ptfe"] - 1.5, 0, "CsI(Tl)", fontsize=8,
             ha="center", va="center", rotation=90, color="#4a3405")
    # стрелка кванта
    ax2.annotate("", xy=(z0 + depth - 0.15, 0), xytext=(zOut + 2.6, 0),
                 arrowprops=dict(arrowstyle="->", color="#b03030", lw=1.6))
    ax2.text(zOut + 2.8, 0.9, "квант", fontsize=8.5, color="#b03030",
             ha="left")
    # размеры
    ax2.annotate("", xy=(z0 + depth, 0.5 * winY + 0.6),
                 xytext=(zOut, 0.5 * winY + 0.6),
                 arrowprops=dict(arrowstyle="<->", color="#1f4e79", lw=0.9))
    ax2.text(z0 + depth + 0.5 * g["wCapWin"], 0.5 * winY + 1.0,
             "%s" % ru(g["wCapWin"]), fontsize=8, ha="center", color="#1f4e79")
    ax2.annotate("", xy=(z0, -yT - 1.2), xytext=(zOut, -yT - 1.2),
                 arrowprops=dict(arrowstyle="<->", color="#1f4e79", lw=0.9))
    ax2.text(z0 + 0.5 * g["wCap"], -yT - 2.6, "крышка %s" % ru(g["wCap"]),
             fontsize=8, ha="center", color="#1f4e79")
    ax2.annotate("наружная плоскость корпуса —\nот неё отсчитаны 10 см",
                 xy=(zOut, -0.5 * winY - 0.5), xytext=(zOut + 3.2, -yT - 2.0),
                 fontsize=7.5, ha="left", va="top", color="#3a3f44",
                 arrowprops=dict(arrowstyle="->", color="#3a3f44", lw=0.7))
    ax2.text(z0 + 0.2, 0.5 * winY - 1.6, "выборка %s мм" % ru(depth),
             fontsize=7.5, ha="left", color="#7a2020")
    ax2.set_xlim(z0 - g["ptfe"] - g["alFoil"] - 5.0, zOut + 17.0)
    ax2.set_ylim(-yT - 8, yT + 3)

    fig.suptitle("AtomSpectra Nano 16 PRO — передняя крышка: алюминий %s мм "
                 "с фрезеровкой до %s мм напротив кристалла"
                 % (ru(g["wCap"]), ru(g["wCapWin"])), fontsize=12, y=0.955)
    fig.text(0.5, 0.015,
             "В пучке стоит стек ОКНА: ПТФЭ %s + фольга %s + Al %s = "
             "%s г/см².   Вне окна крышка %s даёт %s г/см², то есть в %s раза "
             "больше.\n"
             "Сторона выборки источником не задана; принята внутренняя — "
             "наружная плоскость корпуса тогда остаётся на месте и опорная "
             "плоскость замера не меняется."
             % (ru(g["ptfe"]), ru(g["alFoil"]), ru(g["wCapWin"]),
                ru(sd_win, 4), ru(g["wCap"]), ru(sd_rim, 4),
                ru(sd_rim / sd_win, 2)),
             fontsize=8.2, ha="center", color="#555555")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=160)
    print("записано: %s" % OUT)
    print("стек в окне  %s г/см²" % ru(sd_win, 6))
    print("стек вне окна %s г/см²" % ru(sd_rim, 6))
    print("окно %s x %s мм, выборка %s мм, остаток %s мм"
          % (ru(winX), ru(winY), ru(depth), ru(g["wCapWin"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
