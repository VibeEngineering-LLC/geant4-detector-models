# -*- coding: utf-8 -*-
"""Схема замера «прибор лежит на пачке электродов WT-20» — К СОГЛАСОВАНИЮ.

Рисунок выпускается ПЕРЕД прогонами: постановка задана словами оператора и
фотографией пачки, и до счёта её надо увидеть. Размеры ЧИТАЮТСЯ из
`geometry/ASN16Detector.hh` (Nano16Geom), а не дублируются здесь: зеркало
констант в этом дереве уже давало расхождение рисунка с моделью.

Оси как в модели: начало — центр кристалла; +Y к рабочей грани (вверх),
пачка лежит ПОД дном корпуса; +Z к переднему торцу; X — по ширине корпуса.

    python analysis/draw_wt20_setup.py
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

_HERE = os.path.dirname(os.path.abspath(__file__))
HH = os.path.normpath(os.path.join(_HERE, "..", "geometry", "ASN16Detector.hh"))
OUT = os.path.normpath(os.path.join(_HERE, "..", "drawings",
                                    "nano16pro_wt20_setup.png"))


def geom_from_header(path):
    src = open(path, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^\s*(?:double|int)\s+(\w+)\s*=\s*([0-9.]+)\s*;", src,
                         re.MULTILINE):
        out[m.group(1)] = float(m.group(2))
    need = ("cryX", "cryY", "cryZ", "ptfe", "alFoil", "bodyX", "bodyY", "bodyZ",
            "wFront", "wBot", "wSide", "wCap", "pcbT", "sipmT",
            "wt20N", "wt20D", "wt20L", "wt20Pitch", "wt20Wall", "wt20Gap",
            "wt20ThO2", "rhoThO2", "feetH", "feetXY", "feetInset")
    miss = [k for k in need if k not in out]
    if miss:
        raise SystemExit("в %s не найдены поля: %s" % (path, ", ".join(miss)))
    return out


_G = geom_from_header(HH)
CRY_X, CRY_Y, CRY_Z = _G["cryX"], _G["cryY"], _G["cryZ"]
PTFE, ALFOIL = _G["ptfe"], _G["alFoil"]
BODY_X, BODY_Y, BODY_Z = _G["bodyX"], _G["bodyY"], _G["bodyZ"]
W_FRONT, W_BOT, W_SIDE, W_CAP = (_G["wFront"], _G["wBot"], _G["wSide"],
                                 _G["wCap"])
PCB_T, SIPM_T = _G["pcbT"], _G["sipmT"]
N_ROD, D_ROD, L_ROD = int(_G["wt20N"]), _G["wt20D"], _G["wt20L"]
PITCH, WALL, GAP = _G["wt20Pitch"], _G["wt20Wall"], _G["wt20Gap"]
THO2, RHO_THO2 = _G["wt20ThO2"], _G["rhoThO2"]

# --- границы прибора: та же арифметика, что в ASN16Detector.cc ---------------
yCryT, yCryB = +CRY_Y / 2, -CRY_Y / 2
yPtfeT, yPtfeB = yCryT + PTFE, yCryB - PTFE
yFoilT, yFoilB = yPtfeT + ALFOIL, yPtfeB - ALFOIL
yBodyT = yFoilT + W_FRONT
yBodyB = yBodyT - BODY_Y
yCavB = yBodyB + W_BOT
yPcbT, yPcbB = yFoilB, yFoilB - PCB_T

xCry = CRY_X / 2
xPtfe, xFoil = xCry + PTFE, xCry + PTFE + ALFOIL
xBody = BODY_X / 2
xCav = xBody - W_SIDE

zCryF, zCryB = +CRY_Z / 2, -CRY_Z / 2
zPtfeF = zCryF + PTFE
zFoilF = zPtfeF + ALFOIL
zCapFi = zFoilF
zBodyF = zCapFi + W_CAP
zBodyB = zBodyF - BODY_Z
zCapBi = zBodyB + W_CAP
zSipmB = zCryB - SIPM_T
zMid = 0.5 * (zBodyF + zBodyB)

# --- границы пачки -----------------------------------------------------------
R = 0.5 * D_ROD
HALF_SPAN = 0.5 * (N_ROD - 1) * PITCH
FEET_H, FEET_XY, FEET_IN = _G["feetH"], _G["feetXY"], _G["feetInset"]
# Прибор СТОИТ НА НОЖКАХ, а не лежит дном: между дном корпуса и крышкой пенала
# остаётся воздух высотой FEET_H (ОПЕРАТОР, 06.08.2026).
yCaseT = yBodyB - FEET_H
yCaseInT = yCaseT - WALL
yRod = yCaseInT - GAP - R
yCaseInB = yRod - R - GAP
yCaseB = yCaseInB - WALL
xCaseIn = HALF_SPAN + R + GAP
xCaseOut = xCaseIn + WALL
zCaseIn = 0.5 * L_ROD + GAP
zCaseOut = zCaseIn + WALL

# плотность сплава — по правилу смеси, как в ASN16Detector::MakeWT20
RHO_W = 19.30                        # G4_W, база NIST
_f = THO2 / 100.0
RHO_ROD = 1.0 / ((1.0 - _f) / RHO_W + _f / RHO_THO2)
V_PACK = N_ROD * 3.141592653589793 * (R / 10.0) ** 2 * (L_ROD / 10.0)   # см³
M_PACK = V_PACK * RHO_ROD                                              # г

C_AL = "#9aa5ad"
C_CAP = "#7f8a93"
C_CSI = "#d8ad3c"
C_PTFE = "#f4f4f0"
C_FOIL = "#c9cdd1"
C_SIPM = "#4a6fa5"
C_PCB = "#2f7a4a"
C_AIR = "#eef3f7"
C_CASE = "#bcd3ee"
C_ROD = "#8c2f2f"
EC = "#3a3f44"
EC_AIR = "#aab6c2"


def ru(x, nd=2):
    return ("%.*f" % (nd, x)).replace(".", ",")


def rect(ax, x0, x1, y0, y1, fc, label=None, ec=EC, lw=0.7, z=2):
    ax.add_patch(Rectangle((min(x0, x1), min(y0, y1)), abs(x1 - x0),
                           abs(y1 - y0), facecolor=fc, edgecolor=ec,
                           linewidth=lw, label=label, zorder=z))


def dim(ax, p0, p1, text, off, vertical=False, fs=7.2, color="#1f4e79"):
    if vertical:
        ax.annotate("", xy=(off, p1), xytext=(off, p0),
                    arrowprops=dict(arrowstyle="<->", color=color, lw=0.9))
        ax.text(off, 0.5 * (p0 + p1), " " + text, color=color, fontsize=fs,
                ha="left", va="center")
    else:
        ax.annotate("", xy=(p1, off), xytext=(p0, off),
                    arrowprops=dict(arrowstyle="<->", color=color, lw=0.9))
        ax.text(0.5 * (p0 + p1), off, text, color=color, fontsize=fs,
                ha="center", va="bottom")


def draw_pack(ax, axis):
    """Пенал и стержни. axis='x' — поперёк стержней, axis='z' — вдоль."""
    if axis == "x":
        rect(ax, -xCaseOut, xCaseOut, yCaseB, yCaseT, C_CASE,
             label="пенал, акрил %s мм" % ru(WALL), z=1.4)
        rect(ax, -xCaseIn, xCaseIn, yCaseInB, yCaseInT, C_AIR, ec=EC_AIR,
             lw=0.5, z=1.5)
        for i in range(N_ROD):
            x = -HALF_SPAN + i * PITCH
            ax.add_patch(Circle((x, yRod), R, facecolor=C_ROD, edgecolor=EC,
                                linewidth=0.6, zorder=1.7,
                                label=("W + %s %% ThO2, Ø%s мм"
                                       % (ru(THO2, 0), ru(D_ROD, 1)))
                                if i == 0 else None))
    else:
        rect(ax, zMid - zCaseOut, zMid + zCaseOut, yCaseB, yCaseT, C_CASE,
             z=1.4)
        rect(ax, zMid - zCaseIn, zMid + zCaseIn, yCaseInB, yCaseInT, C_AIR,
             ec=EC_AIR, lw=0.5, z=1.5)
        rect(ax, zMid - 0.5 * L_ROD, zMid + 0.5 * L_ROD, yRod - R, yRod + R,
             C_ROD, z=1.7)


# --- рисунок -----------------------------------------------------------------
fig = plt.figure(figsize=(14.2, 9.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1], hspace=0.20,
                      left=0.045, right=0.985, top=0.845, bottom=0.035)

# ===================== поперечный разрез (плоскость Z = 0) ===================
ax = fig.add_subplot(gs[0])
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("Поперёк стержней (плоскость Z = 0): прибор СТОРОНОЙ ПЛАТЫ, "
             "стоит на ножках — излучение входит через дно корпуса и плату",
             fontsize=9.5, pad=6)

rect(ax, -xBody, xBody, yBodyB, yBodyT, C_AL, label="корпус Al")
rect(ax, -xCav, xCav, yCavB, yFoilT, C_AIR, label="воздух полости", ec=EC_AIR,
     lw=0.5, z=2.1)
rect(ax, -xCav, xCav, yPcbB, yPcbT, C_PCB, label="плата %s мм" % ru(PCB_T),
     z=2.8)
rect(ax, -xFoil, xFoil, yFoilB, yFoilT, C_FOIL,
     label="Al-фольга %s" % ru(ALFOIL), z=3)
rect(ax, -xPtfe, xPtfe, yPtfeB, yPtfeT, C_PTFE, label="ПТФЭ %s" % ru(PTFE),
     z=4)
rect(ax, -xCry, xCry, yCryB, yCryT, C_CSI, label="CsI(Tl)", z=5)
draw_pack(ax, "x")
# ножки: в поперечном разрезе они стоят по краям корпуса
for s in (-1, +1):
    rect(ax, s * (xBody - FEET_IN - FEET_XY), s * (xBody - FEET_IN),
         yCaseT, yCaseT + FEET_H, "#3a3a40",
         label="ножки %s мм (силикон)" % ru(FEET_H) if s < 0 else None, z=3)

dim(ax, yRod, yCryB, "%s мм от оси стержня до нижней грани кристалла"
    % ru(yCryB - yRod), off=xCaseOut + 3, vertical=True)
dim(ax, -HALF_SPAN - R, HALF_SPAN + R, "%s мм — вся пачка (%d × Ø%s, шаг %s)"
    % (ru(2 * (HALF_SPAN + R)), N_ROD, ru(D_ROD, 1), ru(PITCH)),
    off=yCaseB - 3.4)
dim(ax, -xBody, xBody, "корпус %s мм" % ru(BODY_X), off=yBodyT + 2.2)
ax.annotate("стек в пучке снизу вверх:\nакрил %s + воздух %s + НОЖКИ %s "
            "(воздух) + Al дна %s + плата %s"
            % (ru(WALL), ru(GAP), ru(FEET_H), ru(W_BOT), ru(PCB_T)),
            xy=(xCry * 0.55, yPcbB), xytext=(xBody + 6, yCaseInB - 2),
            fontsize=7.6, color="#7a2020", ha="left",
            arrowprops=dict(arrowstyle="->", color="#7a2020", lw=0.8))
ax.set_xlim(-xCaseOut - 24, xCaseOut + 30)
ax.set_ylim(yCaseB - 7, yBodyT + 7)
h, lab = ax.get_legend_handles_labels()
ax.legend(h, lab, loc="upper left", fontsize=7.4, frameon=False, ncol=4,
          bbox_to_anchor=(0.0, 1.02), handlelength=1.4)

# ===================== продольный разрез (плоскость X = 0) ==================
ax2 = fig.add_subplot(gs[1])
ax2.set_aspect("equal")
ax2.axis("off")
ax2.set_title("Вдоль стержней (плоскость X = 0): корпус %s мм лежит по "
              "середине стержней %s мм" % (ru(BODY_Z), ru(L_ROD)),
              fontsize=9.5, pad=6)

rect(ax2, zBodyB, zBodyF, yBodyB, yBodyT, C_AL)
rect(ax2, zBodyB, zBodyF, yCavB, yFoilT, C_AIR, ec=EC_AIR, lw=0.5, z=2.1)
rect(ax2, zBodyF, zCapFi, yCavB, yFoilT, C_CAP,
     label="торцевые крышки Al %s" % ru(W_CAP), z=2.6)
rect(ax2, zCapBi, zBodyB, yCavB, yFoilT, C_CAP, z=2.6)
rect(ax2, zCapFi, zCapBi, yPcbB, yPcbT, C_PCB, z=2.8)
rect(ax2, zFoilF, zCryB, yFoilB, yFoilT, C_FOIL, z=3)
rect(ax2, zPtfeF, zCryB, yPtfeB, yPtfeT, C_PTFE, z=4)
rect(ax2, zCryF, zCryB, yCryB, yCryT, C_CSI, z=5)
rect(ax2, zCryB, zSipmB, yCryB, yCryT, C_SIPM, label="SiPM", z=5)
draw_pack(ax2, "z")
for zc in (zMid + 0.5 * BODY_Z - FEET_IN - FEET_XY,
           zMid - 0.5 * BODY_Z + FEET_IN):
    rect(ax2, zc, zc + FEET_XY, yCaseT, yCaseT + FEET_H, "#3a3a40", z=3)

ax2.plot([zMid, zMid], [yCaseB - 4, yBodyT + 4], color="#8a6d3b", lw=0.8,
         ls=":", zorder=6)
ax2.text(zMid, yBodyT + 5, "центр корпуса z = %s" % ru(zMid), fontsize=7.4,
         ha="center", color="#8a6d3b")
ax2.plot([0, 0], [yCaseB - 4, yBodyT + 4], color="#4a3405", lw=0.8, ls="--",
         zorder=6)
ax2.text(0, yCaseB - 5.6, "центр кристалла z = 0", fontsize=7.4, ha="center",
         va="top", color="#4a3405")
dim(ax2, zMid - 0.5 * L_ROD, zMid + 0.5 * L_ROD, "стержни %s мм" % ru(L_ROD),
    off=yCaseB - 3.4)
dim(ax2, zBodyB, zBodyF, "корпус %s мм" % ru(BODY_Z), off=yBodyT + 2.2)
ax2.set_xlim(zMid - zCaseOut - 12, zMid + zCaseOut + 12)
ax2.set_ylim(yCaseB - 12, yBodyT + 9)
h2, lab2 = ax2.get_legend_handles_labels()
ax2.legend(h2, lab2, loc="upper right", fontsize=7.4, frameon=False, ncol=2,
           handlelength=1.4)

fig.suptitle("AtomSpectra Nano 16 PRO на пачке электродов WT-20 — принятая "
             "геометрия замера 01.06.2024 (размеры в мм)",
             fontsize=11.5, y=0.985)
fig.text(0.045, 0.955,
         "ЭТИКЕТКА ПАЧКИ: «Tungsten Electrodes 3.2 mm × 175 mm, 10 Pieces, "
         "2%%Thoriated(WT20), ANSI/AWS A5.12M-98, ISO 6848».\n"
         "ОПЕРАТОР: прибор стороной платы, вдоль электродов, по центру; "
         "крышка пенала акрил %s; ножки-пуговички 3–4 мм по углам; "
         "шаг укладки %s ЗАМЕРЕН ПО ФОТО (±0,3), гнёзда свободные.\n"
         "ДОПУЩЕНИЯ: зазор стержни–крышка %s, плотность ThO2 %s г/см³ "
         "(PNNL-15870). Сплав %s г/см³ по правилу смеси, пачка %s см³ = "
         "%s г, тория %s г при номинале 2 %% масс. (допуск ISO 6848 "
         "1,70…2,20 %%)."
         % (ru(WALL), ru(PITCH), ru(GAP), ru(RHO_THO2), ru(RHO_ROD, 3),
            ru(V_PACK, 3), ru(M_PACK, 1),
            ru(M_PACK * _f * 232.038 / 264.038, 3)),
         fontsize=8.2, va="top", ha="left", linespacing=1.5, color="#7a2020")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=160)
print("записано: %s" % OUT)
print("ось стержней      y = %s мм" % ru(yRod))
print("дно корпуса       y = %s мм" % ru(yBodyB))
print("низ кристалла     y = %s мм (от оси стержня %s мм)"
      % (ru(yCryB), ru(yCryB - yRod)))
print("центр корпуса     z = %s мм" % ru(zMid))
print("пачка: %d × Ø%s × %s мм, шаг %s, масса %s г, плотность %s г/см³"
      % (N_ROD, ru(D_ROD, 1), ru(L_ROD, 0), ru(PITCH), ru(M_PACK, 1),
         ru(RHO_ROD, 3)))
print("розыгрыш GPS: halfx %s, halfy %s, halfz %s, centre 0 %s %s мм"
      % (ru(HALF_SPAN + R, 2), ru(R, 2), ru(0.5 * L_ROD, 1), ru(yRod, 2),
         ru(zMid, 2)))
