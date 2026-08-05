# -*- coding: utf-8 -*-
"""Схема замера: точечный источник по оси кристалла, 10 см от торца.

Рисунок к согласованию ПЕРЕД прогонами. Размеры зеркалят
`geometry/ASN16Detector.hh` (Nano16Geom); при правке .hh сверять здесь —
скрипт своей связи с моделью не имеет, это отдельный источник тех же чисел.

Оси как в модели: начало — центр кристалла, +Z к переднему торцу и к
источнику, +Y к рабочей грани 18 x 57 под стенкой 1,20 мм.

    python analysis/draw_setup.py
"""
import os
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# --- зеркало Nano16Geom (мм) -------------------------------------------------
CRY_X, CRY_Y, CRY_Z = 18.00, 15.00, 57.00
PTFE, ALFOIL = 1.00, 0.10
BODY_X, BODY_Y, BODY_Z = 42.00, 25.00, 86.00
W_FRONT, W_BOT, W_SIDE, W_CAP = 1.20, 2.05, 1.95, 2.00
PCB_T, SIPM_T = 1.60, 1.50

# --- параметры замера --------------------------------------------------------
DIST = 100.00        # ОПЕРАТОР: 10 см от торца
# Конус розыгрыша GPS. Выбран с большим запасом намеренно: усечение конуса
# занижает eps и искажает ФОРМУ кривой (урок точечных сеток Гамма-1С). При 35°
# на плоскости переднего торца конус накрывает радиус 70 мм при 25,9 мм до
# самого дальнего угла корпуса, то есть корпус целиком внутри конуса по всей
# своей длине.
THETA_DEG = 35.0

# --- производные границы (та же арифметика, что в ASN16Detector.cc) ----------
# Порядок слоёв: ПТФЭ НА КРИСТАЛЛЕ, фольга НА ПТФЭ (уточнение оператора
# 05.08.2026). Зеркало обязано повторять ASN16Detector.cc:117-119, а не
# наоборот: до 06.08.2026 здесь стоял до-разворотный порядок по ВСЕМ ТРЁМ
# осям, и отслеживаемая «согласованная схема» рисовала обёртку перевёрнутой
# (найдено аудитом кода). Габариты и все печатаемые числа от порядка не
# зависят — сумма толщин та же, — поэтому расхождение не проявлялось нигде,
# кроме самого рисунка слоёв.
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
zBodyF = zCapFi + W_CAP           # наружный торец корпуса, +31,60
zBodyB = zBodyF - BODY_Z
zCapBi = zBodyB + W_CAP
zSipmB = zCryB - SIPM_T

zSrc = zBodyF + DIST              # источник, +131,60
yBodyC = 0.5 * (yBodyB + yBodyT)  # ось КОРПУСА по высоте: −2,70

FRAC = 0.5 * (1.0 - math.cos(math.radians(THETA_DEG)))

# --- цвета -------------------------------------------------------------------
C_AL = "#9aa5ad"
C_ABS = "#5b5f63"
C_CSI = "#d8ad3c"
C_PTFE = "#f4f4f0"
C_FOIL = "#c9cdd1"
C_SIPM = "#4a6fa5"
C_PCB = "#2f7a4a"
C_AIR = "#eef3f7"
C_CONE = "#c0392b"
EC = "#3a3f44"
EC_AIR = "#aab6c2"


def ru(x, nd=2):
    return ("%.*f" % (nd, x)).replace(".", ",")


def rect(ax, x0, x1, y0, y1, fc, label=None, ec=EC, lw=0.7, z=2):
    ax.add_patch(Rectangle((min(x0, x1), min(y0, y1)),
                           abs(x1 - x0), abs(y1 - y0),
                           facecolor=fc, edgecolor=ec, linewidth=lw,
                           label=label, zorder=z))


def dim(ax, p0, p1, text, off, fs=7.5, color="#1f4e79"):
    ax.annotate("", xy=(p1, off), xytext=(p0, off),
                arrowprops=dict(arrowstyle="<->", color=color, lw=0.9))
    ax.text(0.5 * (p0 + p1), off, text, color=color, fontsize=fs,
            ha="center", va="bottom")


def draw_body(ax, half_lo, half_hi, side):
    """Прибор в разрезе. side='y' — вид сбоку (Y–Z), side='x' — план (X–Z)."""
    if side == "y":
        rect(ax, zBodyB, zBodyF, yBodyB, yBodyT, C_AL, label="корпус Al")
        rect(ax, zBodyB, zBodyF, yCavB, yPtfeT, C_AIR, label="воздух полости",
             ec=EC_AIR, lw=0.5, z=2.1)
        rect(ax, zBodyF, zCapFi, yCavB, yPtfeT, C_ABS, label="крышки ABS", z=2.6)
        rect(ax, zCapBi, zBodyB, yCavB, yPtfeT, C_ABS, z=2.6)
        rect(ax, zCapFi, zCapBi, yPcbB, yPcbT, C_PCB, label="плата", z=2.8)
        rect(ax, zPtfeF, zCryB, yPtfeB, yPtfeT, C_PTFE, label="ПТФЭ 1,00", z=3)
        rect(ax, zFoilF, zCryB, yFoilB, yFoilT, C_FOIL, label="Al-фольга 0,10",
             z=4)
        rect(ax, zCryF, zCryB, yCryB, yCryT, C_CSI, label="CsI(Tl)", z=5)
        rect(ax, zCryB, zSipmB, yCryB, yCryT, C_SIPM, label="SiPM", z=5)
    else:
        rect(ax, zBodyB, zBodyF, -xBody, xBody, C_AL)
        rect(ax, zBodyB, zBodyF, -xCav, xCav, C_AIR, ec=EC_AIR, lw=0.5, z=2.1)
        rect(ax, zBodyF, zCapFi, -xCav, xCav, C_ABS, z=2.6)
        rect(ax, zCapBi, zBodyB, -xCav, xCav, C_ABS, z=2.6)
        rect(ax, zPtfeF, zCryB, -xPtfe, xPtfe, C_PTFE, z=3)
        rect(ax, zFoilF, zCryB, -xFoil, xFoil, C_FOIL, z=4)
        rect(ax, zCryF, zCryB, -xCry, xCry, C_CSI, z=5)
        rect(ax, zCryB, zSipmB, -xCry, xCry, C_SIPM, z=5)


def draw_source(ax, half):
    """Источник, ось розыгрыша и границы конуса."""
    tan = math.tan(math.radians(THETA_DEG))
    # конус до задней плоскости корпуса
    dz = zSrc - zBodyB
    for s in (+1, -1):
        ax.plot([zSrc, zBodyB], [0, s * tan * dz], color=C_CONE, lw=0.9,
                ls="--", zorder=6)
    ax.plot([zSrc, zBodyB], [0, 0], color="#555555", lw=0.7, ls="-.", zorder=6)
    ax.plot([zSrc], [0], marker="o", ms=6, color=C_CONE, zorder=7)
    ax.text(zSrc, half * 0.10, "  точечный источник\n  Cs-137, 661,657 кэВ",
            fontsize=8, ha="left", va="bottom", color=C_CONE, zorder=7)


# --- рисунок -----------------------------------------------------------------
fig = plt.figure(figsize=(13.6, 8.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.18,
                      left=0.035, right=0.985, top=0.855, bottom=0.03)

# ============================ вид сбоку ======================================
ax = fig.add_subplot(gs[0])
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("Вид сбоку (плоскость X = 0). Источник на оси КРИСТАЛЛА (y = 0); "
             "ось корпуса ниже на %s мм" % ru(-yBodyC), fontsize=9.5, pad=6)
draw_body(ax, yBodyB, yBodyT, "y")
draw_source(ax, 45)

# ось корпуса — отдельной линией, чтобы смещение было видно, а не подразумевалось
ax.plot([zBodyB - 4, zBodyF + 12], [yBodyC, yBodyC], color="#8a6d3b", lw=0.8,
        ls=":", zorder=6)
ax.text(zBodyB - 5, yBodyC, "ось корпуса  ", fontsize=7.5, ha="right",
        va="center", color="#8a6d3b")

dim(ax, zBodyF, zSrc, "%s (10 см от НАРУЖНОЙ поверхности крышки)" % ru(DIST),
    off=yBodyT + 6)
dim(ax, zCryF, zSrc, "%s до передней грани кристалла" % ru(zSrc - zCryF),
    off=yBodyT + 15)
ax.annotate("крышка ABS %s + обёртка %s = %s мм\nмежду опорной плоскостью и "
            "кристаллом" % (ru(W_CAP), ru(PTFE + ALFOIL), ru(zBodyF - zCryF)),
            xy=(0.5 * (zCryF + zBodyF), yCryB - 1), xytext=(zBodyF + 16,
                                                           yBodyB - 13),
            fontsize=7.5, color="#7a2020", ha="left",
            arrowprops=dict(arrowstyle="->", color="#7a2020", lw=0.8))
ax.text(zSrc - 12, -34, "конус розыгрыша ±%s°" % ru(THETA_DEG, 0),
        fontsize=7.5, color=C_CONE, ha="right")
ax.set_xlim(zBodyB - 34, zSrc + 46)
ax.set_ylim(-46, 36)
h, lab = ax.get_legend_handles_labels()
ax.legend(h, lab, loc="lower left", fontsize=7.4, frameon=False, ncol=4,
          bbox_to_anchor=(0.0, -0.02), handlelength=1.4)

# ============================ план ===========================================
ax2 = fig.add_subplot(gs[1])
ax2.set_aspect("equal")
ax2.axis("off")
ax2.set_title("План (плоскость Y = 0): источник на оси кристалла и по ширине "
              "корпуса", fontsize=9.5, pad=6)
draw_body(ax2, -xBody, xBody, "x")
draw_source(ax2, 45)
dim(ax2, zBodyF, zSrc, ru(DIST), off=xBody + 6)
ax2.text(zBodyB - 4, xBody + 3, "торец 18 × 15 мм = 2,70 см²; рабочая грань "
                                "18 × 57 мм = 10,26 см², отношение 3,80",
         fontsize=7.5, color="#4a3405", ha="left")
ax2.set_xlim(zBodyB - 34, zSrc + 46)
ax2.set_ylim(-46, 36)

# ============================ шапка ==========================================
fig.suptitle("AtomSpectra Nano 16 PRO — схема опорного замера: точечный "
             "источник по оси, 10 см от торца (размеры в мм)",
             fontsize=11.5, y=0.985)
fig.text(0.035, 0.955,
         "К СОГЛАСОВАНИЮ. Источник поставлен на ось КРИСТАЛЛА (y = 0), а не "
         "корпуса: «по оси детектора» прочитано как ось чувствительного "
         "объёма. Расстояние отсчитано\nот наружной поверхности передней "
         "крышки (z = %s), до передней грани кристалла остаётся ещё %s мм — "
         "на плече 100 мм это ~%s %% по телесному углу.  Конус ±%s°, "
         "доля (1−cos θ)/2 = %s."
         % (ru(zBodyF), ru(zBodyF - zCryF), ru(200 * (zBodyF - zCryF) / DIST, 0),
            ru(THETA_DEG, 0), ru(FRAC, 4)),
         fontsize=8.2, va="top", ha="left", linespacing=1.5, color="#7a2020")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "drawings", "nano16pro_setup_point10cm.png")
out = os.path.normpath(out)
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=160)
print("записано: %s" % out)

print("источник            z = %s мм (ось кристалла, x = y = 0)" % ru(zSrc))
print("наружный торец      z = %s мм" % ru(zBodyF))
print("передняя грань CsI  z = %s мм" % ru(zCryF))
print("ось корпуса         y = %s мм (источник выше неё)" % ru(yBodyC))
print("конус               %s°, доля (1−cos θ)/2 = %s" % (ru(THETA_DEG, 0),
                                                          ru(FRAC, 5)))
print("радиус конуса на плоскости торца  %s мм при %s мм до дальнего угла"
      % (ru(DIST * math.tan(math.radians(THETA_DEG)), 1),
         ru(math.hypot(xBody, max(abs(yBodyB), abs(yBodyT))), 1)))
