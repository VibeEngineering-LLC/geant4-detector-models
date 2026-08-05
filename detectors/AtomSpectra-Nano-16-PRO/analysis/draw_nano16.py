# -*- coding: utf-8 -*-
"""Разрез AtomSpectra Nano 16 PRO — воспроизводимый чертёж из справочных данных.

Опубликованного чертежа прибора нет. Размеры собраны со слов оператора
(29.07.2026, уточнения 05.08.2026), с чертежей профиля экструзии линейки Nano
(чертёж Nano 5 PRO) и с фотографии открытого торца прибора. Всё, что не задано
ни одним из этих источников, помечено звёздочкой (*) и в подписи, и в служебной
печати скрипта. Рисунок предназначен для согласования геометрии ДО постройки
модели Geant4, а не как воспроизведение заводского чертежа.

Уточнения оператора 05.08.2026:
  - свободный объём корпуса — воздух; сплошного «подвала электроники» как
    вещества нет, но печатная плата в модели остаётся: она лежит по дну полости
    во всю длину, кристалл опирается на неё (видно на фото торца);
  - SiPM подключён с торца, на границе кристалла и воздушной зоны, то есть на
    задней грани 18 × 15 мм; обёртки на этой грани нет (оптический контакт);
  - торцевые крышки — ПЛАСТИК типа ABS, не алюминий. Для мягкого края это
    существенно: опорный замер Cs-137 снят через переднюю крышку.

Оси (совпадают с геометрией Geant4), начало — центр кристалла:
  Z — вдоль корпуса (86 мм); +Z к переднему торцу (грань 18 × 15 мм) и к
      источнику, −Z к грани с SiPM. До 05.08.2026 этот чертёж был зеркален
      модели по Z — рисовался до того, как ось развернули; расхождение нашёл
      независимый аудит. Второй чертёж каталога (`draw_setup.py`) с самого
      начала строился по модели, и два чертежа одного прибора смотрели в
      разные стороны;
  X — по ширине корпуса (42 мм);
  Y — по высоте (25 мм); +Y к грани 18 × 57 мм под тонкой стенкой 1,20 мм.

    python analysis/draw_nano16.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# --- параметры геометрии (мм) ------------------------------------------------
# кристалл CsI(Tl): брусок, не цилиндр
CRY_X, CRY_Y, CRY_Z = 18.00, 15.00, 57.00
RHO_CSI = 4.51
# обёртка кристалла — на пяти гранях, кроме задней (там SiPM)
PTFE, RHO_PTFE = 1.00, 2.20
ALFOIL = 0.10
RHO_AL = 2.70               # (*) сплав корпуса не назван
SIPM_T = 1.50               # (*) толщина сборки SiPM данными не задана
PCB_T = 1.60                # (*) толщина платы не измерена, типовая FR4
# корпус (экструзия линейки Nano, габарит Nano 16 PRO)
BODY_X, BODY_Y, BODY_Z = 42.00, 25.00, 86.00
W_FRONT = 1.20              # рабочая стенка над кристаллом (из чертежа профиля)
W_SIDE = 1.95               # (*) перенос с Nano 5 PRO: (39,50 − 35,60) / 2
W_BOT = 2.05                # (*) из выносок профиля
W_CAP = 2.00                # (*) толщина торцевой крышки; материал — ABS

# --- производные границы -----------------------------------------------------
# Y: кристалл упёрт в переднюю стенку — воздушного зазора в лицевом стеке нет
yCryT, yCryB = +CRY_Y / 2, -CRY_Y / 2
# ПТФЭ лежит НА КРИСТАЛЛЕ, фольга — НА ПТФЭ (оператор, 05.08.2026).
yPtfeT, yPtfeB = yCryT + PTFE, yCryB - PTFE
yFoilT, yFoilB = yPtfeT + ALFOIL, yPtfeB - ALFOIL
yBodyT = yFoilT + W_FRONT
yBodyB = yBodyT - BODY_Y
yBodyInB = yBodyB + W_BOT
yPcbT, yPcbB = yFoilB, yFoilB - PCB_T     # кристалл опирается на плату (фото)
AIR_UNDER = yPcbB - yBodyInB              # воздух под платой

# X: кристалл по фото стоит примерно по центру полости; точно не задан (*)
xCry = CRY_X / 2
xPtfe, xFoil = xCry + PTFE, xCry + PTFE + ALFOIL
xBody = BODY_X / 2
xBodyIn = xBody - W_SIDE
# Воздух сбоку меряется от НАРУЖНОЙ поверхности обёртки, то есть от фольги.
# После разворота слоёв (ПТФЭ на кристалле, фольга на ПТФЭ) внешним стал xFoil;
# до 06.08.2026 здесь стоял xPtfe — 9,05 мм против модельных 8,95
# (ASN16Detector.cc: xCav − xFoil), и это число дважды печаталось на разрезе.
AIR_SIDE = xBodyIn - xFoil                # воздух сбоку, с каждой стороны

# Z: обёртка кристалла упёрта в переднюю крышку, сзади SiPM без обёртки
# +Z к переднему торцу — как в ASN16Detector.cc, а не наоборот.
zCryF, zCryB = +CRY_Z / 2, -CRY_Z / 2
zPtfeF = zCryF + PTFE
zFoilF = zPtfeF + ALFOIL
zCapFin = zFoilF                          # внутренняя грань передней крышки
zBodyF = zCapFin + W_CAP
zBodyB = zBodyF - BODY_Z
zCapBin = zBodyB + W_CAP                  # внутренняя грань задней крышки
zSipmB = zCryB - SIPM_T
AIR_BEHIND = abs(zSipmB - zCapBin)        # воздух за SiPM (модуль: −Z вглубь)

# --- цвета -------------------------------------------------------------------
C_AL = "#9aa5ad"
C_ABS = "#5b5f63"
C_CSI = "#d8ad3c"
C_PTFE = "#f4f4f0"
C_FOIL = "#c9cdd1"
C_SIPM = "#4a6fa5"
C_PCB = "#2f7a4a"
C_PCB_PALE = "#d6ead9"      # та же плата на плане — она лежит ниже кристалла
C_AIR = "#eef3f7"
EC = "#3a3f44"
EC_AIR = "#aab6c2"


def ru(x, nd=2):
    """Число в русской записи: десятичная запятая (ГОСТ 8.417)."""
    return ("%.*f" % (nd, x)).replace(".", ",")


def rect(ax, x0, x1, y0, y1, fc, label=None, hatch=None, ec=EC, lw=0.7, z=2):
    """Прямоугольник по двум углам (порядок углов произволен)."""
    ax.add_patch(Rectangle((min(x0, x1), min(y0, y1)),
                           abs(x1 - x0), abs(y1 - y0),
                           facecolor=fc, edgecolor=ec, linewidth=lw,
                           hatch=hatch, label=label, zorder=z))


def dim(ax, p0, p1, text, vertical=False, off=0.0, fs=7.0, color="#1f4e79"):
    """Размерная линия со стрелками и подписью."""
    if vertical:
        ax.annotate("", xy=(off, p1), xytext=(off, p0),
                    arrowprops=dict(arrowstyle="<->", color=color, lw=0.8))
        ax.text(off, 0.5 * (p0 + p1), " " + text, color=color, fontsize=fs,
                ha="left", va="center", rotation=90)
    else:
        ax.annotate("", xy=(p1, off), xytext=(p0, off),
                    arrowprops=dict(arrowstyle="<->", color=color, lw=0.8))
        ax.text(0.5 * (p0 + p1), off, text, color=color, fontsize=fs,
                ha="center", va="bottom")


def frame(ax, title):
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9.5, pad=6)
    ax.axis("off")


# --- рисунок -----------------------------------------------------------------
fig = plt.figure(figsize=(13.6, 9.4))
gs = fig.add_gridspec(2, 2, width_ratios=[86, 46], height_ratios=[30, 46],
                      hspace=0.16, wspace=0.06,
                      left=0.03, right=0.985, top=0.90, bottom=0.03)

# ============================ 1. продольный разрез ===========================
ax = fig.add_subplot(gs[0, 0])
frame(ax, "Продольный разрез (плоскость X = 0): Z — вдоль корпуса, Y — вверх, "
          "к тонкой стенке")

rect(ax, zBodyF, zBodyB, yBodyB, yBodyT, C_AL, label="корпус Al (сплав — *)")
rect(ax, zBodyF, zCapBin + W_CAP, yBodyInB, yFoilT, C_AIR,
     label="воздух полости", ec=EC_AIR, lw=0.5, z=2.1)
rect(ax, zBodyF, zCapFin, yBodyInB, yFoilT, C_ABS,
     label="крышки ABS (толщина — *)", z=2.6)
rect(ax, zCapBin, zBodyB, yBodyInB, yFoilT, C_ABS, z=2.6)
rect(ax, zCapFin, zCapBin, yPcbB, yPcbT, C_PCB, label="плата FR4 (толщина — *)",
     z=2.8)
# обёртка доведена до задней грани кристалла и перекрыта им — сзади её нет
rect(ax, zFoilF, zCryB, yFoilB, yFoilT, C_FOIL, label="Al-фольга 0,10 мм", z=3)
rect(ax, zPtfeF, zCryB, yPtfeB, yPtfeT, C_PTFE, label="ПТФЭ 1,00 мм", z=4)
rect(ax, zCryF, zCryB, yCryB, yCryT, C_CSI, label="CsI(Tl) 18 × 15 × 57 мм",
     z=5)
rect(ax, zCryB, zSipmB, yCryB, yCryT, C_SIPM, label="SiPM (толщина — *)", z=5)

ax.text(0, yCryT - 3.0, "CsI(Tl)", fontsize=8, ha="center", va="top",
        color="#4a3405", zorder=6)
ax.text(0.5 * (zSipmB + zCapBin), 2, "воздух", fontsize=7.5, ha="center",
        va="center", color="#4a5a68", zorder=6)

dim(ax, zBodyF, zBodyB, "86,00", off=yBodyT + 2.2)
dim(ax, zCryF, zCryB, "57,00", off=yCryT + 0.6)
dim(ax, yBodyB, yBodyT, "25,00", vertical=True, off=zBodyB + 2.4)
dim(ax, yCryB, yCryT, "15,00", vertical=True, off=zSipmB + 3.6)
ax.annotate("стенка 1,20 мм\n(рабочая грань 18 × 57 мм = 10,3 см²)",
            xy=(10, 0.5 * (yFoilT + yBodyT)), xytext=(4, yBodyT + 7),
            fontsize=7.5, color="#1f4e79", ha="left",
            arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=0.8))
ax.annotate("передний торец 18 × 15 мм = 2,7 см²;\nопорный замер Cs-137, "
            "10 см — через крышку ABS",
            xy=(zBodyF, -3), xytext=(zCryF - 20, yBodyB - 12),
            fontsize=7.5, color="#7a2020", ha="left",
            arrowprops=dict(arrowstyle="->", color="#7a2020", lw=0.8))
ax.annotate("SiPM на задней грани, оптический\nконтакт (оператор): обёртки нет",
            xy=(zSipmB + 0.5, 0), xytext=(zBodyB + 1, yBodyT + 7),
            fontsize=7.5, color="#24406e", ha="left",
            arrowprops=dict(arrowstyle="->", color="#24406e", lw=0.8))
ax.annotate("плата под обёрткой во всю длину,\nкристалл опирается на неё "
            "(фото торца)",
            xy=(-8, 0.5 * (yPcbT + yPcbB)), xytext=(zBodyB + 1, yBodyB - 12),
            fontsize=7.5, color="#1c5c33", ha="left",
            arrowprops=dict(arrowstyle="->", color="#1c5c33", lw=0.8))
# Ось по возрастанию: +Z (передний торец, источник) — СПРАВА. Без явного
# порядка matplotlib инвертировал бы ось (zBodyF > zBodyB) и подписи легли бы
# зеркально, наезжая друг на друга.
ax.set_xlim(zBodyB - 10, zBodyF + 28)
ax.set_ylim(yBodyB - 16, yBodyT + 13)

# ============================ 2. вид сверху (план) ===========================
ax2 = fig.add_subplot(gs[1, 0])
frame(ax2, "План со стороны рабочей грани (плоскость X–Z): видна грань "
           "18 × 57 мм под стенкой 1,20 мм")

rect(ax2, zBodyF, zBodyB, -xBody, xBody, C_AL)
rect(ax2, zBodyF, zBodyB, -xBodyIn, xBodyIn, C_AIR, ec=EC_AIR, lw=0.5, z=2.1)
rect(ax2, zBodyF, zCapFin, -xBodyIn, xBodyIn, C_ABS, z=2.6)
rect(ax2, zCapBin, zBodyB, -xBodyIn, xBodyIn, C_ABS, z=2.6)
rect(ax2, zCapFin, zCapBin, -xBodyIn, xBodyIn, C_PCB_PALE, ec="#4f9a6a",
     lw=0.6, z=2.8)
rect(ax2, zFoilF, zCryB, -xFoil, xFoil, C_FOIL, z=3)
rect(ax2, zPtfeF, zCryB, -xPtfe, xPtfe, C_PTFE, z=4)
rect(ax2, zCryF, zCryB, -xCry, xCry, C_CSI, z=5)
rect(ax2, zCryB, zSipmB, -xCry, xCry, C_SIPM, z=5)

ax2.text(0, 0, "рабочая грань\n18 × 57 мм", fontsize=8, ha="center",
         va="center", color="#4a3405", zorder=6)
ax2.text(0.5 * (zSipmB + zCapBin), -6, "плата\n(лежит ниже)",
         fontsize=7.5, ha="center", va="center", color="#1c5c33", zorder=6)

dim(ax2, -xBody, xBody, "42,00", vertical=True, off=zBodyB + 2.4)
dim(ax2, -xCry, xCry, "18,00", vertical=True, off=zSipmB + 3.6)
dim(ax2, zBodyF, zBodyB, "86,00", off=xBody + 2.2)
ax2.annotate("боковой зазор полости %s мм — воздух;\nпо фото торца кристалл "
             "стоит примерно\nпо центру, точно не задан (*)" % ru(AIR_SIDE),
             xy=(-10, 0.5 * (xPtfe + xBodyIn)), xytext=(zBodyF - 25, xBody + 8),
             fontsize=7.5, color="#b03030", ha="left",
             arrowprops=dict(arrowstyle="->", color="#b03030", lw=0.8))
ax2.set_xlim(zBodyB - 10, zBodyF + 26)
ax2.set_ylim(-xBody - 8, xBody + 22)

# ============================ 3. поперечный разрез ===========================
ax3 = fig.add_subplot(gs[1, 1])
frame(ax3, "Поперечный разрез через кристалл (плоскость Z = 0)")

rect(ax3, -xBody, xBody, yBodyB, yBodyT, C_AL)
rect(ax3, -xBodyIn, xBodyIn, yBodyInB, yFoilT, C_AIR, ec=EC_AIR, lw=0.5, z=2.1)
rect(ax3, -xBodyIn, xBodyIn, yPcbB, yPcbT, C_PCB, z=2.8)
rect(ax3, -xFoil, xFoil, yFoilB, yFoilT, C_FOIL, z=3)
rect(ax3, -xPtfe, xPtfe, yPtfeB, yPtfeT, C_PTFE, z=4)
rect(ax3, -xCry, xCry, yCryB, yCryT, C_CSI, z=5)

ax3.text(0, 0, "CsI(Tl)\n18 × 15", fontsize=8, ha="center", va="center",
         color="#4a3405", zorder=6)
dim(ax3, -xBody, xBody, "42,00", off=yBodyT + 2.2)
dim(ax3, yBodyB, yBodyT, "25,00", vertical=True, off=xBody + 2.4)
dim(ax3, -xCry, xCry, "18,00", off=yCryT + 0.6)
ax3.annotate("плата %s мм (*), под ней воздух %s мм"
             % (ru(PCB_T), ru(AIR_UNDER)),
             xy=(-4, 0.5 * (yPcbT + yPcbB)), xytext=(-xBody - 4, yBodyB - 6),
             fontsize=7.5, color="#1c5c33", ha="left",
             arrowprops=dict(arrowstyle="->", color="#1c5c33", lw=0.8))
ax3.annotate("боковая стенка %s мм (*)" % ru(W_SIDE),
             xy=(-xBodyIn - 0.5 * W_SIDE, 4), xytext=(-xBody - 4, yBodyB - 10),
             fontsize=7.5, color="#b03030", ha="left",
             arrowprops=dict(arrowstyle="->", color="#b03030", lw=0.8))

# ВЫНОСКА СЛОЁВ лицевого стека, снаружи внутрь. Порядок обёртки задан
# оператором: ПТФЭ прилегает к кристаллу, фольга лежит на ПТФЭ.
stack_layers = [
    ("стенка корпуса Al", W_FRONT, C_AL),
    ("Al-фольга", ALFOIL, C_FOIL),
    ("ПТФЭ (на кристалле)", PTFE, C_PTFE),
    ("CsI(Tl)", None, C_CSI),
]
lx0, lx1 = xBody + 6.0, xBody + 10.0    # столбик-выноска справа
ly = yCryB                              # снизу вверх, вровень с кристаллом
ax3.annotate("", xy=(xFoil, yCryT * 0.6), xytext=(lx0, ly + 4 * 2.4),
             arrowprops=dict(arrowstyle="->", color="#555555", lw=0.8))
ax3.text(lx0, ly + 4 * 2.4 + 1.2, "порядок слоёв рабочей грани,\nснаружи "
         "внутрь:", fontsize=7.0, ha="left", va="bottom", color="#3a3f44")
yc = ly
for name, t, col in reversed(stack_layers):   # снизу — кристалл, вверх наружу
    h = 2.2
    rect(ax3, lx0, lx1, yc, yc + h, col, ec=EC, lw=0.6, z=7)
    lab = name if t is None else "%s %s мм" % (name, ru(t))
    ax3.text(lx1 + 0.8, yc + 0.5 * h, lab, fontsize=7.0, ha="left",
             va="center", color="#3a3f44", zorder=7)
    yc += h + 0.2

ax3.set_xlim(-xBody - 20, xBody + 40)
ax3.set_ylim(yBodyB - 13, yBodyT + 15)

# ============================ 4. легенда и допущения =========================
ax4 = fig.add_subplot(gs[0, 1])
ax4.axis("off")
h, lab = ax.get_legend_handles_labels()
ax4.legend(h, lab, loc="upper left", fontsize=7.4, frameon=False,
           bbox_to_anchor=(-0.02, 1.03), handlelength=1.6)

stack = W_FRONT / 10 * RHO_AL + PTFE / 10 * RHO_PTFE + ALFOIL / 10 * RHO_AL
mass = CRY_X * CRY_Y * CRY_Z / 1000.0 * RHO_CSI
ax4.text(-0.02, 0.40,
         "Лицевой стек рабочей грани: Al %s + ПТФЭ %s + Al %s мм\n"
         "= %s г/см² (воздушного зазора нет — кристалл упёрт в стенку).\n"
         "Кристалл: %s см³, %s г при ρ = %s г/см³.\n"
         "Свободный объём полости — воздух: %s мм под платой,\n"
         "%s мм с каждой стороны кристалла, %s мм за SiPM.\n\n"
         "(*) ДОПУЩЕНИЯ, чертежом Nano 16 PRO не заданы:\n"
         "  • сплав корпуса, принято ρ = %s г/см³;\n"
         "  • боковая стенка %s мм — перенос с Nano 5 PRO;\n"
         "  • толщина торцевой крышки %s мм (материал — ABS, оператор);\n"
         "  • положение кристалла по X и Z внутри полости;\n"
         "  • толщина платы %s мм и сборки SiPM %s мм, состав обоих."
         % (ru(W_FRONT), ru(PTFE), ru(ALFOIL), ru(stack, 3),
            ru(CRY_X * CRY_Y * CRY_Z / 1000.0), ru(mass, 1), ru(RHO_CSI),
            ru(AIR_UNDER), ru(AIR_SIDE), ru(AIR_BEHIND),
            ru(RHO_AL), ru(W_SIDE), ru(W_CAP), ru(PCB_T), ru(SIPM_T)),
         fontsize=7.6, va="top", ha="left", linespacing=1.5)

fig.suptitle("AtomSpectra Nano 16 PRO — геометрия по справочным данным "
             "(размеры в мм, звёздочка — допущение)", fontsize=11.5, y=0.965)

out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "drawings", "nano16pro_section.png")
out = os.path.normpath(out)
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=160)
print("записано: %s" % out)

# --- служебная сверка --------------------------------------------------------
print("лицевой стек            %s г/см²  (в справке 0,571)" % ru(stack, 3))
print("масса кристалла         %s г       (в справке ~69)" % ru(mass, 1))
print("площадь грани 18 × 57   %s см²" % ru(CRY_X * CRY_Z / 100.0))
print("площадь торца 18 × 15   %s см²" % ru(CRY_X * CRY_Y / 100.0))
print("отношение площадей      %s" % ru(CRY_Z / CRY_Y))
print("воздух под платой       %s мм" % ru(AIR_UNDER))
print("воздух сбоку (каждая)   %s мм" % ru(AIR_SIDE))
print("воздух за SiPM          %s мм" % ru(AIR_BEHIND))
# Габарит обёртки — по НАРУЖНОМУ слою (фольга), иначе это поверка не того, что
# меряется по фото: снаружи видна фольга. До 06.08.2026 печаталось 20,00 × 17,00
# по ПТФЭ, тогда как ASN16Detector.cc:15-19 строит довод «два независимых
# размера сошлись» на 20,20 × 17,20.
print("обёртка+кристалл по X   %s мм  (на фото ~20,5 по масштабу полости)"
      % ru(2 * xFoil))
print("обёртка+кристалл по Y   %s мм  (на фото ~17,3)" % ru(yFoilT - yFoilB))
