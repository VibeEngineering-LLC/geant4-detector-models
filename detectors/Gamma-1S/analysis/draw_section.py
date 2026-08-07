# -*- coding: utf-8 -*-
"""Разрез модели Гамма-1С — воспроизводимый рисунок из параметров геометрии.

Заменяет рисованный вручную gamma1s_section.png. Размеры зеркалят
geometry/G1SDetector.hh (HeadGeom, ShieldGeom); при правке .hh — сверять здесь.

Правки по ревью оператора (чертёж ДЦКИ.412131.001):
  - дно кристалла на 5 мм выше верхней грани свинца дна укрытия
    (floorFromCryCentre −45 -> −34,5: Pb-дно сверху на −36,5, кристалл дно −31,5);
  - внешний корпус-труба постоянного Ø и толщины до низа прибора (как в модели);
  - внизу стальная платформа 10 мм; никакого «светло-серого массива» нет.

Оси: начало — центр кристалла, +Z вверх (к торцу и пробе), −Z вниз (ФЭУ, база).
Рисуется правая половина и зеркалится налево (осесимметричный разрез).

    python analysis/draw_section.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# --- параметры геометрии (мм), зеркало G1SDetector.hh ------------------------
cryDia, cryLen = 63.00, 63.00
mgoRad, mgoFace = 3.65, 6.00
alCan, rubber2, alCase = 0.50, 2.00, 1.50
alCaseFace, faceAir = 2.00, 1.00
faceSealW = 5.50
pmtDia, pmtGlass, pmtLen = 71.00, 1.50, 120.00
window, unitLen = 0.50, 315.00
# защита
cavityDia, cavityH = 200.00, 190.00
cu, cd, pb, steel = 1.00, 1.00, 50.00, 3.00
boreDia = 82.00
# ИСПРАВЛЕНО: дно кристалла 5 мм над свинцом дна (было −45)
floorFromCryCentre = -34.50
platformH = 10.00          # стальная платформа снизу (оператор, чертёж поз. 8)

# --- производные границы (как в G1SDetector.cc) ------------------------------
rCry = 0.5 * cryDia                 # 31.5
zCry = 0.5 * cryLen                 # 31.5
rMgo = rCry + mgoRad                # 35.15
rCan = rMgo + alCan                 # 35.65
rRub = rCan + rubber2               # 37.65
rCase = rRub + alCase               # 39.15
zMgoTop = zCry + mgoFace            # 37.5
zCanTop = zMgoTop + alCan           # 38.0
zAirTop = zCanTop + faceAir         # 39.0
zFace = zAirTop + alCaseFace        # 41.0
zWinBot = -zCry - window            # -32.0
zPmtBot = zWinBot - pmtLen          # -152.0
zTail = zFace - unitLen             # -274.0
rPmt = 0.5 * pmtDia                 # 35.5
rSeal = rRub - faceSealW            # 32.15

rCav = 0.5 * cavityDia              # 100
rCu, rCd, rPb, rSt = rCav + cu, rCav + cu + cd, rCav + cu + cd + pb, \
    rCav + cu + cd + pb + steel     # 101,102,152,155
rBore = 0.5 * boreDia               # 41
z0 = floorFromCryCentre             # -34.5 днище полости (верх Cu-дна)
z1 = z0 + cavityH                   # 155.5 потолок полости
zCu0, zCd0 = z0 - cu, z0 - cu - cd  # -35.5, -36.5 (верх свинца дна)
zPb0 = zCd0 - pb                    # -86.5
zSt0 = zPb0 - steel                 # -89.5
zCu1, zCd1, zPb1, zSt1 = z1 + cu, z1 + cu + cd, z1 + cu + cd + pb, \
    z1 + cu + cd + pb + steel

# цвета (близко к vis-палитре модели)
C = dict(nai="#33e04d", mgo="#f2f2f2", al="#b3b3bf", rub="#26262b",
         glass="#9cccff", vac="#d6ebff", el="#7d4d1a", pb="#5a5a66",
         cd="#9a9a80", cu="#cc8033", steel="#8a929c", air="#ffffff",
         pp="#dcdcd6", samp="#c9a978")

# --- сосуд Маринелли 1 л на головке (VesselGeom) ------------------------------
vWall, vWellR, vWellD, vOutR, vH = 2.00, 40.00, 74.00, 75.00, 110.00
vwft = zFace + vWall                # дно колодца (закрытый верх) на торце (43)
vz0 = vwft - vWellD                 # устье колодца — сосуд НАДЕВАЕТСЯ на головку (-31)
vz1 = vz0 + vH                      # верх сосуда (79)

# --- нижняя секция: ФЭУ в свинцовом цилиндре (50) и стальной трубе (3) --------
# Оператор + фото деталей: ниже измерительной камеры ФЭУ окружён СВИНЦОВЫМ
# ЦИЛИНДРОМ 50 мм в СТАЛЬНОЙ ТРУБЕ 3 мм (продолжение защиты вниз), снизу —
# свинцовая заглушка, всё на стальной платформе-тележке.
rPmtBore = rBore                    # 41 — проход под ФЭУ
rPmtPb = rPmtBore + pb              # 91 — свинец 50 мм вокруг ФЭУ
rPmtSt = rPmtPb + steel            # 94 — стальная труба 3 мм
plugH = 50.0                                   # свинцовая заглушка 50 мм (оператор)
zPlug1, zPlug0 = zTail, zTail - plugH          # заглушка -274 .. -324
zLowTop = zSt0                                 # низ главной защиты (-89.5)
zPlat1 = zPlug0
zPlat0 = zPlat1 - platformH
rPlat = 130.0                                  # база тележки шире трубы

# (метка, r0, r1, z0, z1, цвет) — правая половина
parts = [
    # ВНЕШНИЙ СТАЛЬНОЙ КОРПУС — ПОСТОЯННЫЙ Ø сверху донизу (оператор)
    ("St_corpus", rPb, rSt, zPlat1, zPb1, C["steel"]),
    ("St_top", 0, rSt, zPb1, zSt1, C["steel"]),
    # измерительная камера: свинец 50 + облицовка Cd/Cu + крышка
    ("Pb_side", rCd, rPb, zCd0, zCd1, C["pb"]),
    ("Cd_side", rCu, rCd, zCu0, zCu1, C["cd"]),
    ("Cu_side", rCav, rCu, z0, z1, C["cu"]),
    ("Pb_bottom", rBore, rPb, zPb0, zCd0, C["pb"]),
    ("Cd_bottom", rBore, rCd, zCd0, zCu0, C["cd"]),
    ("Cu_bottom", rBore, rCu, zCu0, z0, C["cu"]),
    ("Cu_top", 0, rCu, z1, zCu1, C["cu"]),
    ("Cd_top", 0, rCd, zCu1, zCd1, C["cd"]),
    ("Pb_top", 0, rPb, zCd1, zPb1, C["pb"]),
    # нижняя секция: свинцовый цилиндр вокруг ФЭУ (до дна камеры), воздух до
    # корпуса; свинцовая заглушка 50 мм снизу; стальная платформа
    ("Pb_pmt", rPmtBore, rPmtPb, zPlug1, zPb0, C["pb"]),
    ("Plug_Pb", 0, rPmtPb, zPlug0, zPlug1, C["pb"]),
    ("Platform", 0, rPlat, zPlat0, zPlat1, C["steel"]),
    # сосуд Маринелли 1 л на головке (проба ОИСН + стенки полипропилена)
    ("V_sample_side", vWellR + vWall, vOutR - vWall, vz0 + vWall, vz1 - vWall, C["samp"]),
    ("V_sample_top", 0, vWellR + vWall, vwft, vz1 - vWall, C["samp"]),
    ("V_wall_out", vOutR - vWall, vOutR, vz0, vz1, C["pp"]),
    ("V_wall_top", 0, vOutR, vz1 - vWall, vz1, C["pp"]),
    ("V_wall_well", vWellR, vWellR + vWall, vz0, vwft, C["pp"]),
    ("V_wall_wellfloor", 0, vWellR + vWall, vwft, vwft + vWall, C["pp"]),
    ("V_bottom", vWellR, vOutR, vz0, vz0 + vWall, C["pp"]),
    # устройство детектирования
    ("Electronics", 0, rRub, zTail, zPmtBot, C["el"]),
    ("PMT_vac", 0, rPmt - pmtGlass, zPmtBot, zWinBot, C["vac"]),
    ("PMT_glass", rPmt - pmtGlass, rPmt, zPmtBot, zWinBot, C["glass"]),
    ("AlCase_side", rRub, rCase, zTail, zAirTop, C["al"]),
    ("Rubber2_side", rCan, rRub, zPmtBot, zCanTop, C["rub"]),
    ("Window", 0, rMgo, zWinBot, -zCry, C["glass"]),
    ("MgO_side", rCry, rMgo, -zCry, zCry, C["mgo"]),
    ("AlCan_side", rMgo, rCan, zWinBot, zMgoTop, C["al"]),
    ("MgO_face", 0, rMgo, zCry, zMgoTop, C["mgo"]),
    ("AlCan_face", 0, rCan, zMgoTop, zCanTop, C["al"]),
    ("FaceAir", 0, rSeal, zCanTop, zAirTop, C["air"]),
    ("FaceSeal", rSeal, rRub, zCanTop, zAirTop, C["rub"]),
    ("AlCase_face", 0, rCase, zAirTop, zFace, C["al"]),
    # кристалл — последним, поверх
    ("NaI", 0, rCry, -zCry, zCry, C["nai"]),
]

fig, ax = plt.subplots(figsize=(6.2, 8.4))
for name, r0, r1, za, zb, col in parts:
    for sgn in (+1, -1):                       # зеркало
        x = sgn * r1 if sgn > 0 else sgn * r1
        x0 = min(sgn * r0, sgn * r1)
        ax.add_patch(Rectangle((x0, za), abs(r1 - r0), zb - za,
                                facecolor=col, edgecolor="#33333322",
                                linewidth=0.3, zorder=2))

# осевая линия
ax.axvline(0, color="#0008", lw=0.6, ls=(0, (6, 4)), zorder=3)

# аннотации
def note(txt, xy, xytext, ha="left"):
    ax.annotate(txt, xy=xy, xytext=xytext, ha=ha, va="center", fontsize=9,
                color="#111", zorder=5,
                arrowprops=dict(arrowstyle="-", color="#111", lw=0.7))

note("сосуд Маринелли 1 л\n(проба ОИСН-16)", (vOutR - vWall, 120), (200, 120))
note("кристалл NaI(Tl) Ø63×63", (0, 0), (200, 0))
note("отражатель MgO", (rMgo, 22), (200, 30))
note("оптический гель / ФЭУ", (0, zWinBot - 40), (200, -60))
note("наружный корпус Al\n(труба постоянного Ø)", (rCase, -120), (200, -125))
note("свинец 50 мм", (rPb - 20, 70), (200, 95))
note("сталь 3 мм", (rSt, 130), (200, 155))
note("свинец 50 мм\n(цилиндр вокруг ФЭУ)", (rPmtPb - 22, -170), (200, -165))
note("сталь 3 мм — внешний корпус\n(постоянный Ø)", (rSt, -215), (205, -218))
note("свинцовая заглушка снизу", (0, 0.5 * (zPlug0 + zPlug1)), (200, zPlug1))
note("стальная платформа 10 мм", (0, 0.5 * (zPlat0 + zPlat1)), (200, zPlat0))

# зазор дно кристалла — верх свинца дна (5 мм)
for zz in (-zCry, zCd0):
    ax.plot([-46, 46], [zz, zz], color="#c0392b", lw=0.7, ls=(0, (3, 2)),
            zorder=5)
ax.annotate("зазор 5 мм", xy=(-40, 0.5 * (-zCry + zCd0)), xytext=(-150, -22),
            ha="left", va="center", fontsize=9, color="#c0392b", zorder=7,
            arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.0))

ax.set_xlim(-172, 320)
ax.set_ylim(zPlat0 - 14, zSt1 + 16)
ax.set_aspect("equal")
ax.set_xlabel("радиус, мм")
ax.set_ylabel("ось Z, мм (0 — центр кристалла, +Z к торцу)")
ax.set_title("Разрез расчётной модели Гамма-1С в защите")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "article", "figures", "gamma1s_section_v4.png")
OUT = os.path.abspath(OUT)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print("сохранено:", OUT)
print("зазор дно кристалла — верх свинца дна: %.1f мм" % (-zCry - zCd0))
