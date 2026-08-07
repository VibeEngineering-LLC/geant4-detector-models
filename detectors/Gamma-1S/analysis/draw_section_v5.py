# -*- coding: utf-8 -*-
"""Разрез установки Гамма-1С — перечерчен 1:1 вокруг детектора (v5).

Метод (задание оператора 2026-08-02): детектор УДС-ГЦ-63х63 — единственный
объект с точными размерами (чертёж «Чертеж 63х63.pdf»), он взят линейкой; защита
масштабирована со сборочного разреза (reference/device/drawings/assembly_section.png)
через детектор (k≈0,66 px/мм, проверено И длиной 315 мм = 207 px, И диаметром
Ø78,3 = 52 px), радиальные толщины закрыты паспортной массой свинца.

ЧТО ИЗМЕНИЛОСЬ против v4 (и против модели G1SDetector.cc):
  - полость камеры Ø156 (Маринелли Ø150 + зазор+медь), НЕ Ø200. Прежние Ø200×190
    были подгонкой под массу свинца БЕЗ нижней секции — полость раздута.
  - защита СТУПЕНЧАТАЯ: широкая камера под Маринелли (свинец ~50 мм) переходит в
    узкий канал Ø~84 вокруг тела ФЭУ (свинец там толще, ~85 мм, наружка постоянна).
  - свинец продолжается ВНИЗ вокруг ФЭУ + свинцовая заглушка снизу; масса сходится
    с паспортом (свинец «не менее 165 кг») без раздувания полости.
  - крышка экрана-защиты (свинцовый диск) сверху; тележка — 4-лучевая с колёсами.
  - НЕТ коллиматора (нижний вид чертежа 63х63 с Ø134/90° — коллиматор, у нас его нет).

Оси: начало — центр кристалла, +Z вверх (к торцу/пробе), −Z вниз (ФЭУ, база).
Рисуется полный осесимметричный разрез (правая половина зеркалится).

    python analysis/draw_section_v5.py
"""
import os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

# --- ДЕТЕКТОР (мм), точный, зеркало G1SDetector.hh HeadGeom -------------------
cryDia, cryLen = 63.00, 63.00
mgoRad, mgoFace = 3.65, 6.00
alCan, rubber2, alCase = 0.50, 2.00, 1.50
alCaseFace, faceAir = 2.00, 1.00
faceSealW = 5.50
pmtDia, pmtGlass, pmtLen = 71.00, 1.50, 120.00
window, unitLen = 0.50, 315.00

rCry, zCry = 0.5*cryDia, 0.5*cryLen        # 31.5, 31.5
rMgo = rCry + mgoRad                        # 35.15
rCan = rMgo + alCan                         # 35.65
rRub = rCan + rubber2                       # 37.65
rCase = rRub + alCase                       # 39.15  (Ø78.3)
zMgoTop = zCry + mgoFace                     # 37.5
zCanTop = zMgoTop + alCan                    # 38.0
zAirTop = zCanTop + faceAir                  # 39.0
zFace = zAirTop + alCaseFace                 # 41.0
zWinBot = -zCry - window                     # -32.0
zPmtBot = zWinBot - pmtLen                   # -152.0
zTail = zFace - unitLen                      # -274.0
rPmt = 0.5*pmtDia                            # 35.5
rSeal = rRub - faceSealW                     # 32.15

# --- ЗАЩИТА (мм), вокруг детектора, k=0.66 с чертежа + масса свинца -----------
steel = 3.0
cu, cd = 1.0, 1.0
outerR = 132.0                 # Ø264 наружная сталь (чертёж/детектор)
pbOutR = outerR - steel        # 129 свинец снаружи
cavR = 78.0                    # Ø156 полость камеры (Маринелли Ø150 + зазор + медь)
boreR = 42.0                   # Ø84 канал вокруг тела ФЭУ (детектор Ø78.3)
pbCavInR = cavR + cu + cd      # 80  свинец от полости
pbBoreInR = boreR + cu + cd    # 44  свинец от канала ФЭУ
# вертикаль (z от центра кристалла)
zStep = -36.5                  # верх свинца-дна камеры; дно кристалла −31.5 → зазор 5 мм
zCavTop = 85.0                 # потолок камеры (над Маринелли, верх которой +79)
zLidTop = zCavTop + 50.0       # 135  верх крышки-защиты (свинец 50)
zLowBot = -175.0              # низ нижнего свинца вокруг ФЭУ (подобран под массу)
plugH = 50.0                   # свинцовая заглушка снизу
zPlug1, zPlug0 = zLowBot, zLowBot - plugH
platformH = 10.0
zPlat1 = zPlug0 - 6.0          # небольшой зазор до платформы
zPlat0 = zPlat1 - platformH
rPlat = 165.0                  # 4-лучевая тележка шире корпуса

# --- масса свинца (проверка против паспорта: «не менее 165 кг») ---------------
def vcyl(ri, ro, h):           # см³
    return math.pi*(ro*ro - ri*ri)*abs(h)/1000.0
m_side = vcyl(pbCavInR, pbOutR, zCavTop - zStep)*11.34/1000      # верхняя обечайка
m_lid  = vcyl(0, pbOutR, 50.0)*11.34/1000                         # крышка
m_low  = vcyl(pbBoreInR, pbOutR, zStep - zLowBot)*11.34/1000      # нижний свинец
m_plug = vcyl(0, pbOutR, plugH)*11.34/1000                        # заглушка
m_pb = m_side + m_lid + m_low + m_plug
print("свинец: обечайка %.1f + крышка %.1f + низ %.1f + заглушка %.1f = %.1f кг (паспорт >=165)"
      % (m_side, m_lid, m_low, m_plug, m_pb))

# --- Маринелли 1 л на головке ------------------------------------------------
vWall, vWellR, vWellD, vOutR, vH = 2.00, 40.00, 74.00, 75.00, 110.00
vwft = zFace + vWall            # 43 дно колодца (закрытый верх)
vz0 = vwft - vWellD            # -31 устье колодца — сосуд надет на головку
vz1 = vz0 + vH                 # 79 верх сосуда

C = dict(nai="#33e04d", mgo="#f2f2f2", al="#b3b3bf", rub="#26262b",
         glass="#9cccff", vac="#d6ebff", el="#7d4d1a", pb="#6a6a76",
         cd="#9a9a80", cu="#cc8033", steel="#8a929c", air="#ffffff",
         pp="#dcdcd6", samp="#c9a978")

# (метка, r0, r1, z0, z1, цвет) — правая половина, зеркалится
parts = [
    # --- ЗАЩИТА ---
    # наружная сталь: постоянный Ø сверху донизу (обечайка + верх)
    ("St_wall", pbOutR, outerR, zPlat1, zLidTop, C["steel"]),
    ("St_lidtop", 0, outerR, zLidTop, zLidTop + steel, C["steel"]),
    # крышка экрана-защиты (свинец) + стальная облицовка снизу крышки
    ("Lid_Pb", 0, pbOutR, zCavTop, zLidTop, C["pb"]),
    # верхняя камера: свинец, кадмий, медь (обечайка)
    ("Pb_up", pbCavInR, pbOutR, zStep, zCavTop, C["pb"]),
    ("Cd_up", cavR + cu, pbCavInR, zStep, zCavTop, C["cd"]),
    ("Cu_up", cavR, cavR + cu, zStep, zCavTop, C["cu"]),
    # ступень-дно камеры: свинец от канала до обечайки (верх = zStep)
    ("Pb_shelf", pbBoreInR, pbOutR, zStep - 12, zStep, C["pb"]),
    ("Cu_shelf", boreR, cavR, zStep - 2, zStep, C["cu"]),   # облицовка дна
    # нижняя секция: узкий канал Ø84, свинец толстый вокруг ФЭУ
    ("Pb_low", pbBoreInR, pbOutR, zLowBot, zStep - 12, C["pb"]),
    ("Cd_low", boreR + cu, pbBoreInR, zLowBot, zStep - 12, C["cd"]),
    ("Cu_low", boreR, boreR + cu, zLowBot, zStep - 12, C["cu"]),
    # свинцовая заглушка снизу + стальная платформа-тележка
    ("Plug_Pb", 0, pbOutR, zPlug0, zPlug1, C["pb"]),
    ("Platform", 0, rPlat, zPlat0, zPlat1, C["steel"]),
    # --- Маринелли 1 л на головке ---
    ("V_sample_side", vWellR + vWall, vOutR - vWall, vz0 + vWall, vz1 - vWall, C["samp"]),
    ("V_sample_top", 0, vWellR + vWall, vwft, vz1 - vWall, C["samp"]),
    ("V_wall_out", vOutR - vWall, vOutR, vz0, vz1, C["pp"]),
    ("V_wall_top", 0, vOutR, vz1 - vWall, vz1, C["pp"]),
    ("V_wall_well", vWellR, vWellR + vWall, vz0, vwft, C["pp"]),
    ("V_wall_wellfloor", 0, vWellR + vWall, vwft, vwft + vWall, C["pp"]),
    ("V_bottom", vWellR, vOutR, vz0, vz0 + vWall, C["pp"]),
    # --- ДЕТЕКТОР ---
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
    ("NaI", 0, rCry, -zCry, zCry, C["nai"]),
]

fig, ax = plt.subplots(figsize=(6.6, 9.2))
for name, r0, r1, za, zb, col in parts:
    for sgn in (+1, -1):
        x0 = min(sgn*r0, sgn*r1)
        ax.add_patch(Rectangle((x0, za), abs(r1 - r0), zb - za,
                               facecolor=col, edgecolor="#33333333",
                               linewidth=0.3, zorder=2))
# колёса тележки (2 в разрезе)
for xw in (-rPlat + 22, rPlat - 22):
    ax.add_patch(Circle((xw, zPlat0 - 22), 22, facecolor="#cfcfcf",
                        edgecolor="#333", lw=0.6, zorder=2))
    ax.add_patch(Circle((xw, zPlat0 - 22), 7, facecolor="#888",
                        edgecolor="#333", lw=0.5, zorder=3))

ax.axvline(0, color="#0008", lw=0.6, ls=(0, (6, 4)), zorder=3)

def note(txt, xy, xytext, col="#111", ha="left"):
    ax.annotate(txt, xy=xy, xytext=xytext, ha=ha, va="center", fontsize=8.5,
                color=col, zorder=6,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.7))

note("крышка экрана-защиты\n(свинец 50 мм)", (0, 0.5*(zCavTop + zLidTop)), (250, 120))
note("сосуд Маринелли 1 л\n(проба ОИСН-16)", (vOutR - vWall, 60), (250, 60))
note("кристалл NaI(Tl) Ø63×63", (0, 0), (250, 5))
note("отражатель MgO", (rMgo, 22), (250, 30))
note("свинец камеры ~50 мм", (0.5*(pbCavInR + pbOutR), 20), (250, -20))
note("сталь 3 мм (корпус, пост. Ø)", (outerR, 40), (250, -55))
note("узкий канал вокруг ФЭУ\nсвинец ~85 мм", (0.5*(pbBoreInR + pbOutR), -120),
     (250, -120))
note("оптический гель / ФЭУ", (0, zWinBot - 40), (250, -160))
note("свинцовая заглушка снизу", (0, 0.5*(zPlug0 + zPlug1)), (250, zPlug1))
note("тележка (сталь 10 мм) + колёса", (rPlat - 30, zPlat0), (250, zPlat0 - 10))

# зазор дно кристалла — верх свинца-дна (5 мм)
for zz in (-zCry, zStep):
    ax.plot([-52, 52], [zz, zz], color="#c0392b", lw=0.7, ls=(0, (3, 2)), zorder=5)
ax.annotate("зазор 5 мм", xy=(-46, 0.5*(-zCry + zStep)), xytext=(-190, -30),
            ha="left", va="center", fontsize=8.5, color="#c0392b", zorder=7,
            arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.0))

ax.set_xlim(-210, 380)
ax.set_ylim(zPlat0 - 52, zLidTop + 22)
ax.set_aspect("equal")
ax.set_xlabel("радиус, мм")
ax.set_ylabel("ось Z, мм (0 — центр кристалла, +Z к торцу)")
ax.set_title("Гамма-1С: разрез установки (перечерчен вокруг детектора, свинец %.0f кг)"
             % m_pb, fontsize=10)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "article", "figures", "gamma1s_section_v5.png")
OUT = os.path.abspath(OUT)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print("сохранено:", OUT)
print("полость Ø%.0f, канал Ø%.0f, наружка Ø%.0f, зазор дно-свинец %.1f мм"
      % (2*cavR, 2*boreR, 2*outerR, -zCry - zStep))
