# -*- coding: utf-8 -*-
"""Осевой разрез свинцовой защиты RC-103 с прибором в маринелли 200 мл.

Числа — те же формулы, что в PbShield.cc/RCDetector.cc::BuildVessel() (не
дублируются «на глаз»): слои защиты Cu->Cd->Pb от полости наружу, сосуд m200
из VesselGeom::Preset("m200"), прибор — DeviceGeom по умолчанию. Геометрия уже
прошла /geometry/test/run без пересечений и самопроверку объёма
(shieldrun geom) — рисунок иллюстрирует ПРОВЕРЕННУЮ модель, не черновик.

РАЗБОР ГЕОМЕТРИИ СОСУДА (по коду BuildVessel, не по внешнему виду детали):
    zSlot   — ОТКРЫТЫЙ торец горловины: сюда встаёт колодец, тут ЖЕ
              торцевая стенка (endWall) отделяет пробу от полости колодца.
    zRim    — противоположный конец сосуда, ЗАКРЫВАЕТСЯ отдельной печатной
              крышкой (capDisc + capSkirt, юбка внахлёст).
    Колодец (гильза cOut/cIn) — отдельная плаcтиковая труба ВНУТРИ пробы,
    стенка 1,25 мм, тянется от острия у носа прибора до горловины; проба
    заполняет кольцевой зазор между гильзой и стенкой сосуда, НЕ весь диск.
    Донце туннеля — колодец ЗАКРЫТ с торца (~1,29 мм сплошного пластика
    между наружным и внутренним острием капсулы), не открыт в пробу.
Прежняя версия рисунка это упускала — сосуд выглядел полым внутри, без
торцевой стенки, без гильзы колодца и без донца туннеля (замечания
оператора, второе — со ссылкой на чертёж изготовителя).

Упрощение, заявленное явно: прибор и колодец сосуда НЕ осесимметричны
(скруглённый прямоугольник), а разрез рисует тела как тела вращения — по
БОЛЬШЕЙ поперечной стороне (34 мм из 34×17,5 у прибора, wellOutX/wellInX у
колодца). Реальное поперечное сечение уже, см. подпись на рисунке.

Расположение в полости: сосуд с прибором поставлен на ДНО полости защиты
(зазор 2 мм) — так стоит реальная сборка под действием силы тяжести, а не
подвешена по центру. Это выбор ТОЛЬКО для рисунка (сдвиг всей сборки по z);
размер и центр самой полости (rCav/hzCav/zCav) в этом файле не меняются.

Запуск:  python plot_shield_section.py [pb] [cu] [cd] [out.png]
         python plot_shield_section.py 50 1.5 1.2 rc_shield_section.png
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

MAT = {
    "Pb": ("#7a7f87", "свинец"),
    "Cd": ("#d8c896", "кадмий"),
    "Cu": ("#c1743a", "медь"),
    "air_cav": ("#eef4fb", "полость (воздух)"),
    "vessel": ("#c9b18a", "сосуд, PLA"),
    "sample": ("#4a7c2f", "проба (черника)"),
    "device": ("#2b2b2b", "прибор RC-103"),
    "crystal": ("#ffd43b", "кристалл CsI(Tl)"),
}

# --- геометрия сосуда m200: числа из VesselGeom::Preset (RCDetector.cc) -----
CRYSTAL_Z0 = 12.00
OUTER_R, INNER_R, BARREL_H, END_WALL = 36.05, 33.24, 69.30, 2.80
WELL_TIP, WELL_TIP_OUT, SEAT_GAP = 47.81, 49.10, 0.32   # донце туннеля = разница
WELL_OUT_X, WELL_IN_X = 18.60, 17.35   # больший поперечник гильзы (упрощение)
CAP_T, CAP_H, CAP_SKIRT_R = 2.80, 12.20, 40.25
CASE_X, CASE_Z, CRYSTAL = 34.00, 123.00, 10.00


def build(pb, cu, cd, floor_gap=2.0):
    rows = []

    # --- прибор+сосуд: ИСТИННЫЕ мировые координаты (кристалл в z=0, как в
    # RCDetector.hh) — ось подписана «z от центра кристалла», и эти тела
    # двигать нельзя, иначе подпись врёт. Сдвигается ЗАЩИТА (см. ниже).
    zSlot = -CRYSTAL_Z0 + WELL_TIP - SEAT_GAP           # 35.49, горловина
    zRim = zSlot - BARREL_H                             # -33.81, дно сосуда
    zSmpTop = zSlot - END_WALL                          # 32.69, верх пробы
    zCapBot = zRim - CAP_T                              # -36.61, низ крышки
    zDeviceNose = -CRYSTAL_Z0                            # -12, нос прибора
    zDeviceTail = zDeviceNose + CASE_Z                   # 111, хвост прибора
    zWellTubeTip = zSlot - WELL_TIP_OUT                  # -13.61, дно гильзы (снаружи)
    zWellCavTip = zSlot - WELL_TIP                       # -12.32, дно полости (внутри)
    # донце туннеля — СПЛОШНОЙ пластик между наружным и внутренним острием
    # капсулы (~1,29 мм), закрывает колодец с торца. Пропущено в прежней
    # версии рисунка — по замечанию оператора со ссылкой на чертёж
    # изготовителя (marinelli_device_section.png): тоннель закрыт с торца.

    # --- защита: rCav/hzCav — умолчания ShieldGeom (PbShield.hh); zCav —
    # ПОДБИРАЕТСЯ так, чтобы сборка стояла на дне полости (зазор floor_gap),
    # а не наоборот. Замечание оператора: подпись оси «z от центра кристалла»
    # обязана оставаться верной, значит двигать можно только защиту.
    rCav, hzCav = 50.0, 90.0
    zCav = hzCav + zCapBot - floor_gap   # низ крышки сосуда - зазор = пол полости

    r, h = rCav, hzCav
    for d, tag in ((cu, "Cu"), (cd, "Cd"), (pb, "Pb")):
        if d <= 0:
            continue
        rOut, hOut = r + d, h + d
        rows.append(dict(mat=tag, ri=r, ro=rOut, z0=zCav - h, z1=zCav + h, order=0))
        rows.append(dict(mat=tag, ri=0, ro=rOut, z0=zCav - hOut, z1=zCav - h, order=0))
        rows.append(dict(mat=tag, ri=0, ro=rOut, z0=zCav + h, z1=zCav + hOut, order=0))
        r, h = rOut, hOut
    rOut, hOut = r, h
    rows.append(dict(mat="air_cav", ri=0, ro=rCav, z0=zCav - hzCav, z1=zCav + hzCav, order=0))

    # сосуд: боковая стенка (координаты НАСТОЯЩИЕ, без сдвига)
    # order: явный порядок отрисовки (больше = поверх). НЕ площадь тела —
    # эвристика «больше площадь снизу» рисовала тонкую, но ШИРОКУЮ торцевую
    # стенку (ri=0..36) ПОВЕРХ длинного, но УЗКОГО прибора (ri=0..17,
    # площадь на бумаге больше из-за длины) — прибор пропадал под стенкой.
    # Порядок нужен физический: защита -> сосуд -> проба/колодец -> прибор
    # -> кристалл, каждый следующий перекрывает предыдущий в его радиусе.
    rows.append(dict(mat="vessel", ri=INNER_R, ro=OUTER_R, z0=zRim, z1=zSmpTop, order=1))
    # торцевая стенка у горловины (ОТСУТСТВОВАЛА в прежней версии рисунка)
    rows.append(dict(mat="vessel", ri=0, ro=OUTER_R, z0=zSmpTop, z1=zSlot, order=1))
    # крышка: диск + юбка внахлёст
    rows.append(dict(mat="vessel", ri=0, ro=OUTER_R, z0=zCapBot, z1=zRim, order=1))
    skirtH = CAP_H - CAP_T
    if skirtH > 0.1:
        rows.append(dict(mat="vessel", ri=OUTER_R, ro=CAP_SKIRT_R,
                         z0=zRim, z1=zRim + skirtH, order=1))
    # проба: заполняет кольцо между стенкой сосуда и гильзой колодца
    rows.append(dict(mat="sample", ri=0, ro=INNER_R, z0=zRim, z1=zSmpTop, order=2))
    # гильза колодца (ОТСУТСТВОВАЛА): плаcтиковая труба вокруг прибора,
    # тянется от острия у носа до горловины
    rows.append(dict(mat="vessel", ri=WELL_IN_X, ro=WELL_OUT_X,
                     z0=zWellCavTip, z1=zSlot, order=3))
    # донце туннеля (ОТСУТСТВОВАЛО): закрывает колодец с торца
    rows.append(dict(mat="vessel", ri=0, ro=WELL_OUT_X,
                     z0=zWellTubeTip, z1=zWellCavTip, order=3))
    # прибор: по большей стороне сечения (34 мм), см. упрощение в docstring
    rows.append(dict(mat="device", ri=0, ro=CASE_X / 2,
                     z0=zDeviceNose, z1=zDeviceTail, order=4))
    rows.append(dict(mat="crystal", ri=0, ro=CRYSTAL / 2,
                     z0=-CRYSTAL / 2, z1=CRYSTAL / 2, order=5))

    return rows, (rOut, hOut), (zCav, hzCav, rCav)


def draw(pb, cu, cd, out):
    rows, (rOut, hOut), (zCav, hzCav, rCav) = build(pb, cu, cd)
    fig, ax = plt.subplots(figsize=(7.4, 9.4), dpi=170)

    rows = sorted(rows, key=lambda r: r["order"])
    seen = set()
    for r in rows:
        col, lab = MAT[r["mat"]]
        w = r["ro"] - r["ri"]
        for sign in (+1, -1):
            x0 = r["ri"] if sign > 0 else -r["ro"]
            first = lab not in seen
            ax.add_patch(Rectangle((x0, r["z0"]), w, r["z1"] - r["z0"],
                                   facecolor=col, edgecolor="#212529",
                                   linewidth=0.3, label=lab if first else None))
            seen.add(lab)

    ax.axvline(0, color="#495057", lw=0.6, ls=(0, (12, 4, 2, 4)))
    pad = 15
    ax.set_xlim(-(rOut + pad), rOut + pad)
    ax.set_ylim(zCav - hOut - pad, zCav + hOut + pad)
    ax.set_aspect("equal")
    ax.set_xlabel("радиус, мм")
    ax.set_ylabel("z от центра кристалла, мм")
    ax.set_title("RadiaCode-103 в маринелли m200, свинцовая защита\n"
                 "Pb %.0f + Cd %.1f + Cu %.1f мм, полость r%.0f×h%.0f, "
                 "сосуд на дне (осевой разрез)" % (pb, cd, cu, rCav, hzCav),
                 fontsize=10.2)
    ax.grid(alpha=0.2, lw=0.35)

    ann = [
        (-(rCav + (pb + cu + cd) / 2), zCav + hzCav + (pb + cu + cd) / 2, "стенка Pb+Cd+Cu"),
        (-rCav * 0.5, zCav + hzCav * 0.65, "полость (воздух)"),
        (-36, zCav - hzCav + 20, "крышка сосуда"),
        (-15, zCav - hzCav + 45, "торцевая стенка\n+ гильза колодца"),
        (-8, zCav - hzCav + 60, "проба (черника)"),
        (-17, zCav - hzCav + 100, "прибор RC-103\n(корпус 34×17,5×123)"),
        (5, zCav - hzCav + 25, "кристалл CsI(Tl) 10³"),
    ]
    ax_left = -(rOut + pad)
    for i, (x, y, t) in enumerate(ann):
        ax.annotate(t, xy=(x, y), xytext=(ax_left + 4, zCav + hOut - 14 - 16 * i),
                    fontsize=6.8, color="#212529", ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#868e96",
                                    shrinkA=1, shrinkB=1,
                                    connectionstyle="angle,angleA=0,angleB=90,rad=3"))

    h_, l_ = ax.get_legend_handles_labels()
    ax.legend(h_, l_, fontsize=6.8, loc="upper right", framealpha=0.95)
    fig.text(0.5, 0.010,
             "Упрощение: прибор и колодец сосуда — скруглённый прямоугольник, "
             "на рисунке (тело вращения) взята бо́льшая сторона.\n"
             "Сосуд поставлен на дно полости (зазор 2 мм) — выбор рисунка, "
             "не физической модели (zCav полости не менялся).\n"
             "Числа — из PbShield.cc/RCDetector.cc, геометрия проверена "
             "shieldrun geom (/geometry/test/run без пересечений).",
             ha="center", fontsize=6.4, color="#343a40")
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(out)
    print("записано", out)


if __name__ == "__main__":
    pb = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
    cu = float(sys.argv[2]) if len(sys.argv) > 2 else 1.5
    cd = float(sys.argv[3]) if len(sys.argv) > 3 else 1.2
    out = sys.argv[4] if len(sys.argv) > 4 else "rc_shield_section.png"
    draw(pb, cu, cd, out)
