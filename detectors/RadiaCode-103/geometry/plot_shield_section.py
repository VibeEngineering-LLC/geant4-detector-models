# -*- coding: utf-8 -*-
"""Осевой разрез свинцовой защиты RC-103 с прибором в маринелли 200 мл.

Домик КОРОБЧАТЫЙ: полость 150x150 в плане, 385 мм высотой, верх открыт,
сборка стоит на дне. Прежняя версия этого файла рисовала цилиндр r=50, h=180 —
геометрию, которой больше нет; рисунок при этом выглядел как иллюстрация
проверенной модели.

Числа — те же формулы, что в PbShield.cc/RCDetector.cc::BuildVessel(), не
переписанные на глаз: слои защиты Cu->Cd->Pb от полости наружу (вверх габарит
растёт ТОЛЬКО при крышке), сосуд m200 из VesselGeom::Preset("m200"), прибор —
DeviceGeom по умолчанию.

РАЗБОР ГЕОМЕТРИИ СОСУДА (по коду BuildVessel, не по внешнему виду детали):
    zSlot   — ОТКРЫТЫЙ торец горловины: сюда встаёт колодец, тут ЖЕ
              торцевая стенка (endWall) отделяет пробу от полости колодца.
    zRim    — противоположный конец сосуда, ЗАКРЫВАЕТСЯ отдельной печатной
              крышкой (capDisc + capSkirt, юбка внахлёст).
    Колодец (гильза cOut/cIn) — отдельная пластиковая труба ВНУТРИ пробы,
    стенка 1,25 мм, тянется от острия у носа прибора до горловины; проба
    заполняет кольцевой зазор между гильзой и стенкой сосуда, НЕ весь диск.
    Донце туннеля — колодец ЗАКРЫТ с торца (~1,29 мм сплошного пластика
    между наружным и внутренним острием капсулы), не открыт в пробу.

Упрощение, заявленное явно: прибор и колодец сосуда НЕ осесимметричны
(скруглённый прямоугольник), а разрез рисует их по БОЛЬШЕЙ поперечной стороне
(34 мм из 34x17,5 у прибора, wellOutX/wellInX у колодца). Реальное поперечное
сечение уже, см. подпись на рисунке.

Посадка: zCav подбирается так, чтобы низ сборки лёг на дно полости — ровно как
RCShieldDetector::PlannedZCav() при seatOnFloor, без произвольного зазора.
Подпись оси «z от центра кристалла» обязана оставаться верной, поэтому двигать
можно только защиту.

Запуск:  python plot_shield_section.py [pb] [cu] [cd] [out.png] [lid|nolid]
         python plot_shield_section.py 50 0 0 rc_shield_section.png nolid
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
    "sample": ("#eaf6ff", "внутри сосуда воздух"),
    "device": ("#2b2b2b", "прибор RC-103"),
    "crystal": ("#ffd43b", "кристалл CsI(Tl)"),
}

# --- геометрия сосуда m200: числа из VesselGeom::Preset (RCDetector.cc) -----
CRYSTAL_Z0 = 12.00
OUTER_R, INNER_R, BARREL_H, END_WALL = 36.05, 33.24, 69.30, 2.80
WELL_TIP, WELL_TIP_OUT, SEAT_GAP = 47.81, 49.10, 0.32
WELL_OUT_X, WELL_IN_X = 18.60, 17.35
CAP_T, CAP_H, CAP_SKIRT_R = 2.80, 12.20, 40.25
CASE_X, CASE_Z, CRYSTAL = 34.00, 123.00, 10.00

# --- полость домика: умолчания ShieldGeom (PbShield.hh) ---------------------
HX_CAV, HY_CAV, HZ_CAV = 75.0, 75.0, 192.5
def shield_rows(pb, cu, cd, hx_cav, hy_cav, hz_cav, z_cav, with_lid):
    """
    Прямоугольники осевого разреза коробчатой защиты; вверх габарит растёт только при крышке;
    стенки по z занимают ровно высоту полости своего слоя, иначе перекрываются с крышкой.
    """
    rows = []
    hx = hx_cav
    z_lo = z_cav - hz_cav
    z_hi = z_cav + hz_cav

    for d, mat in [(cu, "Cu"), (cd, "Cd"), (pb, "Pb")]:
        if d <= 0:
            continue
        hx_out = hx + d
        z_lo_out = z_lo - d
        z_hi_out = z_hi + d if with_lid else z_hi

        # Правая стенка
        rows.append(dict(mat=mat, x0=hx, x1=hx_out, z0=z_lo, z1=z_hi, order=0))
        # Левая стенка
        rows.append(dict(mat=mat, x0=-hx_out, x1=-hx, z0=z_lo, z1=z_hi, order=0))
        # Дно
        rows.append(dict(mat=mat, x0=-hx_out, x1=hx_out, z0=z_lo_out, z1=z_lo, order=0))
        # Крышка
        if with_lid:
            rows.append(dict(mat=mat, x0=-hx_out, x1=hx_out, z0=z_hi, z1=z_hi_out, order=0))

        hx = hx_out
        z_lo = z_lo_out
        z_hi = z_hi_out

    return (rows, hx, z_lo, z_hi)

def ring(mat, ri, ro, z0, z1, order):
    """Тело вращения ri..ro в разрезе — два прямоугольника (или один при ri=0)."""
    if ri <= 0:
        return [dict(mat=mat, x0=-ro, x1=ro, z0=z0, z1=z1, order=order)]
    return [dict(mat=mat, x0=ri, x1=ro, z0=z0, z1=z1, order=order),
            dict(mat=mat, x0=-ro, x1=-ri, z0=z0, z1=z1, order=order)]


def build(pb, cu, cd, with_lid):
    rows = []

    # Прибор+сосуд в ИСТИННЫХ мировых координатах (кристалл в z=0, как в
    # RCDetector.hh) — эти тела двигать нельзя, иначе подпись оси врёт.
    zSlot = -CRYSTAL_Z0 + WELL_TIP - SEAT_GAP
    zRim = zSlot - BARREL_H
    zSmpTop = zSlot - END_WALL
    zCapBot = zRim - CAP_T
    zDeviceNose = -CRYSTAL_Z0
    zDeviceTail = zDeviceNose + CASE_Z
    zWellTubeTip = zSlot - WELL_TIP_OUT
    zWellCavTip = zSlot - WELL_TIP

    # Посадка на дно полости — та же арифметика, что PlannedZCav().
    z_cav = zCapBot + HZ_CAV

    shield, hx_out, z_lo_out, z_hi_out = shield_rows(
        pb, cu, cd, HX_CAV, HY_CAV, HZ_CAV, z_cav, with_lid)
    rows += shield
    rows.append(dict(mat="air_cav", x0=-HX_CAV, x1=HX_CAV,
                     z0=z_cav - HZ_CAV, z1=z_cav + HZ_CAV, order=0))

    # order: явный порядок отрисовки (больше = поверх). НЕ площадь тела —
    # эвристика «больше площадь снизу» рисовала тонкую, но ШИРОКУЮ торцевую
    # стенку поверх длинного, но УЗКОГО прибора, и прибор пропадал.
    rows += ring("vessel", INNER_R, OUTER_R, zRim, zSmpTop, 1)
    rows += ring("vessel", 0, OUTER_R, zSmpTop, zSlot, 1)
    rows += ring("vessel", 0, OUTER_R, zCapBot, zRim, 1)
    skirtH = CAP_H - CAP_T
    if skirtH > 0.1:
        rows += ring("vessel", OUTER_R, CAP_SKIRT_R, zRim, zRim + skirtH, 1)
    # Пробы сейчас нет — сосуд ПУСТОЙ, внутри воздух (постановка задачи).
    rows += ring("sample", 0, INNER_R, zRim, zSmpTop, 2)
    rows += ring("vessel", WELL_IN_X, WELL_OUT_X, zWellCavTip, zSlot, 3)
    rows += ring("vessel", 0, WELL_OUT_X, zWellTubeTip, zWellCavTip, 3)
    rows += ring("device", 0, CASE_X / 2, zDeviceNose, zDeviceTail, 4)
    rows += ring("crystal", 0, CRYSTAL / 2, -CRYSTAL / 2, CRYSTAL / 2, 5)

    return rows, (hx_out, z_lo_out, z_hi_out), z_cav


def draw(pb, cu, cd, out, with_lid):
    rows, (hx_out, z_lo_out, z_hi_out), z_cav = build(pb, cu, cd, with_lid)
    fig, ax = plt.subplots(figsize=(7.0, 9.8), dpi=170)

    seen = set()
    for r in sorted(rows, key=lambda q: q["order"]):
        col, lab = MAT[r["mat"]]
        first = lab not in seen
        ax.add_patch(Rectangle((r["x0"], r["z0"]), r["x1"] - r["x0"],
                               r["z1"] - r["z0"], facecolor=col,
                               edgecolor="#212529", linewidth=0.3,
                               label=lab if first else None))
        seen.add(lab)

    ax.axvline(0, color="#495057", lw=0.6, ls=(0, (12, 4, 2, 4)))
    pad = 18
    ax.set_xlim(-(hx_out + pad), hx_out + pad)
    ax.set_ylim(z_lo_out - pad, z_hi_out + pad)
    ax.set_aspect("equal")
    ax.set_xlabel("x, мм")
    ax.set_ylabel("z от центра кристалла, мм")
    ax.set_title("RadiaCode-103 в маринелли m200 (пустой), свинцовый домик\n"
                 "Pb %.0f + Cd %.1f + Cu %.1f мм, полость %.0fx%.0fx%.0f мм, "
                 "верх %s"
                 % (pb, cd, cu, 2 * HX_CAV, 2 * HY_CAV, 2 * HZ_CAV,
                    "закрыт" if with_lid else "ОТКРЫТ"),
                 fontsize=10.0)
    ax.grid(alpha=0.2, lw=0.35)

    h_, l_ = ax.get_legend_handles_labels()
    ax.legend(h_, l_, fontsize=6.8, loc="upper right", framealpha=0.95)
    fig.text(0.5, 0.010,
             "Упрощение: прибор и колодец сосуда — скруглённый прямоугольник, "
             "в разрезе взята бо́льшая сторона.\n"
             "Сборка стоит на дне полости — так же, как в расчёте "
             "(RCShieldDetector::PlannedZCav при seatOnFloor).\n"
             "Числа — из PbShield.cc/RCDetector.cc, геометрия проверена "
             "shieldrun geom (/geometry/test/run без пересечений).",
             ha="center", fontsize=6.4, color="#343a40")
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(out)
    print("записано", out)
    print("наружный габарит: %.0f x %.0f x %.0f мм, z от %.2f до %.2f, zCav=%.2f"
          % (2 * hx_out, 2 * hx_out, z_hi_out - z_lo_out, z_lo_out, z_hi_out,
             z_cav))


if __name__ == "__main__":
    pb = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
    cu = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    cd = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    out = sys.argv[4] if len(sys.argv) > 4 else "rc_shield_section.png"
    with_lid = "lid" in sys.argv[5:] and "nolid" not in sys.argv[5:]
    draw(pb, cu, cd, out, with_lid)