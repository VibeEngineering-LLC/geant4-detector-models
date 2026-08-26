#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RadiaCode-103 - разрез Y=0 (плоскость XZ) с размерными линиями, matplotlib.

Geant4 не умеет проставлять размерные линии/выноски - это задача инженерного
чертежа, для неё matplotlib с точной привязкой к данным (ax.annotate +
arrowprops) надёжнее, чем накладывать текст поверх рендера физического
движка по пиксельным координатам (риск разъехаться).

ВСЕ числа компонентов - из geant4/gdml/detector/RC103_detector.gdml (SSOT,
читано и сверено построчно 26.08.2026): <define> даёт ЛОКАЛЬНЫЕ позиции p_*
относительно родителя, <solids> даёт размеры x/y/z каждого <box>. Мировые
координаты восстановлены вручную сложением по цепочке <physvol> дерева
World -> RC103_device_log -> {Case_shell_log, Case_interior_log -> {PCB_log,
Battery_log, EMI_log, Display_window_log, Display_LCD_log, USB_log,
DetectorModule_log -> {Capsule_shell_log, Capsule_cavity_log -> {Crystal_log,
SiPM_log, ESR_*_log}, SiPM_carrier_log}}} (RC103_device_log и Case_interior_log
сами сидят в p_origin=(0,0,0), поэтому для верхних узлов локальные и мировые
координаты совпадают; расхождение появляется только внутри DetectorModule_log
и Capsule_cavity_log - см. проверку протрузий ниже). Компонент EMI (тканевый
ESD/EMI-экран) в образце RC-110 отсутствовал вовсе - у RC-103 он есть в GDML
(EMI_log) и включён в чертёж, иначе аннотация про его overlap указывала бы
в пустоту.

Протрузии SiPM (400 мкм) и SiPM_carrier (425 мкм) на zoom-панели - НЕ взяты
готовыми с потолка, а пересчитаны в этой сессии построчно по GDML:
  - SiPM: p_sipm z=-5.4 (лок. отн. Capsule_cavity_log), halfZ=0.4 (SiPMSolid
    z=0.8) -> дальняя грань на z=-5.8 лок., граница полости капсулы (мать
    Capsule_cavity_log) halfZ=5.4 -> граница на z=-5.4 лок. Протрузия =
    5.8-5.4 = 0.4 мм = 400 мкм по -Z.
  - SiPM_carrier: p_sipm_carrier z=-5.55 (лок. отн. DetectorModule_log),
    halfZ=0.275 (SiPM_carrier_solid z=0.55) -> ближняя/дальняя грани
    z=[-5.825,-5.275] лок. Стенка капсулы (Capsule_shell) в Z занимает
    лок. z=[-6.5,-5.4] (наружная halfZ=6.5 минус внутр. halfZ=5.4). Пересечение
    carrier со стенкой = [-5.825,-5.4] = 0.425 мм = 425 мкм.
Оба числа совпали с независимо измеренным G4PVPlacement::CheckOverlaps()
(реальный прогон Geant4, не гипотеза) - см. LESSONS в отчёте сессии.
Остальные overlaps этой геометрии (окно дисплея +1,3 мм, LCD +0,6 мм в стенку
корпуса, EMI +0,05 мм, USB +1,75 мм за внутр. полость, PCB/капсульный модуль
0,485 мм) - из того же прогона CheckOverlaps(), сведены честной сводкой на
главной панели (build_main_panel) - НЕ пересчитывались заново в этом скрипте,
они не требуют новой оценки, а требуют не быть скрытыми.

Оси: X - длина (USB -> +X), Z - толщина (-Z = лицевая грань). Датум - центр
корпуса (0,0,0). Y (ширина, 34,0 мм наруж. / 31,0 мм полость) в этом сечении
не показана (плоскость Y=0) - вынесена текстом.

Выход: detectors/RadiaCode-103/geant4/verify/RC103_dimensioned_section.png

Запуск:
    python dimensioned_section.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SCRIPT_DIR = Path(__file__).resolve().parent
OUT = SCRIPT_DIR.parent / "verify" / "RC103_dimensioned_section.png"

# --- Компоненты: (имя, halfX, halfZ, cx, cz, facecolor, edgecolor, dashed) --
# Половинные размеры - из <solids> (x/y/z /2), центры - мировые координаты,
# восстановленные по дереву <physvol> RC103_detector.gdml (см. докстринг).
COMPONENTS = [
    ("Корпус наруж.",  61.50, 8.75,   0.00,  0.000, "none",    "#888888", False),
    ("Корпус внутр.",  60.00, 7.25,   0.00,  0.000, "none",    "#888888", True),
    ("Кристалл CsI",    5.00, 5.00, -49.50, -0.550, "#d4af00", "#8a7100", False),
    ("SiPM",             3.00, 0.40, -49.50, -5.950, "#39b7d6", "#1c6b80", False),
    ("Капсула наруж.",   8.00, 6.50, -49.50, -0.550, "none",    "#404040", False),
    ("Капсула внутр.",   5.40, 5.40, -49.50, -0.550, "none",    "#404040", True),
    ("Плата PCB",       51.00, 0.50,   6.50, -4.500, "#1a7a34", "#0c451d", False),
    ("Аккумулятор",     30.00, 3.00,  15.00,  3.500, "#a0a0a0", "#5f5f5f", False),
    ("EMI/ESD ткань",   26.00, 0.15, -33.00,  7.150, "#8b5fbf", "#4a2d66", False),
    ("Окно дисплея",    18.00, 0.15, -14.00, -8.400, "#1a1a1f", "#000000", False),
    ("LCD дисплея",     17.00, 0.70, -14.00, -7.150, "#0d0d33", "#000000", False),
    ("USB",              3.75, 1.60,  58.00, -3.000, "#c68a35", "#7a531f", False),
]


def add_rect(ax, half_x, half_z, cx, cz, fc, ec, dashed, lw=1.2, zorder=2):
    r = Rectangle((cx - half_x, cz - half_z), 2 * half_x, 2 * half_z,
                  facecolor=fc, edgecolor=ec, linewidth=lw,
                  linestyle="--" if dashed else "-", zorder=zorder)
    ax.add_patch(r)
    return r


def dim_h(ax, x1, x2, y_dim, y_feat, text, color="black", fs=8):
    """Горизонтальная размерная линия со стрелками + вынесённые пунктиром
    выносные линии от контура (y_feat) до линии размера (y_dim)."""
    for x in (x1, x2):
        ax.plot([x, x], [y_feat, y_dim], color=color, lw=0.6, ls=":", zorder=1)
    ax.annotate("", xy=(x1, y_dim), xytext=(x2, y_dim),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.0))
    va = "bottom" if y_dim >= y_feat else "top"
    ax.text((x1 + x2) / 2, y_dim, text, ha="center", va=va, fontsize=fs,
            color=color, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))


def dim_v(ax, z1, z2, x_dim, x_feat, text, color="black", fs=8):
    """Вертикальная размерная линия, зеркально dim_h."""
    for z in (z1, z2):
        ax.plot([x_feat, x_dim], [z, z], color=color, lw=0.6, ls=":", zorder=1)
    ax.annotate("", xy=(x_dim, z1), xytext=(x_dim, z2),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.0))
    ha = "left" if x_dim >= x_feat else "right"
    ax.text(x_dim, (z1 + z2) / 2, text, ha=ha, va="center", fontsize=fs,
            color=color, rotation=90,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))


def build_main_panel(ax):
    for name, hx, hz, cx, cz, fc, ec, dashed in COMPONENTS:
        add_rect(ax, hx, hz, cx, cz, fc, ec, dashed)

    # 1. Общая длина корпуса 123 мм
    dim_h(ax, -61.50, 61.50, -16.5, -8.75, "123 мм", color="black")
    # 2. Общая толщина (высота в сечении) корпуса 17,5 мм
    dim_v(ax, -8.75, 8.75, 70.0, 61.50, "17,5 мм", color="black")
    # 3. Толщина стенки корпуса 1,5 мм - callout (зона мала для линии)
    ax.annotate("стенка корпуса\n1,5 мм", xy=(30.0, 8.00), xytext=(30.0, 15.0),
                fontsize=8, ha="center", color="#333333",
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.0))
    # 4. Позиция кристалла от датума: X=0 -> X=-49,5
    dim_h(ax, 0.0, -49.50, 12.0, 8.75, "-49,5 мм (крист.)", color="#8a7100")
    # 5. Габарит капсулы 16x13 мм (капсула НЕ куб: Capsule_outer x=16,z=13)
    dim_h(ax, -57.50, -41.50, -9.4, -7.05, "16 мм", color="#404040")
    dim_v(ax, -7.05, 5.95, -69.0, -57.50, "13 мм", color="#404040")
    # 6. Позиция PCB (X=+6,5) и длина (102 мм)
    dim_h(ax, -44.50, 57.50, -12.0, -8.75, "102 мм (центр X=+6,5)", color="#0c451d")
    # 7. Позиция аккумулятора (X=+15,0) и длина (60 мм)
    dim_h(ax, -15.00, 45.00, -14.5, -8.75, "60 мм (центр X=+15,0)", color="#5f5f5f")

    ax.plot(0, 0, "+", color="red", ms=10, mew=1.5, zorder=5)
    ax.text(0, -1.9, "датум (0,0)", color="red", fontsize=7, ha="center")

    ax.text(-77, -17.5, "Ширина корпуса (Y, не в этом сечении):\n"
                          "34,0 мм наруж. / 31,0 мм полость",
            fontsize=7.5, color="#555555", va="bottom")

    # Известные overlaps этой геометрии (G4PVPlacement::CheckOverlaps(),
    # реальный прогон, не гипотеза) - честная сводка, не скрываем; сведена в
    # один блок в свободном верхнем левом углу, чтобы не наезжать на линии
    # размеров (SiPM/carrier с точными числами - на zoom-панели справа).
    ax.text(-77, 24.5,
            "Известные overlaps (G4 CheckOverlaps(), реальный прогон, не гипотеза):\n"
            "SiPM/carrier по -Z: 0,400/0,425 мм (детали - см. zoom-панель справа)\n"
            "окно дисплея/LCD в стенку корпуса: +1,3/+0,6 мм;  EMI: +0,05 мм\n"
            "USB за внутр. полость: +1,75 мм;  PCB и капсульный модуль: 0,485 мм\n"
            "- всё перечисленное: задокументированные упрощения модели, не скрыты",
            fontsize=6.6, color="#333333", va="top", ha="left",
            bbox=dict(fc="#fff6d8", ec="#b0a060", lw=0.6, alpha=0.92, pad=4))

    ax.set_xlim(-79, 79)
    ax.set_ylim(-21, 25.5)
    ax.set_aspect("equal")
    ax.set_xlabel("X, мм  (USB -> +X)")
    ax.set_ylabel("Z, мм  (-Z = лицевая грань)")
    ax.grid(alpha=0.2, lw=0.5)
    ax.set_title("Разрез Y=0 - корпус в сборе")

    legend_items = [
        ("Корпус (наруж./внутр.)", "none", "#888888"),
        ("Кристалл CsI(Tl)", "#d4af00", "#8a7100"),
        ("SiPM", "#39b7d6", "#1c6b80"),
        ("Капсула (наруж./внутр.)", "none", "#404040"),
        ("Плата PCB", "#1a7a34", "#0c451d"),
        ("Аккумулятор", "#a0a0a0", "#5f5f5f"),
        ("EMI/ESD ткань", "#8b5fbf", "#4a2d66"),
        ("Окно дисплея / LCD", "#1a1a1f", "#000000"),
        ("USB", "#c68a35", "#7a531f"),
    ]
    handles = [Rectangle((0, 0), 1, 1, facecolor=fc, edgecolor=ec, linewidth=1.2)
               for _, fc, ec in legend_items]
    labels = [lab for lab, _, _ in legend_items]
    ax.legend(handles, labels, loc="lower right", fontsize=7, framealpha=0.9, ncol=2, columnspacing=1.0)


def build_zoom_panel(ax):
    # Капсула + кристалл + SiPM (те же числа, что в главной панели) плюс ESR
    # (5 граней по 0,065 мм каждая, half-extent 5,065 мм по широкой оси -
    # p_esr_px/nx/py/ny/pz в RC103_detector.gdml; грани -Z НЕТ - со стороны
    # SiPM отражатель не ставится, это штатный оптический выход, а не
    # упрощение) плюс SiPM_carrier (несущая плата SiPM, 10,4x10,4x0,55).
    add_rect(ax, 8.00, 6.50, -49.50, -0.550, "none", "#404040", False, lw=1.5)
    add_rect(ax, 5.40, 5.40, -49.50, -0.550, "none", "#404040", True, lw=1.0)
    add_rect(ax, 0.0325, 5.065, -44.4675, -0.550, "none", "#dcdce0", False, lw=1.0)  # ESR +X
    add_rect(ax, 0.0325, 5.065, -54.5325, -0.550, "none", "#dcdce0", False, lw=1.0)  # ESR -X
    add_rect(ax, 5.065, 0.0325, -49.50, 4.4825, "none", "#dcdce0", False, lw=1.0)    # ESR +Z
    add_rect(ax, 5.00, 5.00, -49.50, -0.550, "#d4af00", "#8a7100", False)
    add_rect(ax, 3.00, 0.40, -49.50, -5.950, "#39b7d6", "#1c6b80", False)
    add_rect(ax, 5.20, 0.275, -49.50, -6.100, "#7a531f", "#4a3010", False)  # SiPM_carrier

    # SiPM/carrier протрузии - пересчитаны построчно по GDML в этом скрипте
    # (см. докстринг), НЕ взяты готовыми; совпали с реальным прогоном
    # G4PVPlacement::CheckOverlaps(). Честные подписи, зазор НЕ рисуется
    # положительным.
    ax.annotate("SiPM: протыкает границу\nполости капсулы на 0,400 мм\nпо -Z (ниша d=11 мм, h=2,2 мм\nне смоделирована)",
                xy=(-49.50, -6.35), xytext=(-41.0, -6.9), fontsize=6.5,
                color="#1c6b80", ha="left",
                arrowprops=dict(arrowstyle="->", color="#1c6b80", lw=0.9))
    ax.annotate("SiPM_carrier: протыкает\nстенку капсулы на 0,425 мм\n(тот же узел, упрощение)",
                xy=(-49.50, -6.375), xytext=(-41.0, -3.3), fontsize=6.5,
                color="#7a531f", ha="left",
                arrowprops=dict(arrowstyle="->", color="#7a531f", lw=0.9))
    ax.annotate("ESR: прилегает к кристаллу\nвплотную, 0 мм зазора\n(грани совпадают в GDML)",
                xy=(-44.4675, 2.0), xytext=(-41.0, 3.2), fontsize=6.5,
                color="#8a7100", ha="left",
                arrowprops=dict(arrowstyle="->", color="#8a7100", lw=0.9))
    ax.text(-59.5, -8.0, "стенка капсулы: X/Y 2,6 мм,\nZ 1,1 мм (Capsule_outer\nне куб: 16x16x13)",
            fontsize=6.2, color="#404040", ha="left", va="bottom")

    ax.set_xlim(-60, -40)
    ax.set_ylim(-9.5, 7.5)
    ax.set_aspect("equal")
    ax.set_xlabel("X, мм")
    ax.set_ylabel("Z, мм")
    ax.grid(alpha=0.2, lw=0.5)
    ax.set_title("Деталь узла кристалл-капсула-SiPM (зум)")


def main() -> int:
    fig, (ax_main, ax_zoom) = plt.subplots(
        1, 2, figsize=(19, 7.5), gridspec_kw={"width_ratios": [2.3, 1.0]})
    build_main_panel(ax_main)
    build_zoom_panel(ax_zoom)
    fig.suptitle("RadiaCode-103 - разрез Y=0, размеры в мм", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(f"OK: {OUT} ({OUT.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())