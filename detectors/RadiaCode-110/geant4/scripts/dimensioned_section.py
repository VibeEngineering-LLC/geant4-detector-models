#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RadiaCode-110 - разрез Y=0 (плоскость XZ) с размерными линиями, matplotlib.

Geant4 не умеет проставлять размерные линии/выноски - это задача инженерного
чертежа, для неё matplotlib с точной привязкой к данным (ax.annotate +
arrowprops) надёжнее, чем накладывать текст поверх рендера физического
движка по пиксельным координатам (риск разъехаться).

ВСЕ числа компонентов - из detectors/RadiaCode-110/geant4/geometry/
RC110Detector.cc (SSOT, читано и сверено построчно 26.08.2026, см. функцию
BuildDevice()). Три числа зазоров во втором subplot (SiPM/ESR/капсула) - из
SESSION-STATE.md контура, раздел "Зазоры при Zc=-0,325" (не в RC110Detector.cc
- это исторический расчёт допуска, актуальная геометрия уже флеш-упрощена,
см. GEANT4-MODEL.md "Ограничения модели"). Это ДОПОЛНЕНИЕ к
verify/RC110_align_check.png и к rc110_view_section.png (настоящий Geant4-
рендер разреза) - не замена.

Оси: X - длина (USB -> +X), Z - толщина (-Z = лицевая грань, значок
радиации), см. docs/GEANT4-MODEL.md "Система координат". Датум - центр
корпуса (0,0,0). Y (ширина, 34.1 мм наруж. / 31.1 мм полость) в этом сечении
не показана (плоскость Y=0) - вынесена текстом.

Выход: detectors/RadiaCode-110/geant4/verify/RC110_dimensioned_section.png
(НЕ трогает verify/RC110_align_check.png - тот отдельный, настоящий рендер).

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
OUT = SCRIPT_DIR.parent / "verify" / "RC110_dimensioned_section.png"

# --- Компоненты: (имя, halfX, halfZ, cx, cz, facecolor, edgecolor, dashed) --
# Половинные размеры и центры - дословно из RC110Detector.cc::BuildDevice().
COMPONENTS = [
    ("Корпус наруж.",  63.30, 10.85,   0.00,  0.000, "none",    "#888888", False),
    ("Корпус внутр.",  61.80,  9.35,   0.00,  0.000, "none",    "#888888", True),
    ("Кристалл CsI",    7.00,  7.00, -50.25,  0.000, "#d4af00", "#8a7100", False),
    ("SiPM",             3.00,  0.40, -50.25, -7.400, "#39b7d6", "#1c6b80", False),
    ("Капсула наруж.",   9.00,  9.00, -49.90, -0.325, "none",    "#404040", False),
    ("Капсула внутр.",   7.50,  7.50, -49.90, -0.325, "none",    "#404040", True),
    ("Плата PCB",       50.75,  0.50,   9.85, -5.450, "#1a7a34", "#0c451d", False),
    ("Аккумулятор",     32.50,  4.80,  23.80,  4.550, "#a0a0a0", "#5f5f5f", False),
    ("Окно дисплея",    18.25,  0.50, -19.25, -8.950, "#1a1a1f", "#000000", False),
    ("LCD дисплея",     17.00,  1.10, -19.25, -7.350, "#0d0d33", "#000000", False),
    ("USB",              3.75,  1.60,  58.75, -4.000, "#c68a35", "#7a531f", False),
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

    # 1. Общая длина корпуса 126.6 мм
    dim_h(ax, -63.30, 63.30, -19.0, -10.85, "126,6 мм", color="black")
    # 2. Общая высота корпуса 21.7 мм
    dim_v(ax, -10.85, 10.85, 72.0, 63.30, "21,7 мм", color="black")
    # 3. Толщина стенки корпуса 1.5 мм - callout (зона мала для линии)
    ax.annotate("стенка корпуса\n1,5 мм", xy=(30.0, 10.10), xytext=(30.0, 17.0),
                fontsize=8, ha="center", color="#333333",
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.0))
    # 4. Позиция кристалла от датума: X=0 -> X=-50.25
    dim_h(ax, 0.0, -50.25, 14.0, 10.85, "-50,25 мм (крист.)", color="#8a7100")
    # 5. Габарит капсулы 18x18 мм, по обеим осям
    dim_h(ax, -58.90, -40.90, -11.5, -9.325, "18 мм", color="#404040")
    dim_v(ax, -9.325, 8.675, -70.0, -58.90, "18 мм", color="#404040")
    # 6. Позиция PCB (9.85) и длина (101.5 мм) - одной линией по полному пролёту
    dim_h(ax, -40.90, 60.60, -14.0, -10.85, "101,5 мм (центр X=+9,85)", color="#0c451d")
    # 7. Позиция аккумулятора (23.8) и длина (65 мм)
    dim_h(ax, -8.70, 56.30, -16.5, -10.85, "65 мм (центр X=+23,8)", color="#5f5f5f")

    ax.plot(0, 0, "+", color="red", ms=10, mew=1.5, zorder=5)
    ax.text(0, -2.3, "датум (0,0)", color="red", fontsize=7, ha="center")

    ax.text(-83, -21, "Ширина корпуса (Y, не в этом сечении):\n"
                        "34,1 мм наруж. / 31,1 мм полость",
            fontsize=7.5, color="#555555", va="bottom")

    ax.set_xlim(-85, 80)
    ax.set_ylim(-23, 19)
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
        ("Окно дисплея / LCD", "#1a1a1f", "#000000"),
        ("USB", "#c68a35", "#7a531f"),
    ]
    handles = [Rectangle((0, 0), 1, 1, facecolor=fc, edgecolor=ec, linewidth=1.2)
               for _, fc, ec in legend_items]
    labels = [lab for lab, _, _ in legend_items]
    ax.legend(handles, labels, loc="lower right", fontsize=7, framealpha=0.9)


def build_zoom_panel(ax):
    # Капсула + кристалл + SiPM (те же числа, что в главной панели) плюс ESR
    # (65 мкм плёнка вокруг кристалла, half-extent ~7.065 мм - геометрия
    # RC110Detector.cc, ESR_px/nx/py/ny/pz) - тонкий контур, не заливка.
    add_rect(ax, 9.00, 9.00, -49.90, -0.325, "none", "#404040", False, lw=1.5)
    add_rect(ax, 7.50, 7.50, -49.90, -0.325, "none", "#404040", True, lw=1.0)
    add_rect(ax, 7.065, 7.065, -50.25, 0.0, "none", "#dcdce0", False, lw=1.0)
    add_rect(ax, 7.00, 7.00, -50.25, 0.0, "#d4af00", "#8a7100", False)
    add_rect(ax, 3.00, 0.40, -50.25, -7.40, "#39b7d6", "#1c6b80", False)

    # Три зазора (SESSION-STATE.md, "Зазоры при Zc=-0,325") - callout-текст,
    # НЕ масштабная линия (25-110 мкм физически неразличимы на этом масштабе).
    ax.annotate("SiPM<->капсула\n+0,025 мм", xy=(-50.25, -7.80),
                xytext=(-41.2, -8.6), fontsize=7, color="#1c6b80", ha="left",
                arrowprops=dict(arrowstyle="->", color="#1c6b80", lw=0.9))
    ax.annotate("ESR<->капсула\n+0,110 мм", xy=(-50.25, 7.07),
                xytext=(-41.2, 6.7), fontsize=7, color="#808088", ha="left",
                arrowprops=dict(arrowstyle="->", color="#808088", lw=0.9))
    ax.annotate("ESR-плёнка\n0,065 мм", xy=(-57.25, 0.0),
                xytext=(-41.2, 1.0), fontsize=7, color="#8a7100", ha="left",
                arrowprops=dict(arrowstyle="->", color="#8a7100", lw=0.9))

    ax.set_xlim(-60, -40)
    ax.set_ylim(-10, 8)
    ax.set_aspect("equal")
    ax.set_xlabel("X, мм")
    ax.set_ylabel("Z, мм")
    ax.grid(alpha=0.2, lw=0.5)
    ax.set_title("Деталь узла кристалл-капсула-SiPM (зум)")


def main() -> int:
    fig, (ax_main, ax_zoom) = plt.subplots(
        1, 2, figsize=(19, 6.5), gridspec_kw={"width_ratios": [2.3, 1.0]})
    build_main_panel(ax_main)
    build_zoom_panel(ax_zoom)
    fig.suptitle("RadiaCode-110 - разрез Y=0, размеры в мм", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(f"OK: {OUT} ({OUT.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())