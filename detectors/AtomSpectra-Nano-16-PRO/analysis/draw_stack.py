# -*- coding: utf-8 -*-
"""Слоевой стек («пирог») входного окна AtomSpectra Nano 16 PRO.

Отдельный рисунок заведён потому, что порядок слоёв — место, где уже была
ошибка: зеркало геометрии в `draw_setup.py` до 06.08.2026 рисовало обёртку
перевёрнутой (фольга на кристалле вместо ПТФЭ), и это не ловилось ничем, так
как СУММА толщин от порядка не зависит. Общий разрез прибора слои показывает
полосками в доли миллиметра, различить порядок на нём нельзя.

Толщины НЕ дублируются здесь константами, а читаются из `geometry/
ASN16Detector.hh` — четвёртого зеркала заводить нельзя, дефект именно в этом
классе. Плотности взяты из базы Geant4/NIST.

Замок: сумма лицевого стека сверяется с зашитой здесь величиной
`FACE_STACK_REF`, переписанной из расчёта с шестью знаками. Это сверка `.hh` с
константой, а НЕ с печатью программы: программа печатает три знака, и до
шестого сверять с ней нечего. Торцевой стек не сверяется ни с чем —
независимого числа для него нигде не печатается.

Показаны ДВА стека, потому что они разные, а в пучке стоит второй:
  рабочая грань (Y) — стенка корпуса Al 1,50;
  торец (Z)        — крышка Al 1,50 с фрезеровкой до 0,60 напротив кристалла;
                     в пучке стоит остаток 0,60. Воздушного зазора между
                     обёрткой и крышкой нет: `zCapFi = zFoilF`.

    python analysis/draw_stack.py
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
                                    "nano16pro_stack.png"))

# Поверхностная плотность лицевого стека, посчитанная по принятым толщинам и
# плотностям NIST. Служит замком: если разбор .hh или плотности разъедутся,
# скрипт остановится, а не нарисует красивую неправду. Программа (ReportMasses)
# печатает ту же величину, но с тремя знаками — 0,571, — поэтому шестизначное
# сравнение идёт с этой константой, а не с её печатью.
FACE_STACK_REF = 0.651840      # г/см² (стенка 1,50 + фольга 0,10 + ПТФЭ 1,00)
# Замок СРАБОТАЛ по назначению 06.08.2026: при правке wFront 1,20 -> 1,50 по
# указанию оператора скрипт остановился, вместо того чтобы молча нарисовать
# новую сумму под старой подписью. Константу обновлять СЛЕДОМ за .hh — падение
# здесь и есть напоминание.

# Плотности по базе Geant4/NIST, г/см³. G4_POLYSTYRENE стоит вместо ABS и
# G4_Al вместо сплава — вещества-заменители, помечены на рисунке.
RHO = {"CsI(Tl)": 4.510, "ПТФЭ": 2.200, "Al": 2.699}


def ru(x, nd=2):
    return ("%.*f" % (nd, x)).replace(".", ",")


def geom_from_header(path):
    """Толщины из Nano16Geom. -> dict имя->мм."""
    src = open(path, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^\s*double\s+(\w+)\s*=\s*([0-9.]+)\s*;", src,
                         re.MULTILINE):
        out[m.group(1)] = float(m.group(2))
    need = ("ptfe", "alFoil", "wFront", "wCap", "wCapWin", "cryX", "cryY",
            "cryZ")
    miss = [k for k in need if k not in out]
    if miss:
        raise SystemExit("в %s не найдены поля: %s" % (path, ", ".join(miss)))
    return out


def draw_stack(ax, layers, title, note):
    """layers: список (имя, мм, плотность или None для воздуха, цвет)."""
    y = 0.0
    total = 0.0
    for name, t, rho, col in layers:
        h = t ** 0.55            # сжатие: иначе 0,10 мм не видно рядом с 2,00
        ax.add_patch(Rectangle((0, y), 1.0, h, facecolor=col,
                               edgecolor="#333333", lw=0.9))
        sd = None if rho is None else 0.1 * t * rho
        if sd is not None:
            total += sd
        lab = "%s   %s мм" % (name, ru(t))
        ax.text(1.06, y + h / 2, lab, fontsize=9, va="center", ha="left")
        if sd is not None:
            ax.text(-0.06, y + h / 2, "%s г/см²" % ru(sd, 4), fontsize=8,
                    va="center", ha="right", color="#555555")
        else:
            ax.text(-0.06, y + h / 2, "—", fontsize=8, va="center",
                    ha="right", color="#999999")
        y += h
    ax.add_patch(Rectangle((0, -0.62), 1.0, 0.60, facecolor="#d9a520",
                           edgecolor="#333333", lw=1.1))
    ax.text(1.06, -0.32, "CsI(Tl)  чувствительная область", fontsize=9,
            va="center", ha="left", fontweight="bold")
    ax.annotate("", xy=(0.5, y + 0.30), xytext=(0.5, y + 0.02),
                arrowprops=dict(arrowstyle="<-", lw=1.4, color="#b03030"))
    ax.text(0.5, y + 0.34, "квант снаружи", fontsize=8.5, ha="center",
            color="#b03030")
    ax.set_title(title, fontsize=10.5, pad=16)
    ax.text(0.5, -0.95, note, fontsize=8.5, ha="center", color="#444444")
    ax.text(0.5, -1.30, "суммарно  %s г/см²" % ru(total, 4), fontsize=10,
            ha="center", fontweight="bold")
    ax.set_xlim(-1.35, 2.95)
    ax.set_ylim(-1.45, y + 0.65)
    ax.axis("off")
    return total


def main():
    g = geom_from_header(HH)
    face = [("ПТФЭ (на кристалле)", g["ptfe"], RHO["ПТФЭ"], "#f2f2f2"),
            ("Al-фольга (на ПТФЭ)", g["alFoil"], RHO["Al"], "#c9c9c9"),
            ("стенка корпуса Al*", g["wFront"], RHO["Al"], "#8d9299")]
    # Воздушного слоя в торцевом стеке НЕТ: ASN16Detector.cc:137 задаёт
    # zCapFi = zFoilF тождественно, обёртка упёрта в крышку. Прежняя версия
    # рисовала здесь «воздух полости 0,90 мм» — единственную толщину, не
    # прочитанную из .hh, то есть ровно то, что запрещает докстринг этого же
    # файла; на сумму 0,4590 г/см² она не влияла (аудит кода 06.08.2026).
    end = [("ПТФЭ (на кристалле)", g["ptfe"], RHO["ПТФЭ"], "#f2f2f2"),
           ("Al-фольга (на ПТФЭ)", g["alFoil"], RHO["Al"], "#c9c9c9"),
           ("крышка Al в окне фрезеровки", g["wCapWin"], RHO["Al"], "#8d9299")]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 6.4))
    t_face = draw_stack(
        a1, face, "Рабочая грань %s × %s мм (ось Y)"
        % (ru(g["cryX"] if "cryX" in g else 18.0), ru(g["cryZ"])),
        "в настоящей постановке НЕ в пучке")
    t_end = draw_stack(
        a2, end, "Торец %s × %s мм (ось Z)"
        % (ru(g["cryX"] if "cryX" in g else 18.0), ru(g["cryY"])),
        "В ПУЧКЕ: источник обращён сюда")

    if abs(t_face - FACE_STACK_REF) / FACE_STACK_REF > 1e-3:
        raise SystemExit(
            "стек рабочей грани %s г/см² не сходится с опорной величиной %s — "
            "разошлись разбор .hh или плотности" % (ru(t_face, 6),
                                                    ru(FACE_STACK_REF, 6)))

    fig.suptitle("AtomSpectra Nano 16 PRO — слоевой стек входного окна "
                 "(порядок слоёв: ПТФЭ на кристалле, фольга на ПТФЭ)",
                 fontsize=12, y=0.98)
    fig.text(0.5, 0.035,
             "Толщины прочитаны из geometry/ASN16Detector.hh, плотности — база "
             "Geant4/NIST. Масштаб по толщине сжат (степень 0,55), иначе слой "
             "0,10 мм неразличим рядом с 2,00 мм.\n"
             "* вещество-заменитель: корпус, крышка и фольга — чистый Al при "
             "заявленном сплаве. В торцевом стеке показан ОСТАТОК крышки в "
             "окне фрезеровки (0,60 мм); вне окна там 1,50 мм. Воздушного "
             "зазора между обёрткой и крышкой в модели нет.",
             fontsize=8, ha="center", color="#555555")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    print("записано: %s" % OUT)
    print("стек рабочей грани %s г/см² (опорная величина %s; программа "
          "печатает её же с тремя знаками)"
          % (ru(t_face, 4), ru(FACE_STACK_REF, 4)))
    print("стек торца         %s г/см²  — ЭТОТ стоит в пучке" % ru(t_end, 4))
    return 0


if __name__ == "__main__":
    sys.exit(main())
