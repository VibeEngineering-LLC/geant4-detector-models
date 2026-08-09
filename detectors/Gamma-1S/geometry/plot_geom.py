"""Осевой разрез геометрии ГАММА-1С по выгрузке из построенного дерева Geant4.

Источник данных — geom.csv, который пишет DumpGeometry() в main.cc обходом
G4PhysicalVolumeStore. Рисуется то, что Geant4 СОБРАЛ, а не то, что задумано
в исходнике: расхождение между замыслом и построенным телом здесь видно.

    g1s.exe nop.mac vessel 1.6 OISN16      (при G1S_DUMP_GEOM=...\\geom.csv)
    python plot_geom.py geom.csv g1s_section.png
"""
import sys
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Цвет и подпись по материалу. Ключ — имя материала Geant4.
MAT = {
    "G4_Pb":              ("#7a7f87", "свинец"),
    "G4_Cd":              ("#d8c896", "кадмий"),
    "G4_Cu":              ("#c1743a", "медь"),
    "G4_STAINLESS-STEEL": ("#4a6fa5", "сталь"),
    "G4_SODIUM_IODIDE":   ("#2f9e44", "NaI(Tl)"),
    "MgO_powder":         ("#f1f3f5", "MgO"),
    "G4_Al":              ("#ced4da", "алюминий"),
    "G4_RUBBER_NATURAL":  ("#343a40", "резина, гель"),
    "G4_Pyrex_Glass":     ("#a5d8ff", "стекло ФЭУ"),
    "G4_Galactic":        ("#e7f5ff", "вакуум ФЭУ"),
    "Electronics":        ("#8c6239", "электроника"),
    "G4_POLYPROPYLENE":   ("#ffe066", "полипропилен"),
    "Sample":             ("#a0522d", "проба ОИСН-16"),
    "G4_AIR":             ("#ffffff", "воздух"),
}


def read(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            break
        rd = csv.DictReader([line] + f.readlines())
        for r in rd:
            if r["rin_mm"] == "?":
                print("тело неизвестного класса:", r["name"])
                continue
            rows.append({
                "name": r["name"], "mat": r["material"],
                "ri": float(r["rin_mm"]), "ro": float(r["rout_mm"]),
                "z0": float(r["z0_mm"]), "z1": float(r["z1_mm"]),
            })
    return rows


def draw(rows, out):
    fig, ax = plt.subplots(figsize=(8.2, 11.0), dpi=170)
    # Крупные тела в глубину, мелкие поверх: иначе кристалл закроется свинцом.
    rows = sorted(rows, key=lambda r: -(r["ro"] - r["ri"]) * (r["z1"] - r["z0"]))
    seen = {}
    for r in rows:
        col, lab = MAT.get(r["mat"], ("#ff00ff", r["mat"]))
        for sign in (+1, -1):                     # разрез симметричен по оси
            x = sign * r["ro"] if sign < 0 else r["ri"]
            w = r["ro"] - r["ri"]
            ax.add_patch(Rectangle((x if sign > 0 else -r["ro"], r["z0"]), w,
                                   r["z1"] - r["z0"], facecolor=col,
                                   edgecolor="#212529", linewidth=0.35,
                                   label=lab if lab not in seen else None))
            seen[lab] = col

    ax.axvline(0, color="#495057", lw=0.6, ls=(0, (12, 4, 2, 4)))
    ax.set_xlim(-330, 330)
    ax.set_ylim(-450, 215)
    ax.set_aspect("equal")
    ax.set_xlabel("радиус, мм")
    ax.set_ylabel("z от центра кристалла, мм")
    ax.set_title("ГАММА-1С: осевой разрез построенной модели\n"
                 "УДС-ГЦ-63х63 в экране-защите, сосуд Маринелли 1 л",
                 fontsize=11)
    ax.grid(alpha=0.22, lw=0.4)
    ax.set_xticks([-300, -200, -100, 0, 100, 200, 300])

    # Выноски вынесены за габарит, полка слева, чтобы не лезли на разрез.
    ann = [
        (-146, 150, "крышка экрана-защиты"),
        (-120,  60, "свинец 45 мм"),
        (-91,   30, "вкладыши Cu 1,5 + Cd 1,2"),
        (-70,   85, "полость Ø182 × 128"),
        (-60,   60, "сосуд Маринелли 1 л"),
        (-31,   15, "кристалл NaI(Tl) Ø63 × 63"),
        (-39,   41, "торец головки, z = +41"),
        (-45, -230, "канал Ø90,6 под головку"),
        (-143, -330, "стальной кожух 7,5 мм"),
        (-39, -200, "ФЭУ и электроника"),
        (-40, -390, "свинцовая пробка канала, 52 мм"),
    ]
    for i, (x, y, t) in enumerate(ann):
        ax.annotate(t, xy=(x, y), xytext=(-325, 190 - 40 * i),
                    fontsize=7.6, color="#212529", ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.55, color="#868e96",
                                    shrinkA=2, shrinkB=1,
                                    connectionstyle="angle,angleA=0,angleB=90,rad=3"))

    h, l = ax.get_legend_handles_labels()
    ax.legend(h, l, fontsize=7, loc="lower right", ncol=2, framealpha=0.95)
    fig.text(0.5, 0.012,
             "Форма — рисунок 1.1 руководства по эксплуатации (продольный разрез), "
             "промер сканлиниями; масштаб задан длиной детектора 315 мм.\n"
             "Проверка: свинец 178,4 кг, медь 1,65 кг, кадмий 1,28 кг — "
             "паспорт ДЦКИ.412131.001 ПС, табл. 2.2: «не менее» 165 / 1,6 / 1,2 кг.",
             ha="center", fontsize=7.2, color="#343a40")
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(out)
    print("записано", out)


if __name__ == "__main__":
    draw(read(sys.argv[1]), sys.argv[2])
