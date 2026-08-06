# -*- coding: utf-8 -*-
"""Разложение функции отклика по каналам взаимодействия.

Канал ставится В МОМЕНТ СОБЫТИЯ по истории процессов (`main.cc`, enum `Chan`)
и пишется отдельным файлом `<узел>_chan.csv`. Из готового спектра его
восстановить нельзя: к моменту записи вклады уже сложены, и любое разложение
постфактум было бы подгонкой формы, а не разбором.

Каналы взаимоисключающие и в сумме дают полный спектр — это проверяется самим
прогоном, расхождение печатается как отказ.

    python analysis/draw_channels.py <каталог спектров> [E1 E2 E3]
"""
import io
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(_HERE, "..", "drawings",
                                    "nano16pro_channels.png"))
FWHM_662 = 41.60
# Ширина канала отображения СВОЯ на каждом узле: постоянная ширина либо
# дробит мягкий узел на десяток ступеней, либо превращает жёсткий в частокол
# одиночных отсчётов. Привязка к ПШПВ держит и то, и другое: четыре канала на
# полуширину — предел, за которым размытие уже ничего не сглаживает.
def step_for(e0):
    return max(5.0, round(fwhm(e0) / 4.0 / 5.0) * 5.0)
E_MAX = 3200.0

# Порядок и подписи — по физическому смыслу, а не по алфавиту: сперва то, что
# кормит пик, потом каналы вылета, потом внешние.
ORDER = [
    ("photo",      "фотоэффект, ничего не вылетело",        "#1f5fa8"),
    ("compt_full", "комптон и поглощение, без вылета",      "#4a90d9"),
    ("pair_full",  "пары, оба 511 поглощены",               "#7fb3e8"),
    ("compt_esc1", "однократный комптон, квант ушёл",       "#b5651d"),
    ("compt_escN", "многократный комптон, квант ушёл",      "#e08a3c"),
    ("pair_esc1",  "пары, вылетел один 511",                "#2f7a4a"),
    ("pair_esc2",  "пары, вылетели оба",                    "#5fbf87"),
    ("brems_esc",  "вылет тормозного",                      "#7a2020"),
    ("xray_esc",   "вылет характеристического рентгена",    "#c04040"),
    ("external",   "вторичные из корпуса и обёртки",        "#6b6b6b"),
    ("other",      "остаточный канал (сторож)",             "#b0b0b0"),
]


def fwhm(e):
    return FWHM_662 * math.sqrt(max(e, 1.0) / 661.657)


def read_chan(path):
    head, names, rows = {}, None, []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("# ").split("=", 1)
                head[k.strip()] = v.strip()
            continue
        p = ln.split(",")
        if names is None:
            names = p[1:]
            continue
        rows.append((float(p[0]), [float(x) for x in p[1:]]))
    return head, names, rows


def broaden(pairs, w, nch, STEP):
    """Свёртка с приборным разрешением; вес w переводит отсчёты в 4π."""
    out = [0.0] * nch
    for e, c in pairs:
        if c <= 0:
            continue
        s = fwhm(e) / 2.3548
        lo, hi = max(0.0, e - 4 * s), min(E_MAX, e + 4 * s)
        acc, norm = [], 0.0
        for k in range(int(lo / STEP), min(int(hi / STEP) + 1, nch)):
            g = math.exp(-0.5 * (((k + 0.5) * STEP - e) / s) ** 2)
            acc.append((k, g))
            norm += g
        if norm <= 0:
            out[min(int(e / STEP), nch - 1)] += c * w
            continue
        for k, g in acc:
            out[k] += c * w * g / norm
    return out


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    src = sys.argv[1]
    wanted = [float(x) for x in sys.argv[2:]] or [180.0, 1480.0, 3000.0]

    files = [f for f in sorted(os.listdir(src)) if f.endswith("_chan.csv")]
    if not files:
        raise SystemExit("в %s нет файлов *_chan.csv — прогон старой ревизии"
                         % src)
    index = {}
    for f in files:
        head, names, rows = read_chan(os.path.join(src, f))
        index[float(head["E_prim_keV"])] = (head, names, rows)

    # Общей оси НЕТ намеренно: у каждого узла свой предел энерговыделения,
    # и при общей шкале два верхних узла сжались бы в левую четверть поля.
    fig, axes = plt.subplots(len(wanted), 1, figsize=(12.6, 4.0 * len(wanted)),
                             sharex=False)
    if len(wanted) == 1:
        axes = [axes]

    stamp = "?"
    for ax, e_want in zip(axes, wanted):
        e0 = min(index, key=lambda x: abs(x - e_want))
        head, names, rows = index[e0]
        stamp = head.get("src_sha1", "?")
        w = float(head["solid_angle_frac"]) / float(head["N_primaries"])
        step = step_for(e0)
        nch = int(E_MAX / step)
        cols = [(k + 0.5) * step for k in range(nch)]

        total = broaden([(e, sum(v)) for e, v in rows], w, nch, step)
        ax.step(cols, total, where="mid", lw=1.6, color="#111111",
                label="полный отклик", zorder=5)

        share = []
        for name, label, colour in ORDER:
            if name not in names:
                continue
            j = names.index(name)
            raw = sum(v[j] for _, v in rows)
            if raw == 0:
                continue
            cur = broaden([(e, v[j]) for e, v in rows], w, nch, step)
            ax.step(cols, cur, where="mid", lw=1.0, color=colour,
                    label="%s — %.1f %%" % (label, 100.0 * raw
                                            / sum(sum(v) for _, v in rows)))
            share.append((name, raw))

        ax.set_yscale("log")
        # Шкала по каждому узлу своя: на 180 кэВ отклик кончается на 200, и
        # общая шкала до 3200 оставила бы девять десятых поля пустыми.
        ax.set_xlim(0, min(3200.0, e0 * 1.15))
        top = max(total)
        # Запас сверху — под легенду: иначе она ложится на континуум, который
        # и есть предмет рисунка.
        ax.set_ylim(top / 2.0e3, top * 40)
        ax.set_ylabel("вероятность на квант в 4π")
        ax.xaxis.set_major_locator(MultipleLocator(250))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.grid(True, which="major", alpha=0.26, lw=0.6)
        ax.grid(True, which="minor", axis="x", alpha=0.11, lw=0.4)
        ax.legend(fontsize=7.8, loc="upper right", ncol=2, framealpha=0.94)
        ax.set_title("падающий квант %.0f кэВ, канал отображения %.0f кэВ"
                     % (e0, step), fontsize=10.5, pad=4)

    axes[-1].set_xlabel("Энерговыделение, кэВ")
    fig.suptitle("AtomSpectra Nano 16 PRO: функция отклика, разложенная по "
                 "каналам взаимодействия\n"
                 "канал ставится по истории процессов события; каналы "
                 "взаимоисключающие и в сумме дают полный отклик",
                 fontsize=12, y=0.995)
    fig.text(0.5, 0.004,
             "Свёрнуто с ПШПВ(E) = 41,60·√(E/661,657) кэВ. Штамп исходников "
             "%s. Проценты в легенде — доля канала во ВСЕХ событиях с сигналом "
             "на этом узле." % stamp,
             fontsize=8.4, ha="center", color="#555555")
    fig.subplots_adjust(left=0.075, right=0.985, top=0.925, bottom=0.062,
                        hspace=0.16)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("записано: %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
