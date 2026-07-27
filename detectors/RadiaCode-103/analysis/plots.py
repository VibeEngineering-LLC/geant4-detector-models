# -*- coding: utf-8 -*-
"""Графики к кривым эффективности."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec
RESULTS = rcspec.RESULTS
FIGS = rcspec.rdir("figures")

LABEL = {
    "full_air_0.00": "воздух (предел без самопоглощения)",
    "full_organic_0.50": "органика 0,50 г/см³",
    "full_soil_0.80": "грунт 0,80 г/см³",
    "full_water_1.00": "вода 1,00 г/см³",
    "full_soil_1.20": "грунт 1,20 г/см³",
    "full_soil_1.60": "грунт 1,60 г/см³",
}
ORDER = ["full_air_0.00", "full_organic_0.50", "full_soil_0.80",
         "full_water_1.00", "full_soil_1.20", "full_soil_1.60"]

LINES = [(238.6, "Pb-212"), (351.9, "Pb-214"), (583.2, "Tl-208"),
         (609.3, "Bi-214"), (661.7, "Cs-137"), (911.2, "Ac-228"),
         (1460.8, "K-40"), (1764.5, "Bi-214"), (2614.5, "Tl-208")]


def style(ax):
    ax.grid(True, which="both", alpha=0.25, lw=0.6)
    ax.tick_params(labelsize=9)


def load_eff():
    p = rcspec.rdir("efficiency.csv")
    data = {}
    for line in open(p, encoding="utf-8").readlines()[1:]:
        f = line.strip().split(",")
        data.setdefault(f[0], []).append(
            (float(f[3]), float(f[4]), float(f[5]), float(f[6]), float(f[8])))
    for k in data:
        data[k] = np.array(sorted(data[k]))
    return data


def fig_efficiency(data):
    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.cm.viridis(np.linspace(0, 0.88, len(ORDER)))
    for c, cfg in zip(cmap, ORDER):
        if cfg not in data:
            continue
        a = data[cfg]
        ax.errorbar(a[:, 0], a[:, 1], yerr=a[:, 2], color=c, marker="o", ms=3.5,
                    lw=1.6, capsize=2, label=LABEL[cfg])
    for E, nuc in LINES:
        ax.axvline(E, color="0.75", lw=0.6, ls=":", zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("энергия фотона, кэВ")
    ax.set_ylabel(r"фотопиковая эффективность $\varepsilon_p$, имп/фотон")
    ax.set_title("Сосуд Маринелли 200 мл, RadiaCode 101/102/103\n"
                 "объёмный источник в пробе", fontsize=11)
    ax.legend(fontsize=8.5, framealpha=0.95)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "efficiency.png"), dpi=150)
    plt.close(fig)


def fig_selfabs(data):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    cmap = plt.cm.viridis(np.linspace(0, 0.88, len(ORDER)))
    for c, cfg in zip(cmap, ORDER):
        if cfg not in data or cfg == "full_air_0.00":
            continue
        a = data[cfg]
        ax.plot(a[:, 0], a[:, 4], color=c, marker="o", ms=3.5, lw=1.6,
                label=LABEL[cfg])
    ax.axhline(1.0, color="0.4", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("энергия фотона, кэВ")
    ax.set_ylabel("множитель самопоглощения к пределу")
    ax.set_title("Самопоглощение в пробе: во сколько раз падает пик\n"
                 "относительно случая без вещества", fontsize=11)
    ax.legend(fontsize=8.5)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "selfabsorption.png"), dpi=150)
    plt.close(fig)


def fig_peaktotal(data):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    cmap = plt.cm.viridis(np.linspace(0, 0.88, len(ORDER)))
    for c, cfg in zip(cmap, ORDER):
        if cfg not in data:
            continue
        a = data[cfg]
        ax.plot(a[:, 0], a[:, 1] / a[:, 3], color=c, marker="o", ms=3.5, lw=1.6,
                label=LABEL[cfg])
    ax.set_xscale("log")
    ax.set_xlabel("энергия фотона, кэВ")
    ax.set_ylabel("пик / полная эффективность")
    ax.set_title("Отношение пик-полное: с ростом плотности пробы пик убывает,\n"
                 "а рассеянный континуум под ним прибывает", fontsize=11)
    ax.legend(fontsize=8.5)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "peak_to_total.png"), dpi=150)
    plt.close(fig)


def fig_beta():
    p = rcspec.rdir("beta_transmission.csv")
    if not os.path.exists(p):
        return
    data = {}
    for line in open(p, encoding="utf-8").readlines()[1:]:
        f = line.strip().split(",")
        data.setdefault(f[0], []).append((float(f[1]), float(f[2]), float(f[3])))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    cmap = plt.cm.plasma(np.linspace(0, 0.8, len(data)))
    for c, cfg in zip(cmap, [k for k in ORDER if k in data]):
        a = np.array(sorted(data[cfg]))
        ax.errorbar(a[:, 0], a[:, 1], yerr=a[:, 2], color=c, marker="o", ms=3.5,
                    lw=1.6, capsize=2, label=LABEL.get(cfg, cfg))
    ax.axvline(1150, color="0.35", lw=1.0, ls="--")
    ax.annotate("порог 0,49 г/см²\n(Кац-Пенфолд)", xy=(1150, 1e-4),
                xytext=(1250, 3e-5), fontsize=8.5, color="0.25")
    for nuc, E in [("Cs-137", 514), ("K-40", 1311), ("Tl-208", 1803),
                   ("Pa-234m", 2269), ("Bi-214", 3272)]:
        ax.axvline(E, color="0.8", lw=0.6, ls=":", zorder=0)
        ax.text(E, ax.get_ylim()[1], " " + nuc, rotation=90, fontsize=7.5,
                va="top", color="0.45")
    ax.set_yscale("log")
    ax.set_xlabel("энергия электрона, кэВ")
    ax.set_ylabel("доля, давшая сигнал в кристалле")
    ax.set_title("Проникающая способность беты из пробы", fontsize=11)
    ax.legend(fontsize=8.5, loc="lower right")
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "beta_transmission.png"), dpi=150)
    plt.close(fig)


def fig_field():
    p = os.path.join(RESULTS, "wallfield_spectrum.csv")
    if not os.path.exists(p):
        return
    e, f = [], []
    for line in open(p, encoding="utf-8"):
        if line and line[0].isdigit():
            a, b = line.split(",")
            e.append(float(a))
            f.append(float(b))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.fill_between(e, f, color="#4a7fb5", alpha=0.35)
    ax.plot(e, f, color="#1f4e79", lw=1.2)
    ax.set_yscale("log")
    ax.set_xlabel("энергия фотона, кэВ")
    ax.set_ylabel("флюенс, см⁻²·с⁻¹ на канал 10 кэВ")
    ax.set_title("Поле ЕРН помещения: две трети флюенса — рассеянные фотоны,\n"
                 "а не линии", fontsize=11)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "wallfield.png"), dpi=150)
    plt.close(fig)


def main():
    os.makedirs(FIGS, exist_ok=True)
    data = load_eff()
    fig_efficiency(data)
    fig_selfabs(data)
    fig_peaktotal(data)
    fig_beta()
    fig_field()
    print("графики:", FIGS)
    for f in sorted(os.listdir(FIGS)):
        print("   ", f)


if __name__ == "__main__":
    main()
