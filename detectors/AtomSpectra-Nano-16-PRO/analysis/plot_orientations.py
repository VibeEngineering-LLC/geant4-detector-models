# -*- coding: utf-8 -*-
"""Две ориентации на одном рисунке и их отношение.

Прибор анизотропен: рабочая грань кристалла 18 × 57 мм = 10,26 см² против
торца 18 × 15 мм = 2,70 см², отношение площадей 3,80. Отношение кривых
эффективности с ним не совпадает и от энергии зависит: у торцевой ориентации
кванты проходят вдоль бруска 57 мм и поглощаются лучше, у ориентации на грань
им доступно только 15 мм. Оба обстоятельства действуют навстречу друг другу,
поэтому отношение — самостоятельное проверяемое число, а не пересчёт площадей.

    python analysis/plot_orientations.py
"""
import os
import sys
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(_HERE, "..", "results"))
END = os.path.join(RES, "eff_point_end10cm.csv")
FACE = os.path.join(RES, "eff_point_face10cm.csv")


def ru(x, nd=2):
    return ("%.*f" % (nd, x)).replace(".", ",")


def load(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        rd = csv.DictReader(l for l in f if not l.startswith("#"))
        for r in rd:
            rows.append({k: (v if k == "shelf" else float(v))
                         for k, v in r.items()})
    rows.sort(key=lambda r: r["E_keV"])
    return rows


def main():
    for p in (END, FACE):
        if not os.path.exists(p):
            print("Нет файла кривой: %s\nСначала прогон сетки и "
                  "analysis/export_curve.py" % p)
            return 2
    a, b = load(END), load(FACE)
    if [r["E_keV"] for r in a] != [r["E_keV"] for r in b]:
        print("Сетки энергий двух кривых не совпадают — отношение не строится.")
        return 1
    e = [r["E_keV"] for r in a]
    ye = [r["eps_peak"] for r in a]
    yf = [r["eps_peak"] for r in b]
    de = [r["d_eps_peak"] for r in a]
    df = [r["d_eps_peak"] for r in b]
    k = [f / g for f, g in zip(yf, ye)]
    dk = [x * ((p / q) ** 2 + (s / t) ** 2) ** 0.5
          for x, p, q, s, t in zip(k, df, yf, de, ye)]

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9.6, 8.4), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2.4, 1],
                                                   hspace=0.08))
    ax.errorbar(e, yf, yerr=df, marker="s", ms=4.5, lw=1.4, capsize=2.5,
                color="#b0691f",
                label="на рабочую грань 18 × 57 мм (10,26 см²)")
    ax.errorbar(e, ye, yerr=de, marker="o", ms=4.5, lw=1.4, capsize=2.5,
                color="#1f4e79", label="на торец 18 × 15 мм (2,70 см²)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Эффективность по ППП на 4π")
    ax.grid(True, which="both", lw=0.4, alpha=0.5)
    ax.legend(fontsize=9, frameon=False)
    ax.set_title("AtomSpectra Nano 16 PRO: точечный источник на оси, 10 см от "
                 "наружной поверхности корпуса\nдве ориентации; расчёт Geant4,"
                 " измерением не подтверждён", fontsize=10.5)

    ax2.errorbar(e, k, yerr=dk, marker="o", ms=4, lw=1.3, capsize=2,
                 color="#2f7a4a")
    ax2.axhline(10.26 / 2.70, lw=1.0, ls="--", color="#999999")
    ax2.text(e[0], 10.26 / 2.70 + 0.06, " отношение площадей граней 3,80",
             fontsize=8, color="#777777", va="bottom")
    ax2.set_xscale("log")
    ax2.set_xlabel("Энергия, кэВ")
    ax2.set_ylabel("грань / торец")
    ax2.grid(True, which="both", lw=0.4, alpha=0.5)

    out = os.path.normpath(os.path.join(RES, "..", "drawings",
                                        "nano16pro_eff_orientations.png"))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print("записано: %s" % out)
    print("%9s %12s %12s %8s" % ("E, кэВ", "торец", "грань", "отношение"))
    for x, p, q, r in zip(e, ye, yf, k):
        print("%9.1f %12.4e %12.4e %8s" % (x, p, q, ru(r)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
