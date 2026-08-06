# -*- coding: utf-8 -*-
"""Рисунок кривой эффективности из results/eff_point_end10cm.csv.

Показываются ДВЕ величины, потому что их путают чаще всего:
  eps по ППП   — доля испущенных квантов, дающих отсчёт в пике полного
                 поглощения; именно она входит в расчёт активности;
  eps полная   — доля квантов, давших любой отсчёт (пик плюс континуум).
Отношение первой ко второй — «пик/полная», отдельная кривая внизу: она
показывает, какая доля отклика приходится на пик. Ход НЕМОНОТОННЫЙ: растёт с
0,68 (40 кэВ) до максимума 0,77 (122 кэВ) и падает до 0,06 (3000 кэВ). К
единице не идёт нигде: на мягком крае её держат вылет характеристического
рентгена кристалла и рассеяние в лицевом стеке. Падение на жёстком крае —
утечка вдоль длинной оси бруска: в торцевой ориентации пучок идёт вдоль неё,
стороны, а 15 мм относятся к ориентации на рабочую грань, которая здесь не
считается.

    python analysis/plot_curve.py [<кривая.csv>]
"""
import os
import sys
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter

_HERE = os.path.dirname(os.path.abspath(__file__))
DEF = os.path.normpath(os.path.join(_HERE, "..", "results",
                                    "eff_point_end10cm.csv"))


# Подпись и имя выходного файла берутся ИЗ ШАПКИ кривой, а не из кода. Прежде
# и то и другое было захардкожено под торцевую геометрию: построение кривой
# «на грань» молча перезаписывало торцевой рисунок торцевой же подписью.
# Найдено независимым аудитом 05.08.2026.
FACE = {
    "end10cm": ("торец 18 × 15 мм", "торцу"),
    "face10cm": ("рабочая грань (длинная)", "рабочей грани"),
}


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEF
    if not os.path.exists(path):
        print("Нет файла кривой: %s\nСначала: python analysis/export_curve.py "
              "<каталог спектров>" % path)
        return 2
    head, rows = [], []
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#"):
                head.append(ln[1:].strip())
            else:
                break
    with open(path, encoding="utf-8") as f:
        rd = csv.DictReader(l for l in f if not l.startswith("#"))
        for r in rd:
            rows.append({k: float(v) if k != "shelf" else v
                         for k, v in r.items()})
    rows.sort(key=lambda r: r["E_keV"])
    mac = ""
    for h in head:
        if h.startswith("run_args"):
            mac = h.split("=", 1)[1].strip()
    tag = os.path.splitext(mac)[0].replace("curve_point_", "")
    if tag not in FACE:
        print("В шапке %s нет узнаваемого run_args (найдено %r).\n"
              "Подпись и имя рисунка выводятся из него — без него рисунок "
              "может оказаться подписан не той гранью." % (path, mac))
        return 2
    grain, dative = FACE[tag]
    e = [r["E_keV"] for r in rows]
    ep = [r["eps_peak"] for r in rows]
    de = [r["d_eps_peak"] for r in rows]
    et = [r["eps_total"] for r in rows]
    pt = [a / b for a, b in zip(ep, et)]

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9.6, 8.4), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2.4, 1],
                                                   hspace=0.08))
    ax.errorbar(e, ep, yerr=de, marker="o", ms=4.5, lw=1.4, capsize=2.5,
                color="#1f4e79", label="по ППП (пик полного поглощения)")
    ax.plot(e, et, marker="s", ms=4, lw=1.2, ls="--", color="#b0691f",
            label="полная (любой отсчёт)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Эффективность на 4π")
    ax.grid(True, which="both", lw=0.4, alpha=0.5)
    ax.legend(fontsize=9, frameon=False)
    ax.set_title("AtomSpectra Nano 16 PRO: точечный источник на оси кристалла,"
                 " 10 см от корпуса\nисточник обращён к %s (%s); расчёт "
                 "Geant4, измерением не подтверждён" % (dative, grain),
                 fontsize=10.5)

    ax2.plot(e, pt, marker="o", ms=4, lw=1.3, color="#2f7a4a")
    ax2.set_xscale("log")
    ax2.set_ylim(0, 1.05)
    ax2.set_xlabel("Энергия, кэВ")
    ax2.set_ylabel("пик / полная")
    ax2.grid(True, which="both", lw=0.4, alpha=0.5)

    # Подписи оси энергий — по УЗЛАМ, а не по декадам. На логарифмической шкале
    # matplotlib по умолчанию подписывает только 100 и 1000, и снять с рисунка
    # положение мягкого края нельзя: как раз там узлы стоят гуще всего.
    ticks = [t for t in (20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000) if min(e) * 0.98 <= t <= max(e) * 1.02]
    for a_ in (ax, ax2):
        a_.xaxis.set_major_locator(FixedLocator(ticks))
        a_.xaxis.set_major_formatter(ScalarFormatter())
        a_.xaxis.set_minor_formatter(NullFormatter())
    for lab in ax2.get_xticklabels():
        lab.set_rotation(45)
        lab.set_ha("right")
        lab.set_fontsize(8)

    note = [h for h in head if h.startswith(("src_sha1", "fwhm_662"))]
    fig.text(0.012, 0.012, "  |  ".join(note), fontsize=7.5, color="#666666")
    out = os.path.join(os.path.dirname(path), "..", "drawings",
                       "nano16pro_eff_point_%s.png" % tag)
    out = os.path.normpath(out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print("записано: %s  (%d узлов)" % (out, len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
