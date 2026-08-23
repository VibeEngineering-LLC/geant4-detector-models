# -*- coding: utf-8 -*-
import os
import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import numpy as np
import fit_nuclides as fn
import fit_lines as fl
import read_rcxml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def build_grp(amps, names, A):
    """Имя группы -> кривая на графике (уже цепочки, доп. группировки нет)."""
    label = {"K-40": "K-40", "Ra-226": "цепочка Ra-226",
             "Th-232": "цепочка Th-232", "mu": "мюоны"}
    return {label.get(n, n): a * A[:, names.index(n)]
            for n, a in amps.items() if n in names}


def main():
    # Равновесие Ra-226/Th-232 задано на уровне шаблонов (fn.merge_by_chain) —
    # amps уже содержит ОДНУ амплитуду на цепочку, группировка для картинки
    # больше не нужна: имена в amps совпадают с именами столбцов A напрямую.
    amps, a_mu = fn.main()
    smp = read_rcxml.read(fn.MEASURED)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))[:len(cnt)]
    cps_meas = cnt / smp.live

    names, cols = fn.merge_by_chain(*fn.load_templates())
    mu, pdg = fn.load_muons()
    if mu is not None:
        names.append("mu")
        cols.append(mu)
    A = np.zeros((len(e_meas), len(cols)))
    for k, c in enumerate(cols):
        A[:, k] = fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)

    grp = build_grp(amps, names, A)
    model = sum(grp.values())
    
    # Рисование
    fig, ax = plt.subplots(figsize=(13.0, 7.2))
    
    # Измерение
    ax.step(e_meas, cps_meas, where="mid", lw=1.0, color="#111111", 
            label="измерение: %s, живое %.1f сут" % (os.path.basename(fn.MEASURED), smp.live / 86400))
    
    # Полная модель
    ax.step(e_meas, model, where="mid", lw=1.4, color="#d62728", 
            label="модель, сумма компонентов")
    
    # Компоненты
    order = ["K-40", "цепочка Ra-226", "цепочка Th-232", "мюоны"]
    colors = ["#1f77b4", "#2ca02c", "#9467bd", "#8c564b"]
    mask = (e_meas >= 20) & (e_meas < 2830)
    
    for g, c in zip(order, colors):
        if g in grp:
            ax.step(e_meas, grp[g], where="mid", lw=0.9, color=c,
                    label="%s (%.1f%%)" % (g, grp[g][mask].sum() / model[mask].sum() * 100))
    
    # Оси и настройки
    ax.set_yscale("log")
    ax.set_xlim(0, 3000)
    y_min = max(1e-7, cps_meas[mask][cps_meas[mask] > 0].min()) / 3
    ax.set_ylim(y_min, max(cps_meas) * 3)
    
    ax.set_xlabel("Энергия, кэВ")
    ax.set_ylabel("скорость счёта, отсч/с на канал")
    
    # Сетка
    ax.grid(True, which="major", alpha=0.3, lw=0.6)
    ax.grid(True, which="minor", alpha=0.12, lw=0.4)
    ax.xaxis.set_major_locator(MultipleLocator(250))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    
    # Легенда
    ax.legend(fontsize=9, loc="upper right", framealpha=0.94)
    
    # Заголовок и подзаголовок
    ax.set_title("Открытый фон RadiaCode-103: нуклидное разложение (K-40 + цепочки Ra-226/Th-232 + мюоны)", fontsize=12)
    fig.text(0.5, 0.955, "амплитуды подобраны по нетто-площадям линий (не по полосам); Ra-226 и Th-232 — вековое равновесие ЗАДАНО (одна амплитуда на цепочку), не проверяется постфактум",
             fontsize=9.5, ha="center", color="#444444")
    
    # Подпись внизу. a_mu — поток в мю/с, НЕ множитель; отношение к PDG явно.
    fig.text(0.5, 0.01,
             "K-40 = %.1f Бк/кг; Ra-226 = %.1f (равновесие задано); "
             "Th-232 = %.1f (равновесие задано); мюоны %.0f мю/с = %.2f x PDG. "
             "Шаблоны: %s" % (amps.get("K-40", 0), amps.get("Ra-226", 0),
                              amps.get("Th-232", 0), a_mu,
                              a_mu / pdg if pdg else float("nan"), fn.BG_DIR),
             fontsize=8.2, ha="center", color="#666666")
    
    # Подстройка
    fig.subplots_adjust(left=0.07, right=0.985, top=0.9, bottom=0.09)
    
    # Сохранение. G4MODELS_OUT_SUFFIX - чтобы сравнить два режима шаблонов
    # (таблица линий / полный распад) не перетирая картинку друг друга.
    suffix = os.environ.get("G4MODELS_OUT_SUFFIX", "")
    OUT = os.path.normpath(os.path.join(_HERE, "..", "results", "figures",
                                        "bg_nuclide_decomposition%s.png" % suffix))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("записано: %s" % OUT)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
