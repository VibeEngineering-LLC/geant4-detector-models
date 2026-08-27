# -*- coding: utf-8 -*-
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_chains_material as fcm

matplotlib.rcParams["font.size"] = 9
matplotlib.rcParams["figure.dpi"] = 130
matplotlib.rcParams["axes.grid"] = True
matplotlib.rcParams["grid.alpha"] = 0.3

OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)

CHAIN_COLOR = {"K40": "#1f77b4", "Ra226": "#d62728", "Th232": "#2ca02c"}
MATERIAL_LS = {"brick": "-", "concrete": "--"}

def plot_panel(ax, r, which, subtitle):
    meas_cps = r["cnt"] / r["live"]
    ax.step(r["e_meas"], meas_cps, where="mid", color="0.35", lw=1.0, label="измерение")
    
    pred = r["pred_a"] if which=="A" else r["pred_b"]
    ax.plot(r["e_meas"], pred / r["live"], color="crimson", lw=1.6, label="сумма модели")
    
    amp = r["amp_a"] if which=="A" else r["amp_b"]
    
    zero = []
    for k, name in enumerate(r["names"]):
        if amp[k] <= 0:
            zero.append(name)
            continue
        
        if name == "mu":
            ax.plot(r["e_meas"], amp[k]*r["A_counts"][:,k]/r["live"], color="0.5", lw=0.9, ls=":", alpha=0.9, label="mu %.1f мюон/с" % amp[k])
        else:
            chain, material = name.rsplit("_", 1)  # разбор имени вида "K40_brick"
            ax.plot(r["e_meas"], amp[k]*r["A_counts"][:,k]/r["live"], color=CHAIN_COLOR[chain], lw=1.1, ls=MATERIAL_LS[material], alpha=0.85, label="%s %.1f Бк/кг" % (name, amp[k]))
    
    ax.set_yscale("log")
    ax.set_xlim(fcm.ftc.E_LO, fcm.ftc.E_HI)
    ax.set_ylim(max(meas_cps.max()*1e-5, 1e-7), meas_cps.max()*2)
    ax.set_xlabel("Энергия, кэВ")
    ax.set_ylabel("Скорость счёта, 1/(с·канал)")
    
    chi2ndf = r["chi2ndf_a"] if which=="A" else r["chi2ndf_b"]
    shape = r["shape_a"] if which=="A" else r["shape_b"]
    title = "Критерий %s: chi2/ndf=%.1f, невязка формы=%.4f\n%s" % (which, chi2ndf, shape, subtitle)
    if zero:
        title += "\nобнулены NNLS: %s" % ", ".join(zero)
    ax.set_title(title, fontsize=10)
    
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7.5, framealpha=1.0)

def main():
    r = fcm.prepare()
    
    fig = plt.figure(figsize=(13, 11.5))
    gs = fig.add_gridspec(2, hspace=0.55, height_ratios=[1, 1])
    
    ax1 = fig.add_subplot(gs[0, 0])
    plot_panel(ax1, r, "A", "цепочки K40/Ra226/Th232 x кирпич/бетон (сплошная=кирпич, пунктир=бетон)")
    
    ax2 = fig.add_subplot(gs[1, 0])
    plot_panel(ax2, r, "B", "то же, критерий формы — кривые кирпич/бетон почти совпадают: материал НЕ разделяется по спектру")
    
    fig.suptitle("RadiaCode-103, фон помещения: разложение по цепочкам и РАЗДЕЛЬНЫМ материалам (кирпич/бетон), метод 1", y=0.995)
    
    # Не вызываем tight_layout(), чтобы не съедался hspace
    path = os.path.join(OUT_DIR, "RC103_bg_decomposition_chains_material.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(os.path.abspath(path))

if __name__ == "__main__":
    main()
