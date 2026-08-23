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
import io
import contextlib

def run_for(f):
    fn.RN_LEAK = float(f)
    old_stdout = sys.stdout
    sys.stdout = mystdout = io.StringIO()
    try:
        amps, a_mu = fn.main()
    finally:
        sys.stdout = old_stdout
    smp = read_rcxml.read(fn.MEASURED)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))[:len(cnt)]
    cps = cnt / smp.live
    names, cols = fn.merge_by_chain(*fn.load_templates())
    mu, pdg = fn.load_muons()
    if not None is mu:
        names.append("mu")
        cols.append(mu)
    A = np.zeros((len(e), len(cols)))
    for k, c in enumerate(cols):
        A[:, k] = fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e)
    comp = {n: amps[n] * A[:, names.index(n)] for n in amps if n in names}
    model = sum(comp.values())
    mask = (e >= 20) & (e < 2830)
    ratio = model[mask].sum() / cps[mask].sum()
    return {
        "f": f,
        "e": e,
        "cps": cps,
        "model": model,
        "comp": comp,
        "ratio": ratio,
        "ra": amps.get("Ra-226", 0),
        "k40": amps.get("K-40", 0),
        "th": amps.get("Th-232", 0),
        "amu": a_mu
    }

def main():
    FS = [0.0, 0.10, 0.15, 0.20, 0.25]
    runs = []
    for f in FS:
        sys.stderr.write("Выполняется f = %g\n" % f)
        runs.append(run_for(f))
    
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(16.0, 6.8), gridspec_kw={"width_ratios": [2.15, 1.0]})
    
    # Левая панель — спектр
    smp = read_rcxml.read(fn.MEASURED)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))[:len(cnt)]
    cps = cnt / smp.live
    ax.step(e, cps, color="black", where="post", label="измеренный спектр")
    
    # Модель для f=0 и f=0.25
    ax.step(runs[0]["e"], runs[0]["model"], color="blue", where="post", label="модель, утечка радона 0 %")
    ax.step(runs[-1]["e"], runs[-1]["model"], color="red", where="post", label="модель, утечка радона 25 %")
    
    # Компоненты для f=0.15
    comp = runs[2]["comp"]
    ax.step(runs[2]["e"], comp["K-40"], color="green", where="post", linewidth=0.8, label="K-40")
    ax.step(runs[2]["e"], comp["Ra-226"], color="orange", where="post", linewidth=0.8, label="цепочка Ra-226")
    ax.step(runs[2]["e"], comp["Th-232"], color="purple", where="post", linewidth=0.8, label="цепочка Th-232")
    if "mu" in comp:
        ax.step(runs[2]["e"], comp["mu"], color="brown", where="post", linewidth=0.8, label="мюоны")
    
    ax.set_yscale("log")
    ax.set_xlim(0, 3000)
    ax.set_ylim(min(cps[cps > 0]), max(cps) * 10)
    ax.grid(True, which="major", linestyle="-", linewidth=0.5)
    ax.grid(True, which="minor", linestyle="--", linewidth=0.5)
    ax.xaxis.set_major_locator(MultipleLocator(250))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.legend(loc="upper right")
    ax.set_xlabel("Энергия, кэВ")
    ax.set_ylabel("скорость счёта, отсч/с на канал")
    ax.set_title("открытый фон, проверка гипотезы об утечке радона")
    
    # Правая панель — результат проверки
    ratios = [r["ratio"] for r in runs]
    ra_values = [r["ra"] for r in runs]
    fs = [r["f"] for r in runs]
    
    ax2.plot(fs, ratios, "o-", color="blue", label="модель / измерение")
    ax2.set_ylim(0.0, 1.0)
    ax2.axhline(y=1.0, linestyle="--", color="black", label="полное совпадение")
    ax2.set_xlabel("утечка радона, %")
    ax2.set_ylabel("модель / измерение (полный спектр)")
    
    ax2_twin = ax2.twinx()
    ax2_twin.plot(fs, ra_values, "o-", color="red", label="Ra-226, Бк/кг")
    ax2_twin.set_ylabel("Ra-226, Бк/кг")
    
    # Легенда
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper right")
    
    ax2.set_title("результат проверки гипотезы об утечке радона")
    
    # Подпись снизу
    caption_text = "полное отношение модель/измерение при каждом f: "
    caption_text += ", ".join(["%0.3f" % r["ratio"] for r in runs])
    caption_text += "; активность Ra-226 растёт с %0.3f до %0.3f Бк/кг" % (ra_values[0], ra_values[-1])
    caption_text += " утечка компенсируется амплитудой и по этому спектру НЕ НАБЛЮДАЕМА, потому что все сильные линии цепочки (295, 352, 609, 1120, 1764) стоят ПОСЛЕ радона, а единственная линия ДО радона — Ra-226 186 кэВ — слаба и перекрыта с U-235 185,7"
    
    fig.text(0.5, 0.02, caption_text, ha="center", va="bottom", fontsize=10)
    
    # Сохранение
    save_path = os.path.normpath(os.path.join(_HERE, "..", "results", "figures", "bg_rn_leak_test.png"))
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=150)
    print(save_path)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
