# -*- coding: utf-8 -*-
"""
График сдвига равновесия внутри цепочек K/Ra/Th, интерпретированный через физику эманирования радона (Rn-222, T1/2=3.8 сут)
и торона (Rn-220, T1/2=55.6 с). Три из пяти пар статистически НЕ ОПРЕДЕЛЕНЫ (ошибка больше или сравнима с самим значением).
График честно показывает это широкими усами, а не скрывает. Единственная содержательная пара — Tl208/Ac228 (0.5σ от ожидания,
согласуется с физикой отсутствия утечки торона).
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)

RADON_COLOR = "#1f77b4"
THORON_COLOR = "#2ca02c"

PAIRS = [
    ("Bi214/Pb214",  1.272, 1.885, 1.0000, "радон", "оба ПОСЛЕ Rn-222 — тест внутренней согласованности, не утечки"),
    ("Pb214/Ra226",  1.679, 85.417, 1.0000, "радон", "Ra-226 ДО Rn-222 — прямой тест утечки радона, но Ra-226 сам не определяется"),
    ("Tl208/Ac228",  0.284, 0.168, 0.3594, "торон", "Ac-228 ДО Rn-220, Tl-208 ПОСЛЕ — тест утечки торона"),
    ("Ac228/Pb212",  1.088, 1.912, 1.0000, "торон", "Ac-228 ДО Rn-220, Pb-212 ПОСЛЕ — тест утечки торона"),
    ("Bi212/Pb212",  3.346, 7.542, 1.0000, "торон", "оба ПОСЛЕ Rn-220 — тест внутренней согласованности, не утечки"),
]

def main():
    fig, ax = plt.subplots(figsize=(10, 6), dpi=130)
    XMAX = 6.0

    for i, (label, ratio, error, expected, check_type, phys) in enumerate(PAIRS):
        y = len(PAIRS) - 1 - i
        color = RADON_COLOR if check_type == "радон" else THORON_COLOR
        
        nsig = abs(ratio - expected) / error if error > 0 else float("nan")
        
        # Отрисовка ошибки с обрезкой вниз до нуля
        lower_err = min(error, ratio)
        ax.errorbar([ratio], [y], xerr=[[lower_err], [error]], fmt="o", ms=9, color=color, capsize=4, zorder=3)
        
        # Ожидаемое значение
        ax.plot([expected], [y], marker="|", color="black", ms=20, mew=2, zorder=4)
        
        # Аннотация
        if math.isnan(nsig):
            text = "не определено"
        else:
            text = "%.1f σ — %s" % (nsig, phys)
        # Якорь аннотации прижат к правому краю видимой области (XMAX-0.15):
        # ratio+error у пар с огромной ошибкой (напр. 1.68+85.4) уходит далеко
        # за xlim и текст рисовался бы за пределами холста.
        anchor_x = min(ratio + error, XMAX - 0.15)
        ax.annotate(text, xy=(anchor_x, y), xytext=(8, 0), textcoords="offset points", fontsize=8, va="center")
    
    # Подписи осей. Точки рисуются с y=len(PAIRS)-1-i (первая пара сверху),
    # поэтому подписи тиков идут В ОБРАТНОМ порядке к PAIRS — иначе метки
    # съезжают на чужие строки.
    ax.set_yticks(range(len(PAIRS)))
    ax.set_yticklabels([pair[0] for pair in reversed(PAIRS)])
    ax.set_xlim(-1, XMAX)
    ax.axvline(1.0, color="0.6", ls=":", lw=1, zorder=1)
    
    # Легенда
    ax.scatter([], [], color=RADON_COLOR, label="проверка утечки радона (Rn-222, T½=3.8 сут — успевает уйти)")
    ax.scatter([], [], color=THORON_COLOR, label="проверка утечки торона (Rn-220, T½=55.6 с — НЕ успевает уйти)")
    # Внутри осей легенда перекрывала верхнюю строку (Bi214/Pb214) и её
    # длинную подпись справа — выносим под график.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=8, ncol=1)
    
    ax.set_xlabel("Отношение амплитуд цепочки")
    ax.set_title("RadiaCode-103: сдвиг равновесия внутри цепочек K/Ra/Th\nчёрная чёрточка — ожидание без утечки; далеко от неё вправо = недостаток дочерних (утечка); справа — значимость и что проверяет пара")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    
    path = os.path.join(OUT_DIR, "RC103_chain_equilibrium.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(os.path.abspath(path))

if __name__ == "__main__":
    main()
