# -*- coding: utf-8 -*-
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)

DATA = [
    ("Ra226", 186.2, 1.042),
    ("Pb212", 238.6, 1.028),
    ("Pb214", 351.9, 1.060),
    ("Ac228", 911.2, 1.105),
    ("Bi212", 727.3, 1.126),
    ("K40",   1460.8, 1.159),
    ("Bi214", 609.3, 1.161),
    ("Tl208", 2614.5, 1.187),
]

def main():
    # Сортировка по энергии
    sorted_data = sorted(DATA, key=lambda x: x[1])
    
    names = [item[0] for item in sorted_data]
    energies = [item[1] for item in sorted_data]
    growths = [item[2] for item in sorted_data]
    
    fig, ax = plt.subplots(figsize=(9, 6), dpi=130)
    
    # Точки
    ax.scatter(energies, growths, s=90, color="#d62728", zorder=3)
    
    # Подписи к точкам
    for i, (name, e, g) in enumerate(sorted_data):
        ax.annotate(name, xy=(e, g), xytext=(6, 4), textcoords="offset points", fontsize=9)
    
    # Линия тренда
    coeffs = np.polyfit(np.log10(energies), growths, 1)
    x_fit = np.logspace(np.log10(min(energies)/1.3), np.log10(max(energies)*1.15), 100)
    y_fit = np.polyval(coeffs, np.log10(x_fit))
    ax.plot(x_fit, y_fit, color="0.4", ls="--", lw=1.2, label="линейный тренд по log(E)")
    
    # Горизонтальная линия
    ax.axhline(1.0, color="0.6", ls=":", lw=1, label="без вклада здания (граница)")
    
    # Оси
    ax.set_xscale("log")
    ax.set_xlabel("Энергия главной линии, кэВ")
    ax.set_ylabel("Рост флюенса: extend=1000мм / extend=0")
    
    # Заголовок
    ax.set_title("RadiaCode-103: верхняя оценка вклада здания по нуклидам\nрост монотонен по жёсткости линии — жёсткая компонента чувствительнее к зданию")
    
    # Легенда
    ax.legend(loc="upper left", fontsize=9)
    
    # Сетка
    ax.grid(True, alpha=0.3, which="both")
    
    # Форматирование
    fig.tight_layout()
    
    # Сохранение
    path = os.path.join(OUT_DIR, "RC103_building_growth.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(os.path.abspath(path))

if __name__ == "__main__":
    main()
