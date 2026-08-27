Сгенерируй ОДИН файл Python 3 — `plot_building_growth.py`. Ответ: только код, без пояснений и без markdown-ограждения.

# Назначение

График роста флюенса ЕРН при продолжении здания за ограждениями комнаты (`extend=1000` мм против изолированной комнаты `extend=0`) — верхняя оценка вклада остального дома, по восьми нуклидам, отсортированным по энергии главной линии. Показывает, что рост монотонен по жёсткости — подтверждение гипотезы «жёсткая компонента растёт от здания сильнее мягкой».

# Данные (захардкодить как список кортежей, это уже посчитанный факт этой сессии)

```python
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
```
Каждый кортеж: (имя нуклида, энергия главной диагностической линии в кэВ, коэффициент роста flux_ext/flux_iso при extend=1000мм, для КОМБИНИРОВАННОЙ модели кирпич+бетон равной активности).

# Шапка

```python
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
```

# Логика

1. Отсортировать `DATA` по энергии (второй элемент кортежа) по возрастанию.
2. `fig, ax = plt.subplots(figsize=(9, 6), dpi=130)`.
3. Ось X — энергия линии (кэВ, логарифмическая шкала `ax.set_xscale("log")` — энергии от 186 до 2615 кэВ, больше десяти раз, линейная шкала сожмёт точки). Ось Y — коэффициент роста.
4. Точки: `ax.scatter(energies, growths, s=90, color="#d62728", zorder=3)`. Рядом с каждой точкой подписать имя нуклида текстом (`ax.annotate(name, xy=(e, g), xytext=(6, 4), textcoords="offset points", fontsize=9)`).
5. Линия тренда: полиномиальная регрессия 1-й степени по `log10(energies)` против `growths` (`np.polyfit(np.log10(energies), growths, 1)`), нарисовать сплошной линией по диапазону X (100 плавных точек `np.logspace(np.log10(min(energies))*0.95, np.log10(max(energies))*1.05, 100)` — НЕТ, правильно: `np.logspace(np.log10(min(energies)/1.3), np.log10(max(energies)*1.15), 100)`), цвет `"0.4"`, `ls="--"`, `lw=1.2`, подпись в легенде `"линейный тренд по log(E)"`.
6. Горизонтальная линия `ax.axhline(1.0, color="0.6", ls=":", lw=1, label="без вклада здания (граница)")`.
7. Подписи осей: `ax.set_xlabel("Энергия главной линии, кэВ")`; `ax.set_ylabel("Рост флюенса: extend=1000мм / extend=0")`.
8. Заголовок: `ax.set_title("RadiaCode-103: верхняя оценка вклада здания по нуклидам\\nрост монотонен по жёсткости линии — жёсткая компонента чувствительнее к зданию")`.
9. `ax.legend(loc="upper left", fontsize=9)`.
10. `ax.grid(True, alpha=0.3, which="both")`.
11. `fig.tight_layout()`.
12. Путь `os.path.join(OUT_DIR, "RC103_building_growth.png")`, `fig.savefig(path, bbox_inches="tight")`, `plt.close(fig)`, `print(os.path.abspath(path))`.

Обернуть в `def main():` и вызвать в `if __name__ == "__main__": main()`.

# Требования

Все подписи на русском. Никаких заглушек. Не открывать никаких CSV — данные захардкожены, это уже готовый посчитанный результат.
