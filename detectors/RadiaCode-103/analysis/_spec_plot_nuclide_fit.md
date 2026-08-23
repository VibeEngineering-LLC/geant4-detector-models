# SPEC: analysis/plot_nuclide_fit.py — картинка нуклидного разложения фона

Сгенерируй ОДИН файл Python 3 (matplotlib, Agg). Выведи ТОЛЬКО код, без пояснений.
Первая строка `# -*- coding: utf-8 -*-`. Комментарии — по-русски.

## Назначение

Рисует измеренный спектр открытого фона RadiaCode-103 и его нуклидное
разложение (fit_nuclides.py): полная модель + компоненты, сгруппированные по
цепочкам (K-40, цепочка Ra-226, цепочка Th-232, мюоны).

## Переиспользование (ОБЯЗАТЕЛЬНО импортом, ничего не переписывать)

```python
import os
import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import numpy as np
import fit_nuclides as fn
import fit_lines as fl
import read_rcxml
```

matplotlib импортировать так: `import matplotlib`, `matplotlib.use("Agg")`,
затем `import matplotlib.pyplot as plt`.

## Получение данных

1. Амплитуды берутся вызовом подгонки (она же напечатает свои таблицы в stdout —
   это нормально):

```python
amps, a_mu, filled, CHAIN_OF = fn.main()
```

2. Измерение (точно как в fn.main):

```python
smp = read_rcxml.read(fn.MEASURED)[0]
cnt = smp.counts[:-1].astype(float)
ch = np.arange(len(cnt))
e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))[:len(cnt)]
cps_meas = cnt / smp.live
```

3. Шаблоны и мюоны (те же функции, что использует фит):

```python
names, cols = fn.load_templates()
mu, pdg = fn.load_muons()
if mu is not None:
    names.append("mu")
    cols.append(mu)
A = np.zeros((len(e_meas), len(cols)))
for k, c in enumerate(cols):
    A[:, k] = fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)
```

4. Собрать словарь всех амплитуд: `aa = dict(amps)`, затем `aa.update(filled)`.

5. Кривые компонентов по группам. Для каждого `(n, a)` из `aa.items()`:
   пропустить если `n == "mu"` или `n not in names`; иначе имя группы
   `g = "K-40"` если `n == "K40"`, иначе `"цепочка Ra-226"` если
   `CHAIN_OF[n] == "Ra"`, иначе `"цепочка Th-232"`; накопить
   `grp[g] = grp.get(g, 0) + a * A[:, names.index(n)]` (нулевой вектор при
   первом обращении — использовать `np.zeros(len(e_meas))`).
   Мюоны: если `"mu" in names` и `"mu" in aa` — `grp["мюоны"] = aa["mu"] * A[:, names.index("mu")]`.

6. Полная модель: `model = sum(grp.values())` (поэлементная сумма векторов).

## Рисунок

Один ax, `figsize=(13.0, 7.2)`.

- Измерение: `ax.step(e_meas, cps_meas, where="mid", lw=1.0, color="#111111", label="измерение: <имя файла>, живое N сут")` — имя `os.path.basename(fn.MEASURED)`, живое `smp.live / 86400` с одним знаком после запятой.
- Полная модель: `ax.step(e_meas, model, where="mid", lw=1.4, color="#d62728", label="модель, сумма компонентов")`.
- Компоненты в фиксированном порядке и цветах:
  `"K-40"` — `#1f77b4`, `"цепочка Ra-226"` — `#2ca02c`, `"цепочка Th-232"` —
  `#9467bd`, `"мюоны"` — `#8c564b`; каждый `ax.step(..., lw=0.9)`, в подписи —
  имя группы и её доля в полной модели в диапазоне 20–2830 кэВ в процентах с
  одним знаком: маска `(e_meas >= 20) & (e_meas < 2830)`, доля =
  `grp[g][mask].sum() / model[mask].sum() * 100`.
- Оси: `ax.set_yscale("log")`, `ax.set_xlim(0, 3000)`; нижний предел Y —
  `max(1e-7, минимум положительных cps_meas в маске) / 3`, верхний —
  `max(cps_meas) * 3` через `ax.set_ylim`.
- Подписи: X `"Энергия, кэВ"`, Y `"скорость счёта, отсч/с на канал"`.
- Сетка: major `alpha=0.3, lw=0.6`, minor по X `alpha=0.12, lw=0.4`;
  `MultipleLocator(250)` major и `MultipleLocator(50)` minor по X
  (`from matplotlib.ticker import MultipleLocator`).
- Легенда: `fontsize=9, loc="upper right", framealpha=0.94`.
- Заголовок (`ax.set_title`, fontsize=12):
  `"Открытый фон RadiaCode-103: нуклидное шаблонное разложение (8 звеньев + мюоны)"`.
- Подзаголовок строкой `fig.text(0.5, 0.955, ..., fontsize=9.5, ha="center", color="#444444")`:
  `"амплитуды подобраны по нетто-площадям линий (не по полосам); звенья без своих линий достроены по равновесию цепочки"`.
- Внизу `fig.text(0.5, 0.01, ..., fontsize=8.2, ha="center", color="#666666")`:
  строка вида `"K-40 = X Бк/кг; Ra-226 = Y; Th-232 = Z (по Tl-208); мюоны x M от PDG. Шаблоны: <fn.BG_DIR>"`,
  где X — `amps.get("K40", 0)`, Y — `aa.get("Ra226", 0)`, Z — `aa.get("Tl208", 0)`,
  M — `a_mu` с двумя знаками; числа активностей с одним знаком после запятой.
- `fig.subplots_adjust(left=0.07, right=0.985, top=0.9, bottom=0.09)`.

## Сохранение

```python
OUT = os.path.normpath(os.path.join(_HERE, "..", "results", "figures",
                                    "bg_nuclide_decomposition.png"))
```

`os.makedirs(os.path.dirname(OUT), exist_ok=True)`, `fig.savefig(OUT, dpi=150)`,
`print("записано: %s" % OUT)`.

## Структура

Всё в `def main():` c `if __name__ == "__main__": sys.exit(main())`; main
возвращает 0. Докстрока модуля по-русски: что рисует, что амплитуды берёт из
fit_nuclides.main() импортом (единственный источник, никакого дублирования
подгонки), запуск `python analysis/plot_nuclide_fit.py`.
