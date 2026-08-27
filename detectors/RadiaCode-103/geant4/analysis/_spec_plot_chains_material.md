Сгенерируй ОДИН файл Python 3 — `plot_chains_material.py`. Ответ: только код, без пояснений и без markdown-ограждения.

# Назначение

Рисует разложение фона по цепочкам (K40/Ra226/Th232) и раздельным материалам (кирпич/бетон) + мюоны. Данные и подгонка берутся импортом из модуля `fit_chains_material` (функция `prepare()`), логика подгонки НЕ дублируется.

Docstring: смысл картинки — показать, почему разделение по материалам статистически не определяется (кривые «кирпич» и «бетон» одной цепочки идут почти по одной линии — cond=517.9 против 246.4 у общей модели без разделения материалов, это уже установлено фактом в fit_chains_material).

# Шапка

```python
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
```

# Функция `plot_panel(ax, r, which, subtitle)`

`r` — словарь, возвращённый `fcm.prepare()`. Ключи (используй ровно эти имена): `e_meas` (массив энергий, кэВ), `cnt` (измеренные отсчёты в канале), `live` (живое время, с), `A_counts` (матрица [канал x звено] в отсчётах на 1 Бк/кг или на 1 мюон/с), `names` (список имён столбцов A_counts по порядку — либо вида `"K40_brick"`/`"Ra226_concrete"`/`"Th232_brick"` и т.п., либо `"mu"`), `amp_a`/`amp_b` (массивы амплитуд по критерию A/B, тот же порядок, что `names`), `pred_a`/`pred_b` (предсказанный спектр в отсчётах, сумма модели), `chi2ndf_a`/`chi2ndf_b`, `shape_a`/`shape_b` (числа).

`which` — строка `"A"` или `"B"`, выбирает соответствующий набор (`amp_a`/`pred_a`/... либо `amp_b`/`pred_b`/...).

Логика:
1. `meas_cps = r["cnt"] / r["live"]`.
2. `ax.step(r["e_meas"], meas_cps, where="mid", color="0.35", lw=1.0, label="измерение")`.
3. `pred = r["pred_a"] if which=="A" else r["pred_b"]`; `ax.plot(r["e_meas"], pred / r["live"], color="crimson", lw=1.6, label="сумма модели")`.
4. `amp = r["amp_a"] if which=="A" else r["amp_b"]`.
5. Для каждого `k, name` в `enumerate(r["names"])`: если `amp[k] <= 0` — добавить `name` в список `zero` (обнулённые) и пропустить отрисовку. Иначе:
   - если `name == "mu"`: нарисовать `ax.plot(r["e_meas"], amp[k]*r["A_counts"][:,k]/r["live"], color="0.5", lw=0.9, ls=":", alpha=0.9, label="mu %.1f мюон/с" % amp[k])`;
   - иначе `chain, material = name.rsplit("_", 1)` (имя вида `"K40_brick"` -> `chain="K40"`, `material="brick"`); нарисовать `ax.plot(r["e_meas"], amp[k]*r["A_counts"][:,k]/r["live"], color=CHAIN_COLOR[chain], lw=1.1, ls=MATERIAL_LS[material], alpha=0.85, label="%s %.1f Бк/кг" % (name, amp[k]))`.
6. `ax.set_yscale("log")`.
7. `ax.set_xlim(fcm.ftc.E_LO, fcm.ftc.E_HI)`.
8. `ax.set_ylim(max(meas_cps.max()*1e-5, 1e-7), meas_cps.max()*2)`.
9. `ax.set_xlabel("Энергия, кэВ")`; `ax.set_ylabel("Скорость счёта, 1/(с·канал)")`.
10. Заголовок: `chi2ndf = r["chi2ndf_a"] if which=="A" else r["chi2ndf_b"]`; `shape = r["shape_a"] if which=="A" else r["shape_b"]`; строка `"Критерий %s: chi2/ndf=%.1f, невязка формы=%.4f\n%s" % (which, chi2ndf, shape, subtitle)`; если список `zero` непустой — добавить строку `"\nобнулены NNLS: %s" % ", ".join(zero)`. Установить через `ax.set_title(title, fontsize=10)`.
11. Легенда СТРОГО вне поля данных (в этом проекте легенда внутри осей уже перекрывала данные — не повторять эту ошибку): `ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7.5, framealpha=1.0)`.

# Функция `main()`

1. `r = fcm.prepare()`.
2. `fig = plt.figure(figsize=(13, 11.5))`.
3. `gs = fig.add_gridspec(2, hspace=0.55, height_ratios=[1, 1])` — hspace большой НАМЕРЕННО, иначе заголовок нижней панели наезжает на подпись оси X верхней (уже было такое в этом проекте).
4. `ax1 = fig.add_subplot(gs[0, 0])`; вызвать `plot_panel(ax1, r, "A", "цепочки K40/Ra226/Th232 x кирпич/бетон (сплошная=кирпич, пунктир=бетон)")`.
5. `ax2 = fig.add_subplot(gs[1, 0])`; вызвать `plot_panel(ax2, r, "B", "то же, критерий формы — кривые кирпич/бетон почти совпадают: материал НЕ разделяется по спектру")`.
6. `fig.suptitle("RadiaCode-103, фон помещения: разложение по цепочкам и РАЗДЕЛЬНЫМ материалам (кирпич/бетон), метод 1", y=0.995)`.
7. НЕ вызывать `fig.tight_layout()` — он съедает заданный `hspace` и возвращает наложение заголовков (уже проверено на этом проекте).
8. `path = os.path.join(OUT_DIR, "RC103_bg_decomposition_chains_material.png")`; `fig.savefig(path, bbox_inches="tight")`; `plt.close(fig)`; `print(os.path.abspath(path))`.

В конце — `if __name__ == "__main__": main()`.

# Требования

Комментарии — только там, где логика неочевидна (разбор имени `name.rsplit("_", 1)`, почему `tight_layout()` не вызывается). Никаких заглушек и несуществующих функций. Не переопределять НИЧЕГО, что уже есть в `fit_chains_material.py` — только импортировать и использовать возвращённый словарь.
