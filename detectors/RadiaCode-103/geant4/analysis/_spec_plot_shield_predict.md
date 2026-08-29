Сгенерируй ОДИН файл Python 3 — `plot_shield_predict.py`. Ответ: только код, без пояснений и без markdown-ограждения.

# Назначение

Нарисовать проверку модели ослабления фона свинцовым домиком: измеренный спектр
в домике против предсказания, построенного на активностях из открытого фона БЕЗ
подгонки. Показать раздельно вклад гамма-фона комнаты и космической компоненты —
именно их соотношение и есть содержание рисунка.

# Шапка

`# -*- coding: utf-8 -*-` и docstring: что рисуется, откуда берутся данные
(импорт `predict_shield`, повторной реализации расчёта нет), и предупреждение,
что легенда обязана лежать ВНЕ поля данных — в этой линии уже был дефект, когда
непрозрачная легенда перекрыла спектр и график читался как обрыв данных.

# Импорты и настройки

```python
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc
import fit_physical_chains as fpc
import predict_shield as ps

matplotlib.rcParams["font.size"] = 9
matplotlib.rcParams["figure.dpi"] = 130
matplotlib.rcParams["axes.grid"] = True
matplotlib.rcParams["grid.alpha"] = 0.3

OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PNG = os.path.join(OUT_DIR, "RC103_shield_predict.png")
```

Цвета — светлые и НЕПРОЗРАЧНЫЕ (alpha не использовать нигде, ни в заливках, ни
в линиях; в этой линии оператор прямо потребовал отказаться от прозрачности):

```python
C_MEAS_OPEN = "#9aa0a6"   # открытый фон, измерение — фон-контекст
C_MEAS = "#111111"        # измерение в домике
C_MODEL = "#d62728"       # полная модель
C_GAMMA = "#4e79a7"       # гамма-компонента (заливка)
C_MUON = "#8c8c8c"        # космическая компонента (заливка)
```

# Функция `build()`

Возвращает словарь со всем, что нужно для рисования. Повторно использует функции
модуля `ps` — своей реализации чтения, свёртки и подгонки НЕ писать.

1. Открытый фон: `cnt_o, e_o, live_o, _ = ps.read_meas(ftc.MEAS_NAME, ftc.CAL_ROOM)`.
2. Шаблоны открытые: `cps_o, var_o, hm_o = ps.load_cps(ps.FMT_OPEN, ftc.MUON_CSV)`.
   Если `hm_o` ложно — `sys.exit("нет мюонного шаблона для открытого фона")`.
3. `names, A_o, V_o = ps.columns_on_grid(cps_o, var_o, e_o, hm_o)`;
   `sel = (e_o >= ftc.E_LO) & (e_o < ftc.E_HI)`; применить ко всем.
4. Подгонка: веса `1.0/np.sqrt(np.maximum(y_o, 1.0))`, вызов
   `ftc.fit(A_o*live_o, y_o, w, names, "", "", var_counts=V_o*live_o*live_o)`.
   Печать подгонки подавить: обернуть в
   `io.StringIO()` + `contextlib.redirect_stdout` (импортировать `io`, `contextlib`).
   Взять только `amp` и `pred_o`.
5. Домик: `cnt_s, e_s, live_s, _ = ps.read_meas(ps.MEAS_SHIELD, None)` — калибровка
   приборная, из файла.
6. Шаблоны с домиком: `cps_s, var_s, hm_s = ps.load_cps(ps.FMT_SHIELD, ps.MUON_SHIELD_CSV)`.
   Если `hm_s` ложно — `sys.exit("нет мюонного шаблона С ДОМИКОМ")`: рисовать
   разложение без космики бессмысленно, это и есть предмет рисунка.
7. `names_s, A_s, V_s = ps.columns_on_grid(cps_s, var_s, e_s, hm_s)`, отбор `sel_s`.
8. Перенос амплитуд ПО ИМЕНАМ через `dict(zip(names, amp))`; отсутствие имени —
   `sys.exit` с этим именем.
9. **Раздельные вклады.** Индекс мюонного столбца: `i_mu = names_s.index("mu")`.
   - `pred_mu = A_s[:, i_mu] * amp_s[i_mu]` — космическая компонента, имп/с;
   - `pred_gamma = A_s @ amp_s - pred_mu` — всё остальное, имп/с;
   - `pred_tot = pred_gamma + pred_mu`.
   Все три уже в имп/с, домножать на живое время НЕ надо.
10. Вернуть словарь с ключами: `e_o`, `meas_open` (= `y_o/live_o`), `e_s`,
    `meas_shield` (= `y_s/live_s`), `pred_tot`, `pred_gamma`, `pred_mu`,
    `amp_map`, `live_o`, `live_s`.

# Функция `draw(d)`

Фигура `plt.subplots(2, 1, figsize=(11, 9), gridspec_kw={"height_ratios": [3, 1]})`.

**Верхняя панель — спектры.**
- Открытый фон: `ax.step(d["e_o"], d["meas_open"], where="mid", color=C_MEAS_OPEN, lw=1.0, label="измерено, открытый фон")`.
- Заливки компонент, СТОПКОЙ снизу вверх, без прозрачности и БЕЗ контурных линий
  (`linewidth=0`): сначала `fill_between(e_s, 0, pred_mu, color=C_MUON, label="модель: космическая компонента")`,
  затем `fill_between(e_s, pred_mu, pred_mu + pred_gamma, color=C_GAMMA, label="модель: гамма-фон комнаты")`.
  Космика идёт ПЕРВОЙ (на заднем плане) — так требовал оператор в этой линии.
- Полная модель: линия `plot(e_s, pred_tot, color=C_MODEL, lw=1.4, label="модель, сумма")`.
- Измерение в домике: `step(e_s, meas_shield, where="mid", color=C_MEAS, lw=1.0, label="измерено, домик")`.
- Ось Y: `ax.set_yscale("log")`, пределы — от `max(1e-7, минимум положительных значений измерения в домике)`
  до `1.6 * максимум по открытому фону`. Верх спектра не должен упираться в рамку.
- Ось X: `set_xlim(ftc.E_LO, ftc.E_HI)`, подпись «энергия, кэВ».
- Подпись Y: «скорость счёта, имп/с на канал».
- Заголовок: «RadiaCode-103: ослабление фона свинцовым домиком — предсказание без подгонки».
- Легенда — ВНЕ поля данных: `ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=3, frameon=False)`.
  ⚠ Тогда подпись оси X перекроется легендой — поэтому подпись оси X ставить только
  на НИЖНЕЙ панели, а на верхней вызвать `ax.set_xlabel("")`.

**Нижняя панель — ослабление по полосам.**
Для каждой полосы из `ftc.BANDS` посчитать два числа:
- измеренное ослабление = сумма `meas_open` по полосе / сумма `meas_shield` по полосе;
- модельное ослабление = сумма модели открытого фона по полосе / сумма `pred_tot` по полосе.
  Модель открытого фона взять из `build()` — добавить её в словарь ключом `pred_open`
  (это `pred_o/live_o`, тоже имп/с).
⚠ Полосы у двух измерений лежат на РАЗНЫХ энергетических сетках (`e_o` и `e_s`) —
маску для каждой суммы строить по своей сетке. Это не приближение: сравниваются
интегралы по одинаковым интервалам энергии.

Рисовать столбиками: `ax.bar(x - 0.2, meas_att, width=0.4, color=C_MEAS, label="измерено")`
и `ax.bar(x + 0.2, model_att, width=0.4, color=C_MODEL, label="модель")`, где
`x = np.arange(len(BANDS))`. Подписи делений — «lo–hi» через `set_xticks`/`set_xticklabels`.
Ось Y логарифмическая, подпись «ослабление, раз». Легенда внутри верхнего левого
угла этой панели допустима — данные там столбиками у оси, не перекроются;
`frameon=False`.
Над каждым столбиком подписать его значение форматом `%.1f`, кегль 8.

`fig.tight_layout()`, `fig.savefig(OUT_PNG)`, печать полного пути в stdout.

# main

```python
if __name__ == "__main__":
    draw(build())
```

# Стиль

Комментарии по-русски. Никаких магических чисел без пояснения. Прозрачность не
использовать вообще.
