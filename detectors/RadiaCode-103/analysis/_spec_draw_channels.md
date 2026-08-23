# SPEC: analysis/draw_channels.py для RadiaCode-103

Сгенерируй ОДИН файл Python 3 (matplotlib, Agg). Выведи ТОЛЬКО код, без пояснений.
Комментарии в коде — на русском языке. Кодировка файла UTF-8, первая строка
`# -*- coding: utf-8 -*-`.

## Назначение

Рисует функцию отклика детектора RadiaCode-103, разложенную по каналам
взаимодействия. Данные лежат в файлах `*_chan.csv`, по одному на энергию
падающего кванта.

## Переиспользование чужого модуля (ОБЯЗАТЕЛЬНО, не переписывать своими руками)

Ядро берётся из модуля-донора по абсолютному пути. Вставь ровно так:

```python
_HERE = os.path.dirname(os.path.abspath(__file__))
_DONOR = os.path.normpath(os.path.join(
    _HERE, "..", "..", "AtomSpectra-Nano-16-PRO", "analysis"))
sys.path.insert(0, _DONOR)
import draw_channels as dc
sys.path.insert(0, _HERE)
import rcspec
```

Из донора `dc` используются: `dc.read_chan(path)`, `dc.broaden(pairs, w, nch, step)`,
`dc.ORDER`, и подменяемая функция `dc.fwhm`.

Сразу после импортов подмени ширину линии на ширину ЭТОГО прибора:

```python
dc.fwhm = lambda e: float(rcspec.fwhm(max(e, 1.0)))
```

Комментарием пояснить: подмена зависимости, а не копия формулы — `dc.broaden`
зовёт `dc.fwhm`, и после подмены донорская свёртка работает с разрешением
RadiaCode-103; коэффициенты живут в `rcspec` и уточняются там.

## Константы модуля

```python
SRC = os.path.normpath(os.path.join(_HERE, "..", "results", "m200", "response"))
OUT = os.path.normpath(os.path.join(_HERE, "..", "results", "figures",
                                    "rc103_channels.png"))
E_MAX = 3200.0
```

Функция `step_for(e0)`: возвращает `max(5.0, round(dc.fwhm(e0) / 4.0 / 5.0) * 5.0)`.
Комментарий: канал отображения привязан к ПШПВ, четыре канала на полуширину.

## Формат входного файла `*_chan.csv`

`dc.read_chan(path)` возвращает кортеж `(head, names, rows)`:
- `head` — словарь строк-заголовков вида `# ключ = значение`; нужны ключи
  `E_prim_keV`, `N_primaries`, `src_sha1`;
- `names` — список имён каналов (столбцы кроме первого);
- `rows` — список пар `(энергия_кэВ, [числа по каналам])`.

ВАЖНО: у файлов этого прибора НЕТ ключа `solid_angle_frac`, который есть у
донора. Ничего с этим делать не нужно — вес считается иначе, см. ниже.

## Функция main()

1. `wanted = [float(x) for x in sys.argv[1:]] or [150.0, 662.0, 1461.0, 2614.0]`
2. Собрать список файлов: те из `os.listdir(SRC)`, что заканчиваются на
   `_chan.csv`, отсортированные. Пусто — `raise SystemExit("в %s нет файлов *_chan.csv" % SRC)`.
3. Построить словарь `index`: ключ `float(head["E_prim_keV"])`, значение
   `(head, names, rows)`.
4. `dc.E_MAX = E_MAX`
5. Фигура: `fig, axes = plt.subplots(len(wanted), 1, figsize=(12.6, 4.0 * len(wanted)), sharex=False)`;
   если `len(wanted) == 1`, обернуть `axes` в список.
6. Для каждой пары `(ax, e_want)`:
   - `e0 = min(index, key=lambda x: abs(x - e_want))`, взять `head, names, rows`;
   - `stamp = head.get("src_sha1", "?")`;
   - вес `w = 1.0 / float(head["N_primaries"])` с комментарием: вес на ОДИН
     квант, вошедший внутрь цилиндра поля; телесного угла здесь нет, источник
     уже задан поверхностью, охватывающей прибор;
   - `step = step_for(e0)`, `nch = int(E_MAX / step)`,
     `cols = [(k + 0.5) * step for k in range(nch)]`;
   - `grand = sum(sum(v) for _, v in rows)`;
   - полный отклик: `total = dc.broaden([(e, sum(v)) for e, v in rows], w, nch, step)`,
     нарисовать `ax.step(cols, total, where="mid", lw=1.6, color="#111111",
     label="полный отклик", zorder=5)`;
   - по каждому `(name, label, colour)` из `dc.ORDER`: если `name` есть в
     `names` — взять индекс `j`, посчитать `raw = sum(v[j] for _, v in rows)`,
     при `raw == 0` пропустить, иначе `cur = dc.broaden([(e, v[j]) for e, v in rows], w, nch, step)`
     и `ax.step(cols, cur, where="mid", lw=1.0, color=colour,
     label="%s — %.1f %%" % (label, 100.0 * raw / grand))`;
   - оформление: `ax.set_yscale("log")`; `ax.set_xlim(0, min(E_MAX, e0 * 1.15))`;
     `top = max(total)`; `ax.set_ylim(top / 2.0e3, top * 40)`;
     `ax.set_ylabel("вероятность на квант поля")`;
     major locator 250, minor locator 50 (`MultipleLocator`);
     сетка major `alpha=0.26, lw=0.6`, сетка minor по оси x `alpha=0.11, lw=0.4`;
     легенда `fontsize=7.8, loc="upper right", ncol=2, framealpha=0.94`;
     заголовок `ax.set_title("падающий квант %.0f кэВ, канал отображения %.0f кэВ, событий с сигналом %d" % (e0, step, grand), fontsize=10.5, pad=4)`.
7. `axes[-1].set_xlabel("Энерговыделение, кэВ")`
8. Общий заголовок:
   `fig.suptitle("RadiaCode-103: функция отклика, разложенная по каналам взаимодействия\nгеометрия фона БЕЗ домика (поле комнаты, цилиндрическая поверхность, сосуд m200 с воздухом)", fontsize=12, y=0.995)`
9. Подпись внизу через `fig.text(0.5, 0.004, ..., fontsize=8.4, ha="center", color="#555555")`
   с текстом: `"Свёрнуто с ПШПВ(E) прибора из rcspec.fwhm. Штамп исходников %s. Проценты в легенде — доля канала во ВСЕХ событиях с сигналом на этом узле. Правило приоритета каналов общее с Nano16." % stamp`
10. `fig.subplots_adjust(left=0.075, right=0.985, top=0.925, bottom=0.062, hspace=0.16)`
11. `os.makedirs(os.path.dirname(OUT), exist_ok=True)`, `fig.savefig(OUT, dpi=150)`,
    `print("записано: %s" % OUT)`, `return 0`
12. В конце файла: `if __name__ == "__main__": sys.exit(main())`

## Докстрока модуля

Многострочная, по-русски: что канал ставится в момент события по истории
процессов (`geometry/main.cc`, enum `Chan`) и из готового спектра не
восстанавливается; что чтение, свёртка и порядок каналов заимствованы у донора
Nano16 импортом, потому что правило приоритета у обоих приборов одно; что своё
здесь — ширина линии из `rcspec`, нормировка на квант поля и подписи.
Последняя строка докстроки — способ запуска:
`python analysis/draw_channels.py [E1 E2 ...]`

## Импорты

`io` не нужен. Нужны: `os`, `sys`, `matplotlib` с `matplotlib.use("Agg")` ДО
`import matplotlib.pyplot as plt`, и `from matplotlib.ticker import MultipleLocator`.
Порядок: сначала стандартные, потом matplotlib, потом блок донора и `rcspec`
(с комментарием `# noqa: E402`).
