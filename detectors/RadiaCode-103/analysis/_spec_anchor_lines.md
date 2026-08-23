# SPEC: analysis/anchor_lines.py — опорные пики 2614,5 / 1460,8 из моноотклика

Сгенерируй ОДИН файл Python 3. Выведи ТОЛЬКО код, без пояснений и без markdown-
обрамления (без ```). Первая строка `# -*- coding: utf-8 -*-`. Все комментарии
и докстроки — по-русски. Используй `math`, `os`, `sys`, `numpy as np`.

## Модульная докстрока (напечатать ДОСЛОВНО как есть ниже, это документация решения)

```
"""Опорные пики 2614,5 (Tl-208) и 1460,8 (K-40) — из МОНООТКЛИКА, метод 1.

ЗАЧЕМ. Полный шаблон звена (bg_bare_field_m1_Tl208.csv - весь распад Tl-208
через детектор) даёт нетто пика 2614 всего ~22 отсчёта (1,0 сигма) - шаблон
шумнее измерения (98 отсчётов, 4,4 сигма): статистика МК размазана по всему
каскаду, а не сосредоточена в интересующей линии. Отдельный МОНОэнергетический
прогон (resp_2614.csv, resp_1461.csv - 3e7 первичных на ТОЙ ЖЕ геометрии
источника, что и wallfield) даёт нетто на порядок надёжнее. Чтобы связать его
с активностью [Бк/кг], нужна доля потока ИМЕННО этой линии в полном флюенсе
звена (wf_m1_<N>.csv) - остальное снимается по тому же тождеству Ф=4N/S, что
и в fit_nuclides.load_templates().

ГЕОМЕТРИЯ - ПРОВЕРЕНА ФАКТОМ ПЕРЕД НОРМИРОВКОЙ (21.08, после /compact):
источник сетки откликов (_resp_1.mac, rc_curves): /gps/pos/type Surface,
shape Cylinder, radius 45.000 mm, halfz 82.500 mm, centre z=37.500 mm, ang cos.
Флюенс (fit_nuclides.CYL): r=45.0, z0=-45.0, z1=120.0 ->
  halfz = (120-(-45))/2 = 82.5 mm; centre z = (120+(-45))/2 = 37.5 mm.
СОВПАДАЕТ ЧИСЛО В ЧИСЛО - это одна и та же полость, розыгрыш откликов не
разъехался с розыгрышем флюенса.

ДВЕ РАЗНЫЕ ШИРИНЫ ЛИНИИ, две разные задачи:
- в СЕТКЕ ФЛЮЕНСА (wf_m1_*.csv) линия физически точечная - размыта только
  биннингом wallfield.cc (kBinKeV=2 кэВ), весь вклад в ОДНОМ бине. Окно
  rcspec.fwhm_win (детекторное разрешение) тут физически не к месту - своя
  маленькая функция line_flux_from_field() ниже (готовой в проекте не
  нашлось: line_net_area заточена под гауссов пик детектора).
- в ОТКЛИКЕ ДЕТЕКТОРА (resp_*.csv, свёрнутый rcspec.fold) и в ИЗМЕРЕНИИ -
  тот же line_net_area/fwhm_win, что уже применяется в fit_lines.py и
  fit_nuclides.py: чтобы формы пика в модели и в измерении отбирались
  ОДНИМ окном, а не двумя разными.

Запуск: python anchor_lines.py
"""
```

## Импорты и sys.path (ОБЯЗАТЕЛЬНО этим порядком, дословно)

```python
HERE = os.path.dirname(os.path.abspath(__file__))
for _d in ("analysis", "drivers"):
    _p = os.path.join(HERE, "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))

import rcspec
import read_rcxml
import fit_lines as fl
import fit_room_field as frf
import fit_nuclides as fn
```

## Константы

```python
RESP_DIR = os.path.join(fn.RESULTS, "bare", "response")

ANCHORS = [
    dict(nuc="Tl208", E0=2614.5, resp="resp_2614.csv", chain="Th-232"),
    dict(nuc="K40", E0=1460.8, resp="resp_1461.csv", chain="K-40"),
]
```

## Функция `line_flux_from_field(e, f, E0, gap=2, n_side=4)`

Docstring: "Чистый флюенс ОДНОЙ линии в родном 2-кэВ биннинге wallfield -
континуум оценивается СРЕДНИМ по n_side бинам с каждой стороны, отступив gap
бинов от пика. Возвращает (net, continuum_per_bin, gross, idx)."

Тело:
```python
e = np.asarray(e)
idx = int(np.argmin(np.abs(e - E0)))
lo = f[max(0, idx - gap - n_side):max(0, idx - gap)]
hi = f[idx + gap + 1:idx + gap + 1 + n_side]
side = np.concatenate([lo, hi])
cont = float(side.mean()) if side.size else 0.0
gross = float(f[idx])
return gross - cont, cont, gross, idx
```

## Функция `response_cps_per_bqkg(nuc, E0, resp_csv)`

Docstring: "Отклик детектора на 1 Бк/кг родителя ЧЕРЕЗ моноотклик линии E0,
по тому же тождеству Ф=4N/S, что fit_nuclides.load_templates()."

Тело (площадь поверхности ТОЙ ЖЕ полости, что и в fit_nuclides.load_templates —
переиспользовать fn.CYL, не заводить новых чисел):
```python
r, hz = fn.CYL["r"] / 10.0, 0.5 * (fn.CYL["z1"] - fn.CYL["z0"]) / 10.0
area = 2 * math.pi * r * (r + 2 * hz)

wf_path = os.path.join(fn.BUILD, "%s_%s.csv" % (fn.WF_PREFIX, nuc))
e_f, f_f = frf.read_wallfield(wf_path)
net_flux, cont, gross, idx = line_flux_from_field(e_f, f_f, E0)
if net_flux <= 0:
    raise SystemExit("line_flux_from_field: net<=0 для %s@%.1f (gross=%.3e cont=%.3e)"
                      % (nuc, E0, gross, cont))
rate = net_flux * area / 4.0

resp_path = os.path.join(RESP_DIR, resp_csv)
meta, hist = rcspec.read_spec(resp_path)
n_primaries = float(meta["N_primaries"])
t_run = n_primaries / rate

folded = rcspec.fold(hist, "103")
e_ch = np.arange(len(folded)) + 0.5
r_line = fl.line_net_area(e_ch, folded, E0)
if r_line is None or r_line["area"] <= 0:
    raise SystemExit("line_net_area(response): площадь <=0 для %s@%.1f" % (nuc, E0))
cps = r_line["area"] / t_run
return dict(cps=cps, r_line=r_line, net_flux=net_flux, cont_flux=cont,
            gross_flux=gross, rate=rate, t_run=t_run, n_primaries=n_primaries,
            area_cm2=area)
```

## Функция `measured_cps(E0)`

```python
smp = read_rcxml.read(fn.MEASURED)[0]
ch = np.arange(len(smp.counts))
e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))
r = fl.line_net_area(e_meas, smp.counts.astype(float), E0)
if r is None:
    return None
return dict(cps=r["area"] / smp.live, sd=r["sd"] / smp.live, live=smp.live,
            gross=r["gross"], cont=r["cont"], lo=r["lo"], hi=r["hi"])
```

## Функция `main()`

Печатает заголовок `"=== Опорные пики из моноотклика (метод 1) — Tl-208 2614,5 / K-40 1460,8 ===\n"`.

Заводит `results = {}`. Для каждого `a` из `ANCHORS`:

1. `resp = response_cps_per_bqkg(a["nuc"], a["E0"], a["resp"])`
2. `meas = measured_cps(a["E0"])`
3. Печатает `"--- %s (%s), E0=%.1f кэВ ---" % (a["chain"], a["nuc"], a["E0"])`.
4. Печатает долю линии от полного флюенса звена: прочитать ещё раз
   `frf.read_wallfield(os.path.join(fn.BUILD, "%s_%s.csv" % (fn.WF_PREFIX, a["nuc"])))`,
   взять `total = sum(его второй массив)`, напечатать
   `"  флюенс: gross_bin %.4e  cont/bin %.4e  net %.4e  (доля %.2f%% от полного)"`
   с `resp["gross_flux"], resp["cont_flux"], resp["net_flux"], 100.0*resp["net_flux"]/total`.
5. Печатает `"  моноотклик: N_primaries=%.0f  сеточный t_run=%.4e с/(Бк/кг)"` с `resp["n_primaries"], resp["t_run"]`.
6. Печатает `"  площадь пика в отклике: %.1f ± %.1f (окно %.0f-%.0f кэВ)"` с
   `resp["r_line"]["area"], resp["r_line"]["sd"], resp["r_line"]["lo"], resp["r_line"]["hi"]`.
7. Печатает `"  => модельный отклик %.6e cps на 1 Бк/кг %s" % (resp["cps"], a["chain"])`.
8. Если `meas is None`: печатает `"  ИЗМЕРЕНИЕ: окно вне спектра — линия не оценена\n"` и `continue` к следующему анкору.
9. Печатает `"  измерено: gross=%.1f cont=%.1f net cps=%.6f ± %.6f (окно %.0f-%.0f кэВ, live=%.0f с)"` с
   `meas["gross"], meas["cont"], meas["cps"], meas["sd"], meas["lo"], meas["hi"], meas["live"]`.
10. Если `meas["cps"] <= 0`: печатает `"  измеренная площадь <=0 — активность не определяется этой линией\n"` и `continue`.
11. Иначе: `act = meas["cps"] / resp["cps"]`; `sd_rel = meas["sd"] / meas["cps"]`;
    `sd = abs(act) * sd_rel` (комментарий рядом: "ошибка МК-статистики отклика
    пока не добавлена — см. DECISIONS.md"); печатает
    `"  => %s = %.1f ± %.1f Бк/кг (только по этой линии, ошибка ПОКА без вклада МК)\n" % (a["chain"], act, sd)`;
    кладёт `results[a["chain"]] = dict(act=act, sd=sd)`.

В конце `return results`.

## Точка входа

```python
if __name__ == "__main__":
    main()
```
