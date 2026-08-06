# -*- coding: utf-8 -*-
"""Насколько результат зависит от самопоглощения в вольфраме.

Вольфрам ослабляет сильно на ВСЕХ линиях ряда, а не только на мягких, поэтому
активность, вынутая из спектра, наследует все неточности геометрии стержня.
Здесь это меряется: считаются
  * эффективная глубина съёма 1/mu — на сколько миллиметров прибор вообще видит
    вглубь стержня на каждой линии;
  * выход в окне линии при разном диаметре стержня (3,20 мм по этикетке против
    2,40 мм, заявленных на одной из фотографий пачки);
  * во сколько раз изменится вынутая активность при подмене диаметра — это и
    есть цена ошибки в геометрии.

Перенос считает `wt20_source_scatter` (фотоэффект, Клейн–Нишина, когерентное);
здесь только перебор вариантов и сведение в таблицу.

    python analysis/wt20_source_sensitivity.py [каталог вывода]
"""
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import wt20_source_scatter as S      # noqa: E402

VARIANTS = [
    ("этикетка", 0.320, 0.485),
    ("фото 2,4 мм", 0.240, 0.485),
]
LINES = [238.63, 300.09, 583.19, 911.20, 2614.51]


def run(diam_cm, pitch_cm):
    S.ROD_R = diam_cm / 2.0
    S.PITCH = pitch_cm
    S.CENTRES = (np.arange(S.N_RODS) - (S.N_RODS - 1) / 2.0) * pitch_cm
    S.RNG = np.random.default_rng(20260807)
    return {e: S.trace_line(e, n=120000) for e in LINES}


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(outdir, exist_ok=True)

    print("=== ЭФФЕКТИВНАЯ ГЛУБИНА СЪЁМА В СПЛАВЕ ===")
    print("  радиус стержня по этикетке 1,60 мм")
    print("  E, кэВ    mu, 1/см   1/mu, мм   доля радиуса")
    for e in LINES:
        mu = (S.mu_of(np.array([e / 1000.0]), S.COH)
              + S.mu_of(np.array([e / 1000.0]), S.INC)
              + S.mu_of(np.array([e / 1000.0]), S.PE))[0]
        d = 10.0 / mu
        print("  %7.2f  %9.3f  %9.2f      %5.2f" % (e, mu, d, d / 1.60))

    res = {}
    for name, d, p in VARIANTS:
        res[name] = run(d, p)

    print()
    print("=== ВЫХОД В ОКНЕ ЛИНИИ, ВВЕРХ (сторона прибора) ===")
    print("  E, кэВ   ⌀3,20 мм   ⌀2,40 мм   отношение   вынутая активность")
    rows = []
    for e in LINES:
        a = res["этикетка"][e]["f_peak_up"]
        b = res["фото 2,4 мм"][e]["f_peak_up"]
        # активность обратно пропорциональна выходу и массе: A ~ 1/(f * m),
        # масса пропорциональна квадрату диаметра
        m_ratio = (0.240 / 0.320) ** 2
        act = (a / b) * (1.0 / m_ratio)
        rows.append((e, a, b, b / a, act))
        print("  %7.2f   %8.4f   %8.4f    %7.3f      x %.2f"
              % (e, a, b, b / a, act))

    print()
    print("  Последняя колонка: во сколько раз изменилась бы УДЕЛЬНАЯ активность,")
    print("  если бы стержень оказался ⌀2,40 мм, а не 3,20. Разброс по линиям")
    print("  показывает, что подмена диаметра ломает не только масштаб, но и")
    print("  согласие линий между собой.")

    p = os.path.join(outdir, "wt20_source_sensitivity.csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        f.write("# чувствительность к диаметру стержня, перенос в источнике\n")
        f.write("E_кэВ;выход_окно_вверх_d3.20;выход_окно_вверх_d2.40;"
                "отношение_выходов;множитель_удельной_активности\n")
        for r in rows:
            f.write("%.2f;%.5f;%.5f;%.4f;%.4f\n" % r)
    print()
    print("записано:", p)


if __name__ == "__main__":
    main()
