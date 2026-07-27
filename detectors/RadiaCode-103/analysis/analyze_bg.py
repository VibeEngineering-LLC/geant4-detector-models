# -*- coding: utf-8 -*-
"""Фон и поправка на экранирование фона пробой.

Даёт:
  1) скорость счёта фона в пустом сосуде (сверяется с показаниями прибора);
  2) поканальный множитель k = спектр(проба) / спектр(пустой сосуд);
  3) величину переподчёта -B*(1-k) в областях интереса ключевых линий.
"""
import math
import os
import re

import numpy as np

import sys
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec
import run_bg

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = rcspec.RESULTS
BG = rcspec.rdir("background")

# Площадь охватывающего цилиндра своя у каждого сосуда — берётся из драйвера,
# чтобы нормировка Ф = 4N/S не разъехалась с тем, чем считали.
CYL_S = run_bg.cyl_area(rcspec.vessel())
SPH_S = run_bg.SPH_S

LINES = [(238.6, "Pb-212"), (295.2, "Pb-214"), (351.9, "Pb-214"),
         (583.2, "Tl-208"), (609.3, "Bi-214"), (661.7, "Cs-137"),
         (911.2, "Ac-228"), (1120.3, "Bi-214"), (1460.8, "K-40"),
         (1764.5, "Bi-214"), (2614.5, "Tl-208")]


def field_fluence():
    p = os.path.join(RESULTS, "field_spectrum.mac")
    for line in open(p, encoding="utf-8"):
        m = re.match(r"#\s*FLUENCE_TOTAL_CM2_S\s*=\s*([\d.eE+-]+)", line)
        if m:
            return float(m.group(1))
    raise SystemExit("не нашёл нормировку флюенса в " + p)


def load(shape, name):
    p = os.path.join(BG, "bg_%s_%s.csv" % (shape, name))
    if not os.path.exists(p):
        return None
    meta, hist = rcspec.read_spec(p)
    return meta, hist


def main():
    phi = field_fluence()
    print("флюенс поля: %.3f см⁻²·с⁻¹" % phi)

    files = sorted(f for f in os.listdir(BG) if f.startswith("bg_cyl_")) \
        if os.path.isdir(BG) else []
    if not files:
        raise SystemExit("нет результатов в " + BG)

    specs = {}
    for fn in files:
        name = fn[len("bg_cyl_"):-len(".csv")]
        meta, hist = rcspec.read_spec(os.path.join(BG, fn))
        n = float(meta["N_primaries"])
        # Ф = 4N/S  =>  реальный поток через поверхность = Ф*S/4 частиц/с
        rate = phi * CYL_S / 4.0
        specs[name] = dict(meta=meta, hist=hist, n=n, cps=hist.sum() / n * rate,
                           rho=float(meta.get("density_gcm3", 0)),
                           spec_per_s=hist / n * rate)
        print("%-14s  плотность %.2f  фон %.2f имп/с (>0 кэВ), выше 20 кэВ %.2f"
              % (name, specs[name]["rho"], specs[name]["cps"],
                 specs[name]["spec_per_s"][20:].sum()))

    # проверка нормировки: сфера должна дать тот же счёт
    sph = load("sph", "air_0.00")
    if sph:
        meta, hist = sph
        n = float(meta["N_primaries"])
        cps = hist.sum() / n * (phi * SPH_S / 4.0)
        base = specs.get("air_0.00", {}).get("cps")
        if base:
            print("\nпроверка нормировки поля: цилиндр %.3f имп/с, сфера %.3f "
                  "имп/с, расхождение %.1f %%"
                  % (base, cps, 100 * (cps / base - 1)))

    empty = specs.get("air_0.00")
    if not empty:
        raise SystemExit("нет конфигурации пустого сосуда")

    # поканальный k и переподчёт в ROI
    print("\nПереподчёт фона при вычитании пустого сосуда")
    print("%-10s %-8s" % ("линия", "нуклид"), end="")
    others = [k for k in specs if k != "air_0.00"]
    for k in others:
        print(" %14s" % k, end="")
    print()

    eb = rcspec.fold(empty["spec_per_s"])
    folded = {k: rcspec.fold(specs[k]["spec_per_s"]) for k in others}
    for E, nuc in LINES:
        b0, lo, hi = rcspec.roi(eb, E)
        print("%-10.1f %-8s" % (E, nuc), end="")
        for k in others:
            b1 = folded[k][lo:hi].sum()
            kf = b1 / b0 if b0 > 0 else float("nan")
            print("  k=%.3f  %+.4f" % (kf, -(b0 - b1)), end="")
        print("   имп/с")
    print("\n(последний столбец каждой пары — систематический переподчёт,")
    print(" то есть сколько лишнего вычитается из площади ROI)")

    out = rcspec.rdir("background_k.csv")
    with open(out, "w", encoding="utf-8") as f:
        f.write("E_keV," + ",".join("k_" + k for k in others) + "\n")
        for i in range(20, 3000, 10):
            e0 = eb[i:i + 10].sum()
            if e0 <= 0:
                continue
            vals = [folded[k][i:i + 10].sum() / e0 for k in others]
            f.write("%.1f," % (i + 5) + ",".join("%.4f" % v for v in vals) + "\n")
    print("\nпоканальный k:", out)


if __name__ == "__main__":
    main()
