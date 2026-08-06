# -*- coding: utf-8 -*-
"""Сверка расчётной кривой с библиотечной кривой BecqMoni (слепой тест).

Обе кривые — АБСОЛЮТНЫЕ эффективности по ППП для одной постановки: точечный
источник на оси, 100 мм от наружной плоскости, торец в пучке. Никакого
подгоночного множителя не вводится и вводить нельзя: смысл теста в том, что
уровень сошёлся сам или не сошёлся сам.

Сетки разные (наша 29 узлов, библиотечная 34). Узел считается ОБЩИМ, если
энергии совпадают в пределах 1 кэВ; иначе библиотечная берётся лог-лог
интерполяцией её же кривой, и такая точка помечается — интерполяция вносит
свою погрешность, а объявить её нечем.

Погрешность расхождения комбинируется квадратично: наша d_eps_peak (случайная,
БЕЗ систематики конвенции съёма) и объявленная библиотекой d_eff_pct. Отсюда
следует читать «в пределах k сигма», а не «отличается на столько-то процентов»:
систематика конвенции в сигму не входит и по величине с ней сопоставима.

    python analysis/compare_library.py [<наша кривая.csv>]
"""
import io
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.normpath(os.path.join(_HERE, "..", "reference",
                                    "becqmoni-library-curve.csv"))
OURS = os.path.normpath(os.path.join(_HERE, "..", "results",
                                     "eff_point_end10cm.csv"))
E_TOL_KEV = 1.0


def read_csv(path, cols):
    """-> список dict по именованным столбцам; строки с '#' пропускаются."""
    rows, head = [], None
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split(",")
        if head is None:
            head = parts
            continue
        rec = dict(zip(head, parts))
        rows.append({k: float(rec[k]) for k in cols})
    if not rows:
        raise SystemExit("в %s нет данных" % path)
    return rows


def loglog_interp(xs, ys, x):
    """Интерполяция по логарифмам обеих осей. Вне диапазона -> None."""
    if x < xs[0] or x > xs[-1]:
        return None
    for i in range(1, len(xs)):
        if x <= xs[i]:
            t = ((math.log(x) - math.log(xs[i - 1]))
                 / (math.log(xs[i]) - math.log(xs[i - 1])))
            return math.exp(math.log(ys[i - 1])
                            + t * (math.log(ys[i]) - math.log(ys[i - 1])))
    return ys[-1]


def main():
    ours_path = sys.argv[1] if len(sys.argv) > 1 else OURS
    lib = read_csv(LIB, ("E_keV", "eff_peak", "d_eff_pct"))
    ours = read_csv(ours_path, ("E_keV", "eps_peak", "d_eps_peak"))
    xs = [r["E_keV"] for r in lib]
    ys = [r["eff_peak"] for r in lib]
    ds = [r["d_eff_pct"] for r in lib]

    print("%9s %13s %13s %9s %8s %7s %6s"
          % ("E, кэВ", "наша", "библиотечная", "разн.,%", "сигма,%", "k", "узел"))
    out = []
    for r in ours:
        e = r["E_keV"]
        exact = None
        for i, x in enumerate(xs):
            if abs(x - e) <= E_TOL_KEV:
                exact = i
                break
        if exact is not None:
            v, dpct, kind = ys[exact], ds[exact], "общий"
        else:
            v = loglog_interp(xs, ys, e)
            if v is None:
                print("%9.1f %13.4e %13s   вне диапазона библиотечной сетки"
                      % (e, r["eps_peak"], "—"))
                continue
            # погрешность берётся ближайшего библиотечного узла: своей у
            # интерполированной точки нет, и придумывать её нельзя
            j = min(range(len(xs)), key=lambda i: abs(xs[i] - e))
            dpct, kind = ds[j], "интерп."
        diff = 100.0 * (r["eps_peak"] / v - 1.0)
        s_ours = 100.0 * r["d_eps_peak"] / r["eps_peak"]
        sig = math.hypot(s_ours, dpct)
        print("%9.1f %13.4e %13.4e %+9.2f %8.2f %6.1f  %s"
              % (e, r["eps_peak"], v, diff, sig, abs(diff) / sig, kind))
        out.append((e, diff, sig, kind))

    if out:
        soft = [d for e, d, s, k in out if e <= 150]
        hard = [d for e, d, s, k in out if e >= 1173]
        mid = [d for e, d, s, k in out if 300 <= e <= 1000]
        def rng(v):
            return "нет узлов" if not v else "%+.2f … %+.2f %%" % (min(v), max(v))
        print("\nмягкий край (<=150 кэВ):  %s" % rng(soft))
        print("середина (300-1000):      %s" % rng(mid))
        print("жёсткий край (>=1173):    %s" % rng(hard))
    return 0


if __name__ == "__main__":
    sys.exit(main())
