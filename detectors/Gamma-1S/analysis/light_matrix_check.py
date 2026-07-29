"""Цена f-пересчёта между матрицами: прямой прогон ОИСН-06 против формулы.

Задача #107 (проверка метода, не расследование: разность партий в маринелли
1,32 сигма и сама по себе незначима). Активность 420-17031 (ро=0,64, ОИСН-06)
считалась через eps штатного прогона ро=1,6/ОИСН-16 с пересчётом
f(мю*ро*d_eff), лёгкую среду представляла вода. Прогон run_light_matrix.py
даёт ту же величину прямым транспортом. Здесь — отношение двух способов.

Обе стороны снимаются одним окном (уширение до приборного разрешения, как в
kit_recalc), поэтому сравнивается именно перенос между матрицами, а не
конвенция площади. Геометрия прогонов совпадает с точностью до колодца
40,15/40,00 (снятый конфликт #103, влияние ~1e-4).
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kit_recalc as kr  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
LINES = [583.187, 911.204, 2614.511]
RHO = 0.64

if __name__ == "__main__":
    p = os.path.join(BUILD, "chain_Th232_oisn06.csv")
    if not os.path.exists(p):
        raise SystemExit("нет chain_Th232_oisn06.csv — сперва "
                         "drivers/run_light_matrix.py")
    hist, N = kr.load_hist(p)
    print("Прямой прогон ОИСН-06 ро=0,64 против f-пересчёта с ро=1,6/ОИСН-16\n"
          "(маринелли, обе стороны — уширенное окно ±1 ПШПВ):\n")
    print("%9s %12s %12s %14s" % ("E, кэВ", "eps прямая", "eps формула",
                                  "формула/прямая"))
    rows = []
    for E in LINES:
        fw = kr.FWHM662 * math.sqrt(E / 661.657)
        a = kr.area_sim(hist, E, fwhm=fw, key="oisn06")
        direct = a / N
        # статистика площади: грубо sqrt(a)/a (полки малы против пика)
        da = 1.0 / math.sqrt(max(a, 1.0))
        key = min(kr.MU_W, key=lambda k: abs(k - E))
        pred = kr.eps_per_decay("Marinelli_1L", "Th232chain", E, fw,
                                RHO, kr.MU_W[key])
        r = pred / direct
        dr = r * math.hypot(da, 0.005)
        print("%9.1f %12.4e %12.4e %10.4f±%.4f" % (E, direct, pred, r, dr))
        rows.append((E, direct, pred, r, dr))

    out = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results",
        "light_matrix_check.csv"))
    csvio.write(
        out,
        ["E_keV", "eps_direct", "eps_formula", "formula_over_direct",
         "d_ratio"],
        [("%.3f" % E, "%.6e" % d, "%.6e" % p, "%.4f" % r, "%.4f" % dr)
         for E, d, p, r, dr in rows],
        comments=[
            "Цена f(mu*ро*d_eff)-пересчёта между матрицами: маринелли,"
            " прямой прогон ОИСН-06 ро=0,64 (600k распадов)",
            "против пересчёта со штатного прогона ро=1,6/ОИСН-16 (вода как"
            " лёгкая среда), окна одинаковые.",
            "formula_over_direct < 1 — формула ЗАНИЖАЕТ eps (и завышает"
            " активность) на этой линии.",
        ])
    print("\nтаблица: %s" % out)
