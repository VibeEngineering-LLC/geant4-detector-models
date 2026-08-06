# -*- coding: utf-8 -*-
"""Разброс результата по окну подгонки — он же оценка систематики.

Нижняя граница окна — главный произвол этого расчёта. Ниже примерно 150 кэВ
модель занижает измеренное почти втрое (K-серия дочерних, см.
`docs/wt20-remarks.md`), и этот участок даёт основную долю хи². Пока он не
закрыт, одно число публиковать нельзя: надо показать, насколько результат
зависит от того, где провести границу.

Скрипт прогоняет `wt20_unfold.py` с разными окнами и сводит в таблицу
активности, их отношение и хи². Сам код разложения при этом один и тот же —
окно передаётся переменными окружения, исходник не правится.

    python analysis/wt20_window_scan.py <спектр.xml> <каталог шаблонов> [выход]
"""
import io
import os
import re
import subprocess
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_DET = os.path.dirname(_HERE)

WINDOWS = [(60.0, 3300.0), (100.0, 3300.0), (150.0, 3300.0),
           (200.0, 3300.0), (300.0, 3300.0), (400.0, 3300.0)]


def run(src, tdir, lo, hi, workdir, calib_src):
    env = dict(os.environ)
    env["WT20_FIT_LO"] = "%.1f" % lo
    env["WT20_FIT_HI"] = "%.1f" % hi
    env["PYTHONIOENCODING"] = "utf-8"
    out = os.path.join(workdir, "w%d" % int(lo))
    os.makedirs(out, exist_ok=True)
    # поправка калибровки та же самая во всех прогонах — копируется, а не
    # считается заново: иначе к разбросу по окну примешался бы разброс
    # калибровки, и разделить их было бы уже нельзя
    if calib_src and os.path.exists(calib_src):
        io.open(os.path.join(out, "calibration_fitted.csv"), "w",
                encoding="utf-8").write(
            io.open(calib_src, encoding="utf-8").read())
    r = subprocess.run([sys.executable,
                        os.path.join(_HERE, "wt20_unfold.py"), src, tdir, out],
                       capture_output=True, text=True, encoding="utf-8",
                       env=env)
    txt = r.stdout or ""
    def grab(pat, cast=float):
        m = re.search(pat, txt)
        return cast(m.group(1)) if m else float("nan")
    return dict(
        lo=lo, hi=hi,
        chi2=grab(r"хи²/n = ([\d.]+)"),
        a1=grab(r"A1 \([^)]*\)\s+(\d+)"),
        a2=grab(r"A2 \([^)]*\)\s+(\d+)"),
        fwhm=grab(r"ПШПВ 662 кэВ = ([\d.]+)"),
        ok=(r.returncode == 0),
        tail=(r.stderr or "")[-200:],
    )


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, tdir = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(_DET, "results")
    os.makedirs(outdir, exist_ok=True)
    calib = os.path.join(outdir, "calibration_fitted.csv")
    work = tempfile.mkdtemp(prefix="wt20scan_")

    print("Разброс по окну подгонки. Поправка калибровки одна и та же во всех")
    print("прогонах, меняется только нижняя граница окна.")
    print()
    print("  окно, кэВ      хи²/n       A1, Бк    A2, Бк   A1/A2   ПШПВ 662")
    rows = []
    for lo, hi in WINDOWS:
        r = run(src, tdir, lo, hi, work, calib)
        ratio = r["a1"] / r["a2"] if r["a2"] else float("nan")
        rows.append((r, ratio))
        print("  %4.0f-%4.0f   %9.1f   %8.0f  %8.0f   %5.3f   %6.1f"
              % (lo, hi, r["chi2"], r["a1"], r["a2"], ratio, r["fwhm"]))

    a1 = [r["a1"] for r, _ in rows]
    a2 = [r["a2"] for r, _ in rows]
    rr = [x for _, x in rows]
    print()
    print("  A1: %.0f…%.0f Бк, размах %.2f раза" % (min(a1), max(a1), max(a1) / min(a1)))
    print("  A2: %.0f…%.0f Бк, размах %.2f раза" % (min(a2), max(a2), max(a2) / min(a2)))
    print("  A1/A2: %.3f…%.3f" % (min(rr), max(rr)))
    print()
    print("  Разброс по окну — это систематика метода, а не статистика.")
    print("  Пока участок ниже 150 кэВ не описан моделью, публиковать надо")
    print("  интервал, а не число.")

    p = os.path.join(outdir, "wt20_window_scan.csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        f.write("# разброс результата по нижней границе окна подгонки\n")
        f.write("окно_низ_кэВ;окно_верх_кэВ;хи2_на_канал;A1_Бк;A2_Бк;A1/A2;ПШПВ_662_кэВ\n")
        for r, ratio in rows:
            f.write("%.0f;%.0f;%.4g;%.0f;%.0f;%.4f;%.3g\n"
                    % (r["lo"], r["hi"], r["chi2"], r["a1"], r["a2"],
                       ratio, r["fwhm"]))
    print()
    print("записано:", p)


if __name__ == "__main__":
    main()
