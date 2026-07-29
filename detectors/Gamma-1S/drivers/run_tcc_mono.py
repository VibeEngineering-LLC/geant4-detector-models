"""Моноэнергетические точки под собственную проверку суммирования (TCC).

ЗАЧЕМ ОТДЕЛЬНЫЙ ДРАЙВЕР, А НЕ РАСШИРЕНИЕ ОБЩЕЙ СЕТКИ. Проверке нужны линии
Co-60 (1173,23 и 1332,49 кэВ), которых в `grid_energies.LINES` нет. Дописать их
туда нельзя: файлы `p5cm_E*.csv` разбирают глобом восемь скриптов
(compare_point, export_curves, mda, summing, point_recalc и другие), и новые
узлы молча вошли бы в опубликованные кривые эффективности как полноправные
точки калибровки. Здесь у прогонов свой тег `tcc5cm` — общая сетка не
затрагивается, а `eps_mono_point` параметризован тегом и читает их без правок.

ГЕОМЕТРИЯ ПОВТОРЯЕТ p5cm ТОЧНО: источник на ZFACE+50 мм, конус θmax = 60°,
режим `shield`. Иначе сравнивать декай-прогон с моно-прогоном нельзя — разница
геометрии подменит собой измеряемый эффект.

661,657 кэВ считается вместе с линиями Co-60 намеренно, хотя эта энергия есть и
в общей сетке: совпадение `tcc5cm` с `p5cm` на ней — контроль того, что
отдельный прогон воспроизводит общий, а не живёт своей жизнью.
"""
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
OUT = os.path.join(BUILD, "grid")

if not os.path.exists(os.path.join(BUILD, "g1s.exe")):
    raise SystemExit(
        "Не найдена собранная модель %s.\n"
        "Соберите её или укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % os.path.join(BUILD, "g1s.exe"))

TAG = "tcc5cm"
ZFACE = 41.0        # наружная плоскость = крышка Al 2 мм
DIST = 50.0         # мм от торца — та же точка, что у p5cm
THMAX = 60.0        # градусов: тот же конус, что у p5cm
MODE = "shield"
N = 400000

# Cs-137 — контроль без каскада; обе линии Co-60 — рабочие точки.
LINES = [661.657, 1173.23, 1332.49]


def macro(lines):
    t = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
         "/gps/particle gamma", "/gps/pos/type Point",
         "/gps/pos/centre 0 0 %.1f mm" % (ZFACE + DIST),
         "/gps/ang/type iso", "/gps/ang/maxtheta %.1f deg" % THMAX]
    for e in lines:
        t.append("/gps/energy %.3f keV" % e)
        t.append("/g1s/outFile %s"
                 % os.path.join(OUT, "%s_E%07.1f.csv" % (TAG, e)))
        t.append("/run/beamOn %d" % N)
    return "\n".join(t) + "\n"


if __name__ == "__main__":
    force = "--force" in sys.argv
    os.makedirs(OUT, exist_ok=True)

    frac = (1 - math.cos(math.radians(THMAX))) / 2
    open(os.path.join(OUT, "%s_solidangle.txt" % TAG), "w").write(
        "%.8f\n" % frac)

    left = LINES if force else [
        e for e in LINES
        if not os.path.exists(os.path.join(OUT, "%s_E%07.1f.csv" % (TAG, e)))]
    if not left:
        print("всё посчитано, считать нечего (--force пересчитает)")
        raise SystemExit(0)

    mp = os.path.join(BUILD, "tcc_mono.mac")
    open(mp, "w", encoding="utf-8").write(macro(left))
    print("=== моно под TCC: %d энергий, конус %.0f град, доля угла %.5f ==="
          % (len(left), THMAX, frac), flush=True)

    r = subprocess.run([os.path.join(BUILD, "g1s.exe"), mp, MODE],
                       cwd=BUILD, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    # Наружу пропускаются не только RESULT: фильтр по allow-list уже однажды
    # съел предупреждения о перекрытиях геометрии, которые печатались годами.
    for ln in (r.stdout or "").splitlines():
        s = ln.strip()
        if s.startswith("RESULT") or "WARNING" in s or "ERROR" in s \
                or "G4Exception" in s or "Overlap" in s:
            print("  ", s, flush=True)
    if r.returncode != 0:
        print("!! код возврата", r.returncode)
        print((r.stderr or "")[-1500:])
        sys.exit(1)
    print("готово")
