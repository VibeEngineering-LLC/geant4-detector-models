"""Сетки эффективности для ВСЕХ геометрий комплекта поверки.

Объёмные геометрии (Маринелли, Дента, Петри) считаются при двух плотностях:
этого достаточно, чтобы подгонкой f(mu*ro*d) получить эффективную толщину и
дальше предсказывать любую плотность. У источников комплекта плотности разные
(масса при одном номинальном объёме кюветы), поэтому без этого не обойтись.

Точечные геометрии считаются с ОГРАНИЧЕНИЕМ ТЕЛЕСНОГО УГЛА: изотропный
источник на 25 см тратит 99,4 % событий впустую. Кванты разыгрываются в конус
вокруг направления на детектор, а число «эквивалентных полных» получается
делением на долю телесного угла (1 − cos θmax)/2. В GPS при `iso` угол θ
отсчитывается от −Z, то есть θ = 0 уже смотрит вниз, на детектор.

ОГОВОРКА к конусу: он корректен для ППП, куда попадают практически только
прямые кванты. Полная эффективность при таком розыгрыше занижена — рассеянные
на стенках защиты кванты не разыгрываются.
"""
import math
import os
import subprocess
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402


BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)

EXE = "g1s.exe"

if not os.path.exists(os.path.join(BUILD, EXE)):
    raise SystemExit(
        "Не найдена собранная модель %s.\n"
        "Соберите её (см. common/cmake и README детектора) или укажите\n"
        "G4MODELS_BUILD_GAMMA_1S на каталог, где она уже лежит."
        % os.path.join(BUILD, EXE))
OUT = os.path.join(BUILD, "grid")
os.makedirs(OUT, exist_ok=True)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grid_energies import LINES  # noqa: E402

# сосуд -> (радиус розыгрыша мм, полувысота мм, центр по z мм)
VOL_SRC = {
    "marinelli": (73.0, 45.0, 16.0),
    "denta": (36.0, 18.0, 61.0),    # кювета стоит на торце z=43, высота 35
    "petri": (42.5, 7.0, 50.0),     # z 43..57
}

JOBS = [
    # (метка, сосуд, матрица, плотность, объём мл, режим, N)
    ("water1.00", "marinelli", "water", 1.00, 1000.0, "vessel:marinelli", 400000),
    ("denta0.60", "denta", "OISN16", 0.60, 120.0, "vessel:denta", 400000),
    ("denta1.60", "denta", "OISN16", 1.60, 120.0, "vessel:denta", 400000),
    ("petri0.60", "petri", "OISN16", 0.60, 60.0, "vessel:petri", 400000),
    ("petri1.60", "petri", "OISN16", 1.60, 60.0, "vessel:petri", 400000),
]

# точечные: (метка, расстояние от торца мм, режим экрана, θmax град, N)
POINTS = [("p5cm", 50.0, "shield", 60.0, 400000),
          ("p25cm", 250.0, "open", 30.0, 400000)]

ZFACE = 43.0    # наружная плоскость крышки детектора


# Досчёт вместо пересчёта: энергия, для которой файл уже есть, пропускается.
# Нужно потому, что сетка расширяется по краям (см. grid_energies.py), а гонять
# заново все двадцать четыре точки в девяти геометриях — часы впустую.
# --force считает всё; так надо, если сменилась сборка модели.
FORCE = False


def todo(tag):
    if FORCE:
        return list(LINES)
    left = [e for e in LINES
            if not os.path.exists(os.path.join(OUT, "%s_E%07.1f.csv" % (tag, e)))]
    if len(left) < len(LINES):
        print("   %s: уже посчитано %d из %d, считаю %d"
              % (tag, len(LINES) - len(left), len(LINES), len(left)), flush=True)
    return left


def macro_volume(tag, vessel, n, lines):
    r, hz, zc = VOL_SRC[vessel]
    t = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
         "/gps/particle gamma", "/gps/pos/type Volume",
         "/gps/pos/shape Cylinder", "/gps/pos/centre 0 0 %.1f mm" % zc,
         "/gps/pos/radius %.1f mm" % r, "/gps/pos/halfz %.1f mm" % hz,
         "/gps/pos/confine Sample", "/gps/ang/type iso"]
    for e in lines:
        t.append("/gps/energy %.3f keV" % e)
        t.append("/g1s/outFile %s" % os.path.join(OUT, "%s_E%07.1f.csv" % (tag, e)))
        t.append("/run/beamOn %d" % n)
    return "\n".join(t) + "\n"


def macro_point(tag, dist, thmax, n):
    t = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
         "/gps/particle gamma", "/gps/pos/type Point",
         "/gps/pos/centre 0 0 %.1f mm" % (ZFACE + dist),
         "/gps/ang/type iso", "/gps/ang/maxtheta %.1f deg" % thmax]
    for e in LINES:
        t.append("/gps/energy %.3f keV" % e)
        t.append("/g1s/outFile %s" % os.path.join(OUT, "%s_E%07.1f.csv" % (tag, e)))
        t.append("/run/beamOn %d" % n)
    return "\n".join(t) + "\n"


def run(mpath, args, label):
    print("=== %s ===" % label, flush=True)
    r = subprocess.run([os.path.join(BUILD, "g1s.exe"), mpath] + args,
                       cwd=BUILD, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    ok = 0
    for ln in (r.stdout or "").splitlines():
        if ln.startswith("RESULT"):
            ok += 1
        elif "проба" in ln or "COMMAND NOT" in ln:
            print("   ", ln.strip(), flush=True)
    print("    точек посчитано: %d" % ok, flush=True)
    if r.returncode != 0:
        print("!! код возврата", r.returncode)
        print((r.stderr or "")[-1500:])
        sys.exit(1)


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for tag, ves, mat, rho, vol, mode, n in JOBS:
        if only and only != tag:
            continue
        p = os.path.join(BUILD, "grid_%s.mac" % tag)
        open(p, "w", encoding="utf-8").write(macro_volume(tag, ves, n))
        run(p, [mode, str(rho), mat, str(vol)],
            "%s: %s, %s, ро=%.2f, %.0f мл" % (tag, ves, mat, rho, vol))

    for tag, dist, mode, th, n in POINTS:
        if only and only != tag:
            continue
        p = os.path.join(BUILD, "grid_%s.mac" % tag)
        open(p, "w", encoding="utf-8").write(macro_point(tag, dist, th, n))
        frac = (1 - math.cos(math.radians(th))) / 2
        run(p, [mode], "%s: точечный на %.0f мм, конус %.0f град, доля угла %.5f"
                       % (tag, dist, th, frac))
        # доля телесного угла нужна при разборе — кладём рядом с данными
        open(os.path.join(OUT, "%s_solidangle.txt" % tag), "w").write(
            "%.8f\n" % frac)
    print("готово")
