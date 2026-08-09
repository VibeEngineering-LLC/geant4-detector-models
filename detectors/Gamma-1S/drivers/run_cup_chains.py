"""Цепочки Ra-226 и Th-232 в геометриях «Дента» и Петри + контроль Cs-137.

Нужны для пересчёта источников комплекта в этих кюветах: поправка на
каскадное суммирование зависит от геометрии (телесного угла), и брать её из
маринелльных прогонов нельзя. Выходы линий на распад родителя геометрия не
меняет — они сверяются с маринелльными как контроль целостности.

Плотности — из описи комплекта (масса / номинальный объём):
  Дента: источники 68–192 г на 120 мл -> 0,57–1,60;
  Петри: 34–96 г на 60 мл -> 0,57–1,60.
Прогоны идут при ро = 1,0 (середина): C слабо зависит от плотности, потому
что это отношение эффективностей в ОДНОЙ геометрии; проверено на маринелли.
"""
import os
import subprocess
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gps_region  # noqa: E402


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

NUCS = [
    ("Ra226chain", 88, 226, "214 226 82 88", 150000),
    ("Th232chain", 90, 232, "208 232 81 90", 150000),
    ("Cs137", 55, 137, "137 137 55 56", 300000),   # контроль C = 1
]

GEOMS = [
    ("denta", "vessel:denta", 1.00, 120.0),
    ("petri", "vessel:petri", 1.00, 60.0),
]


def macro(geom, args):
    # Тело розыгрыша — из выгрузки построенной геометрии (gps_region), не
    # константами: константы расходятся с моделью беззвучно (R75).
    t = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
         "/process/had/rdm/verbose 0",
         "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns",
         "/gps/particle ion", "/gps/energy 0 keV"]
    t += gps_region.gps_lines(BUILD, args[0], args[1:])
    t += ["/gps/ang/type iso"]
    for name, z, a, lim, n in NUCS:
        t += ["/process/had/rdm/nucleusLimits " + lim,
              "/gps/ion %d %d 0 0" % (z, a),
              "/g1s/outFile %s" % os.path.join(BUILD, "cup_%s_%s.csv" % (geom, name)),
              "/run/beamOn %d" % n]
    return "\n".join(t) + "\n"


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for geom, mode, rho, vol in GEOMS:
        if only and only != geom:
            continue
        gargs = [mode, str(rho), "OISN16", str(vol)]
        mp = os.path.join(BUILD, "cupch_%s.mac" % geom)
        open(mp, "w", encoding="utf-8").write(macro(geom, gargs))
        print("=== %s: цепочки Ra/Th + Cs, %s ===" % (geom, mode), flush=True)
        res = subprocess.run([os.path.join(BUILD, "g1s.exe"), mp] + gargs,
                             cwd=BUILD, capture_output=True, text=True,
                             encoding="utf-8", errors="replace")
        for ln in (res.stdout or "").splitlines():
            if ln.startswith(("RESULT", "EMIT")) or "проба" in ln:
                print("  ", ln.strip(), flush=True)
        if res.returncode != 0:
            print("!! код возврата", res.returncode)
            print((res.stderr or "")[-1500:])
            sys.exit(1)
    print("готово")
