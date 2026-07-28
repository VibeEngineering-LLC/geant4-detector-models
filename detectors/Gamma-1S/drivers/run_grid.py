"""Сетка моноэнергий для Гамма-1С: Маринелли 1 л, ОИСН-16 (+ плотности для d_eff).

Энергии: 15 линий ЛСРМ (.efr Маринелли) + опорные точки краёв диапазона.
Каждая энергия — отдельный beamOn в одном процессе (геометрия одна).
"""
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

N = 400000

# Плотности для отдельной подгонки эффективной толщины d_eff: та же матрица,
# f(x) = (1-e^-x)/x, x = мю(E)*ро*d — как в radiacode-curves.
DENSITIES = [1.0, 1.6]


def macro(lines, n, tag):
    txt = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
           "/gps/particle gamma", "/gps/pos/type Volume",
           "/gps/pos/shape Cylinder", "/gps/pos/centre 0 0 16 mm",
           "/gps/pos/radius 73 mm", "/gps/pos/halfz 45 mm",
           "/gps/pos/confine Sample", "/gps/ang/type iso"]
    for e in lines:
        txt.append("/gps/energy %.3f keV" % e)
        txt.append("/g1s/outFile %s" % os.path.join(OUT, "%s_E%07.1f.csv" % (tag, e)))
        txt.append("/run/beamOn %d" % n)
    return "\n".join(txt) + "\n"


if __name__ == "__main__":
    # аргументы: [матрица] [плотности через запятую]
    #   run_grid.py                     -> OISN16, 1,00 и 1,60 (сверка с ЛСРМ)
    #   run_grid.py water 1.0           -> вода 1,0 (МИА по паспорту)
    matrix = sys.argv[1] if len(sys.argv) > 1 else "OISN16"
    rhos = ([float(x) for x in sys.argv[2].split(",")]
            if len(sys.argv) > 2 else DENSITIES)
    for rho in rhos:
        tag = "rho%.2f" % rho if matrix == "OISN16" else "%s%.2f" % (matrix, rho)
        mpath = os.path.join(BUILD, "grid_%s.mac" % tag)
        open(mpath, "w", encoding="utf-8").write(macro(LINES, N, tag))
        print("=== %s (%s) : %d энергий x %d событий ==="
              % (tag, matrix, len(LINES), N), flush=True)
        # ВАЖНО: text=True без encoding декодирует вывод в cp1251 (локаль
        # Windows) и падает на первом же байте вне таблицы — Geant4 печатает
        # и UTF-8, и служебные символы. Декодируем сами, с заменой.
        r = subprocess.run([os.path.join(BUILD, "g1s.exe"), mpath, "vessel",
                            str(rho), matrix],
                           cwd=BUILD, capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        for ln in (r.stdout or "").splitlines():
            if ln.startswith("RESULT") or "проба" in ln:
                print("  ", ln, flush=True)
        if r.returncode != 0:
            print("!! код возврата", r.returncode)
            print((r.stderr or "")[-2000:])
            sys.exit(1)
    print("готово")
