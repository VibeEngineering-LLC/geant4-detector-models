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
import gps_region  # noqa: E402
from grid_energies import LINES  # noqa: E402

N = 400000

# Плотности для отдельной подгонки эффективной толщины d_eff: та же матрица,
# f(x) = (1-e^-x)/x, x = мю(E)*ро*d — как в radiacode-curves.
DENSITIES = [1.0, 1.6]


def todo(tag, force=False):
    """Энергии, которых ещё нет: сетка расширяется по краям (grid_energies),
    и гонять заново все двадцать четыре точки — часы впустую. --force считает
    всё, это нужно после смены сборки модели."""
    if force:
        return list(LINES)
    left = [e for e in LINES
            if not os.path.exists(os.path.join(OUT, "%s_E%07.1f.csv"
                                               % (tag, e)))]
    if left and len(left) < len(LINES):
        print("   %s: уже посчитано %d из %d, считаю %d"
              % (tag, len(LINES) - len(left), len(LINES), len(left)),
              flush=True)
    return left


def macro(lines, n, tag, rho, matrix):
    # Тело розыгрыша — из выгрузки ПОСТРОЕННОЙ геометрии, а не константами:
    # прежние 73/45/16 писались под сосуд из таблицы ЛСРМ и после перехода на
    # чертёж изготовителя перестали покрывать пробу (R68, R75).
    txt = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
           "/gps/particle gamma"]
    txt += gps_region.gps_lines(BUILD, "vessel", [str(rho), matrix])
    txt += ["/gps/ang/type iso"]
    for e in lines:
        txt.append("/gps/energy %.3f keV" % e)
        txt.append("/g1s/outFile %s" % os.path.join(OUT, "%s_E%07.1f.csv" % (tag, e)))
        txt.append("/run/beamOn %d" % n)
    return "\n".join(txt) + "\n"


if __name__ == "__main__":
    # аргументы: [матрица] [плотности через запятую]
    #   run_grid.py                     -> OISN16, 1,00 и 1,60 (сверка с ЛСРМ)
    #   run_grid.py water 1.0           -> вода 1,0 (МИА по паспорту)
    args = [a for a in sys.argv[1:] if a != "--force"]
    force = "--force" in sys.argv
    matrix = args[0] if args else "OISN16"
    rhos = ([float(x) for x in args[1].split(",")]
            if len(args) > 1 else DENSITIES)
    for rho in rhos:
        tag = "rho%.2f" % rho if matrix == "OISN16" else "%s%.2f" % (matrix, rho)
        left = todo(tag, force)
        if not left:
            print("=== %s: всё посчитано ===" % tag, flush=True)
            continue
        mpath = os.path.join(BUILD, "grid_%s.mac" % tag)
        open(mpath, "w", encoding="utf-8").write(
            macro(left, N, tag, rho, matrix))
        print("=== %s (%s) : %d энергий x %d событий ==="
              % (tag, matrix, len(left), N), flush=True)
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
