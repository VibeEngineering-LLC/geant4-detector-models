# -*- coding: utf-8 -*-
"""Фон: отклик прибора в сосуде на однородное изотропное поле ЕРН помещения.

Считает то, ради чего всё затевалось: поканальный множитель k(канал) — какая
доля фона доходит до кристалла сквозь пробу по сравнению с пустым сосудом.
При анализе вычитается холостой набор на ПУСТОМ сосуде, а в самом измерении фон
уже искажён пробой, поэтому вычитание даёт систематическую поправку -B*(1-k).
Знак её не универсален: у плотных матриц проба ослабляет фон (k<1), у лёгких
рассеяние в пробе перевешивает ослабление и фон РАСТЁТ (k>1).

ПОЧЕМУ ЦИЛИНДР, А НЕ СФЕРА. Для любой выпуклой охватывающей поверхности с
равномерным розыгрышем точки и косинусным законом внутрь поле внутри однородно
и изотропно, а флюенс равен Ф = 4N/S (средняя хорда 4V/S). Сфера, в которую
влезает хвост прибора, даёт вчетверо большую площадь, то есть вчетверо больше
впустую выпущенных частиц. Проверка нормировки: прогон со сферой обязан дать
тот же отклик.

Цилиндр свой для каждого сосуда: он должен охватывать и сосуд с крышкой, и хвост
прибора, иначе поле внутри уже не однородно (та же ошибка, что была с областью
розыгрыша объёмного источника в run_grid.py).

Запуск:  python run_bg.py m500          — все конфигурации пятисотки
         python run_bg.py m200 --check  — плюс сверка «цилиндр против сферы»
"""
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

# Корни путей — из переменных окружения (common/py/paths.py).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = str(paths.build("RadiaCode-103"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/RadiaCode-103/drivers/run_grid.py\n"
        "Либо укажите G4MODELS_BUILD_RADIACODE_103 на готовый каталог."
        % BUILD)
EXE = os.path.join(BUILD, "rc_curves.exe")
RESULTS = rcspec.RESULTS
SPECTRUM = os.path.join(RESULTS, "field_spectrum.mac")

# Охватывающий цилиндр: R больше максимального радиуса сосуда, z от донышка
# (с крышкой) до конца хвоста прибора (нос в -12, корпус 123 => хвост +111).
#   m200: обвод-цилиндр R 36.05, юбка крышки 40.25, дно z -33.7
#   m500: обвод-полконус до R 55.6, крышка R 55.45 и дно с ней до z -49.7
CYL = {
    "m200": dict(r=45.0, z0=-45.0, z1=120.0),
    "m500": dict(r=58.0, z0=-55.0, z1=120.0),
}

SPH_R = 140.0
SPH_S = 4 * math.pi * (SPH_R / 10) ** 2

# Плотности: пустой сосуд (холостой набор) и сетка проб. organic 0.49 — та самая
# сушёная черника 246 г в пятисотке, на которой сверяется вся модель.
CONFIGS = [
    ("full", "air", 0.0012),    # пустой сосуд — холостой набор
    ("full", "organic", 0.49),
    ("full", "water", 1.00),
    ("full", "soil", 1.20),
    ("full", "soil", 1.60),
]
NEV = 30_000_000


def cyl_area(v):
    """Площадь боковой поверхности с торцами, см²."""
    c = CYL[v]
    r, hz = c["r"] / 10, 0.5 * (c["z1"] - c["z0"]) / 10
    return 2 * math.pi * r * (r + 2 * hz)


def macro(path, outfile, nev, shape, v):
    src = ["/run/verbose 0", "/event/verbose 0", "/tracking/verbose 0",
           "/run/printProgress 0", "",
           "/gps/pos/type Surface"]
    if shape == "cyl":
        c = CYL[v]
        src += ["/gps/pos/shape Cylinder",
                "/gps/pos/radius %.3f mm" % c["r"],
                "/gps/pos/halfz %.3f mm" % (0.5 * (c["z1"] - c["z0"])),
                "/gps/pos/centre 0 0 %.3f mm" % (0.5 * (c["z1"] + c["z0"]))]
    else:
        src += ["/gps/pos/shape Sphere",
                "/gps/pos/radius %.3f mm" % SPH_R,
                "/gps/pos/centre 0 0 0 mm"]
    src += ["/gps/ang/type cos",
            "/control/execute %s" % SPECTRUM.replace("\\", "/"),
            "/rc/outFile %s" % outfile.replace("\\", "/"),
            "/run/beamOn %d" % nev]
    open(path, "w", encoding="utf-8").write("\n".join(src) + "\n")


def run_one(cfg, v, shape="cyl", nev=NEV):
    mode, matrix, rho = cfg
    name = "%s_%.2f" % (matrix, rho)
    outdir = rcspec.rdir("background", v=v)
    os.makedirs(outdir, exist_ok=True)
    tagged = "bg_%s_%s.csv" % (shape, name)
    out = os.path.join(outdir, tagged)
    if os.path.exists(out):
        print("[--] %s/%s уже посчитано" % (v, tagged), flush=True)
        return 0
    mac = os.path.join(BUILD, "bg_%s_%s_%s.mac" % (v, shape, name))
    macro(mac, out, nev, shape, v)

    t0 = time.time()
    with open(os.path.join(outdir, "run_%s_%s.log" % (shape, name)), "w",
              encoding="utf-8") as lf:
        p = subprocess.run([EXE, mac, mode, matrix, "%.6f" % rho, v],
                           cwd=BUILD, stdout=lf, stderr=subprocess.STDOUT)
    print("[%s] %s/%s  %.1f мин" % ("ok" if p.returncode == 0 else "СБОЙ", v,
                                    tagged, (time.time() - t0) / 60), flush=True)
    return p.returncode


def main():
    if not os.path.exists(SPECTRUM):
        raise SystemExit("нет спектра поля: сначала wallfield.exe + analyze_field.py")
    v = rcspec.vessel()
    print("сосуд %s: цилиндр R %.1f мм, z %.1f..%.1f мм"
          % (v, CYL[v]["r"], CYL[v]["z0"], CYL[v]["z1"]))
    print("  S = %.1f см², Ф = 4N/S = N/%.2f" % (cyl_area(v), cyl_area(v) / 4))
    print("сфера S = %.1f см², Ф = 4N/S = N/%.2f" % (SPH_S, SPH_S / 4))

    jobs = [(c, v, "cyl", NEV) for c in CONFIGS]
    if "--check" in sys.argv:
        jobs.append((CONFIGS[0], v, "sph", NEV))
    with ThreadPoolExecutor(max_workers=len(jobs)) as ex:
        rc = list(ex.map(lambda j: run_one(*j), jobs))
    print("сбоев:", sum(1 for r in rc if r != 0))


if __name__ == "__main__":
    main()
