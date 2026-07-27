# -*- coding: utf-8 -*-
"""Драйвер сетки прогонов: генерирует макросы и запускает rc_curves по
конфигурациям (матрица пробы, плотность) параллельно, по процессу на ядро.

Одна конфигурация = один процесс: плотность и состав пробы задаются до
построения геометрии, поэтому в пределах процесса их менять нельзя.

Запуск:
    python run_grid.py gamma          — сетка фотонов (кривые эффективности)
    python run_grid.py beta           — сетка электронов (проникающая способность)
    python run_grid.py gamma --quick  — вчетверо меньше событий, для проверки
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

# Корни путей — из переменных окружения (common/py/paths.py).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

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
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))

# Сетка энергий фотонов, кэВ: логарифмический костяк плюс реальные линии ЕРН
# (Pb-210 46.5, Th-234 92.6, Ra-226 186, Pb-212 238.6, Pb-214 295/352,
#  Tl-208 583/2615, Bi-214 609/1120/1764/2204, Cs-137 662, Ac-228 911/969,
#  Co-60 1173/1333, K-40 1461).
E_GAMMA = [30, 46.5, 59.5, 80, 92.6, 122, 186, 238.6, 295.2, 351.9, 477,
           583.2, 609.3, 661.7, 795, 911.2, 968.9, 1120.3, 1173.2, 1332.5,
           1460.8, 1764.5, 2204, 2614.5, 3000]

# Сетка энергий электронов, кэВ: вокруг порога проникновения ~1.15 МэВ
E_BETA = [200, 400, 600, 800, 1000, 1100, 1200, 1311, 1400, 1600, 1800,
          2000, 2274, 2500, 3000, 3270]

# (режим, матрица, плотность). air — предельный случай «без самопоглощения»,
# organic 0.50 служит проверкой того, что самопоглощение определяется только
# произведением rho*(mu/rho), а не составом матрицы: по этой проверке остальные
# матрицы (зола, песок, продукты) интерполируются без отдельных прогонов.
# Конфигураций ровно столько, сколько рабочих процессов — иначе последняя
# считается в одиночку и удваивает время сетки.
CONFIGS = [
    ("full", "air", 0.0012),
    ("full", "water", 1.00),
    ("full", "soil", 0.80),
    ("full", "soil", 1.20),
    ("full", "soil", 1.60),
    ("full", "organic", 0.50),
]

N_GAMMA = 2_000_000
N_BETA = 400_000


def tag(mode, matrix, rho):
    return "%s_%s_%.2f" % (mode, matrix, rho)


# Область розыгрыша объёмного источника: цилиндр, ОБЯЗАТЕЛЬНО охватывающий всё
# тело пробы, из него отбираются точки внутри объёма «sample». Если взять
# меньше, источник окажется только во внутренней части пробы и эффективность
# выйдет завышенной — на этом я один раз уже обжёгся, сравнивая сосуды.
VESSELS = {
    "m200": dict(radius=33.24, halfz=33.25, centre=-0.56),
    "m500": dict(radius=43.34, halfz=47.10, centre=2.06),
}


def write_macro(path, particle, energies, nev, outdir, vessel):
    """Пропускаем энергии, для которых результат уже посчитан, — прогон сетки
    можно прервать и продолжить."""
    todo = [e for e in energies
            if not os.path.exists(os.path.join(outdir, "E%07.1f.csv" % e))]
    v = VESSELS[vessel]
    src = ["# сгенерировано run_grid.py",
           "/run/verbose 0", "/event/verbose 0", "/tracking/verbose 0",
           "/run/printProgress 0", "",
           "/gps/particle %s" % particle,
           "/gps/ene/type Mono",
           "# объёмный источник, равномерно распределённый по телу пробы",
           "/gps/pos/type Volume",
           "/gps/pos/shape Cylinder",
           "/gps/pos/radius %.2f mm" % v["radius"],
           "/gps/pos/halfz %.2f mm" % v["halfz"],
           "/gps/pos/centre 0 0 %.2f mm" % v["centre"],
           "/gps/pos/confine sample",
           "/gps/ang/type iso", ""]
    for e in todo:
        src.append("/gps/ene/mono %.4f keV" % e)
        src.append("/rc/outFile %s" % os.path.join(outdir, "E%07.1f.csv" % e).replace("\\", "/"))
        src.append("/run/beamOn %d" % nev)
    open(path, "w", encoding="utf-8").write("\n".join(src) + "\n")
    return len(todo)


def run_one(cfg, gps, pname, energies, nev, vessel):
    mode, matrix, rho = cfg
    name = tag(mode, matrix, rho)
    outdir = os.path.join(RESULTS, vessel, pname, name)
    os.makedirs(outdir, exist_ok=True)
    mac = os.path.join(BUILD, "grid_%s_%s_%s.mac" % (vessel, pname, name))
    todo = write_macro(mac, gps, energies, nev, outdir, vessel)
    if todo == 0:
        print("[--] %s  %s  уже посчитано" % (pname, name), flush=True)
        return 0

    t0 = time.time()
    log = os.path.join(outdir, "run.log")
    with open(log, "w", encoding="utf-8") as lf:
        p = subprocess.run([EXE, mac, mode, matrix, "%.6f" % rho, vessel],
                           cwd=BUILD, stdout=lf, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    print("[%s] %s  %s  %.1f мин" % ("ok" if p.returncode == 0 else "СБОЙ",
                                     pname, name, dt / 60), flush=True)
    return p.returncode


def main():
    pname = sys.argv[1] if len(sys.argv) > 1 else "gamma"
    quick = "--quick" in sys.argv
    vessel = next((a for a in sys.argv[1:] if a in VESSELS), "m200")
    configs = CONFIGS
    if pname == "gamma":
        gps, energies, nev = "gamma", E_GAMMA, N_GAMMA
    elif pname == "beta":
        gps, energies, nev = "e-", E_BETA, N_BETA
        # Для беты конфигурация «воздух» бесполезна и неподъёмна: электроны в
        # разреженной среде летят через весь сосуд, давая огромное число шагов
        # многократного рассеяния (первая же точка не считается за 10 минут,
        # тогда как плотные матрицы проходят по две). Предельный случай
        # «без самопоглощения» для беты и не нужен — реальная проба не воздух.
        configs = [c for c in CONFIGS if c[1] != "air"]
    else:
        raise SystemExit("частица: gamma | beta")
    if quick:
        nev //= 4

    total = len(configs) * len(energies) * nev
    print("сосуд %s, конфигураций %d, энергий %d, событий на точку %d => %.1f млн"
          % (vessel, len(configs), len(energies), nev, total / 1e6), flush=True)

    nproc = len(configs)   # все конфигурации одновременно: см. комментарий к CONFIGS
    with ThreadPoolExecutor(max_workers=nproc) as ex:
        rc = list(ex.map(lambda c: run_one(c, gps, pname, energies, nev, vessel),
                         configs))
    print("сбоев:", sum(1 for r in rc if r != 0))


if __name__ == "__main__":
    main()
