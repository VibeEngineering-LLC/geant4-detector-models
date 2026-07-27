"""Прогоны распада нуклидов смеси №SRC-04 в трёх геометриях комплекта.

Для восстановления активностей нужна эффективность НА РАСПАД по каждой
аналитической линии: она включает и выход линии, и каскадное суммирование
(у Eu-152 и Sc-44 оно велико), и ничего не берёт из справочника.

Границы цепочек:
  Am-241: A=241 — дочерний Np-237 (2,1e6 лет) при поднятом пороге времени
          РАСПАЛСЯ БЫ, чего в реальном источнике нет; окно его отсекает.
  Ti-44:  44 44 20 22 — Ti-44 -> Sc-44 -> Ca-44, обе ступени внутри окна,
          позитрон и 1157 кэВ в истинном совпадении.
  Eu-152: A=152 — распадается в стабильные Sm-152/Gd-152.
  Cs-137: A=137.

Матрица смеси: РИСН-379 — состав по массе из ОРИГИНАЛЬНОГО .spe («Поверка
2016», поле MATERIAL): H 4,3 C 33,0 N 1,2 O 34,8 Na 4,1 Mg 2,2 Ca 20,3 %,
ро = 1,0. НЕ вода: кальций на 59,5 кэВ (Am-241) поднимает поглощение.
Первый проход считался на воде — расхождение вода/РИСН-379 по линиям
>120 кэВ в пределах процента, по 59,5 кэВ — проверить сравнением прогонов.
Массы и объёмы — из RAWMASS оригиналов: Маринелли 1000 г / 1000 мл,
Дента 100 г / 100 мл (геометрия «Дента-100»), Петри 60 г / 60 мл.

Запуск: python run_mix_decay.py [геометрия] [матрица]
        матрица по умолчанию risn379; water — для сравнения.
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

NUCS = [
    # (метка, Z, A, окно Amin Amax Zmin Zmax, распадов)
    ("Am241", 95, 241, "241 241 93 95", 300000),
    ("Ti44", 22, 44, "44 44 20 22", 300000),
    ("Eu152", 63, 152, "152 152 62 64", 300000),
    ("Cs137", 55, 137, "137 137 55 56", 300000),
]

# (геометрия, режим, плотность, объём мл, радиус розыгрыша, полу-Z, центр Z)
GEOMS = [
    ("marinelli", "vessel:marinelli", 1.00, 1000.0, 73.0, 45.0, 16.0),
    ("denta", "vessel:denta", 1.00, 100.0, 36.0, 18.0, 61.0),
    ("petri", "vessel:petri", 1.00, 60.0, 42.5, 7.0, 50.0),
]


def macro(geom, r, hz, zc):
    t = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
         "/process/had/rdm/verbose 0",
         "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns",
         "/gps/particle ion", "/gps/energy 0 keV",
         "/gps/pos/type Volume", "/gps/pos/shape Cylinder",
         "/gps/pos/centre 0 0 %.1f mm" % zc,
         "/gps/pos/radius %.1f mm" % r, "/gps/pos/halfz %.1f mm" % hz,
         "/gps/pos/confine Sample", "/gps/ang/type iso"]
    for name, z, a, lim, n in NUCS:
        t += ["/process/had/rdm/nucleusLimits " + lim,
              "/gps/ion %d %d 0 0" % (z, a),
              "/g1s/outFile %s" % os.path.join(BUILD, "mix_%s_%s.csv" % (geom, name)),
              "/run/beamOn %d" % n]
    return "\n".join(t) + "\n"


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    # По умолчанию вода: текущий g1s.exe собран без risn379, а очередь
    # прогонов уже идёт на нём. После пересборки запускать явно:
    #   python run_mix_decay.py "" risn379
    matrix = sys.argv[2] if len(sys.argv) > 2 else "water"
    for geom, mode, rho, vol, r, hz, zc in GEOMS:
        if only and only != geom:
            continue
        mp = os.path.join(BUILD, "mix_%s.mac" % geom)
        open(mp, "w", encoding="utf-8").write(macro(geom, r, hz, zc))
        print("=== %s: 4 нуклида x распады, %s %s ро=%.2f %.0f мл ==="
              % (geom, mode, matrix, rho, vol), flush=True)
        res = subprocess.run([os.path.join(BUILD, "g1s.exe"), mp, mode,
                              str(rho), matrix, str(vol)],
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
