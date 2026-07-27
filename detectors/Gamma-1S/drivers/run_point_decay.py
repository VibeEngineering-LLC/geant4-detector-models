"""Прогоны распада точечных источников ОСГИ на 5 см от крышки.

Зачем: у Eu-152, Ba-133, Co-60, Y-88, Na-22, Bi-207 и цепочки Th-228 на 5 см
каскадное суммирование составляет единицы—десятки процентов, и пересчёт их
записей без него бессмыслен. Эффективность на распад из этих прогонов включает
выход линии, суммирование и аннигиляцию (Na-22, Y-88, Bi-207) сразу.

Конус здесь ЗАПРЕЩЁН: кванты каскада коррелированы через совпадения в
кристалле, и обрезка направлений исказила бы вероятность второму кванту
попасть в детектор. Только полный изотроп, плата — статистика.

На 25 см отдельные прогоны не гоняются: C-1 пропорциональна полной
эффективности партнёра по каскаду, то есть телесному углу; отношение углов
25 см/5 см ~ 0,06, и поправки там ниже процента. Берётся C(25) = 1 +
(C(5)-1) * 0.06 в kit_recalc.

Крышка ЗАКРЫТА (5 см — внутри полости), фон записи — empty_shield_point5cm.
Долгоживущие дочерние отсечены окнами A/Z, чтобы поднятый порог времени не
распадал то, что в реальном источнике не распадается (Np-237 у Am-241!).
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
ZSRC = 93.0     # мм: 43 (крышка) + 50

# (метка, Z, A, окно nucleusLimits, распадов)
# Статистика: сильносуммирующим больше; Cs-137 — контроль C=1.
NUCS = [
    ("Am241", 95, 241, "241 241 93 95", 400000),
    ("Ba133", 56, 133, "133 133 55 56", 400000),
    ("Bi207", 83, 207, "207 207 82 83", 400000),
    ("Cd109", 48, 109, "109 109 47 48", 400000),
    ("Ce139", 58, 139, "139 139 57 58", 400000),
    ("Co57", 27, 57, "57 57 26 27", 400000),
    ("Co60", 27, 60, "60 60 27 28", 400000),
    ("Cs137", 55, 137, "137 137 55 56", 400000),
    ("Eu152", 63, 152, "152 152 62 64", 600000),
    ("Mn54", 25, 54, "54 54 24 25", 300000),
    ("Na22", 11, 22, "22 22 10 11", 400000),
    ("Th228", 90, 228, "208 228 81 90", 300000),
    ("Y88", 39, 88, "88 88 38 39", 400000),
    ("Zn65", 30, 65, "65 65 29 30", 300000),
]


def macro():
    t = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
         "/process/had/rdm/verbose 0",
         "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns",
         "/gps/particle ion", "/gps/energy 0 keV",
         "/gps/pos/type Point", "/gps/pos/centre 0 0 %.1f mm" % ZSRC,
         "/gps/ang/type iso"]
    for name, z, a, lim, n in NUCS:
        t += ["/process/had/rdm/nucleusLimits " + lim,
              "/gps/ion %d %d 0 0" % (z, a),
              "/g1s/outFile %s" % os.path.join(BUILD, "p5_%s.csv" % name),
              "/run/beamOn %d" % n]
    return "\n".join(t) + "\n"


if __name__ == "__main__":
    mp = os.path.join(BUILD, "p5_decay.mac")
    open(mp, "w", encoding="utf-8").write(macro())
    print("=== точечные распады на 5 см: %d нуклидов ===" % len(NUCS), flush=True)
    res = subprocess.run([os.path.join(BUILD, "g1s.exe"), mp, "shield"],
                         cwd=BUILD, capture_output=True, text=True,
                         encoding="utf-8", errors="replace")
    for ln in (res.stdout or "").splitlines():
        if ln.startswith(("RESULT", "EMIT")):
            print("  ", ln.strip(), flush=True)
    if res.returncode != 0:
        print("!! код возврата", res.returncode)
        print((res.stderr or "")[-1500:])
        sys.exit(1)
    print("готово")
