# -*- coding: utf-8 -*-
"""Коэффициенты ослабления матриц на ТОЙ ЖЕ сетке энергий, что и всё остальное.

Зачем отдельный драйвер. Список энергий существовал в трёх копиях: в двух
драйверах на питоне (сведены в `grid_energies.py`) и ТРЕТЬЕЙ — прямо в
`geometry/mucalc.cc`. Копии разъезжаются молча: сетку расширили краями
паспортных зон 45,3 и 3552,5 кэВ, а таблицы mu остались на прежних двадцати
точках. Проявилось это не сразу и не там — подгонка эффективной толщины
упала с «нет mu для E = 45,300 кэВ», то есть через два шага после причины.

Теперь список один: этот драйвер пишет его во временный файл и передаёт
mucalc аргументом.

    python detectors/Gamma-1S/drivers/run_mu.py
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grid_energies import LINES  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
EXE = os.path.join(BUILD, "mucalc.exe")


def main():
    if not os.path.exists(EXE):
        raise SystemExit(
            "Не найден %s — соберите mucalc (см. common/cmake) или укажите\n"
            "G4MODELS_BUILD_GAMMA_1S на каталог, где он лежит." % EXE)
    lst = os.path.join(BUILD, "mu_energies.txt")
    with open(lst, "w", encoding="ascii") as fh:
        fh.write("\n".join("%.3f" % e for e in LINES) + "\n")
    print("сетка: %d энергий, %.1f…%.1f кэВ" % (len(LINES), min(LINES),
                                                max(LINES)), flush=True)
    r = subprocess.run([EXE, lst], cwd=BUILD, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    for ln in (r.stdout or "").splitlines():
        if ln.startswith("RESULT"):
            print("  ", ln.strip())
    if r.returncode != 0:
        print((r.stderr or "")[-1500:])
        raise SystemExit("mucalc вернул %d" % r.returncode)
    # Проверка вслух: таблица обязана накрыть сетку целиком, иначе
    # самопоглощение снова упадёт через два шага после причины.
    for fn in ("mu_oisn16.csv", "mu_water.csv"):
        p = os.path.join(BUILD, fn)
        got = [float(l.split(",")[0]) for l in open(p, encoding="utf-8")
               if l[:1].isdigit()]
        miss = [e for e in LINES if not any(abs(e - g) < 0.001 for g in got)]
        print("   %-14s точек %d%s" % (fn, len(got),
                                       "" if not miss else
                                       "  НЕ ХВАТАЕТ: " + ", ".join(
                                           "%.1f" % m for m in miss)))
        if miss:
            raise SystemExit("таблица mu не накрывает сетку")
    print("готово")


if __name__ == "__main__":
    main()
