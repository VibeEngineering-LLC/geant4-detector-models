# -*- coding: utf-8 -*-
"""Пересчитать ВСЁ, что зависит от сеток моноэнергий, одной командой.

Зачем это отдельным скриптом. Сетка энергий расширяется (края паспортных зон,
см. drivers/grid_energies.py), и после каждого расширения надо перегнать не один
скрипт, а тринадцать: кривые, пересчёт комплекта, самопоглощение, суммирование,
МИА, сводную страницу. Пока это делалось руками, часть чисел в отчёте оставалась
от предыдущей сетки, и заметить это было нечем — все они выглядят одинаково
правдоподобно.

Порядок важен: export_curves готовит results/eff_*.csv, на которые опираются
остальные, а build_web собирает страницу последним, когда все таблицы на месте.

    python detectors/Gamma-1S/analysis/recalc_all.py
    python detectors/Gamma-1S/analysis/recalc_all.py --only kit_recalc
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths  # noqa: E402

# Порядок — зависимостями, а не по алфавиту.
STEPS = [
    ("export_curves.py", [], "кривые ППП в results/eff_*.csv и манифест"),
    # ВНИМАНИЕ: таблицы mu считаются НЕ здесь, а drivers/run_mu.py — им нужен
    # собранный mucalc и окружение Geant4. Если сетка расширялась, запустить
    # его ДО этого скрипта, иначе selfabs_fit падает с «нет mu для E = …».
    ("selfabs_fit.py", [], "эффективная толщина d_eff по парам плотностей"),
    ("summing.py", [], "поправки на каскадное суммирование"),
    ("compare_lsrm.py", [], "сверка с паспортной кривой маринелли"),
    ("compare_cups.py", [], "то же по кюветам"),
    ("compare_point.py", [], "то же по точечным геометриям"),
    ("kit_recalc.py", [], "пересчёт объёмных записей комплекта"),
    ("point_recalc.py", [], "пересчёт точечных записей комплекта"),
    ("deconv.py", [], "связанная деконволюция по группам линий"),
    ("deconv_balance.py", [], "баланс пика и континуума: измерение против модели"),
    ("continuum.py", [], "континуум расчётного спектра вне пиков"),
    ("mda.py", [], "минимальная измеряемая активность"),
    ("build_web.py", [], "сводная страница docs/gamma-1s/"),
]


def main():
    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1]
    build = paths.build("Gamma-1S")
    if not os.path.isdir(os.path.join(str(build), "grid")):
        raise SystemExit(
            "Нет каталога сеток %s/grid.\n"
            "Сначала посчитайте сетки:\n"
            "    python detectors/Gamma-1S/drivers/run_grid.py\n"
            "    python detectors/Gamma-1S/drivers/run_all_grids.py" % build)
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    fails = []
    for name, args, note in STEPS:
        if only and only not in name:
            continue
        print("=" * 70, flush=True)
        print("%-20s %s" % (name, note), flush=True)
        r = subprocess.run([sys.executable, os.path.join(HERE, name)] + args,
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", env=env)
        tail = [ln for ln in (r.stdout or "").splitlines() if ln.strip()][-6:]
        for ln in tail:
            print("   ", ln, flush=True)
        if r.returncode != 0:
            fails.append(name)
            print("!! КОД ВОЗВРАТА %d" % r.returncode, flush=True)
            print((r.stderr or "")[-1200:], flush=True)
    print("=" * 70)
    if fails:
        # Молчаливый провал одного шага из тринадцати — ровно тот случай, из
        # которого берутся числа «от прошлой сетки». Заканчиваем ненулевым кодом.
        raise SystemExit("НЕ ПРОШЛИ: " + ", ".join(fails))
    print("пересчитано всё")


if __name__ == "__main__":
    main()
