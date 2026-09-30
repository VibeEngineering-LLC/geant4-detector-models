# -*- coding: utf-8 -*-
r"""Проверка калибровки калибровочного спектра Th-232 GS2020 PRO MAX по донору
becqmoni.py (common/py) — #CAL-0: своя калибровка проверяется каждый раз.
Источник: spectrum (19).xml (папка загрузок оператора), оператор 2026-09-25,
'референсные спектры ториевого калибровочного источника, положение прибора то же'."""

import sys
sys.stdout.reconfigure(encoding="utf-8")

sys.path.insert(0, r"<REPOS>\geant4-detector-models\common\py")
import becqmoni as bm

PATH = r"<DOWNLOADS>\spectrum (19).xml"

LINES = [238.632, 583.187, 727.330, 860.564, 911.204, 968.971, 1588.19, 2614.511]

NAMES = {
    238.632: "Pb-212",
    583.187: "Tl-208",
    727.330: "Bi-212",
    860.564: "Tl-208",
    911.204: "Ac-228",
    968.971: "Ac-228",
    1588.19: "Ac-228",
    2614.511: "Tl-208"
}

sp, _bg = bm.read(PATH)

print("live=%.1f real=%.1f cal=%r" % (sp.live, sp.real, sp.cal))

print("сумма отсчётов: %d" % sp.n.sum())
print()

print("%-10s %-9s %-9s %-8s %-10s" % ("линия", "ном.кэВ", "центр.кэВ", "d,кэВ", "d/ПШПВ"))

for E0 in LINES:
    try:
        res = bm.peak_find(sp, E0)
        c, f, _h = res
        d = c - E0
        print("%-10s %-9.2f %-9.2f %-8.2f %-10.3f" % (NAMES[E0], E0, c, d, d / f if f else float("nan")))
    except Exception as exc:  # noqa: BLE001
        print("%-10s %-9.2f ОШИБКА: %s" % (NAMES[E0], E0, exc))
