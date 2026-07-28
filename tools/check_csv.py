# -*- coding: utf-8 -*-
"""Проверка выгруженных таблиц штатным парсером CSV.

ЗАЧЕМ. Файлы в results/ — это то, что читает внешний потребитель, и читает он
их обычным csv.DictReader. Значит и проверять их надо тем же: прочитать и
сверить число полей в каждой строке с шапкой.

ОТКУДА ВЗЯЛОСЬ. export_curves.py писал строки склейкой ",".join(...) без
всякого экранирования. Значение geometry у точечных сеток — «точечный, 25 см»,
с запятой внутри. Все 48 строк обеих точечных кривых в efficiency_curves.csv
получили 14 полей при шапке в 13, и csv.DictReader не падал, а ТИХО сдвигал
всё правее на одно поле: E_keV становился словом «ОТКРЫТА», eps_net — числом
энергии. Объёмные кривые при этом читались верно (в их названиях запятых нет),
поэтому дефект и прожил столько времени.

Три строки проверки ловят весь класс. Запускать после каждой выгрузки; в
протокол проверок — рядом с check_paths.py.

    python tools/check_csv.py            # всё дерево detectors/*/results
    python tools/check_csv.py путь.csv   # один файл

Код возврата: 0 — чисто, 1 — есть битые файлы.
"""
import csv
import glob
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(path):
    """[(номер строки, сколько полей)] у строк, не совпавших с шапкой."""
    with open(path, encoding="utf-8", newline="") as fh:
        rows = [(i + 1, r) for i, r in enumerate(csv.reader(fh))
                if r and not r[0].lstrip().startswith("#")]
    if not rows:
        return []
    n = len(rows[0][1])
    return [(i, len(r)) for i, r in rows if len(r) != n]


def main(argv):
    if argv:
        files = argv
    else:
        files = sorted(glob.glob(os.path.join(REPO, "detectors", "*",
                                              "results", "*.csv")))
    bad = 0
    for p in files:
        errs = check(p)
        if errs:
            bad += 1
            rel = os.path.relpath(p, REPO)
            print("БИТ %s: %d строк не совпали с шапкой" % (rel, len(errs)))
            for i, n in errs[:3]:
                print("      строка %d: %d полей" % (i, n))
            print("      причина обычно одна: запятая внутри значения и "
                  "запись без экранирования (csv.writer вместо join)")
    print("проверено файлов: %d, битых: %d" % (len(files), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
