# -*- coding: utf-8 -*-
"""Целостность результатов: файл считается годным, если в шапке есть
N_primaries и E_prim_keV и присутствует хотя бы одна строка данных.
Обрезанные (процесс убит на середине записи) удаляются, чтобы докачка
пересчитала точку заново."""
import os
import sys

RESULTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       "..", "results"))


def check(path):
    n = e = False
    rows = 0
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("# N_primaries"):
                    n = True
                elif line.startswith("# E_prim_keV"):
                    e = True
                elif line and line[0].isdigit():
                    rows += 1
    except Exception:
        return False
    return n and e and rows > 0


def main():
    drop = "--fix" in sys.argv
    bad = []
    total = 0
    for root, _, files in os.walk(RESULTS):
        for fn in files:
            # только спектры прогонов; сводные таблицы постобработки пропускаем
            if not fn.endswith(".csv") or not (
                    fn.startswith("E") or fn.startswith("bg_") or fn.startswith("nuc_")):
                continue
            total += 1
            p = os.path.join(root, fn)
            if not check(p):
                bad.append(p)
    print("файлов %d, битых %d" % (total, len(bad)))
    for p in bad:
        print("   ", os.path.relpath(p, RESULTS))
        if drop:
            os.remove(p)
    if bad and drop:
        print("удалены, докачка пересчитает")


if __name__ == "__main__":
    main()
