# -*- coding: utf-8 -*-
"""Сверка ИСПУЩЕННЫХ квантов прогона с библиотекой МАГАТЭ.

Проверяется не отклик, а сам розыгрыш распада: сколько квантов каждой линии
выпущено на один распад. Если ограничение ряда (Stacking в main.cc) убивает
лишнее или лишнего не убивает, это видно здесь сразу — по чужим линиям или по
пропавшему каскаду. Проверка не зависит от геометрии: считаются кванты в
момент рождения, до какого-либо переноса.

Выходы сравниваются в ОКНЕ ±2 кэВ вокруг библиотечной энергии: гистограмма
испущенных пишется с шагом 1 кэВ, а близкие линии одного нуклида в неё
складываются — поэтому расхождение по слитым линиям ожидаемо и отмечается.

    python analysis/wt20_emit_check.py <каталог прогона>
"""
import csv
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.normpath(os.path.join(_HERE, "..", "reference", "nuclide-lines"))

# имя файла прогона -> имя файла библиотеки
NUC = {"Th232": "232th", "Ra228": "228ra", "Ac228": "228ac", "Th228": "228th",
       "Ra224": "224ra", "Rn220": "220rn", "Po216": "216po", "Pb212": "212pb",
       "Bi212": "212bi", "Tl208": "208tl", "Po212": "212po"}
WIN = 2.0        # окно сравнения, кэВ
MIN_I = 1.0      # сверяются линии интенсивнее этого, %


def read_emit(path):
    head, rows = {}, []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("# ").split("=", 1)
                head[k.strip()] = v.strip()
            continue
        if not ln or ln.startswith("E_keV"):
            continue
        a, b = ln.split(",")
        rows.append((float(a), float(b)))
    return head, rows


def read_lib(nuc):
    p = os.path.join(LIB, "%s_gammas.csv" % nuc)
    out = []
    if not os.path.exists(p):
        return out
    for r in csv.DictReader(io.open(p, encoding="utf-8")):
        try:
            out.append((float(r["energy"]), float(r["intensity"])))
        except (TypeError, ValueError):
            continue
    return out


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    d = sys.argv[1]
    print("%-8s %9s %9s %9s %8s  %s"
          % ("нуклид", "E, кэВ", "МАГАТЭ %", "Geant4 %", "откл. %", "примечание"))
    for key, nuc in NUC.items():
        p = os.path.join(d, "%s_emit.csv" % key)
        if not os.path.exists(p):
            continue
        head, rows = read_emit(p)
        n = float(head["N_primaries"])
        lib = [(e, i) for e, i in read_lib(nuc) if i >= MIN_I and e >= 20.0]
        lib.sort(key=lambda t: -t[1])
        full = read_lib(nuc)
        for e, i in lib[:8]:
            got = sum(c for ee, c in rows if abs(ee - e) <= WIN)
            pc = 100.0 * got / n
            # В окно ±2 кэВ попадают и СЛАБЫЕ соседи — сравнивать надо с их
            # суммой, иначе «превышение выхода» будет означать всего лишь
            # слипшиеся линии. Первый вариант проверки сравнивал с одной
            # линией и давал +17 % на Ac-228 338,32 ровно по этой причине.
            i_sum = sum(y for x, y in full if abs(x - e) <= WIN)
            near = [x for x, y in full if x != e and abs(x - e) <= WIN
                    and y >= 0.05]
            note = ("+%d соседей в окне" % len(near)) if near else ""
            dev = (pc - i_sum) / i_sum * 100.0 if i_sum else float("nan")
            print("%-8s %9.2f %9.2f %9.2f %8.1f  %s"
                  % (key, e, i_sum, pc, dev, note))
    return 0


if __name__ == "__main__":
    sys.exit(main())
