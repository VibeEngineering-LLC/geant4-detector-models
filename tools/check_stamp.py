# -*- coding: utf-8 -*-
"""Сторож штампа провенанса и объявления наблюдаемой в results/.

ЧТО ПРОВЕРЯЕТ. Для каждой таблицы results/*.csv:
  * есть ли штамп вообще (строки «#@»);
  * объявлены ли ОБЯЗАТЕЛЬНЫЕ ключи наблюдаемой;
  * каков вердикт по входам (`src.inputs_verdict`), записанный производителем.

ЧЕГО НЕ ПРОВЕРЯЕТ. Правильность самого объявления: сказать «окно ±1 ПШПВ» и
взять ±6 кэВ штамп не мешает. Сторож ловит МОЛЧАНИЕ, а не ложь — но молчание и
было тем, что стоило вывода: три расходящиеся таблицы одной величины лежали
рядом и внешне ничем не отличались друг от друга.

ПОЧЕМУ НЕ ПАДАЕТ НА ДОЛГЕ. Незаштампованных таблиц сейчас большинство, и
сторож, падающий с первого дня, будет отключён в тот же день. Поэтому: код
возврата 1 — только если заштампованная таблица заявила ПЛОХОЙ вердикт по
входам (stale/mixed) или потеряла обязательные ключи. Долг печатается числом,
чтобы он убывал видимо, а не «когда-нибудь».

Запуск:
    python tools/check_stamp.py [--strict] [детектор ...]
--strict — считать отсутствие штампа тоже ошибкой (для конца миграции).
"""
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common", "py"))
import stamp  # noqa: E402

NEED = ("obs.quantity", "obs.area", "obs.window", "obs.shelf", "obs.blurred")
BAD_VERDICTS = ("stale", "mixed")

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    ".."))


def main(argv):
    strict = "--strict" in argv
    dets = [a for a in argv[1:] if not a.startswith("-")]
    if not dets:
        dets = sorted(os.path.basename(p) for p in
                      glob.glob(os.path.join(REPO, "detectors", "*"))
                      if os.path.isdir(p))
    stamped, plain, errs, warns = [], [], [], []
    for det in dets:
        for p in sorted(glob.glob(os.path.join(
                REPO, "detectors", det, "results", "*.csv"))):
            rel = os.path.relpath(p, REPO).replace("\\", "/")
            st = stamp.read_table_stamp(p)
            if not st:
                plain.append(rel)
                continue
            stamped.append(rel)
            miss = [k for k in NEED if k not in st]
            if miss:
                errs.append("%s: нет ключей наблюдаемой: %s"
                            % (rel, ", ".join(miss)))
            v = st.get("src.inputs_verdict")
            if v in BAD_VERDICTS:
                errs.append("%s: вердикт по входам «%s» (отпечаток %s)"
                            % (rel, v, st.get("src.inputs_sha1", "?")))
            elif v == "unstamped":
                warns.append("%s: входы без штампа (%s из %s) — числа"
                             " непрослежены"
                             % (rel, st.get("src.inputs_unstamped", "?"),
                                st.get("src.inputs_n", "?")))

    total = len(stamped) + len(plain)
    print("Таблиц %d: со штампом %d, без штампа %d."
          % (total, len(stamped), len(plain)))
    for w in warns:
        print("  ?? %s" % w)
    for e in errs:
        print("  !! %s" % e)
    if plain:
        print("  долг миграции (%d) — первые: %s"
              % (len(plain), ", ".join(plain[:6])))
    if errs:
        print("\nОШИБКИ: %d." % len(errs))
        return 1
    if strict and plain:
        print("\n--strict: %d таблиц без штампа." % len(plain))
        return 1
    print("\nОшибок нет." + (" Предупреждений: %d." % len(warns) if warns
                             else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
