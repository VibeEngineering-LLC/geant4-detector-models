#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Сверяет ДВА независимых источника коэффициента конверсии на пересечении:

- data/conversion_coeff_sum_peak_levels.csv (decay_rads, родитель-специфичный,
  основной источник F_B с 09.08.2026, Б1);
- data/adopted_gammas_cross_check.csv (gammas, ADOPTED LEVELS, вне распада,
  fetch_adopted_gammas.py).

Ключ сверки — ЭНЕРГИЯ ГАММА-ПЕРЕХОДА (`E_keV`), допуск ±0,003 кэВ, БЕЗ
привязки к нуклиду. Так, а не по (нуклид, уровень): decay_rads ключует строку
по РОДИТЕЛЮ (в нашем списке — например, Ac228), а gammas — по нуклиду,
которому принадлежит сама схема уровней (для того же перехода это
ДОЧЕРНИЙ — Th228); плюс энергии УРОВНЕЙ в двух выгрузках расходятся на
0,01-0,02 кэВ (разные вводы одной оценки ENSDF), тогда как энергия самого
ПЕРЕХОДА в обеих совпадает практически всегда. Ложное совпадение по случайно
близкой энергии у РАЗНЫХ физических переходов теоретически возможно — набор
всего 13 нуклидов, при подозрении сверять руками по `level`/`start_jp`.

Не подменяет F_B — печатает согласие/расхождение там, где оба источника дают
число.

Запуск:
    python tools/check_cc_cross_source.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PRIMARY_CSV = os.path.join(HERE, "..", "data", "conversion_coeff_sum_peak_levels.csv")
CROSS_CSV = os.path.join(HERE, "..", "data", "adopted_gammas_cross_check.csv")
TOL_KEV = 0.003


def load(path, cc_field):
    out = []
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    for r in csv.DictReader(lines):
        cc_raw = (r.get(cc_field) or "").strip()
        if not cc_raw:
            continue
        try:
            e = float(r["E_keV"])
            cc = float(cc_raw)
        except ValueError:
            continue
        out.append({"nuclide": r["nuclide"], "E_keV": e, "level": r.get("level", ""), "cc": cc})
    return out


def main():
    primary = load(PRIMARY_CSV, "conversion_coeff")
    cross = load(CROSS_CSV, "tot_conv_coeff")

    print("Основной источник (decay_rads): %d строк с CC" % len(primary), file=sys.stderr)
    print("Кросс-источник (gammas, ADOPTED): %d строк с CC" % len(cross), file=sys.stderr)
    print("Допуск по энергии перехода: %.3f кэВ" % TOL_KEV, file=sys.stderr)
    print("", file=sys.stderr)

    matches = []
    for p in primary:
        cands = [c for c in cross if abs(c["E_keV"] - p["E_keV"]) <= TOL_KEV]
        if len(cands) == 1:
            matches.append((p, cands[0]))
        elif len(cands) > 1:
            print("! %s %.3f кэВ: %d кандидатов в кросс-источнике на этой энергии — "
                  "пропущено, неоднозначно (%s)"
                  % (p["nuclide"], p["E_keV"], len(cands),
                     ", ".join("%s/%s" % (c["nuclide"], c["cc"]) for c in cands)),
                  file=sys.stderr)

    print("Однозначных совпадений по энергии: %d" % len(matches), file=sys.stderr)
    print("", file=sys.stderr)

    if not matches:
        print("Совпадений нет — сверять нечего на этом наборе нуклидов.", file=sys.stderr)
        return

    header = "%-8s %-8s %10s %10s %10s %10s" % (
        "нуклид*", "нуклид**", "E, кэВ", "decay_rads", "gammas", "откл.,%")
    print(header)
    print("* по decay_rads (родитель)   ** по gammas (владелец уровня)")
    print("-" * len(header))
    devs = []
    for p, c in sorted(matches, key=lambda pc: pc[0]["E_keV"]):
        dev_pct = 100.0 * abs(p["cc"] - c["cc"]) / p["cc"] if p["cc"] else float("nan")
        devs.append(dev_pct)
        print("%-8s %-8s %10.3f %10.4g %10.4g %9.2f%%"
              % (p["nuclide"], c["nuclide"], p["E_keV"], p["cc"], c["cc"], dev_pct))

    print("")
    devs_sorted = sorted(devs)
    n = len(devs_sorted)
    if n % 2:
        median = devs_sorted[n // 2]
    else:
        median = (devs_sorted[n // 2 - 1] + devs_sorted[n // 2]) / 2.0
    print("Медиана отклонения: %.3f%%" % median, file=sys.stderr)
    print("Максимум отклонения: %.3f%%" % max(devs_sorted), file=sys.stderr)


if __name__ == "__main__":
    main()
