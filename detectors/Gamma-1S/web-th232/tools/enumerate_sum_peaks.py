#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Перебор кандидатов на пары суммирования (SUM_PEAKS), этап 4 задачи #175.

Инвентаризация (задача #176, агент B) отметила: `_sum_peaks_with_fb` в
export_data.py уже программно ПРОВЕРЯЕТ заданную пару (находит переходы,
сшивает по общему уровню, считает F_B) -- но САМИ пары до сих пор
подбираются вручную, просмотром библиотеки. Критерий последовательного
каскада («конец одного перехода = начало другого») в принципе пригоден и
для перебора, а не только для проверки. Здесь -- перебор.

ЧТО ЭТОТ СКРИПТ НЕ ДЕЛАЕТ (сознательно, директива «минимальное
вмешательство ИИ», 09.08.2026): не решает, какую пару включать в
config. Перекрёстная проверка независимым источником (LNHB) и суждение
о надёжности данных ENSDF при неполной оценке -- см. комментарии
export_data.py про Ac-228 214,850/674,750 -- остаются вне этого скрипта:
у него нет данных LNHB. Отбор кандидатов детерминирован (порог на
I1*I2/F_B), дальше -- отчёт для просмотра, не автоматическое решение.

Использование:
    python tools/enumerate_sum_peaks.py [--config PATH] [--csv PATH]
"""
import argparse
import csv
import os

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(HERE, "..", "data", "ensdf_th232_chain_lines.csv")
DEFAULT_CONFIG = os.path.join(HERE, "..", "configs", "th232.yaml")

LEVEL_TOL = 0.02   # кэВ -- совпадение уровней (тот же допуск, что в export_data.py)
MIN_SCORE = 0.5    # порог I1*I2/F_B (%^2/%), ниже -- шум, не показывается вовсе


def load_transitions(csv_path):
    """[(nuclide, E_keV, I_percent, start, end)] только гамма-переходы с
    разобранным уровнем start->end."""
    out = []
    with open(csv_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("line_type") != "gamma":
                continue
            lvl = (row.get("level") or "").strip()
            if not lvl or "->" not in lvl:
                continue
            try:
                start_s, end_s = lvl.split("->")
                start, end = float(start_s), float(end_s)
                e = float(row["E_keV"])
                ip = float(row["I_percent"])
            except ValueError:
                continue
            out.append((row["nuclide"], e, ip, start, end))
    return out


def depopulation(transitions, nuclide, level, tol=LEVEL_TOL):
    """F_B: суммарная I_percent всех переходов, НАЧИНАЮЩИХСЯ на уровне."""
    return sum(ip for n, e, ip, s, end in transitions
               if n == nuclide and abs(s - level) < tol)


def enumerate_pairs(transitions):
    """Все пары (A, B) одного нуклида, где конец A = начало B (в допуске)
    -- то есть A, затем B в одном каскаде через общий промежуточный
    уровень. Пары с A==B и зеркальные дубли (A,B)/(B,A) как разные
    физические переходы разрешены -- это не одно и то же (например
    A: 968,972->0 и B: 968,972->X -- обе валидны как «первая половина»
    разных каскадов), дубль по (E1,E2) отбрасывается на выходе.
    """
    by_nuclide = {}
    for t in transitions:
        by_nuclide.setdefault(t[0], []).append(t)

    seen = set()
    out = []
    for nuc, ts in by_nuclide.items():
        for a in ts:
            _, ea, ia, sa, enda = a
            for b in ts:
                _, eb, ib, sb, endb = b
                if a is b:
                    continue
                if abs(enda - sb) >= LEVEL_TOL:
                    continue
                # общий уровень = конец A = начало B
                level = sb
                key = (nuc, round(min(ea, eb), 3), round(max(ea, eb), 3))
                if key in seen:
                    continue
                seen.add(key)
                fb = depopulation(ts, nuc, level)
                if fb <= 0:
                    continue
                score = ia * ib / fb
                out.append({
                    "nuclide": nuc, "e1_kev": ea, "e2_kev": eb,
                    "i1_pct": ia, "i2_pct": ib, "level_kev": level,
                    "fb_pct": fb, "score": score,
                })
    out.sort(key=lambda r: -r["score"])
    return out


def known_pairs(cfg):
    """Множество (нуклид, E1, E2) уже принятых в конфиг -- для сверки,
    порядок энергий не важен."""
    out = set()
    for s in cfg.get("sum_peaks", []):
        e1, e2 = round(s["e1_kev"], 3), round(s["e2_kev"], 3)
        out.add((s["nuclide"], min(e1, e2), max(e1, e2)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--min-score", type=float, default=MIN_SCORE)
    args = ap.parse_args()

    transitions = load_transitions(args.csv)
    with open(args.config, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    known = known_pairs(cfg)

    # Порог применяется ТОЛЬКО к новым кандидатам (иначе уже принятая в
    # конфиг пара, чей score чуть ниже порога, ложно попадёт в раздел
    # «в конфиге, но перебор не нашёл» -- порог отсекает показ, не факт
    # существования пары).
    all_candidates = enumerate_pairs(transitions)
    matched, new = [], []
    for c in all_candidates:
        key = (c["nuclide"], round(min(c["e1_kev"], c["e2_kev"]), 3),
               round(max(c["e1_kev"], c["e2_kev"]), 3))
        if key in known:
            matched.append(c)
        elif c["score"] >= args.min_score:
            new.append(c)

    print("порог score=I1*I2/F_B >= %.2f (только для новых кандидатов);"
          " в конфиге найдено %d/%d, новых кандидатов выше порога %d\n"
          % (args.min_score, len(matched), len(known), len(new)))

    print("=== УЖЕ В КОНФИГЕ (сверка -- перебор нашёл то же, что и ручной отбор) ===")
    for c in matched:
        print("  %-6s %8.3f + %8.3f  I1=%.3f I2=%.3f F_B=%6.2f  уровень=%.3f  score=%.2f"
              % (c["nuclide"], c["e1_kev"], c["e2_kev"], c["i1_pct"], c["i2_pct"],
                 c["fb_pct"], c["level_kev"], c["score"]))

    print("\n=== НОВЫЕ КАНДИДАТЫ (не в конфиге -- требуют суждения: LNHB, физдопустимость) ===")
    for c in new:
        print("  %-6s %8.3f + %8.3f  I1=%.3f I2=%.3f F_B=%6.2f  уровень=%.3f  score=%.2f"
              % (c["nuclide"], c["e1_kev"], c["e2_kev"], c["i1_pct"], c["i2_pct"],
                 c["fb_pct"], c["level_kev"], c["score"]))

    missing = known - {(c["nuclide"], round(min(c["e1_kev"], c["e2_kev"]), 3),
                         round(max(c["e1_kev"], c["e2_kev"]), 3)) for c in all_candidates}
    if missing:
        print("\n=== В КОНФИГЕ, НО ПЕРЕБОР НЕ НАШЁЛ (расхождение -- разобрать!) ===")
        for nuc, e1, e2 in sorted(missing):
            print("  %-6s %8.3f + %8.3f" % (nuc, e1, e2))


if __name__ == "__main__":
    main()
