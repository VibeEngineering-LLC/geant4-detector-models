#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Выгружает уровни, заселяемые β⁺-ветвью распада (IAEA Live Chart of
Nuclides, `decay_rads`, `rad_types=bp`, поле `daughter_level_energy` +
`intensity_beta`) -- нужно, чтобы `enumerate_sum_peaks.py` мог перебирать
кандидатов суммирования С УЧАСТИЕМ аннигиляционного кванта (511 кэВ).

Зачем отдельный файл (ТЗ внешнего аудита -- Цензор, ретроспектива Теста 3
AmTiCsEu, п.9, 11.08.2026). Критерий прямого перебора `enumerate_pairs()`
(«конец перехода A = начало перехода B») требует уровня у ОБЕИХ линий пары
-- аннигиляционный квант физически НЕ ядерный переход, уровня не несёт,
и структурно НЕ МОЖЕТ быть найден этим критерием ни при каком запуске
(amticseu-remarks.md, разбор Sc-44 1668/2179 кэВ). Позитрон-излучающие
дочерние -- не редкий частный случай (в этом же тесте -- Sc-44), поэтому
дыра постоянная, не разовая.

Физика: аннигиляция происходит практически мгновенно после термализации
позитрона (масштаб пикосекунд) относительно временного разрешения
детектора -- для целей истинного совпадения (TCS) 511 кэВ считается
СОВПАДАЮЩИМ с любым гамма-переходом, НАЧИНАЮЩИМСЯ на уровне, заселённом
именно β⁺-ветвью (не EC -- EC не даёт позитрона и аннигиляционных
квантов). `intensity_ec` того же запроса НЕ используется здесь намеренно.

Запуск:
    python tools/fetch_beta_plus_feeds.py --nuclide Sc44=44SC

Пишет data/beta_plus_feeds_sum_peak_levels.csv (перезаписывает целиком,
детерминировано по входному списку -- коммитить как обычные данные).
"""
import argparse
import csv
import io
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_CSV = os.path.join(HERE, "..", "data", "beta_plus_feeds_sum_peak_levels.csv")
API = "https://nds.iaea.org/relnsd/v1/data?fields=decay_rads&nuclides={nuc}&rad_types=bp"


def fetch(nuc_iaea):
    url = API.format(nuc=nuc_iaea)
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        raw = resp.read().decode("utf-8")
    return list(csv.DictReader(io.StringIO(raw)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nuclide", action="append", default=[], metavar="KEY=IAEA",
                    help="напр. Ti44chain=44SC -- KEY обязан совпадать с "
                         "ключом группы в configs/<источник>.yaml и в "
                         "data/ensdf_<источник>_chain_lines.csv (nuclide-"
                         "колонка), НЕ с именем физического нуклида -- иначе "
                         "enumerate_annihilation_pairs() не сматчит; IAEA-код "
                         "-- дочерний нуклид, который САМ и испускает β+ "
                         "(для цепочки Ti-44->Sc-44 группа называется "
                         "Ti44chain, но β+ даёт Sc-44, поэтому 44SC)")
    args = ap.parse_args()
    if not args.nuclide:
        raise SystemExit("Передать хотя бы один --nuclide KEY=IAEA.")

    out_rows = []
    for spec in args.nuclide:
        nuc_key, nuc_iaea = spec.split("=", 1)
        rows = fetch(nuc_iaea)
        n_kept = 0
        for r in rows:
            ib_raw = (r.get("intensity_beta") or "").strip()
            lvl_raw = (r.get("daughter_level_energy") or "").strip()
            if not ib_raw or not lvl_raw:
                continue  # чистая EC-ветвь (позитрона нет) или уровень не оценён
            try:
                ib = float(ib_raw)
                lvl = float(lvl_raw)
            except ValueError:
                continue
            out_rows.append({
                "nuclide": nuc_key, "daughter_level_kev": lvl,
                "intensity_beta_plus_pct": ib,
                "source": ("IAEA-NDS/ENSDF;extracted=nds.iaea.org/relnsd/v1/data;"
                          "fields=decay_rads;rad_types=bp;date=%s"
                          % (r.get("Extraction_date") or "?")),
            })
            n_kept += 1
        print("%s (%s): %d β+ ветвей с уровнем" % (nuc_key, nuc_iaea, n_kept),
              file=sys.stderr)

    out_rows.sort(key=lambda r: (r["nuclide"], -r["intensity_beta_plus_pct"]))
    fieldnames = ["nuclide", "daughter_level_kev", "intensity_beta_plus_pct", "source"]
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
        f.write("# Уровни, заселяемые бета-плюс-ветвью (не EC) -- для перебора "
                "кандидатов суммирования с аннигиляционным квантом 511 кэВ, "
                "tools/enumerate_sum_peaks.py.enumerate_annihilation_pairs().\n")
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print("Записано %d строк в %s" % (len(out_rows), OUT_CSV), file=sys.stderr)


if __name__ == "__main__":
    main()
