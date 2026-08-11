#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Выгружает полную библиотеку линий эмиссии (гамма И характеристический
рентген СОБСТВЕННОГО распада нуклида, не путать с рентгеном дочернего атома
от внутренней конверсии) для заданного списка нуклидов из IAEA Live Chart of
Nuclides -- ОДИН запрос на нуклид ПОКРЫВАЕТ ОБА `rad_types` (`g` и `x`) по
умолчанию, без отдельного напоминания оператора.

Зачем этот файл (ТЗ внешнего аудита -- Цензор, ретроспектива Теста 3
AmTiCsEu, 11.08.2026, п.1). Первая версия библиотеки этого источника
запрашивала только `rad_types=g` -- сильнейшие K-рентген-линии дочерних
(Eu-152 Kα1 40,117 кэВ, I=37,8 % -- сильнее большинства гамма-линий страницы)
отсутствовали ПОЛНОСТЬЮ, что было главной причиной 55%-го расхождения
методов по Am-241 (amticseu-remarks.md §10, №5). Причина отсутствия --
человеческая: сборка библиотеки для НОВОГО источника делалась ad-hoc
(разовым запросом внутри сессии, без сохранённого инструмента), и запрос
`rad_types=x` для этого источника никто не написал, пока оператор не
заметил дыру по результату. `fetch_conversion_coeff.py`/th232/ra226 эту
же дыру не имеют -- их сборка (более ранняя, до этого скилла) `rad_types=x`
запрашивала. Этот файл делает «оба типа по умолчанию» СВОЙСТВОМ ИНСТРУМЕНТА,
не памятью агента.

Формат вывода -- совместим с data/ensdf_<источник>_chain_lines.csv
(export_data.py.load_full_library()): nuclide,E_keV,I_percent,unc_I_percent,
source,level,line_type. `line_type` — "gamma" для ядерных переходов
(включая аннигиляцию — уровня нет, но это НЕ рентген), "xray" ТОЛЬКО для
характеристического рентгена ДОЧЕРНЕГО атома, встроенного в ответ
`rad_types=g` (внутренняя конверсия) -- его load_full_library() исключает
намеренно (див. докстринг export_data.py, физика уже посчитана Geant4 в
шаблоне метода 1, повторный учёт в методе 2 задвоил бы). Рентген СОБСТВЕННОГО
распада (`rad_types=x`, например Eu-152 K-рентген Sm при EC) размечен
"gamma" ПРЕДНАМЕРЕННО -- он НЕ встроен ни в чей чужой шаблон и обязан
попасть в библиотеку метода 2 как рабочая линия (тот же принцип, что уже
применён к 511 кэВ аннигиляции в этом же файле раньше).

Проверка полноты (НОВОЕ, п.1 ТЗ) -- `check_completeness()`: для каждого
нуклида, чей режим распада (переданный явно в NUCLIDES, не угадывается)
предполагает характеристический рентген (EC/IC), требует хотя бы одну
xray-по-смыслу строку (rad_types=x с I_percent выше порога) -- иначе явная
ошибка `SystemExit`, не тихий пропуск. Разметка режима распада -- РУКАМИ по
конфигу источника (decay_mode: "EC"/"IT"/"beta-"/"beta+"/"alpha" и т.п.),
не выводится из самого запроса -- у скрипта нет физической модели, чтобы
это знать заранее.

Запуск:
    python tools/fetch_line_library.py --source amticseu
(добавить `--source <id>` -- берёт список нуклидов из configs/<id>.yaml,
поле `nuclides[].key`/`nuclides[].iaea_id`/`nuclides[].decay_mode` -- НОВЫЕ
поля схемы конфига, см. SKILL.md v1.4.0 "Config-driven схема".)

Пишет data/ensdf_<source>_chain_lines.csv (перезаписывает целиком,
детерминировано по входному списку нуклидов -- коммитить как обычные данные,
как fetch_conversion_coeff.py).
"""
import argparse
import csv
import io
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "..", "data")

API = ("https://nds.iaea.org/relnsd/v1/data?fields=decay_rads"
       "&nuclides={nuc}&rad_types={rt}")

# Режимы распада, для которых характеристический рентген СОБСТВЕННОГО
# распада (rad_types=x) физически ожидается -- используются в
# check_completeness(). EC/IC/beta+ (позитронный) всегда оставляют
# вакансию во внутренней оболочке дочернего атома; alpha/beta- сами по
# себе рентген не дают (хотя дочерний атом может светить рентгеном от
# ВНУТРЕННЕЙ КОНВЕРСИИ гамма-перехода -- это ДРУГОЙ канал, уже покрыт
# отдельно `rad_types=g`-embedded xray, не этой проверкой).
XRAY_EXPECTED_MODES = {"EC", "EC+beta+", "beta+", "IT"}


def fetch(nuc_iaea, rad_types):
    url = API.format(nuc=nuc_iaea, rt=rad_types)
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        raw = resp.read().decode("utf-8")
    return list(csv.DictReader(io.StringIO(raw)))


def _rows_from_query(nuc_key, nuc_iaea, rad_types, extraction_date_box):
    rows = fetch(nuc_iaea, rad_types)
    out = []
    for r in rows:
        e_raw, i_raw = (r.get("energy") or "").strip(), (r.get("intensity") or "").strip()
        if not e_raw or not i_raw:
            continue  # энергия/интенсивность не оценены ENSDF -- подставить нечем
        try:
            E, I = float(e_raw), float(i_raw)
        except ValueError:
            continue
        sl, el = (r.get("start_level_energy") or "").strip(), (r.get("end_level_energy") or "").strip()
        level = "%s->%s" % (sl, el) if sl and el else ""
        unc_raw = (r.get("unc_i") or "").strip()
        extraction_date_box[0] = r.get("Extraction_date") or extraction_date_box[0]
        # rad_types=g отдаёт вместе с ядерными гамма характеристический
        # рентген ДОЧЕРНЕГО атома (внутренняя конверсия) -- API размечает
        # это отдельным полем `decay` / по отсутствию уровня ПРИ типе "g";
        # эвристика: если запрошено rad_types=g и в ответе явно указан тип
        # строки как рентген (поле `type`, когда есть) -- xray, иначе gamma.
        row_type = (r.get("type") or "").strip().lower()
        if rad_types == "g" and ("x" in row_type or "ray" in row_type):
            line_type = "xray"
        else:
            # rad_types=x (рентген СОБСТВЕННОГО распада) и обычные rad_types=g
            # ядерные переходы -- оба "gamma" по конвенции этого файла (см.
            # докстринг: рентген своего распада обязан войти в библиотеку
            # как рабочая линия, не исключаться).
            line_type = "gamma"
        out.append({
            "nuclide": nuc_key, "E_keV": E, "I_percent": I,
            "unc_I_percent": unc_raw,
            "source": ("IAEA-NDS/ENSDF;extracted=nds.iaea.org/relnsd/v1/data;"
                      "fields=decay_rads;rad_types=%s;date=%s"
                      % (rad_types, extraction_date_box[0] or "?")),
            "level": level, "line_type": line_type,
        })
    return out


def fetch_nuclide_library(nuc_key, nuc_iaea):
    """Обе выгрузки (g И x) для одного нуклида -- всегда обе, без флага."""
    box = [None]
    rows_g = _rows_from_query(nuc_key, nuc_iaea, "g", box)
    rows_x = _rows_from_query(nuc_key, nuc_iaea, "x", box)
    return rows_g + rows_x


def check_completeness(all_rows, nuclide_modes):
    """Для каждого нуклида с ожидаемым рентгеном СОБСТВЕННОГО распада
    (XRAY_EXPECTED_MODES) требует хотя бы одну строку rad_types=x с
    I_percent >= 1 -- иначе SystemExit с точным указанием, чего не хватает.
    `nuclide_modes` -- {nuc_key: decay_mode_str}, из конфига источника."""
    by_nuc = {}
    for r in all_rows:
        by_nuc.setdefault(r["nuclide"], []).append(r)
    problems = []
    for nuc_key, mode in nuclide_modes.items():
        if mode not in XRAY_EXPECTED_MODES:
            continue
        rows = by_nuc.get(nuc_key, [])
        has_own_xray = any(
            "rad_types=x" in r["source"] and r["I_percent"] >= 1.0
            for r in rows)
        if not has_own_xray:
            problems.append(
                "%s (режим распада %s -- рентген СОБСТВЕННОГО распада "
                "ожидается): ни одной строки rad_types=x с I>=1%% не "
                "найдено. Библиотека НЕПОЛНА -- проверить нуклид в IAEA "
                "Live Chart вручную, не публиковать как есть."
                % (nuc_key, mode))
    if problems:
        raise SystemExit(
            "Проверка полноты библиотеки провалена (%d нуклид(ов)):\n  %s"
            % (len(problems), "\n  ".join(problems)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True,
                    help="id источника (совпадает с configs/<id>.yaml)")
    ap.add_argument("--nuclide", action="append", default=[], metavar="KEY=IAEA:MODE",
                    help="переопределить/добавить нуклид явно, напр. "
                         "Eu152=152EU:EC -- для источников без поля "
                         "iaea_id/decay_mode в конфиге (переходный режим, "
                         "пока схема конфига не расширена)")
    args = ap.parse_args()

    if not args.nuclide:
        raise SystemExit(
            "Список нуклидов пуст -- передать --nuclide KEY=IAEA:MODE для "
            "каждого нуклида источника (пока схема configs/<id>.yaml не "
            "несёт iaea_id/decay_mode явно -- см. SKILL.md v1.4.0, задача "
            "довести конфиг до самодостаточного).")

    nuc_iaea, nuc_modes = {}, {}
    for spec in args.nuclide:
        key, rest = spec.split("=", 1)
        iaea_id, mode = rest.split(":", 1)
        nuc_iaea[key] = iaea_id
        nuc_modes[key] = mode

    out_rows = []
    for nuc_key, nuc_id in nuc_iaea.items():
        rows = fetch_nuclide_library(nuc_key, nuc_id)
        n_g = sum(1 for r in rows if "rad_types=g" in r["source"])
        n_x = sum(1 for r in rows if "rad_types=x" in r["source"])
        print("%s (%s): %d строк (g=%d, x=%d)" % (nuc_key, nuc_id, len(rows), n_g, n_x),
              file=sys.stderr)
        out_rows.extend(rows)

    check_completeness(out_rows, nuc_modes)

    out_rows.sort(key=lambda r: (r["nuclide"], r["E_keV"]))
    out_csv = os.path.join(DATA_DIR, "ensdf_%s_chain_lines.csv" % args.source)
    fieldnames = ["nuclide", "E_keV", "I_percent", "unc_I_percent", "source",
                  "level", "line_type"]
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print("Записано %d строк в %s (проверка полноты пройдена)"
          % (len(out_rows), out_csv), file=sys.stderr)


if __name__ == "__main__":
    main()
