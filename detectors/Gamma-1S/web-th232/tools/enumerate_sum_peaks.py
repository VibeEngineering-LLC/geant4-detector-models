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
DEFAULT_CC_CSV = os.path.join(HERE, "..", "data",
                              "conversion_coeff_sum_peak_levels.csv")
DEFAULT_BP_CSV = os.path.join(HERE, "..", "data",
                              "beta_plus_feeds_sum_peak_levels.csv")
ANNIHILATION_KEV = 511.0

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


def load_cc(cc_csv_path):
    """(nuclide, round(E_keV,3)) -> коэффициент полной внутренней конверсии.

    ИСПРАВЛЕНО 09.08.2026 (аудит Б4, задача #2, та же находка, что и Б1 в
    export_data.py._sum_peaks_with_fb, здесь была независимая
    незасинхронизированная копия): без этого файла F_B ниже считался ТОЛЬКО
    по гамма-депопуляции — заведомый недосчёт для низкоэнергетичных
    переходов в тяжёлых ядрах (см. export_data.py, докстринг
    _sum_peaks_with_fb, множитель занижения до ×4,7 на живом примере)."""
    cc_of = {}
    if not os.path.isfile(cc_csv_path):
        return cc_of
    with open(cc_csv_path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(ln for ln in f if not ln.startswith("#")):
            cc_raw = (row.get("conversion_coeff") or "").strip()
            if not cc_raw:
                continue
            try:
                cc_of[(row["nuclide"], round(float(row["E_keV"]), 3))] = \
                    float(cc_raw)
            except ValueError:
                continue
    return cc_of


def depopulation(transitions, nuclide, level, cc_of, tol=LEVEL_TOL):
    """F_B = Σ I_gamma·(1+CC) по всем переходам, НАЧИНАЮЩИХСЯ на уровне
    (полная депопуляция — гамма + внутренняя конверсия, не только гамма).

    Возвращает (fb, all_have_cc): all_have_cc=False, если хоть один
    учтённый переход не нашёл CC в cc_of (тогда для него подставлен 0 —
    поведение как до правки Б4 для этой конкретной строки, F_B может быть
    занижен, вызывающий код обязан это показать, не скрывать)."""
    fb, all_have_cc = 0.0, True
    for n, e, ip, s, end in transitions:
        if n != nuclide or abs(s - level) >= tol:
            continue
        cc = cc_of.get((n, round(e, 3)))
        if cc is None:
            all_have_cc = False
            cc = 0.0
        fb += ip * (1.0 + cc)
    return fb, all_have_cc


def load_beta_plus_feeds(bp_csv_path):
    """[(nuclide, level_kev, intensity_beta_plus_pct)] -- уровни, заселяемые
    β⁺-ветвью (fetch_beta_plus_feeds.py). Пусто, если файла нет (нуклиды
    источника не имеют β⁺-ветви -- законное отсутствие, не ошибка)."""
    out = []
    if not os.path.isfile(bp_csv_path):
        return out
    with open(bp_csv_path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(ln for ln in f if not ln.startswith("#")):
            try:
                out.append((row["nuclide"], float(row["daughter_level_kev"]),
                           float(row["intensity_beta_plus_pct"])))
            except (ValueError, KeyError):
                continue
    return out


def enumerate_annihilation_pairs(transitions, bp_feeds, cc_of):
    """Кандидаты суммирования С УЧАСТИЕМ аннигиляционного кванта (511 кэВ).

    Структурно НЕ находится enumerate_pairs() (см. докстринг модуля и
    fetch_beta_plus_feeds.py) -- аннигиляция не несёт ядерного уровня,
    критерий «конец A = начало B» к ней неприменим. Здесь уровень паре
    задаёт САМА β⁺-ветвь (daughter_level_kev из fetch_beta_plus_feeds.py),
    а не соседний гамма-переход: 511 кэВ считается совпадающим с ЛЮБЫМ
    гамма-переходом, начинающимся на этом уровне (тот же уровень, что
    заселён β⁺-распадом непосредственно перед аннигиляцией).

    score = I_beta+ * I_gamma / F_B(level) -- та же формула и тот же
    смысл, что у enumerate_pairs(), F_B -- полная депопуляция уровня
    (гамма+конверсия), не только рассматриваемый гамма-переход."""
    out = []
    for nuc, level, i_bp in bp_feeds:
        nuc_ts = [t for t in transitions if t[0] == nuc]
        fb, all_have_cc = depopulation(nuc_ts, nuc, level, cc_of)
        if fb <= 0:
            continue
        for n, e, ip, s, end in nuc_ts:
            if abs(s - level) >= LEVEL_TOL:
                continue
            score = i_bp * ip / fb
            suspicious = (not all_have_cc) and abs(fb - ip) < 0.01
            out.append({
                "nuclide": nuc, "e1_kev": ANNIHILATION_KEV, "e2_kev": e,
                "i1_pct": i_bp, "i2_pct": ip, "level_kev": level,
                "fb_pct": fb, "score": score,
                "cc_incomplete": not all_have_cc, "suspicious": suspicious,
                "annihilation": True,
            })
    out.sort(key=lambda r: -r["score"])
    return out


def enumerate_pairs(transitions, cc_of):
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
                fb, all_have_cc = depopulation(ts, nuc, level, cc_of)
                if fb <= 0:
                    continue
                score = ia * ib / fb
                # ИСПРАВЛЕНО 09.08.2026 (Б4, #2): раньше — score==I1/I2 не
                # отслеживался вовсе. Подозрительно точное совпадение F_B с
                # одним из двух I чаще всего значит «уровень депопулируется
                # ЕДИНСТВЕННЫМ учтённым переходом» — тот же класс риска, что
                # и Б1 (186,827 кэВ Ac-228, где F_B без CC был занижен в
                # 4,7 раза). Если для ЭТОГО перехода есть CC — не проблема
                # (F_B==I2 законно при CC=0, единственный радиационный канал
                # без конверсии). Если CC НЕТ (all_have_cc=False) — красный
                # флаг, score не заслуживает доверия без проверки вручную.
                suspicious = (not all_have_cc) and (
                    abs(fb - ia) < 0.01 or abs(fb - ib) < 0.01)
                out.append({
                    "nuclide": nuc, "e1_kev": ea, "e2_kev": eb,
                    "i1_pct": ia, "i2_pct": ib, "level_kev": level,
                    "fb_pct": fb, "score": score,
                    "cc_incomplete": not all_have_cc,
                    "suspicious": suspicious,
                })
    out.sort(key=lambda r: -r["score"])
    return out


def enumerate_two_hop(transitions):
    """Кандидаты «через один промежуточный переход, пропущенный парой»:
    A: X->Y, B: Y->Z, C: Z->W (три перехода одного каскада подряд) --
    пара (A, C) физически валидна как сумм-пик (оба кванта капчены в
    одном распаде, B либо ушёл, либо просто не входит в ЭТУ пару), но
    ПРЯМОЙ перебор enumerate_pairs() её не находит: конец A (Y) не
    совпадает с началом C (Z), совпадение только через посредника B.

    Найдено аудитом Б4 (09.08.2026): пример 277,371+2614,511=2891,9 кэВ
    Tl208 (через посредника 583,187 кэВ, уровень 3197,717->2614,529)
    пропускался.

    НЕ считает score здесь: у пары «через посредника» другая физика,
    чем у прямой пары (в цепочку из вероятности входит ЕЩЁ одна ветвь
    B|A, а какая доля этой ветви приходится именно на путь A->B->C, а не
    на A->(другой B')->... -- отдельный вопрос, которым эта функция не
    занимается). Отчёт -- список троек для РУЧНОЙ оценки, не готовый
    кандидат в SUM_PEAKS (см. докстринг модуля: скрипт не решает, только
    показывает)."""
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
                if a is b or abs(enda - sb) >= LEVEL_TOL:
                    continue
                for c in ts:
                    _, ec, ic, sc, endc = c
                    if c is a or c is b or abs(endb - sc) >= LEVEL_TOL:
                        continue
                    key = (nuc, round(min(ea, ec), 3), round(max(ea, ec), 3),
                           round(eb, 3))
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append({
                        "nuclide": nuc, "e1_kev": ea, "e2_kev": ec,
                        "e_mid_kev": eb, "i1_pct": ia, "i_mid_pct": ib,
                        "i2_pct": ic, "level_start": sa, "level_mid": sb,
                        "level_end": sc,
                    })
    # ИСПРАВЛЕНО 09.08.2026 (самоаудит перед докладом внешнему аудитору):
    # сортировка СНАЧАЛА по нуклиду (алфавит) хоронила флагманский пример
    # (277,371+2614,511 Tl208, ради которого функция и написана) за
    # несколькими сотнями записей Ac228 — вывод main() (топ-40) его
    # физически не показывал НИКОГДА, при том что сам список содержал
    # пару верно. Сортировка только по I1*I2 (без группировки по нуклиду)
    # — самые заметные кандидаты любого нуклида идут первыми.
    out.sort(key=lambda r: -r["i1_pct"] * r["i2_pct"])
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
    ap.add_argument("--cc-csv", default=DEFAULT_CC_CSV)
    ap.add_argument("--bp-csv", default=DEFAULT_BP_CSV)
    ap.add_argument("--min-score", type=float, default=MIN_SCORE)
    args = ap.parse_args()

    transitions = load_transitions(args.csv)
    cc_of = load_cc(args.cc_csv)
    bp_feeds = load_beta_plus_feeds(args.bp_csv)
    with open(args.config, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    known = known_pairs(cfg)

    # Порог применяется ТОЛЬКО к новым кандидатам (иначе уже принятая в
    # конфиг пара, чей score чуть ниже порога, ложно попадёт в раздел
    # «в конфиге, но перебор не нашёл» -- порог отсекает показ, не факт
    # существования пары).
    all_candidates = enumerate_pairs(transitions, cc_of)
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

    def fmt(c):
        flag = ""
        if c.get("suspicious"):
            flag = "  ⚠️ F_B без CC и подозрительно точно совпадает с I1/I2 — не доверять без проверки"
        elif c.get("cc_incomplete"):
            flag = "  (CC неполный для этого уровня — F_B может быть занижен)"
        return ("  %-6s %8.3f + %8.3f  I1=%.3f I2=%.3f F_B=%6.2f  уровень=%.3f  score=%.2f%s"
                % (c["nuclide"], c["e1_kev"], c["e2_kev"], c["i1_pct"], c["i2_pct"],
                   c["fb_pct"], c["level_kev"], c["score"], flag))

    print("=== УЖЕ В КОНФИГЕ (сверка -- перебор нашёл то же, что и ручной отбор) ===")
    for c in matched:
        print(fmt(c))

    print("\n=== НОВЫЕ КАНДИДАТЫ (не в конфиге -- требуют суждения: LNHB, физдопустимость) ===")
    for c in new:
        print(fmt(c))

    missing = known - {(c["nuclide"], round(min(c["e1_kev"], c["e2_kev"]), 3),
                         round(max(c["e1_kev"], c["e2_kev"]), 3)) for c in all_candidates}
    if missing:
        print("\n=== В КОНФИГЕ, НО ПЕРЕБОР НЕ НАШЁЛ (расхождение -- разобрать!) ===")
        for nuc, e1, e2 in sorted(missing):
            print("  %-6s %8.3f + %8.3f" % (nuc, e1, e2))

    # Б4 (#2): каскады через ОДИН пропущенный промежуточный переход —
    # прямой перебор enumerate_pairs() их принципиально не видит (общий
    # уровень A и C совпадает не напрямую, а через посредника B).
    two_hop = enumerate_two_hop(transitions)
    two_hop_new = [t for t in two_hop
                   if (t["nuclide"], round(min(t["e1_kev"], t["e2_kev"]), 3),
                       round(max(t["e1_kev"], t["e2_kev"]), 3)) not in known]
    if two_hop_new:
        print("\n=== ЧЕРЕЗ ОДИН ПРОПУЩЕННЫЙ ПЕРЕХОД (score НЕ считается — см. докстринг "
              "enumerate_two_hop, физика пары через посредника отличается от прямой; "
              "ручная оценка обязательна) ===")
        for t in two_hop_new[:40]:
            print("  %-6s %8.3f + %8.3f  (через %8.3f, I_mid=%.3f)  уровни %.3f->%.3f->%.3f"
                  % (t["nuclide"], t["e1_kev"], t["e2_kev"], t["e_mid_kev"],
                     t["i_mid_pct"], t["level_start"], t["level_mid"], t["level_end"]))
        if len(two_hop_new) > 40:
            print("  ... ещё %d, обрезано (не молча — см. это сообщение)"
                  % (len(two_hop_new) - 40))

    # П.9 ТЗ Цензора (11.08.2026): аннигиляция структурно невидима прямому
    # перебору (нет уровня) -- отдельная функция, отдельный источник данных
    # (fetch_beta_plus_feeds.py, rad_types=bp). Пусто, если у нуклидов
    # источника нет β⁺-ветви -- не ошибка, просто нечего искать.
    if bp_feeds:
        ann = enumerate_annihilation_pairs(transitions, bp_feeds, cc_of)
        ann_matched = [c for c in ann
                       if (c["nuclide"], round(min(c["e1_kev"], c["e2_kev"]), 3),
                           round(max(c["e1_kev"], c["e2_kev"]), 3)) in known]
        ann_new = [c for c in ann
                   if (c["nuclide"], round(min(c["e1_kev"], c["e2_kev"]), 3),
                       round(max(c["e1_kev"], c["e2_kev"]), 3)) not in known
                   and c["score"] >= args.min_score]
        print("\n=== С АННИГИЛЯЦИОННЫМ КВАНТОМ (511 кэВ), %d β+-ветвей "
              "в data/beta_plus_feeds_sum_peak_levels.csv ===" % len(bp_feeds))
        print("--- уже в конфиге ---")
        for c in ann_matched:
            print(fmt(c))
        print("--- новые кандидаты ---")
        for c in ann_new:
            print(fmt(c))
        if not ann_matched and not ann_new:
            print("  (кандидатов выше порога score>=%.2f нет)" % args.min_score)


if __name__ == "__main__":
    main()
