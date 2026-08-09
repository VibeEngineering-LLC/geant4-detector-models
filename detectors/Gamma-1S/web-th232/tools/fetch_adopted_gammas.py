#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Выгружает СХЕМУ УРОВНЕЙ ВНЕ РАСПАДА (ENSDF ADOPTED LEVELS/GAMMAS) для тех же
13 нуклидов, что в `fetch_conversion_coeff.py`, — из IAEA Live Chart of
Nuclides REST, ПОЛЕ `fields=gammas` (не `decay_rads`).

Зачем отдельный источник, если `data/conversion_coeff_sum_peak_levels.csv`
уже есть. Два источника закрывают РАЗНЫЕ вопросы:

- `decay_rads&rad_types=g` (fetch_conversion_coeff.py) — гамма-линии,
  испускаемые ПРИ РАСПАДЕ конкретного родителя (нужен `NUCLIDES[родитель]`).
  Это наш ОСНОВНОЙ источник CC для F_B, используется с 09.08.2026 (Б1).
- `gammas` (этот файл) — принятая схема уровней ДОЧЕРНЕГО нуклида САМОГО ПО
  СЕБЕ, вне привязки к родителю (ADOPTED LEVELS ENSDF). Ключуется по нуклиду
  напрямую — не нужно знать/поддерживать список родителей.

Роль в конвейере — ПЕРЕКРЁСТНАЯ ПРОВЕРКА, не замена. Проверено 09.08.2026
(разбор внешнего документа BecqMoni/tools/nucdb): полнота `tot_conv_coeff` в
`gammas` НИЖЕ, чем у `decay_rads` для реально заселяемых переходов (в выборке
Th-228 — 85 из 591 строк, ≈14%) — не источник для замены Б1-фикса. Ценность —
(1) независимое подтверждение уже посчитанных CC на пересечении, (2) доступ к
схеме нуклида, которого ещё нет в `NUCLIDES` (не нужен родитель), для будущего
расширения (two-hop, новые источники).

Сверка на запуске (см. tools/check_cc_cross_source.py): переход 129,065 кэВ
(186,838→57,773, уровни Th-228) даёт tot_conv_coeff=3,74 — дословно то же
число, что уже используется в F_B через decay_rads Ac-228 (тот же физический
переход, Ac-228 распадается в Th-228). Не новое число — независимое
подтверждение существующего.

Запуск:
    python tools/fetch_adopted_gammas.py

Пишет data/adopted_gammas_cross_check.csv (детерминирован по входным
нуклидам — коммитить как обычные данные).
"""
import csv
import io
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_CSV = os.path.join(HERE, "..", "data", "adopted_gammas_cross_check.csv")

# Тот же список и тот же формат IAEA-имени, что в fetch_conversion_coeff.py —
# но здесь нуклид ключуется САМ СОБОЙ (дочерний), не как родитель распада.
NUCLIDES = {
    # Цепочка Th-232 (configs/th232.yaml)
    "Ac228": "228AC",
    "Bi212": "212BI",
    "Pb212": "212PB",
    "Ra224": "224RA",
    "Rn220": "220RN",
    "Th228": "228TH",
    "Th232": "232TH",
    "Tl208": "208TL",
    # Цепочка Ra-226 (configs/ra226.yaml)
    "Bi214": "214BI",
    "Pb214": "214PB",
    "Po214": "214PO",
    "Ra226": "226RA",
    "Rn222": "222RN",
}

API = "https://nds.iaea.org/relnsd/v1/data?fields=gammas&nuclides={nuc}"


def fetch(nuc_iaea):
    url = API.format(nuc=nuc_iaea)
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        raw = resp.read().decode("utf-8")
    return list(csv.DictReader(io.StringIO(raw)))


def main():
    out_rows = []
    extraction_date = None
    for nuc_key, nuc_iaea in NUCLIDES.items():
        rows = fetch(nuc_iaea)
        n_kept = 0
        n_with_cc = 0
        for r in rows:
            sl = (r.get("start_level_energy") or "").strip()
            el = (r.get("end_level_energy") or "").strip()
            en = (r.get("energy") or "").strip()
            if not sl or not en:
                continue
            try:
                E = float(en)
                s = float(sl)
                e = float(el) if el else 0.0
            except ValueError:
                continue
            cc_raw = (r.get("tot_conv_coeff") or "").strip()
            cc = float(cc_raw) if cc_raw else None
            unc_cc = (r.get("tce_unc") or "").strip()
            extraction_date = r.get("Extraction_date") or extraction_date
            out_rows.append({
                "nuclide": nuc_key,
                "E_keV": E,
                "level": "%s->%s" % (sl, el or "0"),
                "start_jp": (r.get("start_level_jp") or "").strip(),
                "end_jp": (r.get("end_level_jp") or "").strip(),
                "relative_intensity": (r.get("relative_intensity") or "").strip(),
                "multipolarity": (r.get("multipolarity") or "").strip(),
                "mixing_ratio": (r.get("mixing_ratio") or "").strip(),
                "tot_conv_coeff": "" if cc is None else cc,
                "tce_unc": unc_cc,
                "source": "IAEA-NDS/ENSDF via Live Chart of Nuclides "
                          "(gammas, ADOPTED LEVELS — вне распада); "
                          "eval=%s" % (r.get("ensdf_authors") or "?"),
            })
            n_kept += 1
            if cc is not None:
                n_with_cc += 1
        print("%s (%s): %d переходов, %d с tot_conv_coeff (%.0f%%)"
              % (nuc_key, nuc_iaea, n_kept, n_with_cc,
                 100.0 * n_with_cc / n_kept if n_kept else 0.0),
              file=sys.stderr)

    out_rows.sort(key=lambda r: (r["nuclide"], r["E_keV"]))
    fieldnames = ["nuclide", "E_keV", "level", "start_jp", "end_jp",
                  "relative_intensity", "multipolarity", "mixing_ratio",
                  "tot_conv_coeff", "tce_unc", "source"]
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
        f.write("# Извлечено %s из IAEA Live Chart of Nuclides "
                "(nds.iaea.org/relnsd/v1/data?fields=gammas), схема уровней "
                "ВНЕ распада (ADOPTED LEVELS ENSDF), не decay_rads.\n"
                % (extraction_date or "?"))
        f.write("# Роль — ПЕРЕКРЁСТНАЯ ПРОВЕРКА к data/conversion_coeff_sum_"
                "peak_levels.csv, не замена: полнота tot_conv_coeff ниже "
                "(см. докстринг tools/fetch_adopted_gammas.py). Сверка — "
                "tools/check_cc_cross_source.py.\n")
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print("Записано %d строк в %s" % (len(out_rows), OUT_CSV), file=sys.stderr)


if __name__ == "__main__":
    main()
