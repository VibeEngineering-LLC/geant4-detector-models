#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Выгружает коэффициенты полной внутренней конверсии (conversion_coeff) для
γ-переходов конкретных нуклидов из IAEA Live Chart of Nuclides (REST
`decay_rads`, rad_types=g, ПОЛНЫЙ набор полей — не тот подмножественный
экспорт, что лежит в data/ensdf_*_chain_lines.csv).

Зачем отдельный файл. Аудит Б1 (внешний аудитор, 09.08.2026, коммит df5d178):
F_B в _sum_peaks_with_fb() (export_data.py) считался как сумма ТОЛЬКО
гамма-интенсивностей исходящих из уровня переходов — это депопуляция уровня
БЕЗ канала внутренней конверсии, а не истинное заселение уровня. Для
низкоэнергетичных переходов в тяжёлых ядрах (пример: 129,065 кэВ Ac-228,
уровень 186,827 кэВ) конверсия доминирует над гамма-каналом в разы —
численно подтверждено: входящих в уровень (сумма I_percent каскадных гамма
+ прямое бета-питание) 11,01-11,47 %, исходящих гамма-квантов всего 2,42 %
(разница ×4,5-4,7, две независимые проверки — через баланс
population=depopulation по прямому бета-питанию IAEA decay_rads rad_types=bm,
и через I_gamma×(1+CC) по этому файлу — сошлись в пределах 4 %).

Правильная величина F_B = Σ I_gamma·(1+CC) по ВСЕМ исходящим гамма-переходам
уровня (полная депопуляция = гамма + конверсионные электроны). Раньше это
поле (`conversion_coeff`) не извлекалось при построении data/ensdf_*_chain_
lines.csv — не потому, что источник его не даёт (API отдаёт колонку
`conversion_coeff` в том же запросе `rad_types=g`), а потому, что экспорт CSV
брал подмножество колонок и её не сохранил. Здесь колонка добывается заново
из ТОГО ЖЕ источника, тем же запросом.

Область действия: изначально (09.08.2026, аудит Б1) список нуклидов был
закрытый и узкий — только те, что фигурируют в SUM_PEAKS (Ac228, Tl208,
Bi214). РАСШИРЕН тем же вечером (аудит Б4, задача #2) на ВСЕ нуклиды
обеих цепочек: tools/enumerate_sum_peaks.py (перебор НОВЫХ кандидатов,
не только уже принятых пар) должен уметь считать правильный F_B для
любого нуклида библиотеки, не только уже отобранных вручную.

Запуск:
    python tools/fetch_conversion_coeff.py

Пишет data/conversion_coeff_sum_peak_levels.csv (append/overwrite полностью,
файл детерминирован по входным нуклидам — коммитить как обычные данные).
"""
import csv
import io
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_CSV = os.path.join(HERE, "..", "data", "conversion_coeff_sum_peak_levels.csv")

# nuclide-ключ пайплайна (как в data/ensdf_*_chain_lines.csv) -> формат IAEA Live Chart
# Список РАСШИРЕН 09.08.2026 (аудит Б4, задача #2): изначально только три
# нуклида, реально фигурирующих в SUM_PEAKS (Ac228, Tl208, Bi214) -- теперь
# ВСЕ нуклиды обеих цепочек (Th-232 и Ra-226), чтобы tools/enumerate_sum_
# peaks.py мог считать F_B по правильной формуле (гамма+конверсия) для
# ЛЮБОГО кандидата, а не только для уже принятых в конфиг пар.
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

API = "https://nds.iaea.org/relnsd/v1/data?fields=decay_rads&nuclides={nuc}&rad_types=g"


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
        for r in rows:
            sl = (r.get("start_level_energy") or "").strip()
            el = (r.get("end_level_energy") or "").strip()
            if not sl or not el:
                continue  # характеристический рентген / без привязки к уровням
            try:
                E = float(r["energy"])
                I = float(r["intensity"])
                s = float(sl)
                e = float(el)
            except ValueError:
                continue
            cc_raw = (r.get("conversion_coeff") or "").strip()
            cc = float(cc_raw) if cc_raw else None
            unc_cc = (r.get("unc_cc") or "").strip()
            extraction_date = r.get("Extraction_date") or extraction_date
            out_rows.append({
                "nuclide": nuc_key,
                "E_keV": E,
                "I_percent": I,
                "level": "%s->%s" % (r["start_level_energy"], r["end_level_energy"]),
                "conversion_coeff": "" if cc is None else cc,
                "unc_cc": unc_cc,
                "multipolarity": (r.get("multipolarity") or "").strip(),
                "source": "IAEA-NDS/ENSDF via Live Chart of Nuclides "
                          "(decay_rads, rad_types=g, поле conversion_coeff); "
                          "eval=%s" % (r.get("ensdf_authors") or "?"),
            })
            n_kept += 1
        print("%s (%s): %d transitions с уровнями" % (nuc_key, nuc_iaea, n_kept),
              file=sys.stderr)

    out_rows.sort(key=lambda r: (r["nuclide"], r["E_keV"]))
    fieldnames = ["nuclide", "E_keV", "I_percent", "level", "conversion_coeff",
                  "unc_cc", "multipolarity", "source"]
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
        f.write("# Извлечено %s из IAEA Live Chart of Nuclides "
                "(nds.iaea.org/relnsd/v1/data?fields=decay_rads&rad_types=g).\n"
                % (extraction_date or "?"))
        f.write("# Только строки с непустыми start_level_energy/end_level_energy "
                "(нуклидные гамма-переходы, не характеристический рентген).\n")
        f.write("# Назначение: F_B в export_data.py._sum_peaks_with_fb() и в "
                "tools/enumerate_sum_peaks.py (аудит Б1/Б4, коммит df5d178). "
                "Список нуклидов — ОБЕ цепочки целиком (Th-232 и Ra-226), см. "
                "докстринг tools/fetch_conversion_coeff.py.\n")
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print("Записано %d строк в %s" % (len(out_rows), OUT_CSV), file=sys.stderr)


if __name__ == "__main__":
    main()
