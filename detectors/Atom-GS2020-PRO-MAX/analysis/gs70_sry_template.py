# -*- coding: utf-8 -*-
"""#GS-70: шаблон Sr-90+Y-90 в равновесии 1:1 на распад Sr-90 = поканальная сумма шаблонов Sr90 и Y90 с РАВНЫМ числом распадов.
Запуск: python gs70_sry_template.py <папка с mix_Sr90_npsmoff.csv и mix_Y90_npsmoff.csv> → mix_SrY90_npsmoff.csv там же"""
import sys, os; sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_templates_gs2020 as mt
d = sys.argv[1]; a, b = (mt.read_chunk(os.path.join(d, f"mix_{n}_npsmoff.csv")) for n in ("Sr90", "Y90"))
ha, hb = dict(a["header"]), dict(b["header"]); na, nb = int(ha["n_events_processed"]), int(hb["n_events_processed"])
if na != nb or a["bins"] != b["bins"]: raise SystemExit(f"ОТКАЗ: Sr90 n={na}, Y90 n={nb} или разные бины — сумма не равновесие 1:1")
edep = [x + y for x, y in zip(a["edep"], b["edep"])]; light = [x + y for x, y in zip(a["light"], b["light"])]
hdr = [(k, v) for k, v in a["header"]] + [("gs70_sum_of", "mix_Sr90+mix_Y90 (1:1 per Sr-90 decay)")]
mt.write_template(os.path.join(d, "mix_SrY90_npsmoff.csv"), a["comments"], hdr, a["bins"], edep, light)
print(f"SrY90: n={na} распадов Sr-90 (+{nb} Y-90); сумма edep {sum(edep)} = Sr90 {sum(a['edep'])} + Y90 {sum(b['edep'])}")
