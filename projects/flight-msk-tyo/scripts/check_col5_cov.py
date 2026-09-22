import json, sys
sys.path.insert(0, "scripts")
from translate_dict import translate_place
rows = json.load(open("outputs/table9_rows.json", encoding="utf-8"))
vals5 = sorted(set(r[4] for r in rows))
miss = [v for v in vals5 if translate_place(v) is None]
print("col5 total:", len(vals5), "missing:", len(miss))
for m in miss: print(" ", m)
