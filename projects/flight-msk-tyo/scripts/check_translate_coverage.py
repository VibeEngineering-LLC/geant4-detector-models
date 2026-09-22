"""Проверить покрытие словаря translate_dict на всех уникальных reaction/place из table9_rows.json."""
import json, sys
sys.path.insert(0, "scripts")
from translate_dict import translate_reaction, translate_place

rows = json.load(open("outputs/table9_rows.json", encoding="utf-8"))
reactions = sorted(set(r[1] for r in rows))
places = sorted(set(r[3] for r in rows))
miss_r = [r for r in reactions if translate_reaction(r) is None]
miss_p = [p for p in places if translate_place(p) is None]
print("reactions total:", len(reactions), "missing:", len(miss_r))
for m in miss_r: print(" R:", m)
print("places total:", len(places), "missing:", len(miss_p))
for m in miss_p: print(" P:", m)
