"""table9_rows.json (RU) -> research/B-gamma-lines-airframe.en.md (EN), словарный перевод."""
import json, re, sys
sys.path.insert(0, "scripts")
from translate_dict import translate_reaction, translate_yield, translate_place

rows = json.load(open("outputs/table9_rows.json", encoding="utf-8"))
out = ["## 9. БОЛЬШАЯ СВОДНАЯ ТАБЛИЦА ЛИНИЙ", "",
       "| E, keV | Reaction / nuclide | Yield | Where it occurs in the aircraft | Visibility |", "|---|---|---|---|---|"]
fail = []
for r in rows:
    e = re.sub(r"(?<=\d),(?=\d)", ".", r[0])
    reaction = translate_reaction(r[1])
    yld = translate_yield(r[2])
    place = translate_place(r[3])
    vis = translate_place(r[4])
    if None in (reaction, yld, place, vis):
        fail.append(r)
        continue
    out.append(f"| {e} | {reaction} | {yld} | {place} | {vis} |")
open("research/B-gamma-lines-airframe.en.md", "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")
print("rows:", len(rows), "written:", len(out) - 4, "failed:", len(fail))
for f in fail:
    print(" FAIL:", f)
