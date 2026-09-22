"""Проверить, есть ли кириллица в колонках 1(E)/3(Выход)/5(Заметность) таблицы 9."""
import json, re
rows = json.load(open("outputs/table9_rows.json", encoding="utf-8"))
cyr = re.compile(r"[А-Яа-яЁё]")
vals1 = set(r[0] for r in rows if cyr.search(r[0]))
vals3 = set(r[2] for r in rows if cyr.search(r[2]))
vals5 = set(r[4] for r in rows if cyr.search(r[4]))
print("col1 cyr:", len(vals1)); [print(" ", v) for v in sorted(vals1)]
print("col3 cyr:", len(vals3)); [print(" ", v) for v in sorted(vals3)]
print("col5 cyr:", len(vals5)); [print(" ", v) for v in sorted(vals5)]
