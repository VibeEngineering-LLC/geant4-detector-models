"""Извлечь ВСЕ строки таблицы раздела 9 (5 колонок) -> outputs/table9_rows.json,
для перевода колонок 2 (Реакция) и 4 (Материал) на EN. Колонки 1/3/5 (числа, буквы) не трогать."""
import json

lines = open("research/B-gamma-lines-airframe.md", encoding="utf-8").readlines()
start = None
for i, line in enumerate(lines):
    if line.startswith("## 9. БОЛЬШАЯ СВОДНАЯ ТАБЛИЦА ЛИНИЙ"):
        start = i + 1
        break
rows = []
for i in range(start, len(lines)):
    line = lines[i]
    if line.startswith("## "):
        break
    if line.startswith("|") and not line.startswith("| E, кэВ |") and not line.startswith("|---|"):
        cells = [c.strip() for c in line.strip().split("|")[1:-1]]
        if len(cells) == 5:
            rows.append(cells)
print("rows:", len(rows))
json.dump(rows, open("outputs/table9_rows.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
