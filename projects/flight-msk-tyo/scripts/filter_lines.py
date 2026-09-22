"""Убирает из таблицы линий кандидатов от элементов, которых нет в модели (салон + донор). Использование: filter_lines.py <in_prefix> <out_prefix>."""
import sys, csv, re
ABSENT = {"Ga", "Br", "Ca", "Gd", "Fe", "Ti", "Zn", "Ni", "Ac", "Pb", "Bi", "Tl", "Ra"}   # нет в модели: элементы; цепочек Th/U в источниках нет
sup = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")
rows = list(csv.reader(open(sys.argv[1] + "_lines.csv", encoding="utf-8")))
head = rows[0]
keep = [r for r in rows[1:] if not (r[6] == "1" and re.match(r"^[0-9]*([A-Z][a-z]?)", r[1].translate(sup)) and re.match(r"^[0-9]*([A-Z][a-z]?)", r[1].translate(sup)).group(1) in ABSENT)]
csv.writer(open(sys.argv[2] + "_lines.csv", "w", encoding="utf-8", newline="\n")).writerows([head] + keep)
import shutil
shutil.copy(sys.argv[1] + "_unlisted.csv", sys.argv[2] + "_unlisted.csv")
print("строк:", len(rows) - 1, "->", len(keep))
