"""Извлечь уникальные строки реакций/материалов для перевода на EN. -> outputs/neu_strings.json"""
import json, re, sys
sys.path.insert(0, "scripts")
from bline_table import parse_lines

rows = parse_lines("research/B-gamma-lines-airframe.md")
neu = [r for r in rows if re.search(r"\(n[,;]\s*[γg]\)|\(n[,;]\s*n", r["reaction"]) and not re.search(r"GAGG|LaBr", r["place"])]
reactions_all = sorted(set(r["reaction"] for r in rows))   # ВСЕ строки — так же, как candidates в data.js
places_neu = sorted(set(r["place"] for r in neu))
print("reactions_all:", len(reactions_all))
print("places_neu:", len(places_neu))
json.dump({"reactions_all": reactions_all, "places_neu": places_neu}, open("outputs/neu_strings.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
