"""table9_rows.json -> outputs/translate_reactions.jsonl + outputs/translate_places.jsonl (уникальные строки)."""
import json

rows = json.load(open("outputs/table9_rows.json", encoding="utf-8"))
reactions = sorted(set(r[1] for r in rows))
places = sorted(set(r[3] for r in rows))
print("reactions:", len(reactions), "places:", len(places))
with open("outputs/translate_reactions.jsonl", "w", encoding="utf-8") as f:
    for i, r in enumerate(reactions):
        f.write(json.dumps({"id": f"r{i}", "text": r}, ensure_ascii=False) + "\n")
with open("outputs/translate_places.jsonl", "w", encoding="utf-8") as f:
    for i, p in enumerate(places):
        f.write(json.dumps({"id": f"p{i}", "text": p}, ensure_ascii=False) + "\n")
