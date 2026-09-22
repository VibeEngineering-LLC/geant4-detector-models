"""Разбить строки материалов на токены (для словарного перевода) -> печать уникальных токенов."""
import json, re
places = [json.loads(l)["text"] for l in open("outputs/translate_places.jsonl", encoding="utf-8")]
tokens = set()
for p in places:
    p2 = p.replace("**", "")
    for part in re.split(r",\s*|\s*/\s*|\s*—\s*", p2):
        part = part.strip()
        if part:
            tokens.add(part)
print(len(tokens))
for t in sorted(tokens):
    print(repr(t))
