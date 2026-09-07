# -*- coding: utf-8 -*-
"""Режет многофайловый вывод генератора (_gen_all.txt) по маркерам
`=== FILE: <имя> ===` на отдельные файлы рядом. Fence-строки ``` отбрасывает."""
import re, sys, pathlib
sys.stdout.reconfigure(encoding="utf-8")
src = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "_gen_all.txt")
text = src.read_text(encoding="utf-8")
parts = re.split(r"^=== FILE: (\S+) ===\s*$", text, flags=re.M)
if len(parts) < 3:
    print("маркеров не найдено"); sys.exit(1)
for name, body in zip(parts[1::2], parts[2::2]):
    lines = [l for l in body.splitlines() if not l.strip().startswith("```")]
    out = src.parent / name
    out.write_text("\n".join(lines).strip() + "\n", encoding="utf-8")
    print(f"{name:42s} {len(lines):5d} строк")
