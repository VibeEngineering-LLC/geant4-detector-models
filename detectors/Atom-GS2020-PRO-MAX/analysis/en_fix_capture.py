#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""#GS-79: ручные правки src/en/*.html → записи html_fix в patches/en_strings.json (переводчик воспроизводит их из кэша).
Запуск: python en_fix_capture.py [имена]; отказ, если правки не воспроизводят файл побайтно.
Сгенерирован ступенью 2 (qwen3.8:27b) по scripts/specs/SPEC-en_fix_capture.md; индексы слияния исправлены вручную."""
import difflib
import io
import json
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_page_gs2020_en as en
import translate_gs2020_en as tr


def assembled(name):
    ru = en.strip_comments_html(en.read(os.path.join(en.SRC_RU, name)))
    cache = tr.load_cache()
    tr.set_chunk(name)
    chunks = tr.split_chunks(ru)
    for c in chunks:
        key = tr.chunk_key(c)
        if key not in cache:
            raise SystemExit("ОТКАЗ: нет перевода куска в кэше: " + c[:80])
    return "\n".join(cache[tr.chunk_key(c)] for c in chunks)


def make_fixes(name, base, cur):
    a = base.splitlines(keepends=True)
    b = cur.splitlines(keepends=True)
    ops = [op for op in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes() if op[0] != "equal"]
    
    merged_ops = []
    i = 0
    while i < len(ops):
        prev = ops[i]
        j = i + 1
        while j < len(ops):
            next_op = ops[j]
            gap = next_op[1] - prev[2]
            if gap <= 4:
                prev = (prev[0], prev[1], next_op[2], prev[3], next_op[4])   # (tag, i1, i2, j1, j2): исправлено вручную
                j += 1
            else:
                break
        merged_ops.append(prev)
        i = j

    fixes = []
    for op in merged_ops:
        tag, i1, i2, j1, j2 = op
        k = 1
        lo = max(0, i1 - k)
        hi = min(len(a), i2 + k)
        old = "".join(a[lo:hi])
        while base.count(old) != 1 and k < 6:
            k += 1
            lo = max(0, i1 - k)
            hi = min(len(a), i2 + k)
            old = "".join(a[lo:hi])
        
        new = "".join(a[lo:i1]) + "".join(b[j1:j2]) + "".join(a[i2:hi])
        fixes.append({"file": name, "ru": old, "en": new, "n": 1})
        
    return fixes


def main():
    names = sys.argv[1:] if sys.argv[1:] else en.HTML_FILES
    S = en.load_strings()
    S.setdefault("html_fix", [])
    
    for name in names:
        base = assembled(name)
        cur = en.read(os.path.join(en.SRC_EN, name))
        fixes = make_fixes(name, base, cur)
        
        S2 = {"html_fix": fixes}
        out, k = tr.apply_fixes(name, base, S2)
        if out != cur:
            raise SystemExit("ОТКАЗ: правки %s не воспроизводятся (%d из %d применены)" % (name, k, len(fixes)))
        
        other_fixes = [f for f in S["html_fix"] if f["file"] != name]
        S["html_fix"] = other_fixes + fixes
        
        print("%s: правок %d" % (name, len(fixes)))
        
    en.write(en.STRINGS, json.dumps(S, ensure_ascii=False, indent=1) + "\n")


if __name__ == "__main__":
    main()
