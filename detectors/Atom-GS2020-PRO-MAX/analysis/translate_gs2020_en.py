# -*- coding: utf-8 -*-
r"""#GS-79: перевод английской страницы GS2020 локальной моделью (qwen3.8:27b через guarded_generate).
  html [имена]            — src/*.html → src/en/*.html кусками ≤4500 знаков, кэш по sha куска (src/en/.chunks.json),
                            скелет разметки ru=en и 0 кириллицы на каждый кусок; правки после вычитки — en_strings.json "html_fix".
  missing [--accept-counts] — недостающие строки JS/данных/подстановок из patches/en_missing.json (его пишет
                            build_page_gs2020.py --lang en) → patches/en_strings.json.
Словарь терминов в промпте — patches/en_glossary.md. Сгенерирован ступенью 2 по scripts/specs/SPEC-translate_gs2020_en.md."""
import argparse
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.environ.get("GS_VRAM_GUARD_DIR", os.path.join(os.path.expanduser("~"), ".claude", "skills", "workflow", "scripts")))   # папка с vram_guard_reference.py
import build_page_gs2020_en as en
from vram_guard_reference import guarded_generate, wrap_untrusted

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

MODEL = os.environ.get("GS79_MODEL", "qwen3.8:27b")
GLOSSARY = os.path.join(HERE, "patches", "en_glossary.md")
CACHE = os.path.join(en.SRC_EN, ".chunks.json")
PROMPT_VERSION = "gs79-v3"   # v3 (10.10): словарь дополнен вторым раундом терминологии #GS-78 — кэш v1/v2 недействителен
CHUNK_MAX = int(os.environ.get("GS79_CHUNK", "4500"))   # GS79_CHUNK — меньше для файла, где кусок не проходит гейты
CHUNK_FILE = {"berry-pops.html": 2200}   # файл, где кусок 4500 дважды не прошёл скелет (10.10); размер — часть ключа кэша


def set_chunk(name):
    """Размер куска для файла: GS79_CHUNK из окружения, иначе CHUNK_FILE, иначе 4500 (то же у en_fix_capture)."""
    global CHUNK_MAX
    CHUNK_MAX = int(os.environ.get("GS79_CHUNK") or CHUNK_FILE.get(name, 4500))
BATCH_MAX = 5000

RULES = """You translate parts of a scientific web page on gamma-ray spectrometry (NaI(Tl) scintillation spectrometer,
Geant4 Monte Carlo modelling, spectrum fitting) from Russian into English.
Register: scientific English as in Nuclear Instruments and Methods A or Applied Radiation and Isotopes;
impersonal constructions; no colloquialisms; no word-by-word calques of Russian syntax; keep the meaning exactly,
add nothing and omit nothing.
Terminology: use the English terms of the glossary below exactly; never translate a glossary term word by word.
Numbers: keep every value; decimal comma becomes decimal point (0,05 -> 0.05; 1,405·10¹⁰ -> 1.405·10¹⁰);
keep thousands separators as they are (spaces or &nbsp;); units: кэВ keV, МэВ MeV, Бк Bq, Бк/кг Bq/kg, мм mm,
см cm, г g, л L, мл mL, ч h, с s, г/см³ g/cm³, лет years.
Keep unchanged: issue references such as #GS-50, nuclide names (Th-232, K-40, Cs-137), code identifiers,
file names, URLs, the nickname Am6er.
The output must not contain any Cyrillic character.

GLOSSARY:
{GLOSSARY}"""

HTML_TASK = """Translate the HTML fragment given between the DATA tags. Output ONLY the translated fragment, nothing before or
after it, no markdown fences. Keep the HTML exactly: every tag, every attribute and its value, class, id, href, src,
data-* attributes, character entities (&nbsp; &minus; and so on), placeholders in double curly braces such as
{{fw_k}}, and marker comments such as <!--@k40panel--> must stay identical and in the same order. Translate only the
human-readable text, the values of the alt, title and aria-label attributes, and the Russian text inside string
literals of <script> code. Do not add, remove, merge or split elements. Keep the line structure.
The fragment may start or end in the middle of an element; translate it as it is."""

JSON_TASK = """Translate the strings given between the DATA tags. The input is JSON {"items": [{"id": ..., "kind": ..., "ru": ...}]}.
Return JSON {"items": [{"id": ..., "en": ...}]} with the same ids, one item per input item.
kind "js": the value is a raw JavaScript string literal or a template-literal part, INCLUDING its delimiters
(" or ' or ` and, for template parts, the ${ and } boundaries). Keep the delimiters, escapes (\\u00a0, \\n, \\"),
HTML tags with their attributes, and leading and trailing spaces exactly; translate only the Russian words.
kind "data": a text value shown on the page (a table note, a label, a caption). Keep HTML tags, leading and
trailing spaces.
kind "fill": a short value inserted into a sentence (a number with unit, a phrase); keep its grammatical role
so that it fits into an English sentence."""


def glossary():
    return en.read(GLOSSARY)


def ask(prompt, fmt):
    r = guarded_generate(
        MODEL,
        prompt,
        fmt=fmt,
        want_gpu=True,
        priority=50,
        temperature=0.0,
        num_ctx=32768,
        num_predict=12000,
        think=False,
        timeout_s=1800,
        project="GEANT4",
        agent="gs79-translate",
    )
    if r.get("done_reason") != "stop":
        raise RuntimeError("truncated: %r" % r.get("done_reason"))
    text = (r.get("response") or "").strip()
    if not text:
        raise RuntimeError("empty response")
    lines = text.split("\n")
    if lines and lines[0].startswith("```"):
        lines = lines[1:]
    if lines and lines[-1].strip() == "```":
        lines = lines[:-1]
    # модель иногда повторяет обёртку wrap_untrusted (<DATA_xxxxxxxxxxxx>…</DATA_…>) — снять до проверки скелета
    return re.sub(r"</?DATA_[0-9a-f]{12}>\n?", "", "\n".join(lines))


def split_chunks(text):
    lines = text.split("\n")
    chunks = []
    cur = []
    cur_len = 0
    last_boundary = -1

    def is_boundary(i):
        if i >= len(lines) - 1:
            return True
        cur_line = lines[i].rstrip("\r").strip()
        next_line = lines[i + 1].lstrip()
        return cur_line.endswith(">") and next_line.startswith("<")

    for i, line in enumerate(lines):
        line_len = len(line) + (1 if i < len(lines) - 1 else 0)
        if cur and cur_len + line_len > CHUNK_MAX:
            if is_boundary(i - 1):
                chunks.append("\n".join(cur))
                cur = []
                cur_len = 0
                last_boundary = -1
            elif last_boundary >= 0:
                split_at = last_boundary
                chunks.append("\n".join(cur[:split_at + 1]))
                cur = cur[split_at + 1:]
                cur_len = sum(len(l) + 1 for l in cur)
                last_boundary = -1
            else:
                chunks.append("\n".join(cur))
                cur = []
                cur_len = 0
                last_boundary = -1
        cur.append(line)
        cur_len += line_len
        if is_boundary(i):
            last_boundary = len(cur) - 1

    if cur:
        chunks.append("\n".join(cur))
    return chunks


def load_cache():
    if not os.path.exists(CACHE):
        return {}
    with open(CACHE, "r", encoding="utf-8") as f:
        return json.load(f)


def save_cache(c):
    with open(CACHE, "w", encoding="utf-8", newline="") as f:
        json.dump(c, f, ensure_ascii=False, indent=0)


def chunk_key(ru):
    import hashlib
    return hashlib.sha256((PROMPT_VERSION + "\n" + ru).encode("utf-8")).hexdigest()


def translate_chunk(ru, cache):
    key = chunk_key(ru)
    if key in cache:   # кэш перепроверяется текущими гейтами: не прошедший — переводится заново
        c = cache[key]
        if en.skeleton_diff(ru, c) is None and not en.CYR.search(c) and not en.numbers_diff(en.text_only(ru), en.text_only(c)):
            return c
        print("кэш отклонён гейтами, перевожу заново: %s" % ru[:60].replace("\n", " "))

    lead = ru[:len(ru) - len(ru.lstrip())]
    trail = ru[len(ru.rstrip()):]

    base_prompt = RULES.replace("{GLOSSARY}", glossary()) + "\n" + HTML_TASK + "\n" + wrap_untrusted(ru, "DATA")
    out = ask(base_prompt, None)
    out = lead + out.strip() + trail

    reason = None
    d = en.skeleton_diff(ru, out)
    if d is not None:
        reason = "skeleton_diff: %r" % d
    elif en.CYR.search(out):
        reason = "cyrillic found"
    elif en.numbers_diff(en.text_only(ru), en.text_only(out)):
        reason = "numbers changed: " + en.numbers_diff(en.text_only(ru), en.text_only(out))

    if reason:
        retry_prompt = base_prompt + "\n\nYOUR PREVIOUS OUTPUT WAS REJECTED: " + reason + "\nTranslate again, keeping the HTML exactly."
        out = ask(retry_prompt, None)
        out = lead + out.strip() + trail
        d = en.skeleton_diff(ru, out)
        if d is not None:
            raise RuntimeError("skeleton_diff: %r | chunk: %s" % (d, ru[:300]))
        if en.CYR.search(out):
            raise RuntimeError("cyrillic found | chunk: %s" % ru[:300])
        if en.numbers_diff(en.text_only(ru), en.text_only(out)):
            raise RuntimeError("numbers changed: %s | chunk: %s" % (en.numbers_diff(en.text_only(ru), en.text_only(out)), ru[:300]))

    cache[key] = out
    save_cache(cache)
    return out


def apply_fixes(name, text, S):
    applied = 0
    for e in S.get("html_fix", []):
        if e["file"] != name:
            continue
        c = text.count(e["ru"])
        if c != e["n"]:
            print("WARNING: %s: %r expected %d got %d" % (name, e["ru"][:60], e["n"], c))
            continue
        text = text.replace(e["ru"], e["en"])
        applied += 1
    return text, applied


def cmd_html(names):
    if not names:
        names = en.HTML_FILES
    cache = load_cache()
    manifest = {}
    if os.path.exists(en.MANIFEST):
        with open(en.MANIFEST, "r", encoding="utf-8") as f:
            manifest = json.load(f)

    for name in names:
        ru = en.strip_comments_html(en.read(os.path.join(en.SRC_RU, name)))
        set_chunk(name)
        chunks = split_chunks(ru)
        n = len(chunks)
        translated = []
        for i, chunk in enumerate(chunks):
            key = chunk_key(chunk)
            status = "cached" if key in cache else "model"
            out = translate_chunk(chunk, cache)
            print("%s chunk %d/%d %s %d" % (name, i + 1, n, status, len(chunk)))
            translated.append(out)
        out = "\n".join(translated)
        d = en.skeleton_diff(ru, out)
        if d is not None:
            raise RuntimeError("skeleton_diff after join: %r" % d)
        S = en.load_strings()
        out, k = apply_fixes(name, out, S)
        d = en.skeleton_diff(ru, out)
        if d is not None:
            raise RuntimeError("skeleton_diff after fixes: %r" % d)
        en.write(os.path.join(en.SRC_EN, name), out)
        manifest[name] = {
            "ru_sha256": en.sha256_file(os.path.join(en.SRC_RU, name)),
            "model": MODEL,
            "prompt": PROMPT_VERSION,
            "date": time.strftime("%Y-%m-%d %H:%M"),
            "chunks": n,
            "fixes": k,
        }
        with open(en.MANIFEST, "w", encoding="utf-8", newline="") as f:   # после каждого файла: сбой следующего не теряет готовые
            json.dump(manifest, f, ensure_ascii=False, indent=1)


def validate_item(kind, ru, en_text):
    if en.CYR.search(en_text):
        return "cyrillic found"
    nd = en.numbers_diff(ru, en_text)   # #GS-79: числа сохраняются (запятая → точка)
    if nd:
        return "numbers changed: " + nd
    if kind == "js":
        if en_text[0] != ru[0]:
            return "first char mismatch"
        if ru.endswith("${"):
            if not en_text.endswith("${"):
                return "template end mismatch"
        else:
            if en_text[-1] != ru[-1]:
                return "last char mismatch"
        inner_ru = ru[1:-1] if not ru.endswith("${") else ru[1:-2]
        inner_en = en_text[1:-1] if not en_text.endswith("${") else en_text[1:-2]
        lead_ru = inner_ru[:len(inner_ru) - len(inner_ru.lstrip())]
        lead_en = inner_en[:len(inner_en) - len(inner_en.lstrip())]
        if lead_ru != lead_en:
            return "leading whitespace mismatch"
        trail_ru = inner_ru[len(inner_ru.rstrip()):]
        trail_en = inner_en[len(inner_en.rstrip()):]
        if trail_ru != trail_en:
            return "trailing whitespace mismatch"
        if ru.count("${") != en_text.count("${"):
            return "template count mismatch"
    d = en.skeleton_diff(ru, en_text)
    if d is not None:
        return "skeleton_diff: %r" % d
    if kind in ("data", "fill"):
        lead_ru = ru[:len(ru) - len(ru.lstrip())]
        lead_en = en_text[:len(en_text) - len(en_text.lstrip())]
        if lead_ru != lead_en:
            return "leading whitespace mismatch"
        trail_ru = ru[len(ru.rstrip()):]
        trail_en = en_text[len(en_text.rstrip()):]
        if trail_ru != trail_en:
            return "trailing whitespace mismatch"
    return None


def translate_items(items):
    batches = []
    cur = []
    cur_len = 0
    for item in items:
        l = len(item["ru"])
        if cur and cur_len + l > BATCH_MAX:
            batches.append(cur)
            cur = []
            cur_len = 0
        cur.append(item)
        cur_len += l
    if cur:
        batches.append(cur)

    ok = {}
    failed = []
    for batch in batches:
        def mk(b):   # в модель — только id/kind/ru (без служебного n)
            return RULES.replace("{GLOSSARY}", glossary()) + "\n" + JSON_TASK + "\n" + wrap_untrusted(
                json.dumps({"items": [{"id": x["id"], "kind": x["kind"], "ru": x["ru"]} for x in b]}, ensure_ascii=False), "DATA")
        prompt = mk(batch)
        try:
            resp = ask(prompt, "json")
            data = json.loads(resp)
            resp_map = {item["id"]: item["en"] for item in data["items"]}
        except Exception as e:
            resp_map = {}
            batch_error = str(e)

        for item in batch:
            en_text = resp_map.get(item["id"])
            if en_text is None:
                reason = "missing in response"
            else:
                reason = validate_item(item["kind"], item["ru"], en_text)
            if reason:
                retry_prompt = mk([item]) + "\n\nYOUR PREVIOUS OUTPUT WAS REJECTED: " + reason + "\nTranslate again, keeping the HTML exactly."
                try:
                    resp2 = ask(retry_prompt, "json")
                    data2 = json.loads(resp2)
                    resp_map2 = {i["id"]: i["en"] for i in data2["items"]}
                    en_text2 = resp_map2.get(item["id"])
                    if en_text2 is not None:
                        reason2 = validate_item(item["kind"], item["ru"], en_text2)
                        if reason2 is None:
                            ok[item["id"]] = en_text2
                            continue
                        reason = reason2
                except Exception as e2:
                    reason = "retry failed: %s" % e2
                failed.append({"id": item["id"], "ru": item["ru"], "reason": reason})
            else:
                ok[item["id"]] = en_text
    return ok, failed


def cmd_missing(accept_counts):
    miss = json.loads(en.read(en.MISSING))
    S = en.load_strings()

    if miss.get("stale_html") or miss.get("skeleton"):
        if miss.get("stale_html"):
            print("stale_html: %s" % miss["stale_html"])
        if miss.get("skeleton"):
            print("skeleton: %s" % miss["skeleton"])
        print("run: translate_gs2020_en.py html <names>")

    if miss.get("count_mismatch"):
        if accept_counts:
            for m in miss["count_mismatch"]:
                ru_val = m["ru"]
                found = False
                for key in ("js", "js_code"):
                    if key in S:
                        for entry in S[key]:
                            if entry.get("ru") == ru_val:   # n — словарь {файл: число}, правится только свой файл
                                old_n = entry["n"].get(m["file"], 0)
                                if m["actual"]:
                                    entry["n"][m["file"]] = m["actual"]
                                else:
                                    entry["n"].pop(m["file"], None)
                                print("count updated: %s %s %d -> %d" % (m["file"], ru_val[:40], old_n, m["actual"]))
                                found = True
                                break
                    if found:
                        break
                if not found:
                    print("WARNING: count mismatch entry not found: %s" % ru_val[:40])
            for key in ("js", "js_code"):
                if key in S:
                    to_remove = [e for e in S[key] if not e.get("n")]
                    for e in to_remove:
                        S[key].remove(e)
                        print("removed: %s" % e["ru"][:40])
        else:
            for m in miss["count_mismatch"]:
                print("count mismatch: %s expected %d got %d" % (m["ru"][:40], m["expected"], m["actual"]))
            print("--accept-counts")

    items = []
    for i, entry in enumerate(miss.get("js", [])):
        items.append({"id": "js:%d" % i, "kind": "js", "ru": entry["ru"], "n": entry.get("n", 0)})
    for i, entry in enumerate(miss.get("data", [])):   # en_missing.json: data — список строк, не словарей
        items.append({"id": "data:%d" % i, "kind": "data", "ru": entry})
    for i, entry in enumerate(miss.get("fill", [])):
        items.append({"id": "fill:%d" % i, "kind": "fill", "ru": entry["ru"]})

    ok, failed = translate_items(items)

    existing_js = {e["ru"] for e in S.get("js", [])}
    existing_data = {e["ru"] for e in S.get("data", [])}
    existing_fill = {e["ru"] for e in S.get("fill", [])}

    added_js = 0
    added_data = 0
    added_fill = 0

    for item in items:
        en_text = ok.get(item["id"])
        if en_text is None:
            continue
        ru_val = item["ru"]
        if item["kind"] == "js":
            if ru_val not in existing_js:
                S.setdefault("js", []).append({"ru": ru_val, "en": en_text, "n": item.get("n", 0)})
                existing_js.add(ru_val)
                added_js += 1
        elif item["kind"] == "data":
            if ru_val not in existing_data:
                S.setdefault("data", []).append({"ru": ru_val, "en": en_text})
                existing_data.add(ru_val)
                added_data += 1
        elif item["kind"] == "fill":
            if ru_val not in existing_fill:
                S.setdefault("fill", []).append({"ru": ru_val, "en": en_text})
                existing_fill.add(ru_val)
                added_fill += 1

    with open(en.STRINGS, "w", encoding="utf-8", newline="") as f:
        f.write(json.dumps(S, ensure_ascii=False, indent=1) + "\n")

    if failed:
        for f_item in failed:
            print('ollama_failure: {"model": "%s", "error": "%s", "file": "en_strings.json"}' % (MODEL, f_item["reason"]))
        sys.exit(5)
    else:
        print("added: js=%d data=%d fill=%d" % (added_js, added_data, added_fill))
        sys.exit(0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cmd", choices=["html", "missing"])
    parser.add_argument("names", nargs="*")
    parser.add_argument("--accept-counts", action="store_true")
    args = parser.parse_args()

    try:
        if args.cmd == "html":
            cmd_html(args.names)
        elif args.cmd == "missing":
            cmd_missing(args.accept_counts)
    except RuntimeError as e:
        print('ollama_failure: {"model": "%s", "error": "%s", "tier": "translate", "file": "%s"}' % (MODEL, str(e), args.names[0] if args.names else "unknown"), file=sys.stderr)
        sys.exit(4)


if __name__ == "__main__":
    main()
