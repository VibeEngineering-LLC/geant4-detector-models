# -*- coding: utf-8 -*-
r"""#GS-79: английская версия страницы GS2020 → dist/en/ (вызывается из build_page_gs2020.py --lang en; ru_switch — всегда).
Перевод: src/en/*.html (translate_gs2020_en.py html), строки JS/данных/подстановок — patches/en_strings.json
(translate_gs2020_en.py missing). Гейты: свежесть перевода по sha русского исходника (src/en/MANIFEST.json),
скелет разметки ru = en, число вхождений каждой замены в JS, кириллица = 0, запрещённое слово = 0, секрет-скан.
Сгенерирован ступенью 2 (Ollama qwen3.8:27b) по scripts/specs/SPEC-build_page_gs2020_en.md."""
import io
import json
import os
import re
import subprocess
import tempfile
import hashlib
import html.parser
import collections

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = r"<WORKDIR>\GEANT4\web\gs2020-th232-page"
SRC_RU = os.path.join(PAGE, "src")
SRC_EN = os.path.join(PAGE, "src", "en")
MANIFEST = os.path.join(SRC_EN, "MANIFEST.json")
STRINGS = os.path.join(HERE, "patches", "en_strings.json")
MISSING = os.path.join(HERE, "patches", "en_missing.json")
SRC_JSON = os.path.join(PAGE, "gs2020_sources.json")
IMG_EN = os.path.join(PAGE, "src", "img", "en")
NODE_STRIP = os.path.join(HERE, "js_strip_comments.js")
NODE_LITS = os.path.join(HERE, "js_cyr_literals.js")
HTML_FILES = ("index.html", "k40-panel.html", "k40-pops.html", "berry-panel.html", "berry-pops.html")
JS_FILES = ("g1s-th232.js", "g1s-k40.js", "g1s-berry.js", "gs-sources.js")
CYR = re.compile("[А-Яа-яЁё]")
KI = re.compile(r"(?<![А-ЯЁа-яё])\u041a\u0418(?![А-ЯЁа-яё])")
TOKEN = re.compile(r"\{\{[a-z0-9_]+(?:\s+[^\s}]+)*\}\}", re.I)
PCT_EN = (re.compile(r"\d+(?:\.\d+)?\s*%"), re.compile(r"\b\d+\.\d+\s+times\b"))
TRANSLATABLE_ATTRS = ("alt", "title", "aria-label", "placeholder")
SWITCH_RE = re.compile(r'<div class="top">(\s*)<div class="title">')

DIST_RU = None


def refuse(msg):
    raise SystemExit("ОТКАЗ #GS-79: " + msg)


def read(p):
    if not os.path.exists(p):
        refuse("нет файла " + p)
    with io.open(p, encoding="utf-8", newline="") as f:
        return f.read()


def write(p, text):
    d = os.path.dirname(p)
    if d:
        os.makedirs(d, exist_ok=True)
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        f.write(text)


def sha256_file(p):
    h = hashlib.sha256()
    with io.open(p, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def lex_code(code, js):
    parts = [r"/\*[\s\S]*?\*/"]
    if js:
        parts.append(r"//[^\n]*")
    parts.append(r'"(?:\\.|[^"\\\n])*"')
    parts.append(r"'(?:\\.|[^'\\\n])*'")
    if js:
        parts.append(r"`(?:\\.|[^`\\])*`")
    pattern = re.compile("|".join(parts))

    def repl(m):
        s = m.group(0)
        if s.startswith("/*") or s.startswith("//"):
            return ""
        if s.startswith('"'):
            return '""'
        if s.startswith("'"):
            return "''"
        if s.startswith("`"):
            return "``"
        return s

    out = pattern.sub(repl, code)
    out = re.sub(r"\s+", " ", out)
    return out.strip()


def strip_comments_html(text):
    def process_js_block(m):
        content = m.group(1)
        parts = [r'"(?:\\.|[^"\\\n])*"', r"'(?:\\.|[^'\\\n])*'", r"`(?:\\.|[^`\\])*`", r"/\*[\s\S]*?\*/", r"//[^\n]*"]
        pattern = re.compile("|".join(parts))

        def repl(mm):
            s = mm.group(0)
            if s.startswith('"') or s.startswith("'") or s.startswith("`"):
                return s
            return ""

        new_content = pattern.sub(repl, content)
        lines = new_content.split("\n")
        cleaned_lines = [line.rstrip() for line in lines]
        return "<script>" + "\n".join(cleaned_lines) + "</script>"

    text = re.sub(r"<script>([\s\S]*?)</script>", process_js_block, text)

    def process_style_block(m):
        content = m.group(1)
        new_content = re.sub(r"/\*[\s\S]*?\*/", "", content)
        return "<style>" + new_content + "</style>"

    text = re.sub(r"<style>([\s\S]*?)</style>", process_style_block, text)

    def keep_comment(m):
        body = m.group(1)
        if body.startswith("@"):
            return m.group(0)
        return ""

    text = re.sub(r"<!--([\s\S]*?)-->", keep_comment, text)
    return text


class _Skel(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.items = []
        self.cur = None

    def handle_starttag(self, tag, attrs):
        self.items.append(("S", tag, tuple(sorted((k, v or "") for k, v in attrs if k not in TRANSLATABLE_ATTRS))))
        self.cur = tag

    def handle_startendtag(self, tag, attrs):
        self.items.append(("S", tag, tuple(sorted((k, v or "") for k, v in attrs if k not in TRANSLATABLE_ATTRS))))
        self.cur = tag

    def handle_endtag(self, tag):
        self.items.append(("E", tag))
        self.cur = None

    def handle_comment(self, data):
        if data.startswith("@"):
            self.items.append(("C", data))

    def handle_data(self, data):
        if self.cur in ("script", "style"):
            norm = lex_code(data, self.cur == "script")
            if norm:
                self.items.append(("D", self.cur, norm))


def skeleton(text):
    p = _Skel()
    p.feed(text)
    p.close()
    return p.items, sorted(TOKEN.findall(text))


def skeleton_diff(ru_text, en_text):
    ru_items, ru_tokens = skeleton(ru_text)
    en_items, en_tokens = skeleton(en_text)
    if ru_items == en_items and ru_tokens == en_tokens:
        return None
    msgs = []
    if ru_items != en_items:
        idx = 0
        while idx < min(len(ru_items), len(en_items)) and ru_items[idx] == en_items[idx]:
            idx += 1
        if idx < min(len(ru_items), len(en_items)):
            msgs.append("item[%d]: %r vs %r" % (idx, str(ru_items[idx])[:200], str(en_items[idx])[:200]))
        else:
            msgs.append("len: %d vs %d" % (len(ru_items), len(en_items)))
    if ru_tokens != en_tokens:
        c_ru = collections.Counter(ru_tokens)
        c_en = collections.Counter(en_tokens)
        diff = c_ru - c_en
        if diff:
            msgs.append("tokens missing in en: %s" % str(dict(diff)))
        diff2 = c_en - c_ru
        if diff2:
            msgs.append("tokens extra in en: %s" % str(dict(diff2)))
    return "; ".join(msgs)


def load_strings():
    S = json.loads(read(STRINGS))
    for k in ("units", "html", "js_code", "js", "data", "fill"):
        if k not in S:
            S[k] = []
    if "fill_key" not in S:
        S["fill_key"] = {}
    return S


NUM = re.compile(r"\d+(?:\.\d+)?")
MONTHS = ("January February March April May June July August September October November December").split()
DATE_RE = re.compile(r"(?<![\d.])(\d{1,2})(?:([–-])(\d{1,2}))?\.(\d{1,2})\.(\d{4})(?![\d.])")
_M = "|".join(MONTHS)
RANGE_RE = re.compile(r"(\d{1,2}) (%s) (\d{4})(?:[   ]|&nbsp;)*[—–-](?:[   ]|&nbsp;)*(\d{1,2}) (%s) (\d{4})" % (_M, _M))
NO_GROUP = {"126301"}   # номера статей и прочие идентификаторы, разделитель тысяч к ним не применяется
GROUP_RE = re.compile(r"(?<![\d.,\w-])\d{1,3}(?:(?:[   ]|&nbsp;)\d{3})+(?!\d)|(?<![\d.,\w-])\d{5,}(?!\d)")


def en_num(s):
    """#GS-79 (координатор 10.10): даты ДД.ММ.ГГГГ → «9 October 2026»; целые ≥5 знаков — группы через U+202F,
    4-значные без разделителя (ISO 80000-1); десятичный разделитель — точка (ставится раньше)."""
    s = DATE_RE.sub(lambda m: "%d%s %s %s" % (int(m.group(1)), (m.group(2) + str(int(m.group(3)))) if m.group(3) else "",
                                              MONTHS[int(m.group(4)) - 1], m.group(5)) if 1 <= int(m.group(4)) <= 12 else m.group(0), s)
    def rng(m):   # sterile-2 #31: «5 October 2026 — 9 October 2026» -> «5–9 October 2026»
        d1, m1, y1, d2, m2, y2 = m.groups()
        if y1 != y2:
            return m.group(0)
        return "%s–%s %s %s" % (d1, d2, m1, y1) if m1 == m2 else "%s %s – %s %s %s" % (d1, m1, d2, m2, y1)
    s = RANGE_RE.sub(rng, s)
    def grp(m):
        d = re.sub(r"\D", "", m.group(0))
        if d in NO_GROUP:   # sterile-2 #13: идентификаторы (номер статьи) — не количество
            return m.group(0)
        return d if len(d) < 5 else "{:,}".format(int(d)).replace(",", " ")
    return GROUP_RE.sub(grp, s)


def numbers_diff(ru, en_text, ignore=()):
    """#GS-79: числа перевода = числа исходника (десятичная запятая → точка, разделители тысяч и форма дат
    приводятся одним en_num с обеих сторон); пусто — совпали. ignore — осознанные исключения (en_strings num_ignore)."""
    norm = lambda t: [x for x in NUM.findall(en_num(t).replace(" ", "")) if x not in ignore]
    a = collections.Counter(norm(re.sub(r"(?<=\d),(?=\d)", ".", ru)))
    b = collections.Counter(norm(en_text))
    return "" if a == b else "потеряно %s, лишнее %s" % (dict(a - b), dict(b - a))


def text_only(html_text):
    """Видимый текст для сверки чисел: без <style>/<script> (CSS rgba(14,12,8) — не десятичная запятая), без тегов
    и атрибутов (их сверяет скелет); значения alt/title/aria-label переводятся и сохраняются отдельно."""
    t = re.sub(r"<(style|script)\b[\s\S]*?</\1>", " ", html_text)
    alts = " ".join(re.findall(r'(?:alt|title|aria-label)="([^"]*)"', t))
    return re.sub(r"<[^>]*>", " ", t) + " " + alts


def numbers_gate(S):
    """Пары js/data/fill/html_fix, где числа разошлись; num_ok у записи — осознанное исключение (с причиной в why)."""
    bad = []
    for sec in ("js", "data", "fill"):
        for e in S[sec]:
            d = numbers_diff(e["ru"], e["en"])
            if d and not e.get("num_ok"):
                bad.append("%s «%s»: %s" % (sec, e["ru"][:50], d))
            q = re.sub(r"<[^>]*>", "", e["en"][1:-1] if sec == "js" else e["en"])   # текст без тегов и ограничителей
            if "'" in q or (sec != "js" and '"' in q):   # подсказки идут в title='…' без экранирования — только ‘ ’ “ ”
                bad.append("%s «%s»: прямая кавычка/апостроф в тексте" % (sec, e["ru"][:50]))
    return bad


def to_en_value(s, key, S):
    for e in S["fill"]:
        if e.get("ru") == s:
            return e.get("en")
    for pat, repl in S["fill_key"].get(key, []):
        s = re.sub(pat, repl, s)
    for pat, repl in S["units"]:
        s = re.sub(pat, repl, s)
    s = re.sub(r"(?<=\d),(?=\d)", ".", s)
    return en_num(s)


def en_fill(fill, S, missing):
    out = {}
    for k, fn in fill.items():
        def wrapper(*a, _k=k, _fn=fn):
            res = _fn(*a)
            s = str(res)
            conv = to_en_value(s, _k, S)
            if CYR.search(conv):
                if not any(m.get("key") == _k and m.get("ru") == s for m in missing["fill"]):
                    missing["fill"].append({"key": _k, "ru": s})
                return "@@MISSING@@"
            return conv
        out[k] = wrapper
    return out


def node(script, *args):
    r = subprocess.run(["node", script, *args], capture_output=True, text=True, encoding="utf-8")
    if r.returncode != 0:
        refuse("node failed: " + r.stderr)
    return r


def translate_js(name, S, missing, mism):
    src = os.path.join(DIST_RU, "scripts", name)
    tmp_dir = tempfile.mkdtemp()
    tmp_out = os.path.join(tmp_dir, "out.js")
    node(NODE_STRIP, src, tmp_out)
    txt = read(tmp_out)

    entries = list(S["js_code"])
    js_entries = sorted(S["js"], key=lambda e: len(e["ru"]), reverse=True)
    entries.extend(js_entries)

    for e in entries:
        if "files" in e and name not in e["files"]:   # #GS-79: замена только для части файлов (напр. «образец» тория)
            continue
        c = txt.count(e["ru"])
        n = e["n"].get(name, 0)
        if c != n:
            mism.append({"file": name, "ru": e["ru"], "expected": n, "actual": c})
            continue
        if c > 0:
            txt = txt.replace(e["ru"], e["en"])

    tmp_in = os.path.join(tmp_dir, "in.js")
    write(tmp_in, txt)
    out_json = os.path.join(tmp_dir, "lits.json")
    node(NODE_LITS, out_json, tmp_in)
    lits = json.loads(read(out_json))   # {имя_файла: {литерал: число}} — один входной файл
    for lit, count in [(k2, v2) for fd in lits.values() for k2, v2 in fd.items()]:
        found = False
        for m in missing["js"]:
            if m["ru"] == lit:
                m["n"][name] = count
                found = True
                break
        if not found:
            missing["js"].append({"ru": lit, "n": {name: count}})

    return txt


def translate_data(obj, D, missing):
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            if isinstance(k, str) and CYR.search(k):
                refuse("Cyrillic key in data: " + k)
            out[k] = translate_data(v, D, missing)
        return out
    if isinstance(obj, list):
        return [translate_data(v, D, missing) for v in obj]
    if isinstance(obj, str):
        if obj in D:   # #GS-79: словарь — и для строк без кириллицы (напр. «P(γ+ce) 10,34 %» с десятичной запятой)
            return en_num(D[obj])
        if CYR.search(obj):
            if obj not in missing["data"]:
                missing["data"].append(obj)
            return obj
        return obj
    return obj


def data_js(path, var, D, missing, pop_fill=False):
    d = json.loads(read(path))
    if pop_fill:
        d.pop("fill", None)
    d = translate_data(d, D, missing)
    raw = json.dumps(d, ensure_ascii=False, separators=(",", ":"))
    if "</script" in raw.lower():
        refuse("script tag in data")
    return "window.%s=" % var + raw + ";\n"


def pct_guard(text, where):
    matches = []
    for pat in PCT_EN:
        matches.extend(pat.findall(text))
    if matches:
        refuse(where + ": числа расчёта набраны цифрами, нужна подстановка: " + str(matches))


def switch_html(href, label, hreflang, ws):
    return ('<div class="top" style="position:relative"><a href="' + href + '" class="btn lang-switch" hreflang="' + hreflang +
            '" lang="' + hreflang + '" style="position:absolute;top:0;right:0" onclick="location.href=this.getAttribute(\'href\')+location.hash;return false">'
            + label + '</a>' + ws + '<div class="title" style="padding-right:3.6rem">')


def inject_switch(text, href, label, hreflang, where):
    matches = list(SWITCH_RE.finditer(text))
    if len(matches) != 1:
        refuse(where + ": switch count " + str(len(matches)))
    m = matches[0]
    rep = switch_html(href, label, hreflang, m.group(1))
    return text[:m.start()] + rep + text[m.end():]


def ru_switch(bp):
    p = os.path.join(bp.DIST, "index.html")
    write(p, inject_switch(read(p), "en/", "EN", "en", p))
    print("переключатель языка: dist/index.html → en/")


def check_manifest():
    if not os.path.exists(MANIFEST):   # перевода ещё не было — все файлы устарели
        return list(HTML_FILES)
    m = json.loads(read(MANIFEST))
    stale = []
    for name in HTML_FILES:
        p = os.path.join(SRC_RU, name)
        if name not in m:
            stale.append(name)
            continue
        if not os.path.exists(p):
            stale.append(name)
            continue
        if sha256_file(p) != m[name]["ru_sha256"]:
            stale.append(name)
    return stale


def render(path_en, fill, bp, markers):
    txt = read(path_en)
    pct_guard(txt, path_en)
    out, used = bp.substitute(txt, fill, os.path.basename(path_en))
    if not used:
        refuse("no substitutions in " + path_en)
    for mk in markers:
        if mk in out:
            refuse("marker " + mk + " in " + path_en)
    return out, len(used)


def replace_once(t, marker, rep, where):
    c = t.count(marker)
    if c != 1:
        refuse(where + ": marker " + marker + " count " + str(c))
    return t.replace(marker, rep)


def finish(bp, k40, berry):
    global DIST_RU
    DIST_RU = bp.DIST
    OUT = os.path.join(bp.DIST, "en")
    S = load_strings()
    missing = {"js": [], "data": [], "fill": []}
    mism = []

    stale = check_manifest()
    nbad = numbers_gate(S)
    if nbad:
        refuse("числа в переводе строк разошлись с исходником: " + "; ".join(nbad))

    skel_bad = []
    for name in HTML_FILES:
        ru_t, en_t = read(os.path.join(SRC_RU, name)), read(os.path.join(SRC_EN, name))
        diff = skeleton_diff(ru_t, en_t)
        if diff is not None:
            skel_bad.append(name + ": " + diff)
        nd = numbers_diff(text_only(strip_comments_html(ru_t)), text_only(en_t), S.get("num_ignore", ()))   # числа видимого текста сохранены
        if nd:
            skel_bad.append(name + ": числа: " + nd)

    js = {name: translate_js(name, S, missing, mism) for name in JS_FILES}

    D = {e["ru"]: e["en"] for e in S["data"]}
    djs = {
        "data.js": data_js(bp.DATA_JSON, "G1S", D, missing),
        "data-k40.js": data_js(k40.DATA_JSON, "GS_K40", D, missing),
        "data-berry.js": data_js(berry.DATA_JSON, "GS_BERRY", D, missing),
        "sources.js": data_js(SRC_JSON, "GS_SRC", D, missing, pop_fill=True)
    }

    d_th = json.loads(read(bp.DATA_JSON))
    d_k = json.loads(read(k40.DATA_JSON))
    d_b = json.loads(read(berry.DATA_JSON))

    f_th = en_fill(bp.make_fill(d_th), S, missing)
    f_k = en_fill(k40.make_fill(bp, d_k), S, missing)
    f_b = en_fill(berry.make_fill(bp, d_b), S, missing)

    K40M = ("<!--@k40panel-->", "<!--@k40pops-->", "<!--@k40-->")
    BERM = ("<!--@berrypanel-->", "<!--@berrypops-->", "<!--@berry-->")

    th, n_th = render(os.path.join(SRC_EN, "index.html"), f_th, bp, ())
    kp, n_kp = render(os.path.join(SRC_EN, "k40-panel.html"), f_k, bp, K40M)
    kq, n_kq = render(os.path.join(SRC_EN, "k40-pops.html"), f_k, bp, K40M)
    bpn, n_bp = render(os.path.join(SRC_EN, "berry-panel.html"), f_b, bp, BERM)
    bq, n_bq = render(os.path.join(SRC_EN, "berry-pops.html"), f_b, bp, BERM)

    if stale or skel_bad or mism or any(missing.values()):
        miss_out = {
            "stale_html": stale,
            "skeleton": skel_bad,
            "count_mismatch": mism,
            "js": missing["js"],
            "data": missing["data"],
            "fill": missing["fill"]
        }
        write(MISSING, json.dumps(miss_out, ensure_ascii=False, indent=1))
        refuse("stale:%d skel:%d mism:%d js:%d data:%d fill:%d -> %s" % (len(stale), len(skel_bad), len(mism), len(missing["js"]), len(missing["data"]), len(missing["fill"]), MISSING))
    else:
        if os.path.exists(MISSING):
            os.remove(MISSING)

    v = int(max(os.path.getmtime(p) for p in [STRINGS, bp.DATA_JSON, k40.DATA_JSON, berry.DATA_JSON, SRC_JSON] +
                [os.path.join(SRC_EN, n) for n in HTML_FILES] +
                [os.path.join(bp.DIST, "scripts", n) for n in JS_FILES]))

    th = replace_once(th, "<!--@styles-->", '<link rel="stylesheet" href="../styles/g1s-th232.css?v=%d">' % v, "en/index.html")
    th = replace_once(th, "<!--@script-->", '<script src="data.js?v=%d"></script>\n<script src="scripts/g1s-th232.js?v=%d"></script>' % (v, v), "en/index.html")
    th = replace_once(th, "<!--@sources-->", '<script src="sources.js?v=%d"></script>\n<script src="scripts/gs-sources.js?v=%d"></script>' % (v, v), "en/index.html")
    th = replace_once(th, "<!--@k40panel-->", kp, "en/index.html")
    th = replace_once(th, "<!--@k40pops-->", kq, "en/index.html")
    th = replace_once(th, "<!--@k40-->", '<script src="data-k40.js?v=%d"></script>\n<script src="scripts/g1s-k40.js?v=%d"></script>' % (v, v), "en/index.html")
    th = replace_once(th, "<!--@berrypanel-->", bpn, "en/index.html")
    th = replace_once(th, "<!--@berrypops-->", bq, "en/index.html")
    th = replace_once(th, "<!--@berry-->", '<script src="data-berry.js?v=%d"></script>\n<script src="scripts/g1s-berry.js?v=%d"></script>' % (v, v), "en/index.html")

    th = inject_switch(th, "../", "RU", "ru", "en/index.html")

    # #GS-82: рисунок с английской версией в src/img/en/ — копия в en/img/, остальные — общие ../img/
    en_img = set(os.listdir(IMG_EN)) if os.path.isdir(IMG_EN) else set()
    th, n_img = re.subn(r'((?:src|href)=")img/([^"]+)', lambda m: m.group(1) + ("img/" if m.group(2) in en_img else "../img/") + m.group(2), th)
    n_img_en = len(set(re.findall(r'(?:src|href)="img/([^"]+)', th)))
    for name in en_img:  # размеры <img> — по английскому файлу (подписи другой длины меняют обрезку)
        b = open(os.path.join(IMG_EN, name), "rb").read(4096)
        wh = (int.from_bytes(b[16:20], "big"), int.from_bytes(b[20:24], "big")) if name.endswith(".png") else tuple(int(x) for x in re.search(rb'<svg[^>]*?width="(\d+)" height="(\d+)"', b).groups())
        th = re.sub(r'(<img src="img/' + re.escape(name) + r'"[^>]*? width=")\d+(" height=")\d+', lambda m: m.group(1) + str(wh[0]) + m.group(2) + str(wh[1]), th)

    for e in S["html"]:
        c = th.count(e["ru"])
        if c != e["n"]:
            refuse("html replace count mismatch: " + e["ru"] + " expected " + str(e["n"]) + " got " + str(c))
        th = th.replace(e["ru"], e["en"])

    # #GS-79: формат чисел и дат видимого текста (вне <script>/<style> и вне тегов) — тем же en_num
    th = re.sub(r"(<(script|style)\b[\s\S]*?</\2>|<[^>]*>)|([^<]+)", lambda m: m.group(1) or en_num(m.group(3)), th)
    th = '<html lang="en">\n' + th

    files = {"index.html": th}
    files.update(djs)
    for name in JS_FILES:
        files["scripts/" + name] = js[name]

    failures = []
    for rel, text in files.items():
        cyrs = CYR.findall(text)
        if cyrs:
            lines = text.split("\n")
            snippets = []
            for i, line in enumerate(lines):
                if CYR.search(line):
                    snippets.append("%s:%d: %s" % (rel, i + 1, line.strip()[:50]))
                    if len(snippets) >= 5:
                        break
            failures.append("CYR in " + rel + ": " + "; ".join(snippets))
        kis = KI.findall(text)
        if kis:
            failures.append("KI in " + rel)
        sec = bp.secret_scan(text)
        if sec:
            failures.append("secret in " + rel + ": " + str(sec))
        if rel != "index.html" and "</script" in text:
            failures.append("script tag in " + rel)

    for name in sorted(en_img):
        if name.endswith(".svg") and CYR.search(io.open(os.path.join(IMG_EN, name), encoding="utf-8").read()):
            failures.append("CYR in img/en/" + name)
    if failures:
        refuse("; ".join(failures))

    for rel in files:
        print("en/%s  %d КБ  кириллица 0, запрещённое слово 0, секрет-скан чисто" % (rel, len(files[rel].encode("utf-8")) // 1024))

    written = []
    for rel, text in files.items():
        p = os.path.join(OUT, rel)
        write(p, text)
        written.append(rel)

    import shutil
    for name in sorted(en_img):
        os.makedirs(os.path.join(OUT, "img"), exist_ok=True)
        shutil.copy2(os.path.join(IMG_EN, name), os.path.join(OUT, "img", name))
        written.append("img/" + name)

    print("EN: подстановок %d (страница) + %d + %d (K-40) + %d + %d (черника); картинок перенаправлено %d, из них английских %d; v=%d" % (n_th, n_kp, n_kq, n_bp, n_bq, n_img, n_img_en, v))
    return written
