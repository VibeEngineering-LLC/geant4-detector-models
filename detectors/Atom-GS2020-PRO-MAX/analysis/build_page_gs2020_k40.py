# -*- coding: utf-8 -*-
r"""#GS-26 (вариант А): вторая панель страницы GS2020 — K-40 (KCl) — тем же донорским кодом, что панель тория.
JS: src-копия донорского g1s-th232.js (её пишет build_page_gs2020.py) + замены patches/k40_terms.json (число вхождений
каждой сверяется) → src/scripts/g1s-k40.js и dist/scripts/g1s-k40.js; сам донор не правится.
Данные: gs2020_k40_data.json (export_page_gs2020_k40.py) → dist/data-k40.js = window.GS_K40.
Разметка: src/k40-panel.html и src/k40-pops.html, метки {{...}} раскрываются подстановками build_page.py по данным K-40,
вставляются в dist/index.html и в страницу одним файлом по меткам <!--@k40panel-->, <!--@k40pops-->, <!--@k40-->.
Вызывается из build_page_gs2020.py (finish(bp)); отдельный запуск собирает и проверяет только JS и данные.
Спека: scripts\specs\SPEC-build_page_gs2020_k40.md"""

import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = r"D:\cloud-folder\work-folder\GEANT4\web\gs2020-th232-page"
DATA_JSON = os.path.join(PAGE, "gs2020_k40_data.json")
TERMS_JSON = os.path.join(HERE, "patches", "k40_terms.json")
SRC_TH_JS = os.path.join(PAGE, "src", "scripts", "g1s-th232.js")
SRC_K40_JS = os.path.join(PAGE, "src", "scripts", "g1s-k40.js")
PANEL_HTML = os.path.join(PAGE, "src", "k40-panel.html")
POPS_HTML = os.path.join(PAGE, "src", "k40-pops.html")
GATE_WORDS = ("Th-232", "паспорт", "1,405")   # #GS-26: в g1s-k40.js — 0 вхождений (после раскрытия \uXXXX)
INFO_WORDS = ("ветв", "метод 2-2", "window.G1S", "data-pop")
PCT_RE = (re.compile(r"\d+(?:,\d+)?\s*%"), re.compile(r"в\s+\d+,\d+\s+раза"))   # сторож build_page.py: числа расчёта цифрами
MARKERS = ("<!--@k40panel-->", "<!--@k40pops-->", "<!--@k40-->")


def refuse(msg):
    raise SystemExit("ОТКАЗ: " + msg)


def read(p):
    if not os.path.exists(p):
        refuse("нет файла " + p)
    return io.open(p, encoding="utf-8").read()


def unescape_js(txt):
    pat = re.compile(chr(92) * 2 + "u([0-9a-fA-F]{4})")
    return pat.sub(lambda m: chr(int(m.group(1), 16)), txt)


def word_counts(txt, words):
    dec = unescape_js(txt)
    return {w: len(re.findall(re.escape(w), dec, re.I)) for w in words}


def build_js():
    txt = read(SRC_TH_JS)
    terms = json.loads(read(TERMS_JSON))["terms"]
    for i, t in enumerate(terms):
        c = txt.count(t["old"])
        if c != t["n"]:
            refuse("замена %d (%s): «%s» встречается %d раз, ожидалось %d" % (i, t.get("why", ""), t["old"][:70], c, t["n"]))
        txt = txt.replace(t["old"], t["new"])
    g = word_counts(txt, GATE_WORDS)
    if any(v > 0 for v in g.values()):
        refuse("гейт #GS-26 в g1s-k40.js: " + str(g))
    if "</script" in txt.lower():
        refuse("в g1s-k40.js есть закрывающий тег script")
    return txt, len(terms), word_counts(txt, INFO_WORDS)


def data_js():
    raw_text = read(DATA_JSON)
    d = json.loads(raw_text)
    raw = json.dumps(d, ensure_ascii=False, separators=(",", ":"))
    if "</script" in raw.lower():
        refuse("в данных K-40 есть закрывающий тег script")
    return d, "window.GS_K40=" + raw + ";\n"


def make_fill(bp, d):
    f = bp.make_fill(d)
    fw = d["fwhm_cal"]
    m1 = d["method1"]
    m2 = d["method2"]
    meta = d["meta"]

    def _kp_fw_scale():
        return bp.rnum(fw["scale"], 2)

    def _kp_fw_k40_sl():
        return bp.rnum(fw["k40_sl_keV"], 2) + " кэВ"

    def _kp_fw_k40_conv():
        return bp.rnum(fw["k40_conv_keV"], 2) + " кэВ"

    def _kp_fw_k40_law():
        return bp.rnum(fw["k40_law_keV"], 2) + " кэВ"

    def _kp_cal_repr():
        return bp.rnum(meta["cal_sample"]["repr_max_dev_keV"], 1) + " кэВ"

    def _kp_lo():
        return bp.rkev(m1["E_fit_lo"])

    def _kp_hi():
        return bp.rkev(m1["E_fit_hi"])

    def _kp_nnodes():
        return bp.rcnt(m2["n_nodes"])

    def _kp_bg_live():
        return bp.rnum(meta["bg_live_s"], 0) + " с"

    f["kp_fw_scale"] = _kp_fw_scale
    f["kp_fw_k40_sl"] = _kp_fw_k40_sl
    f["kp_fw_k40_conv"] = _kp_fw_k40_conv
    f["kp_fw_k40_law"] = _kp_fw_k40_law
    f["kp_cal_repr"] = _kp_cal_repr
    f["kp_lo"] = _kp_lo
    f["kp_hi"] = _kp_hi
    f["kp_nnodes"] = _kp_nnodes
    f["kp_bg_live"] = _kp_bg_live

    return f


def render(bp, path, fill):
    txt = read(path)
    pcts = [m for r in PCT_RE for m in r.findall(txt)]
    if pcts:
        refuse("%s: числа расчёта набраны цифрами, нужна подстановка: %s" % (os.path.basename(path), pcts))
    out, used = bp.substitute(txt, fill, os.path.basename(path))
    if not used:
        refuse(os.path.basename(path) + ": ни одной подстановки")
    for m in MARKERS:
        if m in out:
            refuse("метка " + m + " внутри " + os.path.basename(path))
    return out, len(used)


def finish(bp):
    js, n_terms, info = build_js()
    with io.open(SRC_K40_JS, "w", encoding="utf-8", newline="") as f:
        f.write(js)

    d, djs = data_js()
    fill = make_fill(bp, d)

    panel, n_p = render(bp, PANEL_HTML, fill)
    pops, n_q = render(bp, POPS_HTML, fill)

    for name, text in (("scripts/g1s-k40.js", js), ("data-k40.js", djs), ("k40-panel.html", panel), ("k40-pops.html", pops)):
        bad = bp.secret_scan(text)
        if bad:
            refuse("секрет-скан в %s: %s" % (name, bad))

    os.makedirs(os.path.join(bp.DIST, "scripts"), exist_ok=True)
    with io.open(os.path.join(bp.DIST, "data-k40.js"), "w", encoding="utf-8") as f:
        f.write(djs)
    with io.open(os.path.join(bp.DIST, "scripts", "g1s-k40.js"), "w", encoding="utf-8") as f:
        f.write(js)

    v = int(max(os.path.getmtime(p) for p in (DATA_JSON, SRC_K40_JS, PANEL_HTML, POPS_HTML, TERMS_JSON)))

    linked = '<script src="data-k40.js?v=%d"></script>\n<script src="scripts/g1s-k40.js?v=%d"></script>' % (v, v)
    inline = "<script>" + djs + "</script>\n<script>\n" + js + "</script>"

    for path, scripts in ((os.path.join(bp.DIST, "index.html"), linked), (bp.SINGLE, inline)):
        t = read(path)
        for marker, rep in zip(MARKERS, (panel, pops, scripts)):
            if t.count(marker) != 1:
                refuse("метка %s в %s встречается %d раз" % (marker, path, t.count(marker)))
            t = t.replace(marker, rep)
        with io.open(path, "w", encoding="utf-8") as f:
            f.write(t)

    print("K-40: замен в JS %d, g1s-k40.js %d КБ, data-k40.js %d КБ, v=%d" % (n_terms, len(js.encode()) // 1024, len(djs.encode()) // 1024, v))
    print("K-40: подстановок в разметке %d (панель) + %d (окна)" % (n_p, n_q))
    print("K-40: гейт %s = 0; справочно (комментарии/недостижимые ветви): %s" % ("/".join(GATE_WORDS), info))

    return ["data-k40.js", os.path.join("scripts", "g1s-k40.js")]


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    js, n_terms, info = build_js()
    d, djs = data_js()
    with io.open(SRC_K40_JS, "w", encoding="utf-8", newline="") as f:
        f.write(js)
    print("K-40 (отдельный запуск): замен %d, g1s-k40.js %d КБ, data-k40.js %d КБ; справочно %s" % (n_terms, len(js.encode()) // 1024, len(djs.encode()) // 1024, info))
    print("разметку и dist собирает build_page_gs2020.py")
