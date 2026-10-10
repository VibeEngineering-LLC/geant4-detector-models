# -*- coding: utf-8 -*-
r"""#GS-72: третья панель страницы GS2020 — «Черника» (Cs-137 + K-40, β Sr-90+Y-90) — тем же донорским кодом, что панели тория и K-40.
JS: src/scripts/g1s-k40.js (его пишет build_page_gs2020_k40.finish — модуль K-40 обязан отработать РАНЬШЕ) + замены patches/berry_terms.json (число вхождений каждой
сверяется; запись с new_file вставляет блок patches/berry_overrides.js) → src/scripts/g1s-berry.js и dist/scripts/g1s-berry.js. Данные: gs2020_berry_data.json
(export_page_gs2020_berry.py) → dist/data-berry.js = window.GS_BERRY. Разметка: src/berry-panel.html и src/berry-pops.html, метки раскрываются подстановками build_page.py,
вставляются в dist/index.html и в страницу одним файлом по меткам <!--@berrypanel-->, <!--@berrypops-->, <!--@berry-->. Спека: scripts\specs\SPEC-build_page_gs2020_berry.md"""

import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = r"<WORKDIR>\GEANT4\web\gs2020-th232-page"
DATA_JSON = os.path.join(PAGE, "gs2020_berry_data.json")
TERMS_JSON = os.path.join(HERE, "patches", "berry_terms.json")
PATCH_DIR = os.path.join(HERE, "patches")
SRC_BASE_JS = os.path.join(PAGE, "src", "scripts", "g1s-k40.js")      # база: готовый JS панели K-40
SRC_BERRY_JS = os.path.join(PAGE, "src", "scripts", "g1s-berry.js")
PANEL_HTML = os.path.join(PAGE, "src", "berry-panel.html")
POPS_HTML = os.path.join(PAGE, "src", "berry-pops.html")
GATE_WORDS = ("Th-232", "паспорт", "1,405")
INFO_WORDS = ("ветв", "метод 2-2", "window.G1S", "data-pop", "window.GS_K40", "k40-")
PCT_RE = (re.compile(r"\d+(?:,\d+)?\s*%"), re.compile(r"в\s+\d+,\d+\s+раза"))
MARKERS = ("<!--@berrypanel-->", "<!--@berrypops-->", "<!--@berry-->")


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
    txt = read(SRC_BASE_JS)
    terms = json.loads(read(TERMS_JSON))["terms"]
    import ru_rules as RU   # #GS-78: база g1s-k40.js уже прошла правила вычитки
    for t in terms:
        t["old"] = RU.js(t["old"])
        if "new" in t:
            t["new"] = RU.js(t["new"])
    for i, t in enumerate(terms):
        c = txt.count(t["old"])
        if c != t["n"]:
            refuse("замена %d (%s): «%s» встречается %d раз, ожидалось %d" % (i, t.get("why", ""), t["old"][:70], c, t["n"]))
        if "new_file" in t:
            new_text = read(os.path.join(PATCH_DIR, t["new_file"])).rstrip("\n")
            if t.get("append_old"):
                new_text = new_text + "\n\n" + t["old"]
        else:
            new_text = t["new"]
        txt = txt.replace(t["old"], new_text)
    g = word_counts(txt, GATE_WORDS)
    if any(v > 0 for v in g.values()):
        refuse("гейт #GS-72 в g1s-berry.js: " + str(g))
    if "</script" in txt.lower():
        refuse("в g1s-berry.js есть закрывающий тег script")
    if "GS_BERRY" not in txt:
        refuse("в g1s-berry.js нет GS_BERRY")
    if 'PFX = "k40-"' in txt:
        refuse('в g1s-berry.js осталась строка PFX = "k40-"')
    return txt, len(terms), word_counts(txt, INFO_WORDS)


def data_js():
    raw_text = read(DATA_JSON)
    d = json.loads(raw_text)
    raw = json.dumps(d, ensure_ascii=False, separators=(",", ":"))
    if "</script" in raw.lower():
        refuse("в данных черники есть закрывающий тег script")
    return d, "window.GS_BERRY=" + raw + ";\n"


def make_fill(bp, d):
    f = bp.make_fill(d)
    fw = d["fwhm_cal"]
    m1 = d["method1"]
    meta = d["meta"]
    fa = d["facts"]
    sl = d["beta1s"]
    pn = m1["per_nuclide"]
    cs = pn["CS137"]
    kk = pn["K40"]

    def sgn(r):
        return ("+" if r >= 1 else "−") + bp.rnum(abs(100 * (r - 1)), 1) + "\u00a0%"

    def pick(group, kind):
        return next(item for item in d["cmp"][group] if item["kind"] == kind)

    f["kp_fw_scale"] = lambda: bp.rnum(fw["scale"], 2)
    f["kp_fw_k40_sl"] = lambda: bp.rnum(fw["k40_sl_keV"], 2) + " кэВ"
    f["kp_fw_k40_conv"] = lambda: bp.rnum(fw["k40_conv_keV"], 2) + " кэВ"
    f["kp_fw_k40_law"] = lambda: bp.rnum(fw["k40_law_keV"], 2) + " кэВ"
    f["kp_cal_repr"] = lambda: bp.rnum(meta["cal_sample"]["repr_max_dev_keV"], 3) + " кэВ"
    f["kp_lo"] = lambda: bp.rkev(m1["E_fit_lo"])
    f["kp_hi"] = lambda: bp.rkev(m1["E_fit_hi"])
    f["kp_bg_live"] = lambda: bp.rnum(meta["bg_live_s"], 0) + " с"

    f["berry_mass"] = lambda: bp.rnum(d["passport"]["mass_g"], 0) + " г"
    f["berry_rho"] = lambda: bp.rnum(meta["matrix_density_g_cm3"], 3) + " г/см³"
    f["berry_live_h"] = lambda: bp.rnum(meta["live_s"] / 3600.0, 1) + " ч"

    def _berry_dates():
        ds = meta["date_start"]
        de = meta["date_end"]
        fmt = lambda s: s[8:10] + "." + s[5:7] + "." + s[0:4]   # #GS-78: единый формат ДД.ММ.ГГГГ
        return fmt(ds) + "–" + fmt(de)

    f["berry_dates"] = _berry_dates

    f["berry_level"] = lambda: "−" + bp.rnum(abs(fa["level_mm"]), 1) + " мм"
    f["berry_cs_a"] = lambda: bp.rcnt(cs["A_Bq"]) + " Бк"
    f["berry_cs_kg"] = lambda: bp.rcnt(cs["per_kg"]) + " Бк/кг"
    f["berry_cs_dstat"] = lambda: "± " + bp.rnum(m1["dA_stat_Bq"], 1) + " Бк"

    f["berry_k_a"] = lambda: bp.rcnt(kk["A_Bq"]) + " Бк"
    f["berry_k_kg"] = lambda: bp.rcnt(kk["per_kg"]) + " Бк/кг"
    f["berry_k_noise"] = lambda: "± " + bp.rcnt(fa["k_noise_Bq"]) + " Бк"
    f["berry_k_sys"] = lambda: "± " + bp.rcnt(fa["k_sys_Bq"]) + " Бк"

    f["berry_cs_sys"] = lambda: bp.rnum(fa["cs_sys_pct"], 0) + "\u00a0%"
    f["berry_cs_var"] = lambda: bp.rnum(fa["variants_cs_pct"], 0) + "\u00a0%"

    f["berry_sl_cs"] = lambda: bp.rnum(sl["cs_now_per_kg"], 1) + " Бк/кг"
    f["berry_sl_cs24"] = lambda: bp.rcnt(sl["cs_2024_per_kg"]) + " Бк/кг"
    f["berry_sl_k"] = lambda: bp.rcnt(sl["k_per_kg"]) + " Бк/кг"
    f["berry_sl_years"] = lambda: bp.rnum(sl["years"], 3) + " года"

    f["berry_ratio_cs"] = lambda: bp.rnum(cs["ref_ratio"], 3)
    f["berry_ratio_k"] = lambda: bp.rnum(kk["ref_ratio"], 3)
    f["berry_dev_cs"] = lambda: sgn(cs["ref_ratio"])
    f["berry_dev_k"] = lambda: sgn(kk["ref_ratio"])

    f["berry_chi2"] = lambda: bp.rnum(m1["chi2_ndof"], 2)
    f["berry_ndof"] = lambda: bp.rcnt(m1["ndof"])
    f["berry_shape662"] = lambda: bp.rnum(m1["shape_662"], 2)
    f["berry_shape1460"] = lambda: bp.rnum(m1["shape_1460"], 2)

    f["berry_edge662"] = lambda: bp.rnum(fa["edge_662_pct"], 1) + "\u00a0%"
    f["berry_edge1461"] = lambda: bp.rnum(fa["edge_1461_pct"], 2) + "\u00a0%"
    f["berry_resid"] = lambda: bp.rnum(fa["resid_pct"], 0) + "\u00a0%"
    f["berry_resid_band"] = lambda: bp.rcnt(fa["resid_lo_keV"]) + "–" + bp.rcnt(fa["resid_hi_keV"]) + "\u00a0кэВ"
    f["berry_gain_w"] = lambda: bp.rnum(fa["gain_w"], 2)
    f["berry_gain_s"] = lambda: bp.rnum(fa["gain_s_pct"], 1) + "\u00a0%"

    f["berry_kfrac_model"] = lambda: bp.rnum(fa["k_frac_model_pct"], 1) + "\u00a0%"
    f["berry_kfrac_found"] = lambda: bp.rnum(fa["k_frac_found_pct"], 1) + "\u00a0%"
    f["berry_beta_img"] = lambda: sl["img"]

    m2 = d["method2"]     # #GS-74: метод 2 посчитан (export_berry_m2.py), числа текста — подстановками
    p2 = m2["per_nuclide"]
    f["berry_m2_cs_a"] = lambda: bp.rcnt(p2["CS137"]["A_Bq"]) + " Бк"
    f["berry_m2_cs_kg"] = lambda: bp.rcnt(p2["CS137"]["per_kg"]) + " Бк/кг"
    f["berry_m2_k_a"] = lambda: bp.rcnt(p2["K40"]["A_Bq"]) + " Бк"
    f["berry_m2_k_kg"] = lambda: bp.rcnt(p2["K40"]["per_kg"]) + " Бк/кг"
    f["berry_m2_dev_cs"] = lambda: sgn(p2["CS137"]["ref_ratio"]) + " (Cs-137)"
    f["berry_m2_dev_k"] = lambda: sgn(p2["K40"]["ref_ratio"]) + " (K-40)"
    f["berry_m2_chi2"] = lambda: bp.rnum(m2["chi2_ndof"], 2)
    f["berry_m2_nodes"] = lambda: bp.rcnt(m2["n_nodes"])
    f["berry_m2_nch"] = lambda: bp.rcnt(m2["n_channels_fit"])
    f["berry_m2_nlines"] = lambda: bp.rcnt(m2["n_lines"])
    f["berry_m2_nxray"] = lambda: bp.rcnt(m2["n_xray_energies"])
    f["berry_m2_shape662"] = lambda: bp.rnum(m2["shape_662"], 2)
    f["berry_m2_shape1460"] = lambda: bp.rnum(m2["shape_1460"], 2)

    return __import__("ru_rules").fill(f)   # #GS-78 r6: тысячи неразрывным пробелом


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
    with io.open(SRC_BERRY_JS, "w", encoding="utf-8", newline="") as f:
        f.write(js)

    d, djs = data_js()
    fill = make_fill(bp, d)

    panel, n_p = render(bp, PANEL_HTML, fill)
    pops, n_q = render(bp, POPS_HTML, fill)

    for name, text in (("scripts/g1s-berry.js", js), ("data-berry.js", djs), ("berry-panel.html", panel), ("berry-pops.html", pops)):
        bad = bp.secret_scan(text)
        if bad:
            refuse("секрет-скан в %s: %s" % (name, bad))

    os.makedirs(os.path.join(bp.DIST, "scripts"), exist_ok=True)
    with io.open(os.path.join(bp.DIST, "data-berry.js"), "w", encoding="utf-8") as f:
        f.write(djs)
    with io.open(os.path.join(bp.DIST, "scripts", "g1s-berry.js"), "w", encoding="utf-8") as f:
        f.write(js)

    v = int(max(os.path.getmtime(p) for p in (DATA_JSON, SRC_BERRY_JS, PANEL_HTML, POPS_HTML, TERMS_JSON, os.path.join(PATCH_DIR, "berry_overrides.js"))))

    linked = '<script src="data-berry.js?v=%d"></script>\n<script src="scripts/g1s-berry.js?v=%d"></script>' % (v, v)
    inline = "<script>" + djs + "</script>\n<script>\n" + js + "</script>"

    for path, scripts in ((os.path.join(bp.DIST, "index.html"), linked), (bp.SINGLE, inline)):
        t = read(path)
        for marker, rep in zip(MARKERS, (panel, pops, scripts)):
            if t.count(marker) != 1:
                refuse("метка %s в %s встречается %d раз" % (marker, path, t.count(marker)))
            t = t.replace(marker, rep)
        with io.open(path, "w", encoding="utf-8") as f:
            f.write(t)

    print("Черника: замен в JS %d, g1s-berry.js %d КБ, data-berry.js %d КБ, v=%d" % (n_terms, len(js.encode()) // 1024, len(djs.encode()) // 1024, v))
    print("Черника: подстановок в разметке %d (панель) + %d (окна)" % (n_p, n_q))
    print("Черника: гейт %s = 0; справочно (комментарии/недостижимые ветви): %s" % ("/".join(GATE_WORDS), info))

    return ["data-berry.js", os.path.join("scripts", "g1s-berry.js")]


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    js, n_terms, info = build_js()
    d, djs = data_js()
    with io.open(SRC_BERRY_JS, "w", encoding="utf-8", newline="") as f:
        f.write(js)
    print("Черника (отдельный запуск): замен %d, g1s-berry.js %d КБ, data-berry.js %d КБ; справочно %s" % (n_terms, len(js.encode()) // 1024, len(djs.encode()) // 1024, info))
    print("разметку и dist собирает build_page_gs2020.py")
