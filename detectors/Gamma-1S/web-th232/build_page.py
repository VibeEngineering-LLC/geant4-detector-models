# -*- coding: utf-8 -*-
"""Сборка страницы «Разложение спектра Th-232 в Маринелли: два метода».

    src/index.html                разметка с метками {{...}}
    src/styles/g1s-th232.css      оформление (адаптация «плаката» ASN16)
    src/scripts/g1s-th232.js      отрисовка и поведение
    g1s_th232_data.json           выгрузка расчёта (export_data.py)

Числа в прозе подставляются ЗДЕСЬ — нераскрытая метка роняет сборку, а не
страницу. Правило то же, что в конвейере ASN16: правки геометрии/анализа
попадают в текст через пересборку, а не через ручное редактирование.

    python export_data.py && python build_page.py

Выход:
    dist/index.html + dist/styles/ + dist/scripts/ + dist/data.js
    g1s_th232.html — то же одним файлом
"""
import io
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "src")
DIST = os.path.join(HERE, "dist")
DATA_JSON = os.path.join(HERE, "g1s_th232_data.json")
SINGLE = os.path.join(HERE, "g1s_th232.html")

TOKEN = re.compile(r"\{\{([a-z0-9_]+)((?:\s+[^\s}]+)*)\}\}", re.I)


class Bad(SystemExit):
    pass


def need(cond, msg):
    if not cond:
        raise Bad("сборка отклонена: " + msg)


# ── формат чисел (RU: запятая, неразрывный пробел разделяет тысячи) ──
def rnum(x, d=1):
    return ("%.*f" % (d, x)).replace(".", ",")


def rcnt(x):
    """Целое с разделением тысяч неразрывным пробелом."""
    s = "%d" % round(float(x))
    out, neg = "", s.startswith("-")
    if neg:
        s = s[1:]
    while len(s) > 3:
        out = "\u00a0" + s[-3:] + out
        s = s[:-3]
    return ("-" if neg else "") + s + out


def rpct(x, d=None):
    if d is None:
        d = 2 if x < 1 else 1
    return rnum(x, d) + "\u00a0%"


def rkev(x):
    return re.sub(r",0$", "", rnum(float(x), 1))


def ract(x):
    """Активность в Бк: округление до целого, разделение тысяч."""
    return rcnt(x) + "\u00a0Бк"


def make_fill(d):
    """Подстановки для ПРОЗЫ страницы.

    Числа, которые на странице переключаются (активность и χ² при учёте
    и без учёта аппаратной ПШПВ, при отобранной и полной библиотеке),
    здесь НЕ подставляются: у них нет одного значения. Их выводит
    скрипт страницы из той же выгрузки — источник тот же, момент
    подстановки другой.
    """
    p = d["passport"]
    m2 = d["method2"]
    m2f = d["method2_full"]
    fw = d["fwhm_cal"]
    lib = d["library"]

    fill = {
        # --- паспорт -----------------------------------------------------
        "apass":      lambda: ract(p["A_Bq"]),
        "apass_dev":  lambda: "± " + rcnt(p["dA_Bq"]) + " Бк",
        "aksp":       lambda: rcnt(p["Bq_per_kg"]) + " Бк/кг",
        "aksp_unc":   lambda: rpct(p["unc_pct"]),
        "mass":       lambda: rnum(p["mass_g"], 0) + " г",
        "date_pass":  lambda: p["date_certified"],
        "date_meas":  lambda: p["date_measured"],
        "decayf":     lambda: rnum(p["decay_factor"], 6),

        # --- измерение и фон --------------------------------------------
        "tlive":      lambda: rnum(d["meta"]["live_s"], 0) + " с",
        "treal":      lambda: rnum(d["meta"]["real_s"], 0) + " с",
        "sysf":       lambda: rpct(d["meta"]["sys_floor_pct"], 0),
        "xr_lo":      lambda: rkev(d["meta"]["k_xray_lo_keV"]),
        "xr_hi":      lambda: rkev(d["meta"]["k_xray_hi_keV"]),

        # --- закон ширины, снятый с этого же спектра ---------------------
        "fw_k":       lambda: rnum(fw["k"], 3),
        "fw_p":       lambda: rnum(fw["p"], 4),
        "fw_rms":     lambda: rpct(fw["rms_dev_pct"], 1),
        "fw_nused":   lambda: rcnt(fw["n_used"]),
        "fw_662":     lambda: rnum(fw["fwhm662_law"], 1) + " кэВ",
        "fw_cs":      lambda: rnum(fw["fwhm662_cs"], 1) + " кэВ",

        # --- состав библиотек метода 2 -----------------------------------
        "m2_nlines":  lambda: rcnt(m2["n_lines"]),
        "m2f_nlines": lambda: rcnt(m2f["n_lines"]),
        "m2_nsum":    lambda: rcnt(m2["n_sum_peaks"]),
        "m2_nsumtot": lambda: rcnt(m2["n_sum_peaks_total"]),
        "m2_nxray":   lambda: rcnt(m2["n_xray_energies"]),
        "m2_ithresh": lambda: rnum(lib["i_threshold_pct"], 0) + " %",
    }
    return fill


def substitute(text, fill, where):
    used = []

    def one(m):
        name, rest = m.group(1), m.group(2).split()
        need(name in fill, "%s: неизвестная подстановка %s" % (where, m.group(0)))
        try:
            out = fill[name](*rest)
        except (TypeError, ValueError) as ex:
            raise Bad("%s: %s — неверные аргументы (%s)"
                      % (where, m.group(0), ex))
        need(out is not None
             and "nan" not in str(out).lower()
             and "inf" not in str(out).lower(),
             "%s: %s не разрешилась" % (where, m.group(0)))
        used.append(m.group(0))
        return str(out)

    out = TOKEN.sub(one, text)
    need("{{" not in out, "%s: остались нераскрытые метки" % where)
    return out, used


# --- секрет-скан (тот же, что в ASN16-конвейере) ------------------------
SECRETS = ((r"[A-Za-z]:[\\/]+Users[\\/]+[^\s\"'<>)]*", "локальный путь"),
           (r"[\w.+-]+@[\w-]+\.[\w.]{2,}", "адрес почты"),
           (r"gh[pousr]_[A-Za-z0-9]{16,}", "токен GitHub"))


def secret_scan(text):
    bad = []
    for pat, why in SECRETS:
        for m in re.findall(pat, text):
            bad.append("%s: %s" % (why, m))
    return bad


def secret_scan_selftest():
    dirty = ("path C:" + chr(92) + "Users" + chr(92) + "someone" + chr(92)
             + "secret.txt and C:/Users/someone/x and a@b.co and "
             + "ghp_" + "0123456789abcdef0123")
    hits = secret_scan(dirty)
    need(len(hits) >= 4, "самопроверка секрет-скана провалена: %s" % hits)
    need(not secret_scan("обычный текст, путь C:/g4work/gamma1s"),
         "секрет-скан срабатывает на чистом тексте")


def main():
    secret_scan_selftest()
    d = json.loads(io.open(DATA_JSON, encoding="utf-8").read())
    html = io.open(os.path.join(SRC, "index.html"), encoding="utf-8").read()
    style = io.open(os.path.join(SRC, "styles", "g1s-th232.css"),
                    encoding="utf-8").read()
    script = io.open(os.path.join(SRC, "scripts", "g1s-th232.js"),
                     encoding="utf-8").read()

    # Сторож против чисел расчёта, набранных цифрами в тексте разметки.
    # Тот же приём: и «12,3 %», и «17 %», и «в 2,0 раза».
    pcts = (re.findall(r"\d+(?:,\d+)?\s*%", html)
            + re.findall(r"в\s+\d+,\d+\s+раза", html))
    need(not pcts,
         "числа расчёта набраны цифрами, нужна подстановка: %s" % pcts)

    fill = make_fill(d)
    html_filled, used = substitute(html, fill, "index.html")
    need(used, "в разметке не нашлось ни одной подстановки")

    raw = json.dumps(d, ensure_ascii=False, separators=(",", ":"))
    need("</script" not in raw.lower(), "в данных есть закрывающий тег script")
    data_js = "window.G1S=" + raw + ";\n"

    need("<!--@styles-->" in html_filled and "<!--@script-->" in html_filled,
         "в разметке нет точек вставки оформления и кода")
    for nm, txt in (("styles/g1s-th232.css", style),
                    ("scripts/g1s-th232.js", script)):
        need("{{" not in txt, "%s: метка подстановки вне разметки" % nm)

    # Метка версии на КАЖДЫЙ файл (data.js/css/js), не только на страницу:
    # без неё браузер после переразвёртывания держит старые модули по их
    # неизменным URL, пока html_filled с меткой ?v= уже новый — страница
    # выглядит «не изменилась», хотя сервер отдаёт свежее (прецедент того
    # же класса на других страницах контура, см. rn-article-style §9).
    # Источник метки — САМЫЙ свежий mtime среди выгрузки И исходников
    # (src/*): одна лишь выгрузка не покрывает правку JS/CSS/HTML без
    # повторного прогона export_data.py — ровно этот пробел поймал на
    # себе в этой же сессии (легенда правилась, данные не пересчитывались,
    # старое v= держало браузер на прежнем скрипте до жёсткой перезагрузки).
    v = int(max(
        os.path.getmtime(DATA_JSON),
        os.path.getmtime(os.path.join(SRC, "index.html")),
        os.path.getmtime(os.path.join(SRC, "scripts", "g1s-th232.js")),
        os.path.getmtime(os.path.join(SRC, "styles", "g1s-th232.css")),
    ))
    linked = (html_filled
              .replace("<!--@styles-->",
                       '<link rel="stylesheet" href="styles/g1s-th232.css?v=%d">' % v)
              .replace("<!--@script-->",
                       '<script src="data.js?v=%d"></script>\n'
                       '<script src="scripts/g1s-th232.js?v=%d"></script>'
                       % (v, v)))
    single = (html_filled
              .replace("<!--@styles-->", "<style>\n" + style + "</style>")
              .replace("<!--@script-->",
                       "<script>" + data_js + "</script>\n"
                       "<script>\n" + script + "</script>"))

    for name, text in (("dist/index.html", linked), ("dist/data.js", data_js),
                       ("dist/styles/g1s-th232.css", style),
                       ("dist/scripts/g1s-th232.js", script),
                       ("g1s_th232.html", single)):
        bad = secret_scan(text)
        if bad:
            print("СЕКРЕТ-СКАН — НАЙДЕНО в %s:" % name)
            for b in bad:
                print("  " + b)
            raise Bad("сборка остановлена, файлы не записаны")

    os.makedirs(os.path.join(DIST, "styles"), exist_ok=True)
    os.makedirs(os.path.join(DIST, "scripts"), exist_ok=True)

    def put(path, text):
        with io.open(path, "w", encoding="utf-8") as g:
            g.write(text)

    put(os.path.join(DIST, "index.html"), linked)
    put(os.path.join(DIST, "data.js"), data_js)
    put(os.path.join(DIST, "styles", "g1s-th232.css"), style)
    put(os.path.join(DIST, "scripts", "g1s-th232.js"), script)
    put(SINGLE, single)

    kb = lambda p: os.path.getsize(p) / 1024.0
    print("подстановок разрешено: %d" % len(used))
    print("dist/index.html          %6.0f КБ" % kb(os.path.join(DIST, "index.html")))
    print("dist/data.js             %6.0f КБ" % kb(os.path.join(DIST, "data.js")))
    print("dist/styles/g1s-th232.css%6.0f КБ" % kb(os.path.join(DIST, "styles", "g1s-th232.css")))
    print("dist/scripts/g1s-th232.js%6.0f КБ" % kb(os.path.join(DIST, "scripts", "g1s-th232.js")))
    print("g1s_th232.html           %6.0f КБ  (одним файлом)" % kb(SINGLE))
    print("секрет-скан: чисто")
    return 0


if __name__ == "__main__":
    sys.exit(main())
