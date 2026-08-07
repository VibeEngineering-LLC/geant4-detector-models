# -*- coding: utf-8 -*-
"""Сборка страницы из отдельных исходников.

    src/index.html          разметка с метками {{...}}
    src/styles/asn16.css    оформление
    src/scripts/asn16.js    отрисовка и поведение
    asn16_data.json         выгрузка расчёта (export_data.py)

Числа в прозе подставляются ЗДЕСЬ, а не в браузере: нераскрытая метка
роняет сборку, а не страницу у читателя. Раньше проверки не было — сборка
молча выпускала битую метку, и на странице от этого падал весь скрипт.

Выход:
    dist/index.html + dist/styles/ + dist/scripts/ + dist/data.js
        — раздельные файлы, для публикации на статическом хостинге;
    asn16_response.html
        — то же одним файлом, для площадок, принимающих один файл.

    python build_page.py
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
DATA_JSON = os.path.join(HERE, "asn16_data.json")
SINGLE = os.path.join(HERE, "asn16_response.html")

TOKEN = re.compile(r"\{\{([a-z0-9_]+)((?:\s+[^\s}]+)*)\}\}", re.I)


class Bad(SystemExit):
    pass


def need(cond, msg):
    if not cond:
        raise Bad("сборка отклонена: " + msg)


# ── формат чисел: те же правила, что и в скрипте страницы ────────────────
def rnum(x, d=1):
    return ("%.*f" % (d, x)).replace(".", ",")


def rcnt(x):
    s = "%d" % round(x)
    out, neg = "", s.startswith("-")
    if neg:
        s = s[1:]
    while len(s) > 3:
        out = "\u00a0" + s[-3:] + out
        s = s[:-3]
    return ("-" if neg else "") + s + out


def rpct(x):
    return rnum(x, 2 if x < 1 else 1) + "\u00a0%"


def rkev(x):
    return re.sub(r",0$", "", rnum(float(x), 1))


def rplural(n, one, few, many):
    a, b = abs(n) % 100, abs(n) % 10
    w = many if (10 < a < 20) or b > 4 or b == 0 else (one if b == 1 else few)
    return "%d %s" % (n, w)


def make_fill(d):
    tabs = {t["e0"]: t for t in d["tabs"]}

    def tab(e0):
        e0 = float(e0)
        need(e0 in tabs, "нет узла %s" % e0)
        return tabs[e0]

    def chan(e0, key):
        t = tab(e0)
        c = [c for c in t["channels"] if c["key"] == key]
        need(c, "на узле %s нет канала %s" % (e0, key))
        return c[0]

    def feat(e0, name):
        t = tab(e0)
        need(name in t["feat"], "на узле %s нет величины %s" % (e0, name))
        return t["feat"][name]

    def run(name):
        need(name in d["run"], "нет параметра прогона %s" % name)
        v = d["run"][name]
        return rkev(v) if isinstance(v, (int, float)) else str(v)

    return {
        "p": lambda key, e0: rpct(chan(e0, key)["pct"]),
        "n": lambda key, e0: rcnt(chan(e0, key)["n"]),
        "f": lambda e0, name: rkev(feat(e0, name)),
        "fn": lambda e0, name: rcnt(feat(e0, name)),
        "run": run,
        "runN": lambda name: rcnt(d["run"][name]),
        "nodes": lambda: rplural(d["run"]["n_nodes"], "узел", "узла", "узлов"),
        "fwhm": lambda: rnum(d["fwhm_662"], 2),
        "mec2": lambda: rnum(d["mec2"], 5),
        "peakhalf": lambda: rkev(d["tabs"][0]["peak_half"]),
        "sang": lambda: rnum(d["run"]["solid_angle_frac"], 4),
        "stamp": lambda: d["stamp"],
        "git": lambda: d["git_describe"],
        "matN": lambda: str(len(d["matrix"]["es"])),
        "matM": lambda: str(len(d["matrix"]["cols"])),
        "psum": lambda e0, *keys: rpct(sum(chan(e0, k)["pct"] for k in keys)),
        "ratio": lambda e0, a, b: rnum(chan(e0, a)["pct"] / chan(e0, b)["pct"], 1),
        # Доля канала СРЕДИ СОБЫТИЙ БЕЗ ВЫЛЕТА — не среди всех. Рядом со
        # словом «пик» доля от всех событий читается неверно.
        "share": lambda key, e0: rpct(
            100.0 * chan(e0, key)["pct"] / tab(e0)["nofly_pct"]),
    }


def substitute(text, fill, where):
    used = []

    def one(m):
        name, rest = m.group(1), m.group(2).split()
        need(name in fill, "%s: неизвестная подстановка %s" % (where, m.group(0)))
        try:
            out = fill[name](*rest)
        except TypeError as ex:
            raise Bad("%s: %s — неверные аргументы (%s)"
                      % (where, m.group(0), ex))
        need(out is not None and "nan" not in str(out).lower(),
             "%s: %s не разрешилась" % (where, m.group(0)))
        used.append(m.group(0))
        return str(out)

    out = TOKEN.sub(one, text)
    need("{{" not in out, "%s: остались нераскрытые метки" % where)
    return out, used


def secret_scan(text):
    bad = []
    for pat, why in ((r"C:\\\\Users\\\\[^\s\"'<]*", "локальный путь"),
                     (r"C:/Users/[^\s\"'<]*", "локальный путь"),
                     (r"[\w.+-]+@[\w-]+\.[\w.]+", "адрес почты"),
                     (r"gh[pousr]_[A-Za-z0-9]{16,}", "токен GitHub")):
        for m in re.findall(pat, text):
            bad.append("%s: %s" % (why, m))
    return bad


def main():
    d = json.loads(io.open(DATA_JSON, encoding="utf-8").read())
    html = io.open(os.path.join(SRC, "index.html"), encoding="utf-8").read()
    style = io.open(os.path.join(SRC, "styles", "asn16.css"), encoding="utf-8").read()
    script = io.open(os.path.join(SRC, "scripts", "asn16.js"), encoding="utf-8").read()

    # Сторож против чисел расчёта, набранных руками. Прецедент в проекте уже
    # был: литералы в шаблоне при заявлении «чисел руками нет».
    pcts = re.findall(r"\d+,\d+\s*%", html)
    need(not pcts, "проценты набраны цифрами, нужна подстановка: %s" % pcts)

    fill = make_fill(d)
    html, used = substitute(html, fill, "index.html")
    need(used, "в разметке не нашлось ни одной подстановки — шаблон не тот?")

    raw = json.dumps(d, ensure_ascii=False, separators=(",", ":"))
    need("</script" not in raw.lower(), "в данных есть закрывающий тег script")
    data_js = "window.ASN16=" + raw + ";\n"

    need("<!--@styles-->" in html and "<!--@script-->" in html,
         "в разметке нет точек вставки оформления и кода")

    linked = (html
              .replace("<!--@styles-->",
                       '<link rel="stylesheet" href="styles/asn16.css">')
              .replace("<!--@script-->",
                       '<script src="data.js"></script>\n'
                       '<script src="scripts/asn16.js"></script>'))
    single = (html
              .replace("<!--@styles-->", "<style>\n" + style + "</style>")
              .replace("<!--@script-->",
                       "<script>" + data_js + "</script>\n"
                       "<script>\n" + script + "</script>"))

    # Секрет-скан ДО записи: раньше файл с утечкой успевал лечь на диск.
    for name, text in (("dist/index.html", linked), ("dist/data.js", data_js),
                       ("dist/styles/asn16.css", style),
                       ("dist/scripts/asn16.js", script),
                       ("asn16_response.html", single)):
        bad = secret_scan(text)
        if bad:
            print("СЕКРЕТ-СКАН — НАЙДЕНО в %s:" % name)
            for b in bad:
                print("  " + b)
            raise Bad("сборка остановлена, файлы не записаны")

    os.makedirs(os.path.join(DIST, "styles"), exist_ok=True)
    os.makedirs(os.path.join(DIST, "scripts"), exist_ok=True)
    io.open(os.path.join(DIST, "index.html"), "w", encoding="utf-8").write(linked)
    io.open(os.path.join(DIST, "data.js"), "w", encoding="utf-8").write(data_js)
    io.open(os.path.join(DIST, "styles", "asn16.css"), "w",
            encoding="utf-8").write(style)
    io.open(os.path.join(DIST, "scripts", "asn16.js"), "w",
            encoding="utf-8").write(script)
    io.open(SINGLE, "w", encoding="utf-8").write(single)

    kb = lambda p: os.path.getsize(p) / 1024.0
    print("подстановок разрешено: %d" % len(used))
    print("dist/index.html      %6.0f КБ" % kb(os.path.join(DIST, "index.html")))
    print("dist/data.js         %6.0f КБ" % kb(os.path.join(DIST, "data.js")))
    print("dist/styles/asn16.css%6.0f КБ" % kb(os.path.join(DIST, "styles", "asn16.css")))
    print("dist/scripts/asn16.js%6.0f КБ" % kb(os.path.join(DIST, "scripts", "asn16.js")))
    print("asn16_response.html  %6.0f КБ  (одним файлом)" % kb(SINGLE))
    print("секрет-скан: чисто")
    return 0


if __name__ == "__main__":
    sys.exit(main())
