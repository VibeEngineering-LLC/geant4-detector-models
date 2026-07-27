# -*- coding: utf-8 -*-
"""Сборка интерактивной статьи из шаблона и посчитанных данных.

Один источник — docs/article.tmpl.html с меткой /*__DATA__*/ вместо данных.
На выходе два файла, отличающихся только обёрткой:

  docs/article-body.html — только содержимое, для публикации артефактом;
  docs/index.html        — полный документ, для GitHub Pages.

Данные (коэффициенты гладкого предела, эффективная толщина, mu/rho, кривые LSRM)
берутся из results/curves.json, который пишет curves.py --json. Никаких чисел
руками: если пересчитана сетка, статья обновляется одной командой.
"""
import base64
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))          # detectors/RadiaCode-103
REPO = os.path.abspath(os.path.join(ROOT, "..", ".."))
# Статья публикуется через GitHub Pages, поэтому живёт в общем docs/ репозитория,
# а не внутри каталога прибора: Pages отдаёт один каталог на весь репозиторий.
DOCS = os.path.join(REPO, "docs", "radiacode-103")
TMPL = os.path.join(DOCS, "article.tmpl.html")
DATA = os.path.join(ROOT, "results", "curves.json")

# Рисунки встраиваются как data-URI: страница обязана открываться и с GitHub
# Pages, и как артефакт, где внешние запросы запрещены политикой безопасности.
FIGS = {
    "__FIG_M500__": os.path.join(ROOT, "results", "m500", "figures", "lsrm_compare.png"),
    "__FIG_M200__": os.path.join(ROOT, "results", "m200", "figures", "lsrm_compare.png"),
}


def data_uri(path):
    with open(path, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode("ascii")

TITLE = "Кривые эффективности для RadiaCode в сосудах Маринелли"
DESC = ("Расчёт в Geant4 фотопиковой эффективности сцинтилляционного гамма-спектрометра "
        "RadiaCode 101/102/103 в авторских сосудах Маринелли 200 и 500 мл: "
        "самопоглощение в пробе, вклад бета-излучения, поправка на искажение фона "
        "пробой, сверка с независимым расчётом LSRM.")

HEAD = """<!doctype html>
<html lang="ru">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s</title>
<meta name="description" content="%s">
<meta property="og:title" content="%s">
<meta property="og:description" content="%s">
<meta property="og:type" content="article">
<style>
*,*::before,*::after{box-sizing:border-box}
html{-webkit-text-size-adjust:100%%}
body{margin:0}
img,canvas{max-width:100%%}
</style>
</head>
<body>
""" % (TITLE, DESC, TITLE, DESC)

TAIL = "\n</body>\n</html>\n"


def main():
    for p in (TMPL, DATA):
        if not os.path.exists(p):
            raise SystemExit("нет файла: " + p
                             + ("\n  сначала: python curves.py --json"
                                if p == DATA else ""))
    tmpl = open(TMPL, encoding="utf-8").read()
    data = json.load(open(DATA, encoding="utf-8"))
    # компактно: статья читает данные, а не человек
    blob = json.dumps(data, ensure_ascii=False, separators=(",", ":"))

    if "/*__DATA__*/" not in tmpl:
        raise SystemExit("в шаблоне нет метки /*__DATA__*/")
    body = tmpl.replace("/*__DATA__*/", blob)

    for mark, path in FIGS.items():
        if mark not in body:
            continue
        if not os.path.exists(path):
            raise SystemExit("нет рисунка: " + path)
        body = body.replace(mark, data_uri(path))

    left = [m for m in FIGS if m in body]
    if left:
        raise SystemExit("не подставлены рисунки: " + ", ".join(left))

    # секрет-скан: в публикуемую страницу не должны попадать локальные пути,
    # имя пользователя и адреса
    bad = [p for p in ("C:\\", "Users\\", "AppData", "@gmail", "USERPROFILE")
           if p in body]
    if bad:
        raise SystemExit("секрет-скан: в тексте найдено " + ", ".join(bad))

    # тема оформления переключается снаружи (артефакт ставит data-theme на :root),
    # для Pages нужен свой переключатель — но не ценой лишнего кода: страница и так
    # следует системной теме через prefers-color-scheme.
    out_body = os.path.join(DOCS, "article-body.html")
    out_page = os.path.join(DOCS, "index.html")
    open(out_body, "w", encoding="utf-8").write(body)
    open(out_page, "w", encoding="utf-8").write(HEAD + body + TAIL)

    for p in (out_body, out_page):
        print("%-40s %6.0f КБ" % (os.path.relpath(p, ROOT),
                                  os.path.getsize(p) / 1024))
    n = len(re.findall(r"<canvas", body))
    print("интерактивных графиков: %d, данных: %.0f КБ" % (n, len(blob) / 1024))


if __name__ == "__main__":
    main()
