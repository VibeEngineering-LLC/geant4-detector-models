# -*- coding: utf-8 -*-
"""Собирает docs/gamma-1s-article/index.html из article/gamma1s-method.md.

Источник — один markdown-файл; здесь только вёрстка (боковое оглавление,
панель сводных чисел). Числа в .stats переносятся ДОСЛОВНО из реферата и
§5.1 текущей редакции статьи — не пересчитываются и не берутся по памяти.
Вкладки .tabs/.panel из article.css не задействованы: в тексте статьи нет
данных по геометриям, готовых для такой раскладки (это отдельная задача).

Запуск: python build_article.py (из этой же папки).
"""
import html
import os
import re

import markdown

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SRC = os.path.join(REPO, "article", "gamma1s-method.md")
OUT_DIR = os.path.join(REPO, "docs", "gamma-1s-article")
OUT = os.path.join(OUT_DIR, "index.html")

TITLE = "Расчётная модель сцинтилляционного гамма-спектрометра: протокол проверки"
DESCRIPTION = ("Двухчастный протокол проверки расчётной модели гамма-спектрометра "
               "методом Монте-Карло: разложение расхождения с аттестованной "
               "кривой на именованные вклады и набор замыканий без эталона.")

# Дословно из article/gamma1s-method.md (реферат и §5.1) — при правке текста
# статьи сверять эти строки заново, не по памяти.
STATS = [
    ("Ø63×63 мм", "кристалл NaI(Tl)", False),
    ("5", "аттестованных геометрий измерения", False),
    ("~5 из 27,6 п.п.", "конвенция съёма площади на жёстком крае", True),
    ("~10 %", "неопределённость модели, средний участок", False),
    ("не закрыт", "баланс вкладов на жёстком крае", False),
    ("7", "замыканий модели без эталона", False),
]


def slugify(text):
    text = re.sub(r"<[^>]+>", "", text)
    text = text.lower()
    text = re.sub(r"[^\w\s-]", "", text, flags=re.UNICODE)
    text = re.sub(r"[\s_]+", "-", text.strip())
    return text or "sec"


def build_toc(md_text):
    """Строит id/оглавление по заголовкам ## и ###, слаг — детерминированный
    (не зависит от порядка вызова, в отличие от toc-расширения markdown)."""
    entries = []
    seen = {}
    for line in md_text.splitlines():
        m = re.match(r"^(#{2,3})\s+(.*)$", line)
        if not m:
            continue
        level = len(m.group(1))
        text = m.group(2).strip()
        slug = slugify(text)
        n = seen.get(slug, 0)
        seen[slug] = n + 1
        if n:
            slug = "%s-%d" % (slug, n)
        entries.append((level, text, slug))
    return entries


def inject_ids(html_body, entries):
    """Проставляет id= в <h2>/<h3> в порядке их появления — тот же порядок,
    что дал build_toc по исходному markdown."""
    it = iter(entries)
    positions = {2: [], 3: []}
    for lvl, text, slug in entries:
        positions[lvl].append(slug)
    counters = {2: 0, 3: 0}

    def repl(m):
        lvl = int(m.group(1))
        slugs = positions.get(lvl, [])
        idx = counters[lvl]
        counters[lvl] = idx + 1
        slug = slugs[idx] if idx < len(slugs) else slugify(m.group(2))
        return '<h%d id="%s">%s</h%d>' % (lvl, slug, m.group(2), lvl)

    return re.sub(r"<h([23])>(.*?)</h\1>", repl, html_body)


def render_stats():
    cells = []
    for v, k, acc in STATS:
        vcls = ' class="v acc"' if acc else ' class="v"'
        cells.append('<div><div%s>%s</div><div class="k">%s</div></div>'
                     % (vcls, html.escape(v), html.escape(k)))
    return '<div class="stats">%s</div>' % "".join(cells)


def render_toc(entries):
    out = ['<span class="toc-kicker">Гамма-1С</span>',
           '<span class="toc-title">Протокол проверки модели</span>']
    for lvl, text, slug in entries:
        if lvl != 2:
            continue
        text_plain = re.sub(r"[*`]", "", text)
        out.append('<a href="#%s">%s</a>' % (slug, html.escape(text_plain)))
    out.append('<div class="toc-foot">Черновик; текст не прошёл повторный '
                'состязательный аудит после отзыва 02.08.2026. '
                '<a href="https://github.com/VibeEngineering-LLC/'
                'geant4-detector-models/blob/main/article/gamma1s-method.md">'
                'исходник на GitHub</a></div>')
    return "\n".join(out)


TMPL = """<meta charset="utf-8">
<title>{title}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="{description}">
<meta name="theme-color" content="#ffffff" media="(prefers-color-scheme: light)">
<meta name="theme-color" content="#1a1a1a" media="(prefers-color-scheme: dark)">
<link rel="stylesheet" href="../assets/page.css">
<link rel="stylesheet" href="../assets/article.css">

<div class="shell">
<nav class="toc">
{toc}
</nav>
<main class="page" lang="ru">
{stats}
{body}
</main>
</div>
"""


def main():
    md_text = open(SRC, encoding="utf-8").read()
    entries = build_toc(md_text)
    body = markdown.markdown(md_text, extensions=["tables", "fenced_code",
                                                    "sane_lists"])
    body = inject_ids(body, entries)
    toc_html = render_toc(entries)
    stats_html = render_stats()
    out = TMPL.format(title=html.escape(TITLE),
                       description=html.escape(DESCRIPTION),
                       toc=toc_html, stats=stats_html, body=body)
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(out)
    print("написано:", OUT, "(%d байт)" % len(out.encode("utf-8")))
    print("разделов в оглавлении:", sum(1 for e in entries if e[0] == 2))


if __name__ == "__main__":
    main()
