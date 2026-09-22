"""Склеивает index.html и data.js в один самодостаточный файл (data.js вставляется вместо <script src="data.js"></script>).
Использование: python inline_page.py <каталог страницы> <выходной .html>"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8")
d, out = sys.argv[1], sys.argv[2]
html = open(os.path.join(d, "index.html"), encoding="utf-8").read()
data = open(os.path.join(d, "data.js"), encoding="utf-8").read()
tag = '<script src="data.js"></script>'
if html.count(tag) != 1:
    sys.exit("в index.html нет ровно одного тега " + tag)
open(out, "w", encoding="utf-8", newline="\n").write(html.replace(tag, "<script>\n" + data + "\n</script>"))
print("записано", out, os.path.getsize(out) // 1024, "КБ")
