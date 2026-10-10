# #GS-82: EN-версии рисунков GS2020 без правки русских генераторов: перевод Text через patches/en_figures.json, SVG — копией.
import sys, io, json, os, re
sys.stdout.reconfigure(encoding="utf-8"); os.chdir(os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use("Agg"); import matplotlib.text as mt, matplotlib.pyplot as plt
from PIL import Image
C = json.load(io.open("patches/en_figures.json", encoding="utf-8")); CYR = re.compile("[А-Яа-яЁё]"); LANG = "ru" if "--ru" in sys.argv else "en"
def tr(s, d):
    if LANG == "ru" or not isinstance(s, str) or not CYR.search(s): return s
    if s not in d: raise SystemExit("ОТКАЗ: нет перевода: " + s)
    return d[s]
_set = mt.Text.set_text; mt.Text.set_text = lambda self, s: _set(self, tr(s, C["text"]))
for gen, out, width in C["figures"]:
    out = out if LANG == "en" else os.path.join(os.environ.get("TEMP", "."), "ru_" + os.path.basename(out)); os.makedirs(os.path.dirname(out), exist_ok=True)
    src = re.sub(r"plt\.savefig\('[^']+'", lambda m: "plt.savefig(%r" % (out + ".tmp.png"), io.open(gen, encoding="utf-8").read())
    exec(compile(src, gen, "exec"), {"__name__": "__main__"}); plt.close("all")
    im = Image.open(out + ".tmp.png").convert("RGB"); im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
    im.quantize(256).save(out, optimize=True); os.remove(out + ".tmp.png"); print(out, im.size, os.path.getsize(out))
for src, out in C["svg"] if LANG == "en" else []:
    t = re.sub(r"<!--[\s\S]*?-->\s*", "", io.open(src, encoding="utf-8").read())
    t = re.sub(r"(<(?:text|title|tspan)\b[^>]*>)([^<]*)(<)", lambda m: m.group(1) + tr(m.group(2), C["svg_text"]) + m.group(3), t)
    if CYR.search(t): raise SystemExit("ОТКАЗ: кириллица в " + out)
    io.open(out, "w", encoding="utf-8", newline="").write(t); print(out, len(t))
