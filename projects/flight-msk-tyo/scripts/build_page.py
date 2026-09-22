"""Сборка web/flight-dad-sgn: text.md + bands.md + page.js + data.js -> index.html и склеенный spectrum-dad-sgn.html.
Использование: build_page.py <каталог страницы> <bands.md> [--lang en]"""
import sys, os, markdown, base64
d, bands = sys.argv[1], sys.argv[2]
LANG = "en" if "--lang" in sys.argv and sys.argv[sys.argv.index("--lang") + 1] == "en" else "ru"
rd = lambda n: open(os.path.join(d, n), encoding="utf-8").read()
import re
RENDER_PNG = os.path.join(os.path.dirname(bands), "..", "..", "outputs", "asn16_geant4_overview.png")
if os.path.exists(RENDER_PNG):
    _b64 = base64.b64encode(open(RENDER_PNG, "rb").read()).decode("ascii")
    _render_txt = {
        "ru": ('Рендер геометрии прибора АтомНано 16 в Geant4', 'Геометрия прибора в расчёте (Geant4, вид «три четверти», разрез корпуса): '
               'кристалл CsI(Tl) — золотистый параллелепипед, плата электроники — снизу, корпус и отражатель показаны как каркас/полупрозрачные грани.'),
        "en": ('Render of the АтомНано 16 instrument geometry in Geant4', 'Instrument geometry in the calculation (Geant4, three-quarter view, housing cutaway): '
               'the CsI(Tl) crystal is the golden block, the electronics board is below it, the housing and reflector are shown as a wireframe/translucent faces.')
    }[LANG]
    RENDER = ('<div class="fig" style="text-align:center;max-width:44rem;margin:0 auto">'
              '<img src="data:image/png;base64,%s" alt="%s" '
              'style="max-width:100%%;background:#0b0c10;border-radius:4px">'
              '<p style="font-size:.85rem;color:var(--dim);margin-top:.4rem">%s</p></div>') % (_b64, _render_txt[0], _render_txt[1])
else:
    RENDER = ""
md = lambda t: markdown.markdown(t, extensions=["tables"])
LANGBTN = (('<a href="@@LANGHREF@@" class="btn lang-switch" style="position:absolute;top:1rem;right:1rem">EN</a>') if LANG == "ru"
           else ('<a href="@@LANGHREF@@" class="btn lang-switch" style="position:absolute;top:1rem;right:1rem">RU</a>'))
if LANG == "en":
    HEAD = ('<div class="top" style="position:relative">%s<div class="title"><p class="eyebrow">АтомНано 16 · CsI 16.2 cm³ · Airbus A320/A321 · Geant4 11.4.2 + PARMA/EXPACS</p><h1>Spectrum at cruise altitude: flight DAD–SGN</h1><p class="stand">Model spectrum of cosmic radiation in the cabin at 9800 m, 20.09.2026, window seat, economy class, range 20 keV – 10 MeV. Units — counts per hour per keV; convolved with the instrument response function (FWHM 41.6 keV at 661.7 keV, ∝√E).</p></div></div>') % LANGBTN
else:
    HEAD = ('<div class="top" style="position:relative">%s<div class="title"><p class="eyebrow">АтомНано 16 · CsI 16,2 см³ · Airbus A320/A321 · Geant4 11.4.2 + PARMA/EXPACS</p><h1>Спектр на эшелоне: рейс DAD–SGN</h1><p class="stand">Модельный спектр космического излучения в салоне на высоте 9800 м, 20.09.2026, место у окна, эконом-класс, диапазон 20 кэВ – 10 МэВ. Единицы — отсчёты в час на кэВ; свёртка с аппаратной функцией (ПШПВ 41,6 кэВ на 661,7 кэВ, ∝√E).</p></div></div>') % LANGBTN
if LANG == "en":
    bt = "".join('<button type="button" class="btn" data-lo="%d" data-hi="%d">%s</button>' % x for x in [(20,100,"20–100"),(100,511,"100–511"),(511,1500,"511–1500"),(1500,3000,"1.5–3 MeV"),(3000,10000,"3–10 MeV"),(20,3000,"20 keV – 3 MeV"),(20,10000,"all")])
else:
    bt = "".join('<button type="button" class="btn" data-lo="%d" data-hi="%d">%s</button>' % x for x in [(20,100,"20–100"),(100,511,"100–511"),(511,1500,"511–1500"),(1500,3000,"1,5–3 МэВ"),(3000,10000,"3–10 МэВ"),(20,3000,"20 кэВ – 3 МэВ"),(20,10000,"весь")])
if LANG == "en":
    CHART = '<p class="method-lede">Drag on the chart to zoom in, double-click for the full range, click a table row to zoom to that line.</p><div class="toolbar" id="toolbar"><span class="grp"><span class="seg-label">Energy</span><div class="seg">%s</div></span><span class="sp"></span><label><input type="checkbox" id="logy-toggle" checked> logarithmic axis</label></div><div class="toolbar" style="margin-top:.6rem"><span class="grp"><span class="seg-label">Breakdown</span><div class="seg" id="mode-seg"><button type="button" class="btn" data-mode="comp">By field particle</button><button type="button" class="btn" data-mode="orig">By origin</button></div></span><span class="sp"></span><button type="button" class="btn" id="reset">Reset</button></div><p class="method-lede" id="mode-note"></p><div id="series-box"></div><canvas id="cv" tabindex="0" role="img" aria-label="Model spectrum"></canvas><div id="readout"></div>' % bt
else:
    CHART = '<p class="method-lede">Протяжка мышью по графику — приближение, двойной щелчок — весь диапазон, щелчок по строке таблицы — приближение к линии.</p><div class="toolbar" id="toolbar"><span class="grp"><span class="seg-label">Энергия</span><div class="seg">%s</div></span><span class="sp"></span><label><input type="checkbox" id="logy-toggle" checked> логарифмическая ось</label></div><div class="toolbar" style="margin-top:.6rem"><span class="grp"><span class="seg-label">Разложение</span><div class="seg" id="mode-seg"><button type="button" class="btn" data-mode="comp">По частицам поля</button><button type="button" class="btn" data-mode="orig">По происхождению</button></div></span><span class="sp"></span><button type="button" class="btn" id="reset">Сброс</button></div><p class="method-lede" id="mode-note"></p><div id="series-box"></div><canvas id="cv" tabindex="0" role="img" aria-label="Модельный спектр"></canvas><div id="readout"></div>' % bt
if LANG == "en":
    LINES = ('<div class="tabs" id="lines-bar" role="tablist"></div><div id="lines-detail" style="margin-top:.8rem"></div>'
             '<h2>Peaks outside the list</h2><div class="tabs" id="unlisted-bar" role="tablist"></div><div id="unlisted-detail" style="margin-top:.8rem"></div>')
else:
    LINES = ('<div class="tabs" id="lines-bar" role="tablist"></div><div id="lines-detail" style="margin-top:.8rem"></div>'
             '<h2>Пики вне перечня</h2><div class="tabs" id="unlisted-bar" role="tablist"></div><div id="unlisted-detail" style="margin-top:.8rem"></div>')
_scheme_lbl = {"aria": "A320 fuselage cross-section in the calculation", "people": "people (0.11 g/cm³ layer)",
    "cargo": "cargo hold", "cabin": "cabin, air", "inst": "instrument (R 14 cm tube)", "floor": "floor",
    "skin": "skin, framing: aluminum alloy 1.85 g/cm²"} if LANG == "en" else {
    "aria": "Сечение фюзеляжа A320 в расчёте", "people": "люди (слой 0,11 г/см³)", "cargo": "багажный отсек",
    "cabin": "салон, воздух", "inst": "прибор (трубка R 14 см)", "floor": "пол",
    "skin": "обшивка, набор: алюминиевый сплав 1,85 г/см²"}
SCHEME = ('<svg viewBox="-215 -215 430 430" role="img" aria-label="%(aria)s" style="max-width:26rem;width:100%%;display:block;margin:0 auto">'
  '<circle r="197.5" fill="none" stroke="currentColor" stroke-width="3"/><circle r="185" fill="none" stroke="currentColor" stroke-width="1"/>'
  '<polygon points="-75,165 75,165 150,90 -150,90" fill="#c9a86a" opacity=".55"/><rect x="-170" y="-60" width="340" height="120" fill="#e0868a" opacity=".45"/>'
  '<line x1="-170" y1="60" x2="170" y2="60" stroke="currentColor" stroke-width="3"/><circle cx="-150" cy="-15" r="14" fill="#f6d31c" stroke="currentColor" stroke-width="2"/>'
  '<g font-size="13" fill="currentColor" text-anchor="middle"><text y="-20">%(people)s</text><text y="112">%(cargo)s</text><text y="-168">%(cabin)s</text><text x="-150" y="-36" font-size="11">%(inst)s</text><text y="50" font-size="11">%(floor)s</text><text y="212" font-size="11">%(skin)s</text></g></svg>') % _scheme_lbl
NUMRE = re.compile(r"^[\s0-9.,%±×·eE\-−+()кмгµмкМэВ°]+$")
def mark_num(html):   # ставит class="num" числовым <td> (центрирование в CSS; текстовые ячейки не трогает)
    def repl(m):
        text = re.sub(r"<[^>]+>", "", m.group(1)).strip()
        return f'<td class="num">{m.group(1)}</td>' if text and NUMRE.match(text) else m.group(0)
    return re.sub(r"<td>(.*?)</td>", repl, html, flags=re.S)
_sfx = ".en.md" if LANG == "en" else ".md"
text_src = "text.en.md" if LANG == "en" else "text.md"
bands_src = bands.replace(".md", _sfx) if LANG == "en" else bands
field_src = os.path.join(os.path.dirname(bands), "field" + _sfx)
body = md(rd(text_src)).replace("<table>", '<div style="overflow-x:auto"><table>').replace("</table>", "</table></div>")
for key, val in (("HEAD", HEAD), ("CHART", CHART), ("BANDS", '<div style="overflow-x:auto">' + md(open(bands_src, encoding="utf-8").read()) + "</div>"), ("LINES", LINES), ("FIELD", '<div style="overflow-x:auto">' + md(open(field_src, encoding="utf-8").read()) + "</div>"), ("SCHEME", SCHEME)):
    body = body.replace("<p>@@%s@@</p>" % key, val)
svg = re.sub(r"<\?xml.*?\?>|<!DOCTYPE.*?>", "", open("results/eff/eff_curves.svg", encoding="utf-8").read(), flags=re.S)
for k, f in (("EFF", "results/final/eff" + _sfx), ("FUEL", "results/final/fuel" + _sfx), ("DOSE", "results/final/dose" + _sfx)):
    body = body.replace("<p>@@%s@@</p>" % k, '<div style="overflow-x:auto">' + md(open(f, encoding="utf-8").read()) + "</div>" if os.path.exists(f) else "")
nlf = "research/G-neutron-lines-summary" + _sfx
if os.path.exists(nlf):
    nl_text = "\n".join(open(nlf, encoding="utf-8").read().split("\n")[4:])   # без заголовка и вводного абзаца (уже в text.md)
    nl_table = md(nl_text).replace("<table>", '<table id="neu-lines">', 1)
    _nl_sum = "Show table (150 rows, click a column header to sort)" if LANG == "en" else "Показать таблицу (150 строк, клик по заголовку столбца сортирует)"
    _nl_locale = "en" if LANG == "en" else "ru"
    nl_html = ('<details><summary style="cursor:pointer">' + _nl_sum + '</summary>'
               '<div style="overflow-x:auto">' + nl_table + "</div></details>"
               '<script>(function(){const t=document.getElementById("neu-lines");if(!t)return;'
               'const tb=t.tBodies[0];[...t.tHead.rows[0].cells].forEach((th,i)=>{th.style.cursor="pointer";let asc=true;'
               'th.addEventListener("click",()=>{const rows=[...tb.rows];const numeric=i===0||i===5||i===6;'
               'rows.sort((a,b)=>{let x=a.cells[i].textContent.trim(),y=b.cells[i].textContent.trim();'
               'if(numeric){x=parseFloat(x.replace(",","."))||-Infinity;y=parseFloat(y.replace(",","."))||-Infinity;'
               'return asc?x-y:y-x;}return asc?x.localeCompare(y,"' + _nl_locale + '"):y.localeCompare(x,"' + _nl_locale + '");});'
               'asc=!asc;rows.forEach(r=>tb.appendChild(r));});});})();</script>')
    body = body.replace("<p>@@NEUTRONLINES@@</p>", nl_html)
body = body.replace("<p>@@EFFFIG@@</p>", '<div class="fig" style="overflow-x:auto;max-width:56rem">' + svg + "</div>")
body = body.replace("<p>@@RENDER@@</p>", RENDER)
body = mark_num(body)
_title = "Spectrum at cruise altitude" if LANG == "en" else "Спектр на эшелоне"
_lang_set = '<script>window.PAGE_LANG="en";</script>' if LANG == "en" else ""
_sp = "../" if LANG == "en" else ""   # styles/data.js на уровень выше из web/flight-dad-sgn/en/
_pagehead = ('<!doctype html><html lang="%s"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
        '<title>%s</title><link rel="stylesheet" href="%sstyles/g1s-th232.css"><link rel="stylesheet" href="%sstyles/flight.css">%s@@DATA@@</head>'
        '<body><div class="app doc">') % (LANG, _title, _sp, _sp, _lang_set)
page = _pagehead + body + "<script>\n" + rd("stack.js") + "\n" + rd("stackui.js") + "\n" + rd("page.js") + "\n</script></div></body></html>"
out_dir = os.path.join(d, "en") if LANG == "en" else d
os.makedirs(out_dir, exist_ok=True)
_data_name = "data.en.js" if LANG == "en" else "data.js"
_href_index = "en/index.html" if LANG == "ru" else "../index.html"
_href_glued = "en/spectrum-dad-sgn-en.html" if LANG == "ru" else "../spectrum-dad-sgn.html"
open(os.path.join(out_dir, "index.html"), "w", encoding="utf-8", newline="\n").write(
    page.replace("@@DATA@@", '<script src="%s%s"></script>' % (_sp, _data_name)).replace("@@LANGHREF@@", _href_index))
css = "".join("<style>\n" + rd("styles/" + n) + "\n</style>" for n in ("g1s-th232.css", "flight.css"))
one = page.replace('<link rel="stylesheet" href="%sstyles/g1s-th232.css"><link rel="stylesheet" href="%sstyles/flight.css">%s@@DATA@@' % (_sp, _sp, _lang_set),
                    _lang_set + css + "<script>\n" + rd(_data_name) + "\n</script>").replace("@@LANGHREF@@", _href_glued)
out_html = "spectrum-dad-sgn-en.html" if LANG == "en" else "spectrum-dad-sgn.html"
open(os.path.join(out_dir, out_html), "w", encoding="utf-8", newline="\n").write(one)
print("ok", len(page) // 1024, len(one) // 1024, "KB")
