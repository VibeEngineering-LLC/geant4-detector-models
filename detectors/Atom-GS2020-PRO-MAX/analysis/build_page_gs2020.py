# -*- coding: utf-8 -*-
r"""Сборка страницы GS2020 Th-232 донорским build_page.py (Gamma-1S/web-th232) без правок: подменяются только пути.
CSS — побайтно донорский; в JS заменены только видимые подписи «паспорт» (список TERMS, число вхождений сверяется). Запуск: python build_page_gs2020.py"""
import hashlib, json, os, shutil, sys
sys.stdout.reconfigure(encoding="utf-8")
DONOR = r"D:\repos-folder\repos\geant4-detector-models\detectors\Gamma-1S\web-th232"
PAGE = r"D:\cloud-folder\work-folder\GEANT4\web\gs2020-th232-page"
sys.path.insert(0, DONOR)
# #PUB-1 / #GS-39 (29.09): вкладка K-40 показала числа прежней шкалы — gs2020_sources.json не перевыгрузили после
# подгонок. Гейт: JSON страницы новее каждой подгонки, из которой он собран, иначе ОТКАЗ со списком устаревших.
_TH, _KC = r"C:\g4work\gs2020\run_marinelli\out_v5_oisn10", r"C:\g4work\gs2020\kcl_v4w85_83"   # KCl: Маринелли 2 v4, конус колодца (29.09)
_FRESH = {"gs2020_th232_data.json": [_TH], "gs2020_k40_data.json": [_KC], "gs2020_sources.json": [_KC]}
_stale = []
for _pj, _dirs in _FRESH.items():
    _src = [os.path.join(d, f) for d in _dirs for f in os.listdir(d) if f.startswith("fit_") and f.endswith(".json")]
    _new = [s for s in _src if os.path.getmtime(s) > os.path.getmtime(os.path.join(PAGE, _pj))]
    if _new: _stale.append("%s старше %s" % (_pj, ", ".join(os.path.basename(s) for s in _new)))
if _stale: raise SystemExit("ОТКАЗ #PUB-1: перевыгрузите данные страницы — " + "; ".join(_stale))
import build_page as bp
# #CHART-1 (оператор 26.09 «сделай увеличение по выделению мышью»): донор не имел drag-zoom на графиках
# метода 1/2 (только на вкладке "калибровка") — пробел донора, не дефект переноса (сверено grep-ом по
# mousedown в g1s-th232.js). Патч сгенерирован ступенью 2 (Ollama qwen3.6:27b, SPEC-gs2020-zoom-patch.md).
ZOOM_JS = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "patches", "gs2020_zoom.js"), encoding="utf-8").read().rstrip("\n")
TERMS = {"g1s-th232.css": [("text-align:justify; text-wrap:pretty}", "text-align:left; text-wrap:pretty}", 2)],  # #GS-10
         "g1s-th232.js": [("    D.nuclides.forEach(function (nuc) {\n      // SECOND (вторичные пики)",
                           "    // #GS-45 (30.09, оператор «да»): легенда по вкладу слоёв (интеграл стека), как порядок отрисовки на графике\n"
                           "    var stkL = (elId === 'legendM2') ? STACK2() : STACK1();\n"
                           "    var arrL = D.nuclides.map(function (n, i) { var v = stkL && stkL[n.key], s = 0;\n"
                           "      if (v) for (var q = 0; q < v.length; q++) s += v[q]; return { n: n, s: s, i: i }; });\n"
                           "    arrL.sort(function (a, b) { return (b.s - a.s) || (a.i - b.i); });\n"
                           "    arrL.map(function (o) { return o.n; }).forEach(function (nuc) {\n      // SECOND (вторичные пики)", 1),
                          ('cell("против паспорта"', 'cell("к известной активности"', 2),
                          ("<th class='num'>к паспорту</th>", "<th class='num'>к известной активности</th>", 1),
                          ('lab: "паспорт"', 'lab: "известная"', 1), ('row("cmp-pass", "паспорт"', 'row("cmp-pass", "известная активность"', 1),
                          ('" паспорта, "', '" известной, "', 2),
                          ('" сумм-пиков + K-рентген")', '" сумм-пиков (суммы 860+2614 = 3475 кэВ в библиотеке нет); рентген K/L учтён в числе линий (#XR-1)")', 1),
                          ('"<td>по цезию комплекта " + num(fw.fwhm662_cs, 1) + " кэВ</td></tr>"', '"<td>таблица пиков прибора</td></tr>"', 1),
                          ('"корневой закон по записи цезия"', '"корневой закон через ту же точку 662 кэВ"', 1),
                          ('"ПШПВ по модели (цезий)"', '"фон: Маринелли+вода"', 1),
                          ('"ПШПВ по линиям спектра"', '"фон как снят (без сосуда)"', 1),
                          ('"все известные линии"', '"метод 2-2: все известные линии ENSDF"', 1),
                          ('"отобранная библиотека"', '"донорская библиотека"', 1),
                          ('"<tr><td>коэффициенты</td><td>"', '"<tr><td>коэффициенты своей шкалы (по реперам этого спектра)</td><td>"', 1),
                          # оператор 26.09 «фон без сосуда убери»: кнопка "lines" снята из разметки (index.html),
                          # дефолт состояния синхронизирован на единственный оставшийся вариант "cs".
                          ('fwhmLaw: "lines", // закон ширины линии: lines (по спектру) | cs (цезий)',
                           'fwhmLaw: "cs", // #PAGE-фон-без-сосуда (26.09): единственный вариант — фон Маринелли+вода', 1),
                          # #PAGE-3 (оператор 25.09 «убери просвечивающие линии»): контур слоя — сразу после ЕГО заливки,
                          # следующий (меньший) слой его перекрывает; отдельный проход контуров поверх всех заливок снят.
                          ("      fillRuns(vec, [[0, e.length - 1]], order[oi].nuc.color);\n    }",
                           "      fillRuns(vec, [[0, e.length - 1]], order[oi].nuc.color);\n      strokeLayer(oi);\n    }", 1),
                          ("    for (var oj = 0; oj < order.length; oj++) {\n      var ks = order[oj].nuc.key;",
                           "    function strokeLayer(oj) {\n      var ks = order[oj].nuc.key;", 1),
                          # #GS-8 (оператор 27.09 «часть шаблонов уехали за график»): заливки слоёв (fillRuns/
                          # segPath) рисуются по ПОЛНОМУ диапазону каналов без обрезки по xLo/xHi — при xLo=0
                          # (донорский дефолт) это не было заметно (e[0]~0 совпадало с левым краем), но #GS-5
                          # ниже подвинул xLo на 25 кэВ, и каналы e<25 стали рисоваться ЛЕВЕЕ границы графика.
                          # Обрезка (clip) области построения — до заливок и до финального strokeRect в конце.
                          ('    g.fillText("энергия, кэВ", (m.l + W - m.r) / 2, H - 2);\n\n    g.strokeStyle = p.rule; g.lineWidth = 2;\n    g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);\n\n    // Заливки',
                           '    g.fillText("энергия, кэВ", (m.l + W - m.r) / 2, H - 2);\n\n    g.strokeStyle = p.rule; g.lineWidth = 2;\n    g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);\n    g.save(); g.beginPath(); g.rect(m.l, m.t, W - m.r - m.l, H - m.b - m.t); g.clip();\n\n    // Заливки', 1),
                          # #CHART-1: zoom-aware диапазон отображения (drawSpectrum) вместо жёстких 0..e[last]
                          # #GS-5 (оператор 27.09 «нижний порог отображения сделай от 25 кэВ»): дефолт (без zoom) —
                          # 25 кэВ, порог прибора, а не 0 (окно ПОДГОНКИ 150-3600 кэВ #SUM-1 не меняется, это чисто вид).
                          ("    var xLo = 0, xHi = e[e.length - 1];",
                           "    var xLo = ST.zoom ? ST.zoom.xLo : 25, xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];", 1),
                          # #CHART-1: курсор-подсказка (attachCursor) — тоже с учётом zoom
                          ("      var e = D.spectrum.e_of_ch;\n      var xHi = e[e.length - 1];\n      if (x < m.l || x > r.width - m.r) ST.cursorE = null;\n      else ST.cursorE = ((x - m.l) / (r.width - m.r - m.l)) * xHi;",
                           "      var e = D.spectrum.e_of_ch;\n      var xLo = ST.zoom ? ST.zoom.xLo : 25, xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];\n      if (x < m.l || x > r.width - m.r) ST.cursorE = null;\n      else ST.cursorE = xLo + ((x - m.l) / (r.width - m.r - m.l)) * (xHi - xLo);", 1),
                          # #GS-4 (оператор 27.09 «не подсвечивается зона выбора мышью»): рамка выделения во время
                          # протяжки — эталон донора CAL.drag (ra226.js:791-800), портировано на ST.drag.
                          ("    if (ST.cursorE !== null) {\n      var xC = mapX(ST.cursorE, xLo, xHi, m.l, W - m.r);\n      g.strokeStyle = p.rule; g.lineWidth = 1; g.setLineDash([4, 4]);\n      g.beginPath(); g.moveTo(xC, m.t); g.lineTo(xC, H - m.b); g.stroke();\n      g.setLineDash([]);\n    }\n    g.strokeStyle = p.rule; g.lineWidth = 2;",
                           "    if (ST.cursorE !== null) {\n      var xC = mapX(ST.cursorE, xLo, xHi, m.l, W - m.r);\n      g.strokeStyle = p.rule; g.lineWidth = 1; g.setLineDash([4, 4]);\n      g.beginPath(); g.moveTo(xC, m.t); g.lineTo(xC, H - m.b); g.stroke();\n      g.setLineDash([]);\n    }\n    if (ST.drag) {\n      var xa3 = mapX(Math.min(ST.drag.e0, ST.drag.e1), xLo, xHi, m.l, W - m.r);\n      var xb3 = mapX(Math.max(ST.drag.e0, ST.drag.e1), xLo, xHi, m.l, W - m.r);\n      g.fillStyle = \"rgba(246,211,28,.22)\";\n      g.fillRect(xa3, m.t, xb3 - xa3, H - m.b - m.t);\n      g.strokeStyle = \"#16140f\"; g.lineWidth = 1.5; g.setLineDash([4, 4]);\n      g.strokeRect(xa3, m.t, xb3 - xa3, H - m.b - m.t);\n      g.setLineDash([]);\n    }\n    g.restore();\n    g.strokeStyle = p.rule; g.lineWidth = 2;", 1),
                          # #CHART-1: сами функции zEfromX/wireZoom — вставлены перед комментарием "перерисовка"
                          ("  /* ── перерисовка активной вкладки ───────────────────────────── */",
                           ZOOM_JS + "\n\n  /* ── перерисовка активной вкладки ───────────────────────────── */", 1),
                          # #GS-19 (оператор 27.09 «тут нужно добавить метод 2-2»; «2-2 это все известные линии»): в сравнении
                          # метод 2 (донорская библиотека) и метод 2-2 (все линии ENSDF) — отдельными строками, не через переключатель.
                          ('      { lab: "метод 2", A: M2().A_Bq, dA: M2().dA_Bq, col: "#c8541c" }\n    ];',
                           '      { lab: "метод 2", A: SRC().method2.A_Bq, dA: SRC().method2.dA_Bq, col: "#c8541c" },\n'
                           '      { lab: "метод 2-2", A: SRC().method2_full.A_Bq, dA: SRC().method2_full.dA_Bq, col: "#7a3b12" }\n    ];', 1),
                          ('    var pass = D.passport, m1 = M1(), m2 = M2();',
                           '    var pass = D.passport, m1 = M1(), m2 = SRC().method2, m2f = SRC().method2_full;', 1),
                          ('          + libTxt + ", " + modeTxt)\n',
                           '          + "донорская библиотека, " + modeTxt)\n'
                           '    + row("cmp-m2", "метод 2-2: функция ПП + все известные линии ENSDF, " + cnt(m2f.n_lines) + " линий",\n'
                           '          m2f.A_Bq, m2f.dA_Bq,\n'
                           '          num(m2f.A_Bq / pass.A_Bq, 3) + " известной, " + signedPct(m2f.A_Bq / pass.A_Bq)\n'
                           '          + "; χ²/ν = " + num(m2f.chi2_ndof, 2) + " на " + cnt(m2f.n_channels_fit) + " каналах окон пиков; " + modeTxt)\n', 1),
                          ('метод 1 относительно метода 2 при "', 'метод 1 относительно метода 2 (донорская библиотека) при "', 1),
                          # #CHART-1: подключение wireZoom к обеим канвам (метод 1, метод 2) — синхронный zoom (ST.zoom общий)
                          ('    attachCursor("cvM1", "m1-tip", function () {',
                           '    wireZoom("cvM1"); wireZoom("cvM2");\n    attachCursor("cvM1", "m1-tip", function () {', 1)]}
# #GS-10 (оператор 27.09 «опять с шириной текста проблемы. исправь и запомни»; скилл web-publish §3): донорский CSS
# снимает предел строки (`max-width:none` у .stand/.ai-note/.method-lede) — абзацы шли во всю ширину. Предел ставит
# СБОРЩИК на ВСЕ абзацы и пункты страницы, а не автор раздела; проверка — MEASURE_CHECK ниже и замер в браузере.
# #GS-16 (оператор 27.09 «текст должен быть по ширине страницы»): #GS-10 было прочитано НЕВЕРНО (W-152) — оператор
# хотел текст на всю ширину, а не предел 78 знаков. Сборщик снимает любые пределы ширины у текстовых блоков.
MEASURE_CSS = ("\n/* #GS-16: текст на всю ширину страницы, ставится сборщиком */\n"
               "body,p,li,figcaption,caption,dd{text-align:left}\n"
               ".app p,.app li,.app figcaption,.pop p,.pop li,.pop dd{max-width:none !important; text-wrap:pretty}\n")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
for rel in (("styles", "g1s-th232.css"), ("scripts", "g1s-th232.js")):
    src, dst = os.path.join(DONOR, "src", *rel), os.path.join(PAGE, "src", *rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    txt = open(src, encoding="utf-8").read()
    for old, new, n in TERMS.get(rel[1], ()):  # оператор 25.09: не «паспорт», а «образец с известной активностью»
        if txt.count(old) != n:
            raise SystemExit("ОТКАЗ: в донорском %s «%s» встречается %d раз, ожидалось %d" % (rel[1], old, txt.count(old), n))
        txt = txt.replace(old, new)
    if rel[1].endswith(".css"):
        txt += MEASURE_CSS
        if "justify" in txt.replace("justify-content", "").replace("justify-items", "").replace("justify-self", ""):
            raise SystemExit("ОТКАЗ: в CSS есть выключка по ширине (web-publish §3)")
    open(dst, "w", encoding="utf-8", newline="").write(txt)
    print("донорский %s sha256 %s, замен терминов %d" % ("/".join(rel), sha(src)[:16], len(TERMS.get(rel[1], ()))))
# #GS-7: паспортные величины прибора (не расчёт) — источник README-референсы.md §8, сайт Gammaspectacular 27.09.
DEV = {"dev_res_mfr": "менее 7,5 % (паспорт; вариант с CsI(Tl) — менее 7,0 %)",
       # #GS-22: состав материала ОИСН-10, документация ЛСРМ (README-референсы, «Материал ОИСН-10»)
       "dev_oisn_fe": "49 %", "dev_oisn_comp": "Fe 49 %, C 31 %, O 12 %, H 4,4 %, N 3,6 %"}
# #GS-24: вкладки источников — данные export_sources_gs2020.py (метки {{k40_*}} + dist/sources.js), код src/scripts/gs-sources.js
SRC_JSON, SRC_APP = os.path.join(PAGE, "gs2020_sources.json"), os.path.join(PAGE, "src", "scripts", "gs-sources.js")
SRCD = json.load(open(SRC_JSON, encoding="utf-8")); DEV.update(SRCD.pop("fill"))
_make_fill = bp.make_fill
def _fill_dev(d):
    f = _make_fill(d)
    f.update({k: (lambda v=v: v) for k, v in DEV.items()})
    return f
bp.make_fill = _fill_dev
bp.SRC, bp.DIST = os.path.join(PAGE, "src"), os.path.join(PAGE, "dist")
bp.DATA_JSON, bp.SINGLE = os.path.join(PAGE, "gs2020_th232_data.json"), os.path.join(PAGE, "gs2020_th232.html")
bp.main()
src_js = "window.GS_SRC=" + json.dumps(SRCD, ensure_ascii=False, separators=(",", ":")) + ";\n"
app_js = open(SRC_APP, encoding="utf-8").read()
for nm, txt in (("sources.js", src_js), ("gs-sources.js", app_js)):
    if bp.secret_scan(txt) or "</script" in txt.lower(): raise SystemExit("ОТКАЗ: секрет-скан или тег script в " + nm)
v = int(max(os.path.getmtime(SRC_JSON), os.path.getmtime(SRC_APP)))
open(os.path.join(bp.DIST, "sources.js"), "w", encoding="utf-8").write(src_js)
open(os.path.join(bp.DIST, "scripts", "gs-sources.js"), "w", encoding="utf-8").write(app_js)
for p, rep in ((os.path.join(bp.DIST, "index.html"), '<script src="sources.js?v=%d"></script>\n<script src="scripts/gs-sources.js?v=%d"></script>' % (v, v)),
               (bp.SINGLE, "<script>" + src_js + "</script>\n<script>\n" + app_js + "</script>")):
    t = open(p, encoding="utf-8").read()
    if t.count("<!--@sources-->") != 1: raise SystemExit("ОТКАЗ: метка <!--@sources--> в %s встречается %d раз" % (p, t.count("<!--@sources-->")))
    open(p, "w", encoding="utf-8").write(t.replace("<!--@sources-->", rep))
# #GS-26 (вариант А): панель K-40 «как у тория» — тот же донорский JS вторым экземпляром (src-копия выше + patches/k40_terms.json),
# данные gs2020_k40_data.json (export_page_gs2020_k40.py), разметка src/k40-panel.html и src/k40-pops.html
import build_page_gs2020_k40 as k40
K40_FILES = k40.finish(bp)
# #GS-29 (оператор 28.09 «слово КИ для данного спектрометра нигде не пиши»): гейт по собранным файлам, слово целиком
import re
KI = re.compile(r"(?<![А-ЯЁа-яё])КИ(?![А-ЯЁа-яё])")
for fn in ("index.html", "sources.js", "data.js") + tuple(K40_FILES):
    if KI.search(open(os.path.join(bp.DIST, fn), encoding="utf-8").read()): raise SystemExit("ОТКАЗ #GS-29: слово «КИ» в dist/" + fn)
print("источники: dist/sources.js %d КБ, dist/scripts/gs-sources.js %d КБ, v=%d" % (len(src_js.encode()) // 1024, len(app_js.encode()) // 1024, v))
# #GS-7: рендеры раздела «Спектрометр» (подготовлены из model/…png, JPG ≤1000 px) — рядом со страницей в dist/img
os.makedirs(os.path.join(bp.DIST, "img"), exist_ok=True)
for fn in sorted(os.listdir(os.path.join(PAGE, "img"))):
    shutil.copy2(os.path.join(PAGE, "img", fn), os.path.join(bp.DIST, "img", fn))
    print("картинка dist/img/%s %d КБ" % (fn, os.path.getsize(os.path.join(PAGE, "img", fn)) // 1024))
