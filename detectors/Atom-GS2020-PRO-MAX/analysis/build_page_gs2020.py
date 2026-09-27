# -*- coding: utf-8 -*-
r"""Сборка страницы GS2020 Th-232 донорским build_page.py (Gamma-1S/web-th232) без правок: подменяются только пути.
CSS — побайтно донорский; в JS заменены только видимые подписи «паспорт» (список TERMS, число вхождений сверяется). Запуск: python build_page_gs2020.py"""
import hashlib, os, shutil, sys
sys.stdout.reconfigure(encoding="utf-8")
DONOR = r"D:\Claude_files\repos\geant4-detector-models\detectors\Gamma-1S\web-th232"
PAGE = r"D:\GoogleDrive\Рабочая папка ИИ\GEANT4\web\gs2020-th232-page"
sys.path.insert(0, DONOR)
import build_page as bp
# #CHART-1 (оператор 26.09 «сделай увеличение по выделению мышью»): донор не имел drag-zoom на графиках
# метода 1/2 (только на вкладке "калибровка") — пробел донора, не дефект переноса (сверено grep-ом по
# mousedown в g1s-th232.js). Патч сгенерирован ступенью 2 (Ollama qwen3.6:27b, SPEC-gs2020-zoom-patch.md).
ZOOM_JS = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "patches", "gs2020_zoom.js"), encoding="utf-8").read().rstrip("\n")
TERMS = {"g1s-th232.css": [("text-align:justify; text-wrap:pretty}", "text-align:left; text-wrap:pretty}", 2)],  # #GS-10
         "g1s-th232.js": [('cell("против паспорта"', 'cell("к известной активности"', 2),
                          ("<th class='num'>к паспорту</th>", "<th class='num'>к известной активности</th>", 1),
                          ('lab: "паспорт"', 'lab: "известная"', 1), ('row("cmp-pass", "паспорт"', 'row("cmp-pass", "известная активность"', 1),
                          ('" паспорта, "', '" известной, "', 2),
                          ('" сумм-пиков + K-рентген")', '" сумм-пиков (суммы 860+2614 = 3475 кэВ в библиотеке нет); рентген K/L учтён в числе линий (#XR-1)")', 1),
                          ('"<td>по цезию комплекта " + num(fw.fwhm662_cs, 1) + " кэВ</td></tr>"', '"<td>таблица пиков прибора</td></tr>"', 1),
                          ('"корневой закон по записи цезия"', '"корневой закон через ту же точку 662 кэВ"', 1),
                          ('"ПШПВ по модели (цезий)"', '"фон: Маринелли+вода, промежуточный замер"', 1),
                          ('"ПШПВ по линиям спектра"', '"фон как снят (без сосуда)"', 1),
                          ('"все известные линии"', '"ENSDF, порог в полпроцента"', 1),
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
DEV = {"dev_res_mfr": "менее 7,5 % (паспорт; вариант с CsI(Tl) — менее 7,0 %)"}
_make_fill = bp.make_fill
def _fill_dev(d):
    f = _make_fill(d)
    f.update({k: (lambda v=v: v) for k, v in DEV.items()})
    return f
bp.make_fill = _fill_dev
bp.SRC, bp.DIST = os.path.join(PAGE, "src"), os.path.join(PAGE, "dist")
bp.DATA_JSON, bp.SINGLE = os.path.join(PAGE, "gs2020_th232_data.json"), os.path.join(PAGE, "gs2020_th232.html")
bp.main()
# #GS-7: рендеры раздела «Спектрометр» (подготовлены из model/…png, JPG ≤1000 px) — рядом со страницей в dist/img
os.makedirs(os.path.join(bp.DIST, "img"), exist_ok=True)
for fn in sorted(os.listdir(os.path.join(PAGE, "img"))):
    shutil.copy2(os.path.join(PAGE, "img", fn), os.path.join(bp.DIST, "img", fn))
    print("картинка dist/img/%s %d КБ" % (fn, os.path.getsize(os.path.join(PAGE, "img", fn)) // 1024))
