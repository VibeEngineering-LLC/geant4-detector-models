# -*- coding: utf-8 -*-
r"""Сборка страницы GS2020 Th-232 донорским build_page.py (Gamma-1S/web-th232) без правок: подменяются только пути.
CSS — побайтно донорский; в JS заменены только видимые подписи «паспорт» (список TERMS, число вхождений сверяется). Запуск: python build_page_gs2020.py"""
import hashlib, os, shutil, sys
sys.stdout.reconfigure(encoding="utf-8")
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DONOR = os.path.join(REPO_ROOT, "detectors", "Gamma-1S", "web-th232")
PAGE = os.environ.get("GS2020_PAGE_WORK", os.path.join(REPO_ROOT, ".work", "gs2020-th232-page"))
sys.path.insert(0, DONOR)
import build_page as bp
TERMS = {"g1s-th232.js": [('cell("против паспорта"', 'cell("к известной активности"', 2),
                          ("<th class='num'>к паспорту</th>", "<th class='num'>к известной активности</th>", 1),
                          ('lab: "паспорт"', 'lab: "известная"', 1), ('row("cmp-pass", "паспорт"', 'row("cmp-pass", "известная активность"', 1),
                          ('" паспорта, "', '" известной, "', 2),
                          ('" сумм-пиков + K-рентген")', '" сумм-пиков (суммы 860+2614 = 3475 кэВ в библиотеке нет); рентген K/L учтён в числе линий (#XR-1)")', 1),
                          ('"<td>по цезию комплекта " + num(fw.fwhm662_cs, 1) + " кэВ</td></tr>"', '"<td>таблица пиков прибора</td></tr>"', 1),
                          ('"корневой закон по записи цезия"', '"корневой закон через ту же точку 662 кэВ"', 1),
                          ('"ПШПВ по модели (цезий)"', '"фон, ослабленный сосудом (K-40)"', 1),
                          ('"ПШПВ по линиям спектра"', '"фон как снят (без сосуда)"', 1),
                          ('"все известные линии"', '"ENSDF, порог в полпроцента"', 1),
                          ('"отобранная библиотека"', '"донорская библиотека"', 1),
                          ('"<tr><td>коэффициенты</td><td>"', '"<tr><td>коэффициенты своей шкалы (по реперам этого спектра)</td><td>"', 1),
                          # #PAGE-3 (оператор 25.09 «убери просвечивающие линии»): контур слоя — сразу после ЕГО заливки,
                          # следующий (меньший) слой его перекрывает; отдельный проход контуров поверх всех заливок снят.
                          ("      fillRuns(vec, [[0, e.length - 1]], order[oi].nuc.color);\n    }",
                           "      fillRuns(vec, [[0, e.length - 1]], order[oi].nuc.color);\n      strokeLayer(oi);\n    }", 1),
                          ("    for (var oj = 0; oj < order.length; oj++) {\n      var ks = order[oj].nuc.key;",
                           "    function strokeLayer(oj) {\n      var ks = order[oj].nuc.key;", 1)]}
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
for rel in (("styles", "g1s-th232.css"), ("scripts", "g1s-th232.js")):
    src, dst = os.path.join(DONOR, "src", *rel), os.path.join(PAGE, "src", *rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    txt = open(src, encoding="utf-8").read()
    for old, new, n in TERMS.get(rel[1], ()):  # оператор 25.09: не «паспорт», а «образец с известной активностью»
        if txt.count(old) != n:
            raise SystemExit("ОТКАЗ: в донорском %s «%s» встречается %d раз, ожидалось %d" % (rel[1], old, txt.count(old), n))
        txt = txt.replace(old, new)
    open(dst, "w", encoding="utf-8", newline="").write(txt)
    print("донорский %s sha256 %s, замен терминов %d" % ("/".join(rel), sha(src)[:16], len(TERMS.get(rel[1], ()))))
bp.SRC, bp.DIST = os.path.join(PAGE, "src"), os.path.join(PAGE, "dist")
bp.DATA_JSON, bp.SINGLE = os.path.join(PAGE, "gs2020_th232_data.json"), os.path.join(PAGE, "gs2020_th232.html")
bp.main()
