"""Отчёт в markdown из текста страницы: text.md + таблицы (поле, полосы, линии). Использование: build_report.py <text.md> <каталог results/final> <out.md>"""
import sys, csv, os
t, d, out = sys.argv[1], sys.argv[2], sys.argv[3]
rd = lambda p: open(p, encoding="utf-8").read()
head = ("# Спектр на эшелоне: рейс DAD–SGN (отчёт)\n\nАтомНано 16 (CsI 16,2 см³), Airbus A320/A321, 20.09.2026, высота 9800 м, место у окна. "
        "Модельный спектр космического излучения 20 кэВ – 10 МэВ: Geant4 11.4.2 + PARMA/EXPACS 4.10. Интерактивная версия с графиком — `web/flight-dad-sgn/index.html`.\n")
rows = [r for r in csv.DictReader(open(os.path.join(d, "prod3fin_lines.csv"), encoding="utf-8")) if r["found"] == "1"]
lines = "| E, кэВ | реакция | площадь ± σ, отсч./ч | значимость |\n|---|---|---|---|\n" + "".join(
    f"| {float(r['E_keV']):.1f} | {r['reaction']} | {float(r['net_cph']):.0f} ± {float(r['sigma_cph']):.0f} | {float(r['signif']):.1f} |\n" for r in rows)
sub = {"HEAD": head, "CHART": "*(график — в интерактивной версии)*", "SCHEME": "*(схема сечения — в интерактивной версии)*",
       "RENDER": "*(рендер геометрии прибора — в интерактивной версии)*",
       "FIELD": rd(os.path.join(d, "field.md")), "BANDS": rd(os.path.join(d, "bands.md")), "LINES": lines}
s = rd(t)
sub["EFFFIG"] = "*(график кривых — в интерактивной версии; значения — в таблице ниже)*"
for k, f in (("EFF", "eff.md"), ("FUEL", "fuel.md"), ("DOSE", "dose.md")):
    sub[k] = rd(os.path.join(d, f)) if os.path.exists(os.path.join(d, f)) else ""
nlf = "research/G-neutron-lines-summary.md"
sub["NEUTRONLINES"] = "\n".join(rd(nlf).split("\n")[4:]) if os.path.exists(nlf) else ""
for k, v in sub.items():
    s = s.replace("@@%s@@" % k, v)
open(out, "w", encoding="utf-8", newline="\n").write(s)
print("отчёт:", out, len(s), "символов; осталось меток:", s.count("@@"))
