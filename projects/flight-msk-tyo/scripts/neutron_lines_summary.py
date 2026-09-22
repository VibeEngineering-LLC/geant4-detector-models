"""Сводка нейтрон-индуцированных линий (n,g) и (n,n') из B-gamma-lines-airframe.md со статусом по расчёту.
-> research/G-neutron-lines-summary.md. Использование: neutron_lines_summary.py [ru|en]"""
import csv, re, sys
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, "scripts")
from bline_table import parse_lines

LANG = sys.argv[1] if len(sys.argv) > 1 else "ru"
TABLE_MD = "research/B-gamma-lines-airframe.en.md" if LANG == "en" else "research/B-gamma-lines-airframe.md"
PFX = "en_" if LANG == "en" else ""
OUT_MD = "research/G-neutron-lines-summary.en.md" if LANG == "en" else "research/G-neutron-lines-summary.md"

rows = parse_lines(TABLE_MD)
NEUTRON_CAPTURE_RE = re.compile(r"\(n[,;]\s*[γg]\)")   # #REM-6: «Выход на 100 захватов» осмыслен ТОЛЬКО для
# (n,γ) — поле "yield" у (n,n′)-строк в источнике держит СЕЧЕНИЕ (мб при энергии), а не выход в процентах;
# без фильтра в колонку попадал текст вида «σ 868 мб @2 МэВ» — нечитаемо и по смыслу не то же самое число
neu = [r for r in rows if re.search(r"\(n[,;]\s*[γg]\)|\(n[,;]\s*n", r["reaction"]) and not re.search(r"GAGG|LaBr", r["place"])]
# GAGG/LaBr3 — справочные линии ДРУГИХ приборов из донорской таблицы; в АтомНано 16 (кристалл CsI) этих материалов нет
found = {round(float(r["E_keV"]), 2): r for r in csv.DictReader(open(f"results/final/{PFX}prod3fin_lines.csv", encoding="utf-8"))}
unl = [(float(r["E_keV"]), float(r["signif"])) for r in csv.DictReader(open(f"results/final/{PFX}prod3fin_unlisted.csv", encoding="utf-8"))]
raw = {round(float(r["E_keV"]), 2): float(r["signif"]) for r in csv.DictReader(open(f"results/final/{PFX}prod3raw_lines.csv", encoding="utf-8"))}
# #REM-4 (оператор, 22.09): значимость нужно проверять и в НЕЙТРОННОМ компоненте, а не только в сумме — реакция
# запущена нейтроном, и в сумме её может маскировать гамма-континуум других компонент (57,6/58,1 кэВ: -2,2σ в
# сумме, но 18,0σ в самом компоненте). Источник: results/final/prod3_total_neutron.csv -> raw_to_cph.py ->
# lines_report.py (тот же порог 4σ) -> results/final/prod3rawneu_lines.csv.
raw_neu = {round(float(r["E_keV"]), 2): float(r["signif"]) for r in csv.DictReader(open(f"results/final/{PFX}prod3rawneu_lines.csv", encoding="utf-8"))}
def nearest_peak(E, tol):   # ближайший значимый пик из unl в пределах tol, или None
    best, bd = None, tol
    for ue, us in unl:
        d = abs(ue - E)
        if d <= bd:
            bd, best = d, (ue, us)
    return best


ST = {"ru": {"found": "найдена", "seen_neu": "видна в компоненте нейтронов, не выделяется в сумме",
             "near": "рядом значимый пик, не отождествлён", "cand": "кандидат", "oow": "вне окна поиска"},
      "en": {"found": "found", "seen_neu": "seen in the neutron component, not resolved in the sum",
             "near": "significant peak nearby, unidentified", "cand": "candidate", "oow": "outside search window"}}[LANG]

def status(E):   # -> (статус, значимость в СУММЕ, (E,sigma) ближайшего пика, значимость в компоненте нейтронов)
    key = round(E, 2)
    tol = max(3, E * 0.01)
    sig_neu = raw_neu.get(key)
    if key in found and found[key]["found"] == "1":
        return ST["found"], float(found[key]["signif"]), None, sig_neu
    sig0 = raw.get(key)
    if sig_neu is not None and sig_neu >= 4:
        return ST["seen_neu"], sig0, None, sig_neu
    peak = nearest_peak(E, tol)
    if peak and (sig0 is None or sig0 < 4):
        return ST["near"], sig0, peak, sig_neu
    if sig0 is not None:
        return ST["cand"], sig0, None, sig_neu
    return ST["oow"], None, None, sig_neu

def elem_of(reaction):
    m = re.match(r"^[⁰-⁹¹²³\s]*([A-Za-z]+)", reaction)
    return m.group(1) if m else "?"

HEADER_RU = ["# G. Нейтрон-индуцированные линии (n,g) и (n,n') — сводка (22.09.2026, ред. 4)", "",
       "Источник реакций: `research/B-gamma-lines-airframe.md`, раздел 9. Значимость КАНДИДАТА считается строго в табличной энергии "
       "уровня (`scripts/lines_report.py`, окно 3 канала, порог 4σ) — это не то же самое, что «есть ли рядом реальный пик»: табличная "
       "энергия и центр наблюдаемого пика могут не совпадать на несколько кэВ. Приводятся ДВЕ значимости: в СУММАРНОМ сыром спектре "
       "(что реально нарисовано на графике страницы) и отдельно — в компоненте, порождённой нейтронами поля (`prod3_total_neutron.csv`), "
       "потому что гамма-континуум других компонент может маскировать линию в сумме, даже когда сама реакция статистически надёжна. "
       "«Найдена» — кандидат сам прошёл порог значимости в сумме (это и показано на графике страницы). «Видна в компоненте нейтронов, "
       "не выделяется в сумме» — в сумме порог не пройден, но в нейтронной компоненте — пройден с запасом (типичный случай — узкая "
       "линия на высоком континууме других частиц). «Рядом значимый пик, не отождествлён» — в табличной точке кандидат не значим ни в "
       "сумме, ни в компоненте, но независимый поиск пиков нашёл значимый максимум в сумме в пределах допуска; программа их не связала "
       "автоматически. «Кандидат» — ни в точке, ни рядом ничего значимого ни в сумме, ни в компоненте. «Вне окна поиска» — соседняя "
       "энергия занята другой линией. Колонка «Выход на 100 захватов» — вероятность того, что нейтрон, захваченный ИМЕННО ЭТИМ нуклидом, "
       "даст гамма-квант именно этой энергии (не вероятность самого захвата); источник — IAEA PGAA adopted database, "
       "`research/data/nglist_a.dat`/`promptgammas.xls` (см. `research/B-gamma-lines-airframe.md` строка 53). Для (n,n′) — прочерк: это "
       "неупругое рассеяние, захвата не происходит.", "",
       "| E, кэВ | Элемент | Реакция | Материал | Статус | Значимость в сумме, σ | Значимость в компоненте нейтронов, σ | Ближайший пик | Выход на 100 захватов |",
       "|---|---|---|---|---|---|---|---|---|"]
HEADER_EN = ["# G. Neutron-induced lines (n,g) and (n,n') — summary (22.09.2026, rev. 4)", "",
   "Reaction source: `research/B-gamma-lines-airframe.en.md`, section 9. CANDIDATE significance is computed strictly at the "
   "tabulated level energy (`scripts/lines_report.py`, 3-channel window, 4σ threshold) — not the same as \"is there a real peak "
   "nearby\": the tabulated energy and the centre of the observed peak may not coincide by a few keV. TWO significances are "
   "given: in the SUMMED raw spectrum (what is actually plotted on the page) and separately in the component produced by field "
   "neutrons (`prod3_total_neutron.csv`), because the gamma continuum of other components can mask the line in the sum even when "
   "the reaction itself is statistically robust.",
   "\"Found\" — the candidate passed the significance threshold in the sum (that is what is shown on the page). \"Seen in the "
   "neutron component, not resolved in the sum\" — the threshold is not passed in the sum but is passed with margin in the "
   "neutron component (typical case — a narrow line on a high continuum of other particles). \"Significant peak nearby, "
   "unidentified\" — not significant at the tabulated point in either the sum or the component, but an independent peak search "
   "found a significant maximum in the sum within tolerance; the program did not link them automatically. \"Candidate\" — "
   "nothing significant at the point or nearby, in either the sum or the component. \"Outside search window\" — the adjacent "
   "energy is occupied by another line.",
   "The \"Yield per 100 captures\" column is the probability that a neutron captured by THIS SPECIFIC nuclide emits a gamma ray "
   "of exactly this energy (not the probability of capture itself); source — the IAEA PGAA adopted database, "
   "`research/data/nglist_a.dat`/`promptgammas.xls` (see `research/B-gamma-lines-airframe.md` line 53). For (n,n′) — a dash: "
   "this is inelastic scattering, no capture occurs.", "",
   "| E, keV | Element | Reaction | Material | Status | Significance in the sum, σ | Significance in neutron component, σ | Nearest peak | Yield per 100 captures |",
   "|---|---|---|---|---|---|---|---|---|"]
out = HEADER_RU if LANG == "ru" else HEADER_EN
for r in sorted(neu, key=lambda x: x["E_keV"]):
    st, sg, peak, sgn = status(r["E_keV"])
    sgtxt = "%.1f" % sg if sg is not None else "—"
    sgntxt = "%.1f" % sgn if sgn is not None else "—"
    keV = "кэВ" if LANG == "ru" else "keV"
    ptxt = "%.1f %s, %.1fσ (Δ=%.1f)" % (peak[0], keV, peak[1], peak[0] - r["E_keV"]) if peak else "—"
    is_capture = bool(NEUTRON_CAPTURE_RE.search(r["reaction"]))
    yld = ((r.get("yield") or "").strip() or "—") if is_capture else "—"
    out.append(f"| {r['E_keV']:.2f} | {elem_of(r['reaction'])} | {r['reaction']} | {r['place']} | {st} | {sgtxt} | {sgntxt} | {ptxt} | {yld} |")
open(OUT_MD, "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")
print(len(neu), "rows" if LANG == "en" else "строк")
