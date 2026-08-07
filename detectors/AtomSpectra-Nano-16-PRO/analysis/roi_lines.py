
# -*- coding: utf-8 -*-
"""Каталог линий из веб-конструктора ROI (Becqmoni ROI wizard).

Источник: https://vibeengineering-llc.github.io/becqmoni-roi-wizard/. Оператор
подбирает линии по спектру, конструктор экспортирует их в XML формата BecqMoni
ROIConfigData. Здесь этот XML читается и переводится в единый список dict-ов, с
которого дальше работают `wt20_lines.py` (оконный съём), `wt20_unfold.py`
(синтетические шаблоны рентгена) и `wt20_report.py` (таблица линий в отчёте).

Каждой линии в конструкторе присвоено имя вида «Ac-228 (Th-232) (129–141)», где
скобка с диапазоном — окно в кэВ, назначенное вручную по форме измеренного пика.
Если окна в имени нет, берётся окно ±1 ПШПВ вокруг энергии пика (ПШПВ по
подгонке в `wt20_unfold`).

Формат XML: <ROIDefinitionData><Name>…</Name><PeakEnergy>…</PeakEnergy>
<Intencity>…</Intencity><HalfLife>годы</HalfLife>…</ROIDefinitionData>. Поле
Intencity — выход в процентах на один распад для γ-линий нуклидов; для линий
характеристического рентгена (XI, XW) там нормированные интенсивности серии, не
выход на распад цепочки.

    python analysis/roi_lines.py <xml> [csv]
"""
import io
import os
import re
import sys
import xml.etree.ElementTree as ET

_HERE = os.path.dirname(os.path.abspath(__file__))
_DET = os.path.dirname(_HERE)

# Ключ нуклида — совпадает с именем файла шаблона в wt20_templates/. XI и XW —
# синтетические компоненты для характеристического рентгена (иод в кристалле,
# вольфрам в самом источнике); шаблоны для них строятся отдельно
# (`build_xray_templates.py`), в цепочке ENSDF их нет.
NUC_MATCH = [
    ("Ac228", re.compile(r"^Ac-228\b")),
    ("Bi212", re.compile(r"^Bi-212\b")),
    ("Pb212", re.compile(r"^Pb-212\b")),
    ("Ra224", re.compile(r"^Ra-224\b")),
    ("Rn220", re.compile(r"^Rn-220\b")),
    ("Th228", re.compile(r"^Th-228\b")),
    ("Th232", re.compile(r"^Th-232\b")),
    ("Tl208", re.compile(r"^Tl-208\b")),
    ("XI",    re.compile(r"^ХРИ\s*I\b")),
    ("XW",    re.compile(r"^ХРИ\s*W\b")),
]

# Окно в имени: (число–число), допускаем и en-dash, и hyphen-minus.
WIN_RE = re.compile(r"\(\s*(\d+(?:[.,]\d+)?)\s*[–—\-]\s*"
                    r"(\d+(?:[.,]\d+)?)\s*\)")

DEFAULT_XML = os.path.join(_DET, "reference", "roi", "wizard_lines_iaea.xml")


def parse_xml(path):
    """Читает XML, отдаёт список dict-ов с полями:
       key, name, E, yield_pct, half_life_y, window (или None).
    """
    root = ET.parse(path).getroot()
    out = []
    for r in root.findall(".//ROIDefinitionData"):
        if not (r.findtext("Enabled", "") or "").lower().startswith("true"):
            continue
        name = (r.findtext("Name", "") or "").strip()
        e = float(r.findtext("PeakEnergy", "0"))
        y = float(r.findtext("Intencity", "0"))
        t_half_y = float(r.findtext("HalfLife", "0"))
        key = None
        for k, pat in NUC_MATCH:
            if pat.search(name):
                key = k
                break
        m = WIN_RE.search(name)
        win = None
        if m:
            a = float(m.group(1).replace(",", "."))
            b = float(m.group(2).replace(",", "."))
            win = (a, b)
        out.append(dict(key=key, name=name, E=e, yield_pct=y,
                        half_life_y=t_half_y, window=win))
    out.sort(key=lambda r: r["E"])
    return out


def by_nuclide(lines):
    """Группирует по key нуклида -> список линий."""
    out = {}
    for r in lines:
        out.setdefault(r["key"], []).append(r)
    return out


def to_csv(lines, path):
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        f.write("# каталог линий из веб-конструктора ROI (Becqmoni)\n")
        f.write("# https://vibeengineering-llc.github.io/becqmoni-roi-wizard/\n")
        f.write("нуклид;E_кэВ;выход_%;окно_низ_кэВ;окно_верх_кэВ;T_1_2_лет;имя\n")
        for r in lines:
            lo, hi = (r["window"] or ("", ""))
            f.write("%s;%.3f;%.4g;%s;%s;%.6g;%s\n" %
                    (r["key"] or "-", r["E"], r["yield_pct"],
                     "%.2f" % lo if lo != "" else "",
                     "%.2f" % hi if hi != "" else "",
                     r["half_life_y"], r["name"]))


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_XML
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        _DET, "reference", "roi", "wizard_lines.csv")
    lines = parse_xml(src)
    to_csv(lines, out)
    print("линий: %d (из %d уникальных нуклидов)" %
          (len(lines), len(set(r["key"] for r in lines))))
    for k, rs in by_nuclide(lines).items():
        print("  %-6s %2d линий" % (k or "-", len(rs)))
    print("записано:", out)


if __name__ == "__main__":
    main()
