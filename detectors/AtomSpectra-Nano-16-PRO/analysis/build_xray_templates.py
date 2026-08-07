# -*- coding: utf-8 -*-
"""Синтетические шаблоны характеристического рентгена (иод и вольфрам).

В цепочке Th-232 таких линий нет — их вносит взаимодействие γ с материалом
кристалла (K-серия иода в CsI, 28-33 кэВ) и с материалом самого источника
(K-серия вольфрама электрода, 58-69 кэВ). Оба процесса рождают квант с
энергией, характерной для оболочки атома-мишени, и появляются в спектре как
самостоятельные пики. В шаблонах Geant4 они частично выходят и без отдельной
обработки (Livermore-физика знает флуоресценцию иода), но затухают на грани
чувствительности и в подгонке уходят в подложку.

Здесь шаблоны строятся из каталога линий конструктора ROI: одна δ-линия для XI
(Kα1 иода) и две δ-линии для XW (Kα1 + Kβ2 вольфрама) с долями, взятыми из
конструктора. Сумма counts нормирована на 1 — «одна испущенная эмиссия»; тогда
коэффициент подгонки в `wt20_unfold` — темп эмиссий, эмиссий/с. Это не
активность цепочки: XI темпом эмиссии описывает поток γ на кристалл, XW — поток
γ на пачку. Соотнесение с активностью — отдельная задача.

Штамп src_sha1 и масса пачки копируются из существующего шаблона Ac228, чтобы
подгонка не спотыкалась о разные ревизии геометрии.

    python analysis/build_xray_templates.py [<каталог шаблонов>] [<xml>]
"""
import csv
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_DET = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)

import roi_lines as R  # noqa: E402

DEFAULT_TDIR = "C:/g4work/asn16/build/wt20_templates"
LIB = os.path.join(_DET, "reference", "nuclide-lines")


def branch_to_tl208():
    """Доля распадов Bi-212 на Tl-208 — из библиотеки МАГАТЭ, не числом."""
    p = os.path.join(LIB, "212bi_gammas.csv")
    for r in csv.DictReader(io.open(p, encoding="utf-8")):
        d, pc = r.get("decay"), r.get("decay_%")
        if d and d.strip().upper().startswith("A") and pc:
            try:
                return float(pc) / 100.0
            except ValueError:
                pass
    raise SystemExit("в %s нет альфа-ветви Bi-212" % p)


def read_head(path):
    head = []
    for ln in io.open(path, encoding="utf-8"):
        if not ln.startswith("#"):
            break
        head.append(ln.rstrip("\n"))
    return head


def write_template(path, head_source, lines, particle):
    """lines: [(E_keV, weight)], weight должно суммироваться в 1."""
    tot = sum(w for _, w in lines)
    lines = [(e, w / tot) for e, w in lines]
    head = read_head(head_source)
    # правим N_primaries: единица эмиссии; particle — на своё имя, чтобы шапка
    # не врала; чистим/переписываем поля, которые для линейчатого шаблона
    # теряют смысл
    out = []
    seen = set()
    for ln in head:
        if ln.startswith("# N_primaries"):
            out.append("# N_primaries = 1")
            seen.add("N_primaries")
        elif ln.startswith("# run_args"):
            out.append("# run_args = synthetic_xray_lines")
        elif ln.startswith("# particle"):
            out.append("# particle = " + particle)
            seen.add("particle")
        elif ln.startswith("# E_prim_keV"):
            e_str = ", ".join("%.2f" % e for e, _ in lines)
            out.append("# E_prim_keV = " + e_str)
        elif ln.startswith("# N_with_signal"):
            continue
        elif ln.startswith("# note"):
            continue
        else:
            out.append(ln)
    if "N_primaries" not in seen:
        out.append("# N_primaries = 1")
    if "particle" not in seen:
        out.append("# particle = " + particle)
    out.append("# note = synthetic X-ray template, lines and weights from "
               "reference/roi/wizard_lines_iaea.xml")
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        for ln in out:
            f.write(ln + "\n")
        f.write("E_keV,counts\n")
        for e, w in lines:
            f.write("%.3f,%.6g\n" % (e, w))


def main():
    tdir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_TDIR
    xml = sys.argv[2] if len(sys.argv) > 2 else R.DEFAULT_XML
    if not os.path.exists(tdir):
        raise SystemExit("нет каталога шаблонов: " + tdir)
    src_head = os.path.join(tdir, "Ac228.csv")
    if not os.path.exists(src_head):
        raise SystemExit("нет шаблона Ac228.csv в " + tdir)

    lines = R.parse_xml(xml)
    xi = [(r["E"], r["yield_pct"]) for r in lines if r["key"] == "XI"]
    xw = [(r["E"], r["yield_pct"]) for r in lines if r["key"] == "XW"]
    if not xi:
        raise SystemExit("в каталоге нет XI-линий")
    if not xw:
        raise SystemExit("в каталоге нет XW-линий")

    # K-серия дочерних из внутренней конверсии: в каталоге это строки с « X K»
    # в имени. Выходы — на распад РОДИТЕЛЯ строки; линии Tl-208 приводятся к
    # распаду подцепочки A2 множителем ветвления Bi-212.
    br = branch_to_tl208()
    xray_rows = [r for r in lines if " X K" in r["name"]]
    xd1 = [(r["E"], r["yield_pct"]) for r in xray_rows if r["key"] == "Ac228"]
    xd2 = ([(r["E"], r["yield_pct"]) for r in xray_rows if r["key"] == "Pb212"]
           + [(r["E"], r["yield_pct"] * br) for r in xray_rows
              if r["key"] == "Tl208"])

    for name, data in (("XI", xi), ("XW", xw), ("XD1", xd1), ("XD2", xd2)):
        if not data:
            print("%s: строк в каталоге нет — шаблон не построен" % name)
            continue
        p = os.path.join(tdir, "%s.csv" % name)
        write_template(p, src_head, data, name)
        print("%s.csv: %s" % (name,
              ", ".join("%.2f (вес %.4g)" % (e, y) for e, y in data)))


if __name__ == "__main__":
    main()
