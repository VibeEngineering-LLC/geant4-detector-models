# -*- coding: utf-8 -*-
"""Сверка эмиссии шаблона Ac-228 (Geant4 МК) с библиотекой ENSDF.

Первый пункт плана задач #170/#151 («сверить выходы линий Ac-228 с
ENSDF напрямую»): iso_Ac228_emit.csv — фотоны, реально испущенные
Geant4 при распаде (гистограмма 1 кэВ, N_primaries розыгрышей); файл
web-th232/data/ensdf_th232_chain_lines.csv — та же ветвь по данным
IAEA-NDS/ENSDF, использованная как библиотека метода 2 (провенанс —
построчно в колонке source самого файла, сверка ENSDF/LNHB уже сделана
при постройке библиотеки).

Оба источника бинуются в те же интервалы 1 кэВ (центр X.5 — так, как
записан iso_*_emit.csv), внутри бина суммируются близкие линии
(детектор их всё равно не разрешит). Строки library с line_type=xray
исключены: это атомная релаксация (K/L-рентген), считается отдельной
сущностью в пайплайне, а не частью ядерной гамма-эмиссии, и сравнивать
её с "гамма"-гистограммой Geant4 напрямую нельзя (см. R76/R69).

Запуск (нужны собранные шаблоны в build/):
    cd C:\\g4work\\gamma1s\\build
    python C:\\g4work\\geant4-detector-models\\detectors\\Gamma-1S\\analysis\\ac228_line_check.py

Не часть основного конвейера export_data.py — диагностика по требованию,
как phantom_check.py.
"""
import csv
import io
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LIB_CSV = os.path.join(HERE, "..", "web-th232", "data",
                        "ensdf_th232_chain_lines.csv")


def load_g4_emit(path):
    g4 = {}
    ntot = None
    with io.open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("# N_primaries"):
                ntot = float(ln.split("=")[1].strip())
            if ln.startswith("#") or not ln.strip() or ln.startswith("E_keV"):
                continue
            parts = ln.strip().split(",")
            if len(parts) < 2:
                continue
            e, n = float(parts[0]), float(parts[1])
            g4[e] = g4.get(e, 0.0) + n
    if ntot is None:
        raise SystemExit("в %s нет строки # N_primaries" % path)
    return g4, ntot


def load_lib_binned(path, nuclide):
    bins = {}
    with io.open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["nuclide"] != nuclide or row["line_type"] != "gamma":
                continue
            e = float(row["E_keV"])
            ip = float(row["I_percent"])
            b = math.floor(e) + 0.5
            bins.setdefault(b, []).append((e, ip))
    return bins


def compare(nuclide, emit_csv, tol_flag=0.30):
    """-> (rows, summary). rows отсортированы по убыванию I_lib.
    tol_flag — порог |ratio-1| для пометки строки как расходящейся,
    применяется только к бинам с I_lib >= 0.3 % (статистика МК на
    слабых бинах сама по себе даёт разброс такого порядка)."""
    g4, ntot = load_g4_emit(emit_csv)
    lib_bins = load_lib_binned(LIB_CSV, nuclide)

    rows = []
    for b, items in lib_bins.items():
        ip_sum = sum(ip for _, ip in items)
        g4pct = g4.get(b, 0.0) / ntot * 100.0
        ratio = (g4pct / ip_sum) if ip_sum > 0 else float("inf")
        label = "+".join("%.2f" % e for e, _ in items)
        flagged = ip_sum >= 0.3 and abs(ratio - 1.0) > tol_flag
        rows.append({"E_bin": b, "label": label, "I_lib_pct": ip_sum,
                      "I_g4_pct": g4pct, "ratio": ratio, "flagged": flagged})
    rows.sort(key=lambda r: -r["I_lib_pct"])

    lib_tot = sum(r["I_lib_pct"] for r in rows)
    g4_tot = sum(r["I_g4_pct"] for r in rows)
    summary = {
        "n_bins": len(rows),
        "lib_total_pct": lib_tot,
        "g4_total_pct": g4_tot,
        "overall_ratio": g4_tot / lib_tot if lib_tot else float("nan"),
        "n_flagged": sum(1 for r in rows if r["flagged"]),
        "flagged": [r for r in rows if r["flagged"]],
    }
    return rows, summary


def main():
    build = sys.argv[1] if len(sys.argv) > 1 else "."
    nuclide = sys.argv[2] if len(sys.argv) > 2 else "Ac228"
    emit_csv = os.path.join(build, "iso_%s_emit.csv" % nuclide)
    if not os.path.isfile(emit_csv):
        raise SystemExit("нет файла %s — запускать из каталога build "
                          "или передать его первым аргументом" % emit_csv)
    rows, summary = compare(nuclide, emit_csv)

    print("%-18s %8s %8s %8s" % ("E_lib(бин)", "I_lib%", "I_g4%", "ratio"))
    for r in rows:
        mark = "  <<<" if r["flagged"] else ""
        print("%-18s %8.3f %8.3f %8.3f%s" % (
            r["label"], r["I_lib_pct"], r["I_g4_pct"], r["ratio"], mark))
    print("---")
    print("бинов сверено: %d, из них расходятся (|ratio-1|>30%%, "
          "I_lib>=0.3%%): %d" % (summary["n_bins"], summary["n_flagged"]))
    print("сумма по всем сверенным бинам: lib=%.2f%% g4=%.2f%% "
          "(общее отношение %.3f)" % (
              summary["lib_total_pct"], summary["g4_total_pct"],
              summary["overall_ratio"]))
    if summary["flagged"]:
        print("\nрасходящиеся линии/дублеты (I_lib >= 0.3 %):")
        for r in summary["flagged"]:
            print("  %-18s I_lib=%.3f%% I_g4=%.3f%% ratio=%.3f" % (
                r["label"], r["I_lib_pct"], r["I_g4_pct"], r["ratio"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
