# -*- coding: utf-8 -*-
"""Манифест прогона: отпечаток каждого расчётного спектра, из которого построена
кривая.

Сами спектры в репозиторий не коммитятся (`.gitignore`: «сотни файлов, десятки
мегабайт»), поэтому из клона нельзя проверить, ТЕ ЛИ файлы дали опубликованную
кривую. Манифест закрывает эту дыру: на каждый узел пишется sha1 файла,
энергия, число первичных и число событий с откликом. Пересчитать кривую из
других спектров и выдать её за эту больше нельзя — расхождение видно построчно.

Это тот же класс дефекта, из-за которого отозвана сверка с записью Cs-137 и
заархивирована библиотечная кривая: вход, лежащий вне дерева, не прослеживается.

    python analysis/make_manifest.py <каталог со спектрами> [<выход.csv>]
"""
import glob
import hashlib
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(_HERE, "..", "results",
                                    "runs_manifest.csv"))


def head_of(path):
    """Шапка спектра -> dict. Читается только до строки данных."""
    out = {}
    for ln in io.open(path, encoding="utf-8"):
        if not ln.startswith("#"):
            break
        if "=" in ln:
            k, v = ln[1:].split("=", 1)
            out[k.strip()] = v.strip()
    return out


def sha1_of(path):
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    src = sys.argv[1]
    out = os.path.normpath(sys.argv[2]) if len(sys.argv) > 2 else OUT
    files = sorted(glob.glob(os.path.join(src, "*.csv")))
    if not files:
        raise SystemExit("в %s нет спектров" % src)

    rows, stamps, macs = [], set(), set()
    for f in files:
        h = head_of(f)
        stamps.add(h.get("src_sha1", "БЕЗ-ШТАМПА"))
        macs.add(h.get("run_args", "?"))
        rows.append((os.path.basename(f), sha1_of(f),
                     h.get("E_prim_keV", "?"), h.get("N_primaries", "?"),
                     h.get("N_with_signal", "?"),
                     h.get("solid_angle_frac", "?")))
    if len(stamps) > 1:
        raise SystemExit("спектры от РАЗНЫХ сборок: %s — манифест по ним "
                         "недействителен" % ", ".join(sorted(stamps)))
    if "БЕЗ-ШТАМПА" in stamps:
        raise SystemExit("спектры без штампа исходников: манифест не строится")

    stamp = stamps.pop()
    with io.open(out, "w", encoding="utf-8", newline="") as g:
        g.write("#@ stamp.version = 1\n")
        g.write("#@ src.script = detectors/AtomSpectra-Nano-16-PRO/analysis/"
                "make_manifest.py\n")
        g.write("#@ src.spectra_sha1 = %s\n" % stamp)
        g.write("#@ src.inputs_n = %d\n" % len(rows))
        g.write("#@ src.inputs_verdict = stamped\n")
        g.write("#@ obs.quantity = отпечатки расчётных спектров сетки; "
                "величина НЕ числовая\n")
        g.write("#@ obs.area = не применимо — площади здесь не снимаются\n")
        g.write("#@ obs.window = не применимо\n")
        g.write("#@ obs.shelf = не применимо\n")
        g.write("#@ obs.blurred = нет — спектры не размывались\n")
        g.write("#@ run.macro = %s\n" % ", ".join(sorted(macs)))
        g.write("file,sha1_12,E_prim_keV,N_primaries,N_with_signal,"
                "solid_angle_frac\n")
        for r in rows:
            g.write("%s\n" % ",".join(r))
    print("записано: %s  (%d узлов, штамп %s)" % (out, len(rows), stamp))
    return 0


if __name__ == "__main__":
    sys.exit(main())
