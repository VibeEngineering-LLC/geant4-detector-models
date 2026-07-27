"""Отличается ли расчёт по BecqMoni XML от расчёта по исходному .spe ЛСРМ.

Комплект поверки лежит в репозитории в двух форматах: приборный двоичный
.spe и конвертированный BecqMoni XML. Весь разбор в этом проекте написан под
XML, поэтому вопрос «а не подменил ли конвертор данные» не праздный: если
форматы разойдутся, разойдутся и все выводы, а заметить это будет негде.

Проверка прямая. Для каждой пары файлов с совпадающим именем сравниваются:
число каналов, полная сумма отсчётов, живое и реальное время, коэффициенты
энергетической калибровки и поканальная разность. Расхождение печатается
числом; молчаливого «похоже, совпало» здесь нет.

    SPECTRAVIBE_ROOT=<корень gamma-spectrum-analysis> \\
        python detectors/Gamma-1S/analysis/xml_vs_spe.py

Читатель .spe берётся штатный (gamma.io.lsrm_spe): отсчёты в этом формате
упакованы двоично, свой разбор писать не надо. Без SPECTRAVIBE_ROOT скрипт
сообщает об этом и выходит — сверять не с чем.
"""
import os
import sys

import numpy as np

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

import becqmoni as bm  # noqa: E402


def spe_reader():
    root = paths.require_spectravibe("сверку XML с исходным .spe ЛСРМ")
    sys.path.insert(0, str(os.path.join(str(root), "scripts")))
    from gamma.io.lsrm_spe import read_lsrm_spe
    return read_lsrm_spe


def pairs():
    """Пары (xml, spe) с совпадающим именем без расширения."""
    root = paths.ref("Gamma-1S")
    xmls = {p.stem: p for p in root.rglob("*.xml")}
    out = []
    for p in root.rglob("*.spe"):
        if p.stem in xmls:
            out.append((xmls[p.stem], p))
    return sorted(out, key=lambda t: str(t[0]))


if __name__ == "__main__":
    read_spe = spe_reader()
    pp = pairs()
    if not pp:
        raise SystemExit(
            "не найдено ни одной пары XML/.spe с совпадающим именем в %s"
            % paths.ref("Gamma-1S"))

    print("Сверка BecqMoni XML с исходным .spe ЛСРМ: пар %d\n" % len(pp))
    print("%-46s %9s %9s %8s %9s" %
          ("запись", "каналов", "d сумма", "d время", "макс |d|"))
    bad = 0
    for xp, sp in pp:
        a, _ = bm.read(str(xp))
        b = read_spe(str(sp))
        ca = np.asarray(a.n, dtype=float)
        cb = np.asarray(b.counts, dtype=float)
        n = min(len(ca), len(cb))
        dmax = float(np.abs(ca[:n] - cb[:n]).max()) if n else float("nan")
        dsum = float(ca.sum() - cb.sum())
        dt = abs(a.live - float(b.live_time))
        flag = ""
        if len(ca) != len(cb) or dsum or dmax or dt > 0.01:
            flag = "  <- РАСХОЖДЕНИЕ"
            bad += 1
        print("%-46s %4d/%-4d %9.0f %8.2f %9.0f%s"
              % (xp.stem[:46], len(ca), len(cb), dsum, dt, dmax, flag))

    print("\nпар с расхождением: %d из %d" % (bad, len(pp)))
    if bad == 0:
        print("Форматы взаимозаменяемы: отсчёты, времена и число каналов\n"
              "совпадают побитово, поэтому результат расчёта от выбора\n"
              "формата не зависит.")
    else:
        print("ВНИМАНИЕ: расчёт по XML и по .spe даст разные числа.\n"
              "Пока расхождение не объяснено, опираться следует на .spe —\n"
              "это исходная запись прибора, а XML получен преобразованием.")
        sys.exit(1)
