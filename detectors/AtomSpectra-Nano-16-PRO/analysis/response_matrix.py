# -*- coding: utf-8 -*-
"""Функция полного отклика детектора по сетке моноэнергетических прогонов.

Отклик НЕ считается отдельно: спектр энерговыделения пишется на каждом узле
самим прогоном (`main.cc`, гистограмма 3700 каналов по 1 кэВ), поэтому матрица
собирается из уже посчитанного. Единственное, что здесь делается заново, —
нормировка, свёртка с приборным разрешением и укладка в матрицу.

Что на выходе:

* `response_matrix_raw_10keV.csv` — отклик БЕЗ размытия, только сложенный в
  каналы по 10 кэВ. Это чистая физика энерговыделения;
* `response_matrix_10keV.csv` — он же, свёрнутый с гауссовым разрешением
  ПШПВ(E) = 41,6 · sqrt(E / 661,657) кэВ (собственная запись Cs-137 прибора).
  Это и есть приборная функция отклика.

Строка матрицы — энергия падающего кванта, столбец — канал энерговыделения,
значение — вероятность на ОДИН испущенный квант в 4π. Сумма строки равна
эффективности регистрации на этом узле; сумма по каналам вокруг E равна
эффективности по ППП с точностью до ширины окна.

Разрешение навешивается ЗДЕСЬ, а не в прогоне: физика от него не зависит, и
переснимать сетку при уточнении ПШПВ не требуется.

    python analysis/response_matrix.py <каталог спектров> [шаг канала, кэВ]
"""
import io
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(_HERE, "..", "results"))
FWHM_662 = 41.60          # кэВ, собственная запись Cs-137
E_MAX = 3200.0            # верх шкалы энерговыделения


def read_spectrum(path):
    head, rows = {}, []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("# ").split("=", 1)
                head[k.strip()] = v.strip()
            continue
        if not ln or ln.startswith("E_keV"):
            continue
        e, c = ln.split(",")
        rows.append((float(e), float(c)))
    return head, rows


def fwhm(e):
    return FWHM_662 * math.sqrt(max(e, 1.0) / 661.657)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    src = sys.argv[1]
    step = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
    nch = int(E_MAX / step)

    # Спутники прогона (`_emit.csv`, `_chan.csv`) имеют ту же шапку и другой
    # состав колонок: при чтении как спектра они ломаются на первой же строке
    # данных, то есть уже после того, как выглядели годными.
    files = sorted(f for f in os.listdir(src)
                   if f.endswith(".csv")
                   and not f.endswith("_emit.csv")
                   and not f.endswith("_chan.csv"))
    if not files:
        raise SystemExit("в %s нет спектров" % src)

    stamps, rows_raw, rows_bro, es = set(), [], [], []
    diag = []
    for fn in files:
        head, spec = read_spectrum(os.path.join(src, fn))
        e0 = float(head["E_prim_keV"])
        n = float(head["N_primaries"])
        frac = float(head["solid_angle_frac"])
        stamps.add(head.get("src_sha1", "?"))
        # Вес одного отсчёта: доля телесного угла делит, число первичных — тоже.
        # Тот же переход к 4pi, что и в export_curve.py: eps = C * f / N.
        w = frac / n

        raw = [0.0] * nch
        bro = [0.0] * nch
        tot = 0.0
        for e, c in spec:
            tot += c
            k = int(e / step)
            if k >= nch:
                k = nch - 1
            raw[k] += c * w
            # Свёртка с разрешением: отсчёт размазывается гауссианой шириной
            # ПШПВ(e). Ядро нестационарное, поэтому строится на каждый отсчёт,
            # а не один раз на всю строку.
            s = fwhm(e) / 2.3548
            lo = max(0.0, e - 4.0 * s)
            hi = min(E_MAX, e + 4.0 * s)
            k0, k1 = int(lo / step), int(hi / step) + 1
            acc = []
            norm = 0.0
            for k in range(k0, min(k1, nch)):
                x = (k + 0.5) * step
                g = math.exp(-0.5 * ((x - e) / s) ** 2)
                acc.append((k, g))
                norm += g
            if norm <= 0:
                bro[min(int(e / step), nch - 1)] += c * w
                continue
            for k, g in acc:
                bro[k] += c * w * g / norm

        es.append(e0)
        rows_raw.append(raw)
        rows_bro.append(bro)
        # Диагностика: отсчёты в окне шириной ПШПВ у самой линии и в среднем
        # на элемент разрешения по всей шкале. Это и есть предел применимости
        # матрицы: при единицах отсчётов на элемент она шумит сильнее, чем
        # отличается от соседней строки.
        fw = fwhm(e0)
        at_line = sum(c for e, c in spec if abs(e - e0) <= fw / 2.0)
        mean_res = tot * fw / E_MAX if tot else 0.0
        diag.append((e0, tot, fw, at_line, mean_res))

    if len(stamps) > 1:
        raise SystemExit("спектры разных ревизий: %s" % ", ".join(sorted(stamps)))
    stamp = stamps.pop()

    def dump(name, rows, title):
        path = os.path.join(RES, name)
        with io.open(path, "w", encoding="utf-8", newline="") as g:
            g.write("#@ stamp.version = 1\n")
            g.write("#@ src.script = detectors/AtomSpectra-Nano-16-PRO/"
                    "analysis/response_matrix.py\n")
            g.write("#@ src.spectra_sha1 = %s\n" % stamp)
            g.write("#@ src.inputs_n = %d\n" % len(rows))
            g.write("#@ src.inputs_verdict = stamped\n")
            g.write("#@ obs.quantity = %s\n" % title)
            g.write("#@ obs.area = не применимо — матрица, не площадь\n")
            g.write("#@ obs.window = канал %g кэВ\n" % step)
            g.write("#@ obs.blurred = %s\n"
                    % ("да, гауссиана ПШПВ(E) = 41,60*sqrt(E/661,657) кэВ"
                       if "приборн" in title else "нет"))
            g.write("#@ obs.reference_plane = наружная поверхность торцевой "
                    "крышки корпуса\n")
            g.write("# строка — энергия падающего кванта, кэВ; столбцы — "
                    "каналы энерговыделения по %g кэВ\n" % step)
            g.write("# значение — вероятность на ОДИН испущенный квант в 4pi\n")
            g.write("E_inc_keV," + ",".join("%g" % ((k + 0.5) * step)
                                            for k in range(nch)) + "\n")
            for e0, row in zip(es, rows):
                g.write("%.3f," % e0 + ",".join("%.6e" % v for v in row) + "\n")
        return path

    p1 = dump("response_matrix_raw_%gkeV.csv" % step, rows_raw,
              "функция отклика без размытия, вероятность на квант в 4pi")
    p2 = dump("response_matrix_%gkeV.csv" % step, rows_bro,
              "приборная функция отклика, вероятность на квант в 4pi")

    print("узлов %d, каналов %d по %g кэВ, штамп %s" % (len(es), nch, step, stamp))
    print("записано: %s" % p1)
    print("записано: %s" % p2)
    print()
    print("  E, кэВ  событий  ПШПВ, кэВ  в пике на ПШПВ  в среднем на элемент")
    for e0, tot, fw, at_line, mean_res in diag:
        if round(e0) not in (30, 180, 680, 1480, 2480, 3000):
            continue
        print("%8.0f %8.0f %10.1f %15.0f %21.1f"
              % (e0, tot, fw, at_line, mean_res))
    return 0


if __name__ == "__main__":
    sys.exit(main())
