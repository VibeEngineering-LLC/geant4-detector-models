"""Сверка расчётных кривых ППП в точечных геометриях с кривыми ЛСРМ.

Точечная геометрия — самый чистый репер модели ДЕТЕКТОРА: сосуда нет,
самопоглощения в пробе нет, остаются только кристалл, отражатель, корпус и
защита. Диапазон у ЛСРМ здесь самый широкий: 59,5–2614 кэВ (24 точки на 5 см,
20 на 25 см), то есть проверяются и мягкий край, где всё решают корпус и
мёртвые слои, и жёсткий, где решает объём кристалла.

Розыгрыш в сетке — конус вокруг направления на детектор; эффективность
восстанавливается делением на долю телесного угла, записанную рядом с данными.
Для ППП это корректно (в пик идут практически только прямые кванты).

Крышка защиты: 5 см — ЗАКРЫТА (источник внутри полости), 25 см — ОТКРЫТА
(источник над защитой). Это разные прогоны, см. run_all_grids.py.
"""
import glob
import math
import os
import re
import sys

import numpy as np

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

# parse_efr живёт в инструментах репозитория, а не среди данных
sys.path.insert(0, str(paths.tools()))
from fetch_efr import parse_efr  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)
REF = str(paths.ref("Gamma-1S"))

GEOM = [("p5cm", "Точечная-5см", "5 см, крышка закрыта"),
        ("p25cm", "Точечная-25см", "25 см, крышка открыта")]
WIN = 6.0


from compare_lsrm import marinelli_k as _marinelli_k  # noqa: E402


def load(path):
    hist, N = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            continue
        if line and line[0].isdigit():
            e, c = line.split(",")
            hist[float(e)] = int(c)
    return hist, N


def mc_curve(tag):
    saf = os.path.join(BUILD, "grid", "%s_solidangle.txt" % tag)
    if not os.path.exists(saf):
        return {}
    frac = float(open(saf).read().strip())
    out = {}
    for p in glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv")):
        m = re.search(r"_E(\d+\.\d)\.csv$", p)
        if not m:
            continue
        E = float(m.group(1))
        hist, N = load(p)
        # сетка моноэнергий: блендов нет, узкое окно корректно
        peak = sum(c for e, c in hist.items() if abs(e - E) <= WIN)
        if peak > 0 and N:
            out[round(E, 1)] = ((peak / N) * frac, math.sqrt(peak) / N * frac)
    return out


if __name__ == "__main__":
    for tag, efr, title in GEOM:
        mc = mc_curve(tag)
        if not mc:
            print("%-6s сетка не готова" % tag)
            continue
        path = paths.efficiency_curve(efr)
        if not os.path.exists(path):
            print("%-6s нет %s" % (tag, efr))
            continue
        pts = [p for s in parse_efr(paths.read_text(path))
               for p in s["points"]]
        print("\n===== %s (%s): точек ЛСРМ %d" % (tag, title, len(pts)))
        # Гладкая кривая по расчётным точкам: у ЛСРМ 24 линии, а в сетке
        # 20 узлов, и совпадают далеко не все. Интерполяция полиномом в log-log
        # по СВОИМ точкам законна — та же процедура, что строит рабочую кривую;
        # χ²/dof подгонки печатается, чтобы интерполяция не была слепой.
        Eg = np.array(sorted(mc))
        yg = np.array([mc[e][0] for e in Eg])
        dyg = np.array([mc[e][1] for e in Eg])
        deg = 5
        cf = np.polyfit(np.log(Eg), np.log(yg), deg, w=yg / np.maximum(dyg, 1e-30))
        rr = (np.log(yg) - np.polyval(cf, np.log(Eg))) * yg / np.maximum(dyg, 1e-30)
        print("   интерполяция: полином %d-й степени по %d узлам, "
              "chi2/dof = %.2f" % (deg, len(Eg), (rr ** 2).sum() / (len(Eg) - deg - 1)))
        print("%9s %-8s %12s %12s %8s %8s %s" %
              ("E, кэВ", "нуклид", "эксп", "МК", "МК/эксп", "±", "источник МК"))
        logs, ws = [], []
        for E, eff, dp, nuc in sorted(pts):
            key = min(mc, key=lambda k: abs(k - E)) if mc else None
            if key is not None and abs(key - E) <= 1.0:
                m, dm = mc[key]
                src = "узел"
            elif Eg[0] <= E <= Eg[-1]:
                m = math.exp(np.polyval(cf, math.log(E)))
                dm = m * 0.02      # погрешность интерполяции, оценка
                src = "интерп."
            else:
                print("%9.1f %-8s %12.3e   вне диапазона сетки" % (E, nuc, eff))
                continue
            r = m / eff
            dr = r * math.hypot(dm / m, dp / 100)
            logs.append(math.log(r))
            ws.append(1.0 / (dr / r) ** 2)
            print("%9.1f %-8s %12.3e %12.3e %8.3f %8.3f %s"
                  % (E, nuc, eff, m, r, dr, src))
        if logs:
            lw = sum(l * w for l, w in zip(logs, ws)) / sum(ws)
            k = math.exp(lw)
            dev = [math.exp(l - lw) - 1 for l in logs]
            rms = math.sqrt(sum(d * d for d in dev) / len(dev))
            chi2 = (sum(w * (l - lw) ** 2 for l, w in zip(logs, ws))
                    / max(1, len(logs) - 1))
            print("\n   средневзвешенное МК/эксп = %.3f по %d точкам"
                  % (k, len(logs)))
            print("   разброс формы RMS = %.1f %%, chi2/dof = %.2f"
                  % (100 * rms, chi2))
            print("   для сравнения: маринелли %s; RadiaCode K_NORM 0,833"
                  % _marinelli_k())
