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


# --- Зонная аппроксимация по образцу ЛСРМ --------------------------------
#
# EffCalcMC/Efficiency строит кривую НЕ одним полиномом по всему диапазону,
# а ЗОНАМИ С ПЕРЕКРЫТИЕМ и сшивкой (оператор показал на EffReg.efa: зоны
# 20–245 / 70–859 / 459–3010 кэВ, степени 5/3/3, χ² = 4,47). Один полином
# 5-й степени по нашим 24 узлам давал χ²/ν = 11,5 — он гнёт середину, чтобы
# дотянуться до краёв. Зоны ниже — те же ЛСРМ-овские, масштабированные на
# диапазон сетки; в перекрытии соседние ветви смешиваются линейно по log E,
# что и есть «сшивка».
ZONES = [(None, 245.4, 5), (70.0, 859.4, 3), (459.2, None, 3)]


def zoned_fit(Eg, yg, dyg):
    """[(lo, hi, coeffs)] по зонам + функция интерполяции log-log."""
    lE, ly = np.log(Eg), np.log(yg)
    w = yg / np.maximum(dyg, 1e-30)
    fits = []
    for lo, hi, deg in ZONES:
        lo = Eg[0] if lo is None else lo
        hi = Eg[-1] if hi is None else hi
        m = (Eg >= lo * 0.999) & (Eg <= hi * 1.001)
        if m.sum() < deg + 2:      # зоне нужен запас узлов над степенью
            deg = max(1, m.sum() - 2)
        cf = np.polyfit(lE[m], ly[m], deg, w=w[m])
        rr = (ly[m] - np.polyval(cf, lE[m])) * w[m]
        chi2 = (rr ** 2).sum() / max(1, m.sum() - deg - 1)
        fits.append((lo, hi, deg, cf, chi2, int(m.sum())))

    def ev(E):
        x = math.log(E)
        # ветви, чья зона накрывает E; в перекрытии — линейная сшивка по logE
        hit = [(lo, hi, cf) for lo, hi, _d, cf, _c, _n in fits
               if lo * 0.999 <= E <= hi * 1.001]
        if not hit:
            lo, hi, _d, cf, _c, _n = fits[0] if E < fits[0][1] else fits[-1]
            return math.exp(np.polyval(cf, x))
        if len(hit) == 1:
            return math.exp(np.polyval(hit[0][2], x))
        (l1, h1, c1), (l2, h2, c2) = hit[0], hit[1]
        a, b = math.log(max(l1, l2)), math.log(min(h1, h2))
        t = 0.5 if b <= a else min(1.0, max(0.0, (x - a) / (b - a)))
        return math.exp((1 - t) * np.polyval(c1, x) + t * np.polyval(c2, x))

    return fits, ev


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
        # 20 узлов, и совпадают далеко не все. Интерполяция — ЗОННАЯ, по
        # образцу ЛСРМ (см. ZONES выше): один полином на весь диапазон гнул
        # середину ради краёв. χ²/ν печатается по каждой зоне.
        Eg = np.array(sorted(mc))
        yg = np.array([mc[e][0] for e in Eg])
        dyg = np.array([mc[e][1] for e in Eg])
        fits, ev = zoned_fit(Eg, yg, dyg)
        print("   интерполяция зонами (образец ЛСРМ), %d узлов:" % len(Eg))
        for lo, hi, deg, _cf, chi2, n in fits:
            print("      %6.1f–%6.1f кэВ  степень %d, узлов %2d, chi2/nu = %.2f"
                  % (lo, hi, deg, n, chi2))
        print("%9s %-8s %12s %12s %8s %8s %s" %
              ("E, кэВ", "нуклид", "эксп", "МК", "МК/эксп", "±", "источник МК"))
        logs, ws = [], []
        for E, eff, dp, nuc in sorted(pts):
            key = min(mc, key=lambda k: abs(k - E)) if mc else None
            if key is not None and abs(key - E) <= 1.0:
                m, dm = mc[key]
                src = "узел"
            elif Eg[0] <= E <= Eg[-1]:
                m = ev(E)
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
