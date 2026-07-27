"""Сверка расчётных кривых ППП в кюветах «Дента» и Петри с кривыми ЛСРМ.

Зачем отдельно. Пересчёт записей в кюветах дал систематическое превышение
активности над паспортом (1,24–1,84), причём ВНУТРИ одной геометрии цезий
даёт 1,07, а торий 1,84. Так геометрия не ошибается — ошибка так себя не
ведёт. Значит надо развести два подозрения:
  (а) геометрия кюветы в модели — тогда поедет ВСЯ кривая целиком;
  (б) прогоны цепочек Ra/Th в кюветах — тогда кривая по моноэнергиям будет
      в порядке, а поедут только нуклиды с цепочками.
Здесь проверяется (а): моно-сетка против кривой ЛСРМ, цепочки не участвуют.

Кюветы стоят на торце детектора (Distance = 0 в .efa), сетки считались при
ОИСН-16 ρ = 1,6 — ровно та засыпка, для которой построены кривые (192 г
в Денте и 96 г в Петри, см. описи комплекта).
"""
import glob
import math
import os
import re
import sys

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
WIN = 6.0

# Геометрия называется так же, как в именах файлов ЛСРМ; сам файл ищется
# по набору эталонных данных (paths.efficiency_curve), потому что раскладка
# у скачанного и у закоммиченного набора разная.
CASES = [("denta1.60", "Дента", "«Дента» 120 мл"),
         ("petri1.60", "Петри", "Петри 60 мл"),
         ("rho1.60", "Маринелли", "Маринелли 1 л (для сравнения)")]


def mc_curve(tag):
    out = {}
    for p in glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv")):
        m = re.search(r"_E(\d+\.\d)\.csv$", p)
        if not m:
            continue
        E = float(m.group(1))
        hist, N = {}, None
        for line in open(p, encoding="utf-8"):
            if line.startswith("#"):
                if "N_primaries" in line:
                    N = int(line.split("=")[1])
                continue
            if line and line[0].isdigit():
                e, c = line.split(",")
                hist[float(e)] = int(c)
        peak = sum(c for e, c in hist.items() if abs(e - E) <= WIN)
        if peak > 0 and N:
            out[round(E, 1)] = (peak / N, math.sqrt(peak) / N)
    return out


if __name__ == "__main__":
    for tag, efr, title in CASES:
        mc = mc_curve(tag)
        path = paths.efficiency_curve(efr)
        if not mc or not os.path.exists(path):
            print("%-12s данных нет" % tag)
            continue
        pts = [p for s in parse_efr(paths.read_text(path))
               for p in s["points"]]
        print("\n===== %s: точек ЛСРМ %d" % (title, len(pts)))
        print("%9s %-8s %12s %12s %8s" %
              ("E, кэВ", "нуклид", "эксп", "МК", "МК/эксп"))
        logs, ws = [], []
        for E, eff, dp, nuc in sorted(pts):
            key = min(mc, key=lambda k: abs(k - E))
            if abs(key - E) > 1.0:
                continue
            m, dm = mc[key]
            r = m / eff
            dr = r * math.hypot(dm / m, dp / 100)
            logs.append(math.log(r))
            ws.append(1.0 / (dr / r) ** 2)
            print("%9.1f %-8s %12.3e %12.3e %8.3f" % (E, nuc, eff, m, r))
        if logs:
            lw = sum(l * w for l, w in zip(logs, ws)) / sum(ws)
            k = math.exp(lw)
            dev = [math.exp(l - lw) - 1 for l in logs]
            rms = math.sqrt(sum(d * d for d in dev) / len(dev))
            print("   средневзвешенное МК/эксп = %.3f по %d точкам, "
                  "RMS формы %.1f %%" % (k, len(logs), 100 * rms))
    print("\nЕсли у кювет отношение около 1, геометрия кювет НЕ виновата —")
    print("тогда искать в прогонах цепочек. Если поехало — виновата геометрия.")
