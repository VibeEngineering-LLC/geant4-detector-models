"""Введена ли поправка на каскадное суммирование в саму кривую ЛСРМ?

Вопрос решающий: если НЕ введена — к расчёту надо применять поправку C, если
введена — применять её будет ошибкой, вносящей систематику в 10–20 % там, где
её нет. Ответ ищется в данных самого ЛСРМ, без обращения к нашей модели.

Три независимых теста.

ТЕСТ 1. Точка Cs-137 (каскада нет) против интерполяции между соседними
точками Bi-214 (каскад есть). Без поправки соседи занижены, и цезий обязан
торчать над ними на величину C.

ТЕСТ 2. Гладкость кривой. Кривая эффективности физически гладкая. Сравниваем
хи-квадрат подгонки полиномом в log-log для кривой КАК ЕСТЬ и для кривой,
УМНОЖЕННОЙ на C (то есть «исправленной задним числом»). Если умножение портит
гладкость, поправка там уже была.

ТЕСТ 3. Куда ложатся некаскадные точки относительно кривой, проведённой
ТОЛЬКО по каскадным. Без поправки они лягут выше.
"""
import sys
import math
import os

import numpy as np

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402


# parse_efr живёт в инструментах репозитория
sys.path.insert(0, str(paths.tools()))
from fetch_efr import parse_efr  # noqa: E402

# Поправки из summing.py (прогоны Geant4, контроль на Cs-137 и K-40 пройден).
# 1764,491 исключена: там сумм-пик садится на саму линию, и C неустойчива
# к ширине окна (0,95…1,19) — в тесты такую величину брать нельзя.
C = {583.187: 1.146, 609.320: 1.115, 768.360: 1.149,
     1120.294: 1.156, 2614.511: 1.185}
NOCASCADE = {661.657: "Cs-137", 1460.822: "K-40"}


def fit_chi2(E, y, dy, deg):
    x, ly = np.log(E), np.log(y)
    w = y / np.maximum(dy, 1e-30)
    c = np.polyfit(x, ly, deg, w=w)
    r = (ly - np.polyval(c, x)) * w
    return float((r ** 2).sum()) / max(1, len(E) - deg - 1), c


if __name__ == "__main__":
    # Точки берём прямо из измеренной кривой ЛСРМ, а не из промежуточного
    # efr_points.json: тот файл — продукт загрузчика и в репозиторий не
    # коммитится, а .efr лежит в эталонном наборе.
    curve = paths.efficiency_curve("Маринелли", "efr")
    if curve is None:
        raise SystemExit("не найдена измеренная кривая .efr для Маринелли "
                         "в %s" % paths.ref("Gamma-1S"))
    pts = {}
    for sec in parse_efr(paths.read_text(curve)):
        for E, eff, dpct, nuc in sec["points"]:
            pts[round(E, 3)] = (eff, dpct, nuc)

    # ---- ТЕСТ 1 -----------------------------------------------------------
    LO, HI, TEST = 609.32, 768.36, 661.657
    e1, d1, _ = pts[LO]
    e2, d2, _ = pts[HI]
    ec, dc, _ = pts[TEST]
    s = (math.log(e2) - math.log(e1)) / (math.log(HI) - math.log(LO))
    pred = math.exp(math.log(e1) + s * (math.log(TEST) - math.log(LO)))
    wgt = (math.log(TEST) - math.log(LO)) / (math.log(HI) - math.log(LO))
    dpred = pred * math.hypot((1 - wgt) * d1 / 100, wgt * d2 / 100)
    r = ec / pred
    dr = r * math.hypot(dc / 100, dpred / pred)
    print("ТЕСТ 1. Cs-137 против интерполяции 609,3 <-> 768,4 (обе Bi-214)")
    print("   точка Cs-137 выше интерполяции в %.3f ± %.3f раза" % (r, dr))
    print("   ждём 1,00 если исправлено;  %.3f если нет" % C[609.320])
    print("   отклонение от «исправлено»    %.1f сигма" % (abs(r - 1) / dr))
    print("   отклонение от «не исправлено» %.1f сигма"
          % (abs(r - C[609.320]) / dr))

    # ---- ТЕСТ 2 -----------------------------------------------------------
    have = sorted(set(list(C) + list(NOCASCADE)) & set(pts))
    E = np.array(have)
    y = np.array([pts[e][0] for e in have])
    dy = np.array([pts[e][0] * pts[e][1] / 100 for e in have])
    k = np.array([C.get(e, 1.0) for e in have])
    print("\nТЕСТ 2. Гладкость кривой по %d точкам (полином 3-й степени в log-log)"
          % len(have))
    for deg in (2, 3):
        c0, _ = fit_chi2(E, y, dy, deg)
        c1, _ = fit_chi2(E, y * k, dy * k, deg)
        print("   степень %d:  как есть chi2/dof = %5.2f   после умножения на C = %5.2f  -> %s"
              % (deg, c0, c1, "исправлено" if c0 < c1 else "НЕ исправлено"))

    # ---- ТЕСТ 3 -----------------------------------------------------------
    casc = [e for e in have if e in C]
    Ec = np.array(casc)
    yc = np.array([pts[e][0] for e in casc])
    dyc = np.array([pts[e][0] * pts[e][1] / 100 for e in casc])
    _, coef = fit_chi2(Ec, yc, dyc, 2)
    print("\nТЕСТ 3. Некаскадные точки против кривой ТОЛЬКО по каскадным")
    for e, nm in NOCASCADE.items():
        if e not in pts:
            continue
        model = math.exp(np.polyval(coef, math.log(e)))
        print("   %-7s %8.3f кэВ:  точка / кривая = %.3f   (ждём 1,00 если "
              "исправлено, %.2f если нет)"
              % (nm, e, pts[e][0] / model, np.interp(e, Ec, [C[x] for x in casc])))
