"""Эффективная толщина d_eff сосуда Маринелли 1 л Гамма-1С.

Метод тот же, что в radiacode-curves: отношение эффективностей при двух
плотностях той же матрицы описывается отношением поправок самопоглощения
    eff(ро2)/eff(ро1) = f(мю*ро2*d) / f(мю*ро1*d),  f(x) = (1-e^-x)/x,
и d — единственный параметр. ЛСРМ в .efa записал Thick=31, DThick=2 мм.

Массовый коэффициент ослабления ОИСН-16 берётся из `build/mu_oisn16.csv` —
его считает mucalc.exe ТОЙ ЖЕ физикой (EmStandardPhysics_option4), что и
транспорт, по той же смеси (G1SDetector::MakeMatrix). Так исключён разнобой
между сечениями расчёта и сечениями поправки.

Вписанная руками таблица XCOM оставлена только как независимая сверка:
расхождение двух источников печатается. Числа, введённые вручную, в этом
проекте не используются в расчёте — на таком уже обжигались.
"""
import glob
import math
import os
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402


BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)

# NIST XCOM, мю/ро (см²/г), total attenuation WITHOUT coherent. ТОЛЬКО для
# сверки с mucalc, в подгонку не идёт. Узлы: кэВ.
XCOM = {
    "H": [(50, 0.3355), (60, 0.3260), (80, 0.3091), (100, 0.2944), (150, 0.2651),
          (200, 0.2429), (300, 0.2112), (400, 0.1893), (500, 0.1729),
          (600, 0.1599), (800, 0.1405), (1000, 0.1263), (1250, 0.1129),
          (1500, 0.1027), (2000, 0.0876), (3000, 0.0691)],
    "C": [(50, 0.1871), (60, 0.1753), (80, 0.1610), (100, 0.1514), (150, 0.1347),
          (200, 0.1229), (300, 0.1066), (400, 0.0954), (500, 0.0870),
          (600, 0.0805), (800, 0.0707), (1000, 0.0636), (1250, 0.0569),
          (1500, 0.0518), (2000, 0.0444), (3000, 0.0356)],
    "N": [(50, 0.1889), (60, 0.1751), (80, 0.1594), (100, 0.1493), (150, 0.1324),
          (200, 0.1207), (300, 0.1046), (400, 0.0936), (500, 0.0853),
          (600, 0.0790), (800, 0.0693), (1000, 0.0624), (1250, 0.0558),
          (1500, 0.0508), (2000, 0.0436), (3000, 0.0350)],
    "O": [(50, 0.1958), (60, 0.1789), (80, 0.1610), (100, 0.1501), (150, 0.1325),
          (200, 0.1206), (300, 0.1044), (400, 0.0934), (500, 0.0851),
          (600, 0.0788), (800, 0.0692), (1000, 0.0623), (1250, 0.0557),
          (1500, 0.0507), (2000, 0.0436), (3000, 0.0351)],
    "Fe": [(50, 1.744), (60, 1.100), (80, 0.5306), (100, 0.3327), (150, 0.1744),
           (200, 0.1329), (300, 0.1043), (400, 0.0909), (500, 0.0821),
           (600, 0.0755), (800, 0.0659), (1000, 0.0591), (1250, 0.0529),
           (1500, 0.0482), (2000, 0.0418), (3000, 0.0344)],
}
W = {"H": 0.022, "C": 0.206, "N": 0.009, "O": 0.049, "Fe": 0.714}


def mu_xcom(E):
    """Сверочное значение по вписанной таблице XCOM (в подгонку НЕ идёт)."""
    tot = 0.0
    for el, w in W.items():
        xs, ys = zip(*XCOM[el])
        e = min(max(E, xs[0]), xs[-1])
        i = 1
        while i < len(xs) - 1 and xs[i] < e:
            i += 1
        t = (math.log(e) - math.log(xs[i - 1])) / (math.log(xs[i]) - math.log(xs[i - 1]))
        tot += w * math.exp(math.log(ys[i - 1]) + t * (math.log(ys[i]) - math.log(ys[i - 1])))
    return tot


def load_mu():
    """мю/ро ОИСН-16 из mucalc: та же физика и та же смесь, что в транспорте."""
    path = os.path.join(BUILD, "mu_oisn16.csv")
    tab = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or line.startswith("E_keV"):
            continue
        e, m = line.split(",")
        tab[round(float(e), 3)] = float(m)
    return tab


MU = load_mu()


def mu_mass(E):
    key = min(MU, key=lambda k: abs(k - E))
    if abs(key - E) > 0.5:
        raise KeyError("нет mu для E = %.3f кэВ" % E)
    return MU[key]


def f(x):
    return (1 - math.exp(-x)) / x if x > 1e-9 else 1.0


def read_run(path):
    N, E0 = None, None
    peak = 0
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            elif "E_prim_keV" in line:
                E0 = float(line.split("=")[1])
            continue
        if line and line[0].isdigit():
            e, c = line.split(",")
            if abs(float(e) - E0) <= 6.0:
                peak += int(c)
    return E0, peak, N


def curve(tag):
    out = {}
    for fn in sorted(glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv"))):
        E0, p, N = read_run(fn)
        if p > 0:
            out[round(E0, 3)] = (p / N, math.sqrt(p) / N)
    return out


# Пары сеток по геометриям. Опорные значения ЛСРМ: .efa этого экземпляра
# детектора и таблица кювет «Прецизионные измерения», с. 11.
VESSELS = {
    # ключ: (метка сетки 1, ро1, метка сетки 2, ро2, .efa мм, dEfa, табл. мм)
    "marinelli": ("rho1.00", 1.00, "rho1.60", 1.60, 31.0, 2.0, 26.0),
    "denta":     ("denta0.60", 0.60, "denta1.60", 1.60, 33.0, 3.0, 36.0),
    "petri":     ("petri0.60", 0.60, "petri1.60", 1.60, 10.0, 1.0, 15.0),
}


def fit_d(c1, c2, rho1, rho2, Es, dmin=2.0, dmax=60.0, step=0.25):
    """Скан по d: возвращает (d, lo, hi по dchi2=1, chi2/dof)."""
    def chi2_of(d_mm):
        d = d_mm / 10.0
        s = 0.0
        for E in Es:
            m1, dm1 = c1[E]
            m2, dm2 = c2[E]
            r = m2 / m1
            dr = r * math.sqrt((dm1 / m1) ** 2 + (dm2 / m2) ** 2)
            mu = mu_mass(E)
            model = f(mu * rho2 * d) / f(mu * rho1 * d)
            s += ((r - model) / dr) ** 2
        return s

    grid = [dmin + step * i for i in range(int((dmax - dmin) / step) + 1)]
    chis = [(chi2_of(d), d) for d in grid]
    best_chi, best = min(chis)
    band = [d for c, d in chis if c <= best_chi + 1.0]
    return best, min(band), max(band), best_chi / max(1, len(Es) - 1)


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    print("Эффективная толщина по парам плотностей. Модель f(x)=(1-e^-x)/x,\n"
          "d — единственный параметр; совпадение точек между собой — проверка.\n")
    for name, (t1, r1, t2, r2, efa, defa, tab) in VESSELS.items():
        if only and only != name:
            continue
        c1, c2 = curve(t1), curve(t2)
        Es = sorted(set(c1) & set(c2))
        if len(Es) < 5:
            print("%-10s сеток нет или мало точек (%d) — прогоны не готовы"
                  % (name, len(Es)))
            continue
        d, lo, hi, chi = fit_d(c1, c2, r1, r2, Es)
        dev = [100 * (mu_mass(E) - mu_xcom(E)) / mu_xcom(E) for E in Es]
        print("%-10s ро %.2f/%.2f, точек %d" % (name, r1, r2, len(Es)))
        print("   d_eff = %.1f мм  (dchi2=1: %.1f..%.1f)  chi2/dof = %.2f"
              % (d, lo, hi, chi))
        print("   ЛСРМ .efa: %.0f ± %.0f мм;  таблица кювет: %.0f мм"
              % (efa, defa, tab))
        # Строгий критерий: складываем свою полосу и заявленную ЛСРМ.
        sig = math.hypot(0.5 * (hi - lo), defa)
        z = abs(d - efa) / sig if sig > 0 else 99
        verdict = ("СОВПАДАЕТ" if z < 1 else
                   "на границе" if z < 2 else "РАСХОДИТСЯ")
        print("   -> с .efa %s: %+.1f мм = %.1f сигма" % (verdict, d - efa, z))
        print("   мю: mucalc против вписанной XCOM %+.1f..%+.1f %%\n"
              % (min(dev), max(dev)))
