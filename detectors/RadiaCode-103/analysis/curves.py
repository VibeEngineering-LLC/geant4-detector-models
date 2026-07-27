# -*- coding: utf-8 -*-
"""Рабочая кривая эффективности: гладкий предел × самопоглощение × нормировка.

ЗАЧЕМ НУЖЕН ЭТОТ СЛОЙ. Точечные значения eps_p из сетки на высоких энергиях
шумные (до 7 % на точку), и брать их напрямую нельзя: у отношений вылезает
f_abs > 1, то есть «самопоглощение увеличивает эффективность». Это не физика, а
статистика. Физика же гладкая, и её удаётся выразить тремя множителями:

    eps_p(E, матрица, rho) = eps_lim(E) * f(mu(E) * rho * d) * k_norm

    eps_lim(E)  — предел без самопоглощения (конфигурация «воздух»), гладкая
                  функция энергии: подгонка полиномом по log E с весами 1/сигма²;
    f(x)        — (1-e^-x)/x, средняя вероятность пройти пробу без взаимодействия
                  при равномерном распределении длин пути; d — чисто
                  геометрическая эффективная толщина, одна на сосуд;
    k_norm      — 0,80 ± 0,02, нормировка абсолютной шкалы по аттестованной пробе.

Проверка законности сглаживания — в конце: невязка по ВСЕМ 150 точкам сетки
сравнивается с их собственной статистической погрешностью. Если совпадает,
сглаживание не выбросило ничего, кроме шума.

Запуск:  python curves.py m500          — отчёт о качестве сведения
         python curves.py m500 --json   — выгрузка для статьи
"""
import csv
import json
import os
import sys

import numpy as np

# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec
from selfabs_fit import load_mu, model as fabs

HERE = os.path.dirname(os.path.abspath(__file__))

# Нормировка абсолютной шкалы. Множитель определён так, чтобы ИМЕННО ЭТА
# сглаженная кривая воспроизводила аттестованную привязку на реперной пробе
# (Cs-137 662 кэВ, сушёная черника 246 г в пятисотке): 3,2095e-4 измеренное
# против 3,8541e-4 расчётного до нормировки. Пересчитывается normalization.py;
# если пересчитана сетка, сверить оттуда.
#
# ВАЖНО: раньше здесь стояло 0,80 — множитель к ТОЧКЕ СЕТКИ, у которой на 662 кэВ
# самопоглощение вышло 1,0012, то есть заведомо шумное. Кривая и нормировка должны
# быть согласованы между собой, иначе калькулятор не воспроизводит репер.
K_NORM = 0.833
D_K_NORM = 0.03        # аттестация 2 % + подгонка пика 0,8 % + сведение кривой 3 %
DEG = 5                # степень полинома по log E; выбор обоснован в report()

# Нижняя граница гладкого описания. Ниже неё лежат K-края иода (33,2 кэВ) и
# цезия (36,0 кэВ): на них eps_lim(E) имеет РАЗРЫВ, и никакой полином его не
# воспроизведёт. Точка 30 кэВ поэтому из подгонки исключена и остаётся только
# табличной.
E_SMOOTH_MIN = 46.0

# Рабочая область кривой. Снизу — 50 кэВ: ниже прибор всё равно не годится для
# количественного анализа (K-края кристалла, поглощение в корпусе и стенке
# колодца, разрешение), да и формула самопоглощения там неверна.
E_USE_MIN, E_USE_MAX = 50.0, 3000.0

# Точность сведения, % (ско невязки по точкам сетки; см. report()).
ACCURACY = {"50..80": 4.0, "80..800": 3.0, "800..3000": "в пределах статистики"}

# Область, на которой определяется эффективная толщина d. Снизу — там, где
# формула (1-e^-x)/x ещё верна: она предполагает, что в пик попадают только
# НЕрассеянные фотоны, а на 30-60 кэВ комптоновское рассеяние отнимает меньше
# ширины пика (на 30 кэВ рассеяние на 90° — это 1,7 кэВ при ширине 12 кэВ), и
# рассеянные фотоны остаются в пике. Сверху — там, где самопоглощение ещё
# измеримо: выше 800 кэВ f_abs = 0,97 при любом разумном d, и d не определяется.
E_FIT_LO, E_FIT_HI = 80.0, 800.0


def load_eff(vessel):
    """-> {config: (E, eps_p, d_eps_p, eps_t, d_eps_t, matrix, rho)}"""
    out = {}
    with open(rcspec.rdir("efficiency.csv", v=vessel), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out.setdefault(r["config"], []).append(
                (float(r["E_keV"]), float(r["eps_p"]), float(r["d_eps_p"]),
                 float(r["eps_t"]), float(r["d_eps_t"]),
                 r["matrix"], float(r["rho_gcm3"])))
    for k in out:
        out[k].sort()
    return out


def fit_limit(E, y, dy, deg=DEG, cov=False):
    """Гладкий предел: полином по log E, веса 1/сигма² в логарифмах.

    При cov=True возвращает ещё и ковариационную матрицу коэффициентов.
    numpy масштабирует её на chi2/dof, то есть в неё уже входит превышение
    разброса точек над чистой статистикой. Из этой матрицы получается
    энергетический ход погрешности самой кривой:
        sigma_ln eps(E) = sqrt( v(E)^T C v(E) ),  v = [x^deg, ..., x, 1].
    Так же считают погрешность кривой у ЛСРМ (через ортогональные полиномы)
    и у ORTEC [GammaVision, разд. 6.12.7].
    """
    x, ly = np.log(E), np.log(y)
    w = y / np.maximum(dy, 1e-30)          # d(ln y) = dy/y
    if cov:
        return np.polyfit(x, ly, deg, w=w, cov=True)
    return np.polyfit(x, ly, deg, w=w)


def sigma_limit(cov_m, E, deg=DEG):
    """Относительная погрешность гладкого предела в точке E."""
    v = np.array([np.log(float(E)) ** k for k in range(deg, -1, -1)])
    return float(np.sqrt(v @ np.asarray(cov_m) @ v))


def limit_at(c, E):
    return np.exp(np.polyval(c, np.log(np.asarray(E, dtype=float))))


def mu_of(mu_tab, matrix, E):
    """mu/rho(E) для матрицы, интерполяция в логарифмах. Зависит только от
    состава, поэтому годится для любой плотности."""
    pts = sorted((e, v) for (m, _r, e), v in mu_tab.items() if m == matrix)
    xs = np.array([p[0] for p in pts])
    ys = np.array([p[1] for p in pts])
    return np.exp(np.interp(np.log(np.asarray(E, dtype=float)),
                            np.log(xs), np.log(ys)))


def collect_f(eff, mu_tab, cp, lo=0.0, hi=1e9):
    """Точки самопоглощения: (E, mu*rho, F измеренное, вес, стат. ошибка)."""
    E, MU, F, W, S = [], [], [], [], []
    for _cfg, rows in eff.items():
        for Ei, ep, dep, _et, _det, matrix, rho in rows:
            if matrix == "air" or ep <= 0 or not (lo <= Ei <= hi):
                continue
            E.append(Ei)
            F.append(ep / float(limit_at(cp, Ei)))
            MU.append(float(mu_of(mu_tab, matrix, Ei)) * rho)
            W.append((ep / dep) ** 2)
            S.append(dep / ep)
    return map(np.array, (E, MU, F, W, S))


def fit_d(vessel, eff, mu_tab, cp):
    """Эффективная толщина d — одна на сосуд, по области E_FIT_LO..E_FIT_HI.

    Два отличия от selfabs_fit.py: отношение берётся не к ШУМНОЙ точке «воздух»,
    а к сглаженному пределу, и подгонка ограничена областью, где формула вообще
    применима (см. пояснение к E_FIT_LO/E_FIT_HI).
    """
    _E, MU, F, W, _S = collect_f(eff, mu_tab, cp, E_FIT_LO, E_FIT_HI)
    ds = np.linspace(0.2, 8.0, 2000)
    chi = [(W * (fabs(MU * d) - F) ** 2).sum() for d in ds]
    d = ds[int(np.argmin(chi))]
    resid = (fabs(MU * d) - F) / F
    return d, resid, len(F)


class Curve:
    """Готовая кривая для одного сосуда."""

    def __init__(self, vessel):
        self.vessel = vessel
        self.eff = load_eff(vessel)
        self.mu = load_mu()
        air = [r for r in self.eff["full_air_0.00"] if r[0] >= E_SMOOTH_MIN]
        E = np.array([r[0] for r in air])
        self.cp, self.cov_p = fit_limit(E, np.array([r[1] for r in air]),
                                        np.array([r[2] for r in air]), cov=True)
        self.ct = fit_limit(E, np.array([r[3] for r in air]),
                            np.array([r[4] for r in air]))
        self.d, self.resid, self.npts = fit_d(vessel, self.eff, self.mu, self.cp)
        self.Emin, self.Emax = float(E[0]), float(E[-1])

    def sigma_p(self, E):
        """Относительная погрешность подгонки кривой ППП в точке E."""
        return sigma_limit(self.cov_p, E)

    def f_abs(self, E, matrix, rho):
        return fabs(mu_of(self.mu, matrix, E) * rho * self.d)

    def eps_p(self, E, matrix="water", rho=1.0, norm=True):
        v = limit_at(self.cp, E) * self.f_abs(E, matrix, rho)
        return v * (K_NORM if norm else 1.0)

    def eps_t(self, E, matrix="water", rho=1.0, norm=True):
        v = limit_at(self.ct, E) * self.f_abs(E, matrix, rho)
        return v * (K_NORM if norm else 1.0)


def report(vessel):
    c = Curve(vessel)
    print("сосуд %s: сведение сетки к формуле" % vessel)
    print("  предел без самопоглощения: полином степени %d по log E" % DEG)

    # обоснование степени: точки выше K-краёв
    air = [r for r in c.eff["full_air_0.00"] if r[0] >= E_SMOOTH_MIN]
    E = np.array([r[0] for r in air])
    y = np.array([r[1] for r in air])
    dy = np.array([r[2] for r in air])
    print("  подгонка по %d точкам от %.1f кэВ (ниже — K-края иода и цезия)"
          % (len(E), E[0]))
    print("\n  %-8s %10s %14s" % ("степень", "хи2/спс", "макс.невязка"))
    for deg in (3, 4, 5, 6):
        cc = fit_limit(E, y, dy, deg)
        r = (limit_at(cc, E) - y) / dy
        chi = (r ** 2).sum() / max(len(E) - deg - 1, 1)
        print("  %-8d %10.2f %13.1f %%"
              % (deg, chi, 100 * np.abs((limit_at(cc, E) - y) / y).max()))

    print("\n  эффективная толщина d = %.3f см (по %d точкам в %.0f..%.0f кэВ)"
          % (c.d, c.npts, E_FIT_LO, E_FIT_HI))
    print("  невязка формулы там: ско %.1f %%, максимум %.1f %%"
          % (100 * c.resid.std(), 100 * np.abs(c.resid).max()))

    # главная проверка: невязка против собственной статистики точек, по полосам
    print("\n  %-16s %5s %12s %12s %7s" %
          ("полоса, кэВ", "точек", "невязка,ско", "статистика", "отн."))
    bands = [(29, 46, "ниже K-краёв"), (46, 80, ""), (80, 800, "область d"),
             (800, 3100, "самопогл. мало")]
    for lo, hi, note in bands:
        res, sig = [], []
        for _cfg, rows in c.eff.items():
            for Ei, ep, dep, _t, _dt, matrix, rho in rows:
                if not (lo <= Ei < hi):
                    continue
                pred = float(c.eps_p(Ei, matrix, rho, norm=False))
                res.append(ep / pred - 1.0)
                sig.append(dep / ep)
        if not res:
            continue
        res, sig = np.array(res), np.array(sig)
        print("  %-16s %5d %11.1f %% %11.1f %% %6.1f  %s"
              % ("%d..%d" % (lo, hi), len(res), 100 * res.std(),
                 100 * sig.mean(), res.std() / sig.mean(), note))
    return c


def dump_json(vessels):
    """Выгрузка для интерактивной статьи: коэффициенты, mu/rho, кривые LSRM."""
    from compare_lsrm import LSRM_FILE, REF, read_lsrm
    out = {"k_norm": K_NORM, "d_k_norm": D_K_NORM, "vessels": {}, "mu": {},
           "lsrm": {}, "E_use": [E_USE_MIN, E_USE_MAX], "accuracy": ACCURACY,
           "deg": DEG, "E_fit": [E_FIT_LO, E_FIT_HI]}
    mats = set()
    for v in vessels:
        c = Curve(v)
        out["vessels"][v] = {
            "poly_p": list(c.cp), "poly_t": list(c.ct), "d_cm": c.d,
            "cov_p": [list(row) for row in c.cov_p],
            "E_min": c.Emin, "E_max": c.Emax,
            "resid_pct": round(100 * float(c.resid.std()), 2),
            "points": {cfg: [[r[0], r[1], r[2], r[3], r[4]] for r in rows]
                       for cfg, rows in c.eff.items()},
        }
        for cfg, rows in c.eff.items():
            mats.add(rows[0][5])
        p = os.path.join(REF, LSRM_FILE[v])
        if os.path.exists(p):
            El, epl, unc = read_lsrm(p)
            good = unc < 0.30
            out["lsrm"][v] = {"E": list(El[good]), "eps": list(epl[good]),
                              "unc": list(unc[good])}
    # проникновение электронов: доля первичных, давших сигнал в кристалле
    for v in vessels:
        p = rcspec.rdir("beta_transmission.csv", v=v)
        if not os.path.exists(p):
            continue
        b = {}
        with open(p, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                b.setdefault(r["config"], []).append(
                    [float(r["E_keV"]), float(r["frac"]), float(r["d_frac"]),
                     float(r["mean_edep_keV"])])
        for k in b:
            b[k].sort()
        out["vessels"][v]["beta"] = b

    mu_tab = load_mu()
    Egrid = sorted(set(e for (_m, _r, e) in mu_tab))
    out["mu"]["E"] = Egrid
    for m in sorted(mats):
        out["mu"][m] = [float(mu_of(mu_tab, m, e)) for e in Egrid]
    p = os.path.join(rcspec.RESULTS, "curves.json")
    json.dump(out, open(p, "w", encoding="utf-8"), ensure_ascii=False)
    print("выгрузка:", p, "%.0f КБ" % (os.path.getsize(p) / 1024))


if __name__ == "__main__":
    vs = [a for a in sys.argv[1:] if a in rcspec.VESSELS] or ["m200", "m500"]
    for v in vs:
        report(v)
        print()
    if "--json" in sys.argv:
        dump_json(["m200", "m500"])
