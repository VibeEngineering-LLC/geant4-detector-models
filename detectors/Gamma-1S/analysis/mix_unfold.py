"""Разложение спектров смеси №SRC-04 на составляющие (задание оператора).

Метод — полноспектральная подгонка: измеренный спектр описывается линейной
комбинацией МОДЕЛЬНЫХ спектров распада компонентов (Am-241, Ti-44/Sc-44,
Eu-152, Cs-137 в той же геометрии, прогоны mix_<геом>_<нукл>.csv), уширенных
разрешением прибора, плюс ИЗМЕРЕННЫЙ фон, приведённый по живому времени.
Амплитуды неотрицательны (NNLS), веса пуассоновские.

Почему шаблоны, а не набор гауссиан. В области 40–110 кэВ налагаются
Am-241 59,5, Ti-44 67,9 и 78,3 И РЕНТГЕН СВИНЦА защиты (Kα 72,8/75,0,
Kβ 84,9/87,3, ослабленный вкладышем Cd+Cu примерно в пять раз), а ПШПВ там
около 15 кэВ — все линии в пределах полутора ПШПВ. Руками это не разделить.
Зато модельный шаблон содержит и линии, и флуоресценцию свинца, и
самопоглощение в матрице — потому что защита и проба стоят в той же геометрии.
Поэтому мягкая область не выносится в отдельную процедуру, а входит в общую
подгонку.

Энергошкала записей смесей сползла (Cs читается как 657 вместо 661,7).
Поэтому шкала — ПАРАМЕТР подгонки: внешний цикл по (a, b, ПШПВ662) поверх
E' = a + b*E_файла, внутри — NNLS по амплитудам. Так соблюдён порядок: сперва
шкала, потом состав, а не наоборот.

Что это даёт сверх попикового пересчёта:
- работают континуумы и СУММ-ПИКИ (511+1157 у Sc-44) — они в шаблонах есть;
- перекрытия разделяются формой, а не назначением «пик -> нуклид»;
- активности Ti-44 и Eu-152 получаются ВПЕРВЫЕ: в файлах комплекта их нет,
  паспорт называет только Am-241 (Маринелли) и Cs-137 (Дента, Петри).
"""
import glob
import math
import os
import sys

import numpy as np
from scipy.optimize import nnls, minimize

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402
from kit_mixture import find_peaks, selfcalib  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)
KIT = str(paths.ref("Gamma-1S"))

NUCS = ["Am241", "Ti44", "Eu152", "Cs137"]
GEOM_FILES = {
    "marinelli": ("Marinelli_1L", "*AmTiCsEu*.xml", 1.000),
    "denta": ("Denta_120mL", "*Am-Ti*.xml", 0.100),
    "petri": ("Petri_60mL", "*Am-Ti*.xml", 0.060),
}
EFIT0, EFIT1 = 45.0, 1900.0      # диапазон подгонки, кэВ
FINE = 0.25                       # шаг тонкой сетки шаблонов, кэВ

# Систематический пол весов. При миллионе отсчётов чисто пуассоновские веса
# отдают подгонку нескольким самым сильным пикам, а любая трёхпроцентная
# неточность формы (разрешение, хвост пика, состав матрицы) раздувает
# chi2/dof до сотен и делает его бессмысленным как мера качества.
# sigma^2 = y + фон + (SYS * модель)^2 — стандартный приём.
SYS = 0.03


def load_hist(path):
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


def fine_template(geom, nuc, egrid, w662):
    """Уширенный модельный спектр на тонкой сетке, отсчётов на ОДИН распад."""
    p = os.path.join(BUILD, "mix_%s_%s.csv" % (geom, nuc))
    if not os.path.exists(p):
        return None
    hist, N = load_hist(p)
    out = np.zeros_like(egrid)
    erf = np.vectorize(math.erf)
    for E0, c in hist.items():
        if c == 0:
            continue
        sig = w662 * math.sqrt(max(E0, 8.0) / 661.657) / 2.3548
        lo = np.searchsorted(egrid, E0 - 5 * sig)
        hi = np.searchsorted(egrid, E0 + 5 * sig)
        if hi <= lo:
            continue
        z = (egrid[lo:hi] - E0) / (sig * math.sqrt(2.0))
        g = np.exp(-0.5 * ((egrid[lo:hi] - E0) / sig) ** 2)
        s = g.sum()
        if s > 0:
            out[lo:hi] += c * g / s
    return out / N


def rebin(cum, egrid, edges):
    """Интеграл тонкого шаблона по каналам записи (границы в кэВ)."""
    v = np.interp(edges, egrid, cum, left=cum[0], right=cum[-1])
    return np.diff(v)


def unfold(geom, verbose=True):
    kdir, mask, mass_kg = GEOM_FILES[geom]
    files = glob.glob(os.path.join(KIT, kdir, mask))
    if not files:
        print("нет записи смеси для", geom)
        return None
    s, b = bm.read(files[0])

    # начальное приближение шкалы — по опорным линиям
    pk = find_peaks(s, b)
    a0, k0, anchors = selfcalib(pk)
    w0 = 0.075 * 661.657

    ch = np.arange(len(s.n) + 1, dtype=float) - 0.5
    Efile = s.energy(ch)                      # границы каналов по файлу
    y = s.n.astype(float)
    bgE = b.energy(np.arange(len(b.n), dtype=float))
    tscale = s.live / b.live

    egrid = np.arange(0.0, 3100.0, FINE)
    cache = {}

    def build(w662):
        key = round(w662, 3)
        if key not in cache:
            cols = []
            for nuc in NUCS:
                t = fine_template(geom, nuc, egrid, w662)
                if t is not None:
                    cols.append((nuc, np.cumsum(t)))
            cache[key] = cols
        return cache[key]

    def objective(par):
        a, k, w = par
        if not (0.9 < k < 1.15 and abs(a) < 60 and 20 < w < 90):
            return 1e9
        edges = a + k * Efile
        mids = 0.5 * (edges[:-1] + edges[1:])
        sel = (mids >= EFIT0) & (mids <= EFIT1)
        if sel.sum() < 50:
            return 1e9
        bg = np.interp(mids, bgE, b.n) * tscale
        cols = build(w)
        if not cols:
            return 1e9
        A = np.stack([rebin(c, egrid, edges)[sel] for _, c in cols], axis=1)
        net = y[sel] - bg[sel]
        # два шага: сперва статистические веса, потом с систематическим полом
        # по полученной модели — иначе пол пришлось бы брать от данных, а он
        # должен отражать неточность МОДЕЛИ
        sig = np.sqrt(np.maximum(y[sel] + bg[sel], 1.0))
        coef, _ = nnls(A / sig[:, None], net / sig)
        mod = A @ coef
        sig = np.sqrt(np.maximum(y[sel] + bg[sel], 1.0) + (SYS * mod) ** 2)
        coef, _ = nnls(A / sig[:, None], net / sig)
        r = (A @ coef - net) / sig
        return float((r ** 2).sum()) / max(1, sel.sum() - len(coef) - 3)

    # Полнота шаблонов — проверять ДО подгонки: при отсутствии части прогонов
    # оптимизатор уходит на границу параметра и обнуляет компоненты, выдавая
    # правдоподобные с виду, но бессмысленные амплитуды.
    have = [n for n in NUCS
            if os.path.exists(os.path.join(BUILD, "mix_%s_%s.csv" % (geom, n)))]
    if len(have) < len(NUCS):
        print("\n===== %s: ШАБЛОНЫ НЕПОЛНЫ (%s из %s) — подгонка не имеет "
              "смысла, прогоны распада не готовы"
              % (geom, ", ".join(have) or "нет", ", ".join(NUCS)))
        return None

    best = minimize(objective, [a0, k0, w0], method="Nelder-Mead",
                    options={"maxiter": 400, "xatol": 0.02, "fatol": 0.01})
    a, k, w = best.x
    chi = best.fun
    # Уход на границу — признак негодной подгонки, а не результата
    at_bound = (abs(k - 0.9) < 1e-3 or abs(k - 1.15) < 1e-3
                or abs(w - 20) < 1e-2 or abs(w - 90) < 1e-2)
    if at_bound:
        print("\n===== %s: подгонка ушла на ГРАНИЦУ параметра "
              "(k=%.5f, ПШПВ=%.1f) — результат отброшен" % (geom, k, w))
        return None

    # финальные амплитуды при найденной шкале
    edges = a + k * Efile
    mids = 0.5 * (edges[:-1] + edges[1:])
    sel = (mids >= EFIT0) & (mids <= EFIT1)
    bg = np.interp(mids, bgE, b.n) * tscale
    cols = build(w)
    A = np.stack([rebin(c, egrid, edges)[sel] for _, c in cols], axis=1)
    net = y[sel] - bg[sel]
    sig = np.sqrt(np.maximum(y[sel] + bg[sel], 1.0))
    coef, _ = nnls(A / sig[:, None], net / sig)
    mod = A @ coef
    sig = np.sqrt(np.maximum(y[sel] + bg[sel], 1.0) + (SYS * mod) ** 2)
    coef, _ = nnls(A / sig[:, None], net / sig)

    res = {}
    if verbose:
        print("\n===== %s : %s" % (geom, os.path.basename(files[0])))
        print("шкала E' = %+.2f %+.5f*E ; ПШПВ(662) = %.1f кэВ (%.1f %%) ;"
              " chi2/dof = %.2f ; каналов в подгонке %d"
              % (a, k, w, 100 * w / 661.657, chi, sel.sum()))
        print("%-8s %14s %14s" % ("нуклид", "распадов/с", "Бк/кг"))
    for (nuc, _), c in zip(cols, coef):
        act = c / s.live
        res[nuc] = act / mass_kg
        if verbose:
            print("%-8s %14.1f %14.0f" % (nuc, act, act / mass_kg))
    return res


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    total = {}
    for geom in GEOM_FILES:
        if only and only != geom:
            continue
        r = unfold(geom)
        if r:
            total[geom] = r

    if total:
        print("\n--- сверка с паспортом (дата измерения 2016-05-30/31) ---")
        # Am-241 4200 Бк/кг (10 %) на 31-05-2002, T1/2 432,6 г -> x0,9778
        # Cs-137 2210 Бк/кг (5 %) на 31-05-2002, T1/2 30,08 г  -> x0,7242
        if "marinelli" in total and "Am241" in total["marinelli"]:
            ref = 4200 * 0.9778
            got = total["marinelli"]["Am241"]
            print("Am-241 маринелли: %5.0f против %5.0f Бк/кг  -> %.3f"
                  % (got, ref, got / ref))
        for g in ("denta", "petri"):
            if g in total and "Cs137" in total[g]:
                ref = 2210 * 0.7242
                got = total[g]["Cs137"]
                print("Cs-137 %-9s: %5.0f против %5.0f Бк/кг  -> %.3f"
                      % (g, got, ref, got / ref))
        print("\nTi-44 и Eu-152 паспортом не заданы — это новый результат;")
        print("их согласованность между тремя геометриями и есть проверка.")
        for nuc in ("Ti44", "Eu152"):
            v = [(g, total[g][nuc]) for g in total if nuc in total[g]]
            if len(v) > 1:
                m = sum(x for _, x in v) / len(v)
                sp = ", ".join("%s %.0f" % (g, x) for g, x in v)
                print("   %-6s %s Бк/кг ; среднее %.0f, разброс %.0f %%"
                      % (nuc, sp, m, 100 * (max(x for _, x in v)
                                            - min(x for _, x in v)) / m))
