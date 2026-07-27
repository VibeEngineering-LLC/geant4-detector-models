"""Пересчёт всех записей комплекта поверки: активность из спектра против паспорта.

Для каждой записи: A_изм = R_пика / eps_на_распад, сравнение с паспортной
активностью, пересчитанной на дату измерения. Отношение A_изм/A_пасп по всем
сорока записям — итоговая таблица качества модели.

Откуда eps_на_распад:
- объёмные одиночные нуклиды: сетка моноэнергий своей геометрии, умноженная
  на выход линии на распад родителя (из прогонов цепочек) и делённая на
  поправку суммирования C своей геометрии (из прогонов распада);
  плотность источника учитывается пересчётом f(mu*ро*d) с подогнанной d_eff;
- смеси: НАПРЯМУЮ из прогонов распада в своей геометрии (mix_*.csv) — там
  выход, суммирование и матрица уже внутри;
- точечные: сетка с конусом (деление на долю телесного угла), суммирование
  из прогонов распада на 5 см; на 25 см C масштабируется отношением телесных
  углов (C-1 пропорциональна полной эффективности партнёра по каскаду).

Поправка на наложения exp(2*тау*R) — по loading.py, тау_форм = 3 мкс.

Плато Комптона и вторичные структуры (замечание оператора): площади из
модельного спектра распада снимаются ТЕМ ЖЕ алгоритмом (окно/полки в долях
ПШПВ той же ширины в кэВ), что и из измеренного, поэтому гладкий континуум
сокращается в отношении. Сумм-пики опознаны отдельно (kit_mixture.py).

Фон: у каждой записи свой вложенный фон; имя опорного файла проверяется на
соответствие геометрии (замечание оператора: у 25 см крышка открыта и фон
другой).
"""
import glob
import math
import os
import re
import sys
from datetime import date

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402

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
TAU_SHAPE = 3.0e-6

# ---------------------------------------------------------------------------
# Справочник записей. Активность/погрешность/дата — из Note (kit_inventory).
# T12 в годах. Линии: аналитические, без перекрытий с другими линиями записи.
# density: масса из описи / номинальный объём заполнения.
T12 = {"Cs-137": 30.08, "K-40": 1.248e9, "Ra-226": 1600.0, "Th-232": 1.405e10,
       "Am-241": 432.6, "Eu-152": 13.517, "Ti-44": 59.1, "Ba-133": 10.551,
       "Cd-109": 461.9 / 365.25, "Ce-139": 137.64 / 365.25,
       "Co-57": 271.74 / 365.25, "Co-60": 5.2712, "Mn-54": 312.2 / 365.25,
       "Na-22": 2.6018, "Y-88": 106.63 / 365.25, "Zn-65": 243.93 / 365.25,
       "Bi-207": 31.55, "Th-228": 1.9116}

# Аналитические линии объёмных источников. Прогон распада — один на нуклид,
# он же даёт эффективность НА РАСПАД со всем, что внутри: выходом линии,
# каскадным суммированием, блендами и континуумом.
VLINES = {
    "Cs-137": ([661.657], "Cs137"),
    "K-40": ([1460.822], "K40"),
    "Ra-226": ([351.932, 609.32, 1120.294], "Ra226chain"),
    "Th-232": ([583.187, 911.204, 2614.511], "Th232chain"),
}

# Где лежит прогон распада: (геометрия, ключ) -> база имени файла
RUNBASE = {
    ("Marinelli_1L", "Cs137"): "decay_Cs137",
    ("Marinelli_1L", "K40"): "decay_K40",
    ("Marinelli_1L", "Ra226chain"): "chain_Ra226",
    ("Marinelli_1L", "Th232chain"): "chain_Th232",
    ("Denta_120mL", "Cs137"): "cup_denta_Cs137",
    ("Denta_120mL", "Ra226chain"): "cup_denta_Ra226chain",
    ("Denta_120mL", "Th232chain"): "cup_denta_Th232chain",
    ("Petri_60mL", "Cs137"): "cup_petri_Cs137",
    ("Petri_60mL", "Ra226chain"): "cup_petri_Ra226chain",
    ("Petri_60mL", "Th232chain"): "cup_petri_Th232chain",
}
# Плотность и матрица, при которых СЧИТАЛСЯ прогон распада
RUNRHO = {"Marinelli_1L": (1.60, "OISN16"), "Denta_120mL": (1.00, "OISN16"),
          "Petri_60mL": (1.00, "OISN16")}

# объёмные записи: (геометрия, маска файла, нуклид, A Бк/кг, dA %, дата пасп.,
#                   масса г, объём заполнения мл)
VOLUME_RECORDS = [
    ("Marinelli_1L", "*M_cs_*", "Cs-137", 1890, 5, "1997-05-30", 570, 1000),
    ("Marinelli_1L", "*M_k_*", "K-40", 2540, 10, None, 665, 1000),
    ("Marinelli_1L", "*M_ra_*", "Ra-226", 1850, 10, None, 622, 1000),
    ("Marinelli_1L", "*Th232*", "Th-232", 1940, 6, "2007-09-17", 1600, 1000),
    ("Denta_120mL", "*Cs137*", "Cs-137", 1760, 5, "2002-05-24", 68, 120),
    ("Denta_120mL", "*K40*", "K-40", 2530, 6, "2002-05-24", 79, 120),
    ("Denta_120mL", "*Ra226*", "Ra-226", 1780, 6, "2002-05-24", 74, 120),
    ("Denta_120mL", "*Th232*", "Th-232", 1940, 6, "2007-09-17", 192, 120),
    ("Petri_60mL", "*Cs137*", "Cs-137", 1760, 5, "2002-05-24", 34, 60),
    ("Petri_60mL", "*K40*", "K-40", 2530, 6, "2002-05-24", 40, 60),
    ("Petri_60mL", "*Ra226*", "Ra-226", 1780, 6, "2002-05-24", 37, 60),
    ("Petri_60mL", "*Th232*", "Th-232", 1940, 6, "2007-09-17", 96, 60),
]

# сетки моноэнергий: геометрия -> (метка сетки, её плотность, d_eff мм)
# d_eff денты/петри подгоняются отдельно (пары 0,60/1,60), пока — геометрическая
# толщина слоя как начальное приближение, уточняется в selfabs по кюветам.
# d_eff — ПОДОГНАННЫЕ значения (selfabs_fit.py), не начальные догадки:
# маринелли 31,5; дента 27,0; петри 12,5 мм.
GRIDS = {"Marinelli_1L": ("rho1.60", 1.60, 31.5),
         "Denta_120mL": ("denta1.60", 1.60, 27.0),
         "Petri_60mL": ("petri1.60", 1.60, 12.5)}

# поправки суммирования по геометриям: файл прогона -> спектр распада
CUP_TAG = {"Marinelli_1L": "decay_%s.csv", "Denta_120mL": "cup_denta_%s.csv",
           "Petri_60mL": "cup_petri_%s.csv"}
CHAIN_ALIAS = {"Cs137": "Cs137", "K40": None,          # K-40 без каскада: C=1
               "Ra226chain": "Ra226chain", "Th232chain": "Th232chain"}


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


def mu_of(fn):
    tab = {}
    for line in open(os.path.join(BUILD, fn), encoding="utf-8"):
        if line.startswith("#") or line.startswith("E_keV"):
            continue
        e, m = line.split(",")
        tab[round(float(e), 3)] = float(m)
    return tab


MU_O = mu_of("mu_oisn16.csv")
MU_W = mu_of("mu_water.csv")

# Разрешение прибора: измерено по пику 662 записи цезия (7,5 %); паспорт ≤8 %.
FWHM662 = 49.9


def fx(x):
    return (1 - math.exp(-x)) / x if x > 1e-9 else 1.0


_BROAD = {}


def area_sim(hist, E, fwhm=None, key=None, win=6.0, bg0=30.0, bg1=10.0):
    """Площадь пика в модельном спектре.

    Если задана ПШПВ — спектр уширяется до разрешения прибора и площадь
    берётся ТЕМ ЖЕ окном (±1 ПШПВ + полки), что и в измерении. Это
    обязательно там, где линии сливаются: Ac-228 911,2 + 968,97 в NaI —
    один пик, и узкое окно по модели давало завышение активности тория
    по этой линии в полтора раза.
    Без ПШПВ — прежнее узкое окно (для сеток моноэнергий, где блендов нет).
    """
    if fwhm:
        if key not in _BROAD:
            _BROAD[key] = bm.broaden(hist)
        a, _ = bm.area_broadened(_BROAD[key], E, fwhm)
        return a
    gross = sum(c for e, c in hist.items() if abs(e - E) <= win)
    side = sum(c for e, c in hist.items() if E - bg0 <= e <= E - bg1)
    return gross - side / (bg0 - bg1) * (2 * win + 1)


def eps_mono(tag, E):
    p = os.path.join(BUILD, "grid", "%s_E%07.1f.csv" % (tag, round(E, 1)))
    if not os.path.exists(p):
        cand = glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv"))
        best = None
        for c in cand:
            m = re.search(r"_E(\d+\.\d)\.csv$", c)
            if m and abs(float(m.group(1)) - E) < 1.0:
                best = c
        if not best:
            return None
        p = best
    hist, N = load_hist(p)
    return area_sim(hist, E) / N


def yield_line(fn, E):
    p = os.path.join(BUILD, fn)
    if not os.path.exists(p):
        return None
    emit, N = load_hist(p)
    tot = sum(c for e, c in emit.items() if abs(e - E) <= 2.0)
    return tot / N if N else None


def C_of(geom, chain_key, E):
    """Поправка суммирования: eps_моно(сетка 1,6) / eps_на_распад(1,0...)

    ВАЖНО: у объёмных геометрий прогоны распада сделаны при ро=1,0, сетка при
    1,6 — из отношения плотностной множитель НЕ сокращается. Поэтому C здесь
    строится как отношение eps_расп(контроль Cs)/eps_расп(линии), приведённое
    к контролю: C = [eps_моно(E)/eps_моно(662)] / [eps_расп(E)/eps_расп(662)]
    — двойное отношение, из которого выпадают и плотность, и нормировка.
    Затем абсолют восстанавливается контролем C(Cs)=1.
    """
    alias = CHAIN_ALIAS.get(chain_key)
    if alias is None:
        return 1.0
    tagf = CUP_TAG[geom]
    sp = os.path.join(BUILD, tagf % alias)
    csp = os.path.join(BUILD, tagf % "Cs137")
    if not (os.path.exists(sp) and os.path.exists(csp)):
        return None
    hist, N = load_hist(sp)
    emitf = ("decay_%s_emit.csv" if geom == "Marinelli_1L"
             else tagf.replace(".csv", "_emit.csv"))
    em, Ne = load_hist(os.path.join(BUILD, (tagf % alias).replace(".csv", "_emit.csv")))
    ch, Nc = load_hist(csp)
    ec, _ = load_hist(os.path.join(BUILD, (tagf % "Cs137").replace(".csv", "_emit.csv")))

    def eps_dec(h, e, n, E0, tagkey):
        nem = sum(c for x, c in e.items() if abs(x - E0) <= 2.0)
        if nem <= 200:
            return None
        # ВАЖНО: площадь из спектра распада берётся уширенным окном — иначе
        # бленды (911+969 у Ac-228) считаются иначе, чем в измерении
        return area_sim(h, E0, fwhm=FWHM662 * math.sqrt(E0 / 661.657),
                        key=tagkey) / nem

    ed = eps_dec(hist, em, N, E, "dec:%s:%s" % (geom, alias))
    edc = eps_dec(ch, ec, Nc, 661.657, "dec:%s:Cs137" % geom)
    m = eps_mono(GRIDS[geom][0], E)
    mc = eps_mono(GRIDS[geom][0], 661.657)
    if None in (ed, edc, m, mc):
        return None
    return (m / mc) / (ed / edc)


def decay_factor(nuc, d0, d1):
    if not d0 or not d1:
        return 1.0
    y0, m0, dd0 = (int(x) for x in d0.split("-"))
    y1, m1, dd1 = (int(x) for x in d1.split("-"))
    dt = (date(y1, m1, dd1) - date(y0, m0, dd0)).days / 365.25
    return 0.5 ** (dt / T12[nuc])


def eps_per_decay(geom, ckey, E, fwhm, rho_src, mu_src):
    """Эффективность НА РАСПАД РОДИТЕЛЯ в окне ±1 ПШПВ с полками.

    Берётся из прогона полного распада: внутри уже выход линии, каскадное
    суммирование, бленды соседних линий и континуум под пиком. Ничего из
    справочника и никаких перемножений руками.

    Плотность прогона и источника разные, поэтому вводится отношение
    поправок самопоглощения f(mu*ro*d). Матрица «лёгких» источников
    (68–79 г вместо 192) в файлах НЕ записана — берётся вода как
    представитель лёгкой среды; ОИСН-16 здесь неприменима, в ней 71 % железа.
    """
    base = RUNBASE.get((geom, ckey))
    if not base:
        return None
    p = os.path.join(BUILD, base + ".csv")
    if not os.path.exists(p):
        return None
    hist, N = load_hist(p)
    a = area_sim(hist, E, fwhm=fwhm, key="dec:" + base)
    if a <= 0 or not N:
        return None
    rho_run, mat_run = RUNRHO[geom]
    dmm = GRIDS[geom][2]
    key = min(MU_O, key=lambda k: abs(k - E))
    mu_run = MU_O[key] if mat_run == "OISN16" else MU_W[key]
    corr = fx(mu_src * rho_src * dmm / 10) / fx(mu_run * rho_run * dmm / 10)
    return (a / N) * corr


if __name__ == "__main__":
    print("Пересчёт объёмных записей комплекта.\n"
          "eps на распад берётся из прогона распада УШИРЕННЫМ окном ±1 ПШПВ —\n"
          "тем же, что в измерении: так бленды (911+969 у Ac-228), плато\n"
          "Комптона и сумм-пики учитываются одинаково с двух сторон.\n")
    print("%-13s %-8s %9s %8s %11s %9s %8s" %
          ("геометрия", "нуклид", "E, кэВ", "имп/с", "eps/распад",
           "A, Бк/кг", "A/пасп"))
    missing, rows = [], []
    for geom, mask, nuc, aspec, dpct, d0, mass, vol in VOLUME_RECORDS:
        # Раскладка комплекта — reference_spectra/reference_kits*/<геометрия>/
        # <нуклид>/sample_*.xml, а не <геометрия>/* прямо в корне эталонов.
        # Стоял плоский glob по корню, он не находил НИЧЕГО, и скрипт молча
        # печатал пустую таблицу: ни одной записи, ни одной жалобы. Каталог
        # ищет paths.kit_dir(), он же выбирает нужный формат.
        kd = paths.kit_dir(geom)
        files = sorted(str(p) for p in kd.rglob(mask)) if kd else []
        if not files:
            missing.append("%s / %s (%s)" % (geom, mask, nuc))
            continue
        s, b = bm.read(files[0])
        # дата измерения и опорный фон
        txt = open(files[0], encoding="utf-8", errors="replace").read()
        mdate = re.search(r"<StartTime>(\d{4}-\d{2}-\d{2})", txt)
        mdate = mdate.group(1) if mdate else None
        bgref = re.search(r"<BackgroundSpectrumFile>([^<]*)</", txt)
        bgref = os.path.basename(bgref.group(1)) if bgref else "?"
        rho = mass / vol
        gtag, grho, dmm = GRIDS[geom]
        A0 = aspec * mass / 1000.0 * decay_factor(nuc, d0, mdate)
        R = float(s.n.sum()) / s.live
        pile = math.exp(2 * TAU_SHAPE * R)
        lines, ckey = VLINES[nuc]
        # ПШПВ прибора — единый закон по всем записям, чтобы окна модели и
        # измерения совпадали в точности (по своему пику мерить нельзя:
        # у слабых линий он не находится)
        for E in lines:
            fw = FWHM662 * math.sqrt(E / 661.657)
            r = bm.net_rate(s, b, E, fw, roi=1.0, side=1.0)
            if r is None or r[0] <= 0:
                continue
            rate = r[0] * pile
            key = min(MU_O, key=lambda k: abs(k - E))
            # матрица источника: при ро≈1,6 это ОИСН-16 (192/96/1600 г —
            # ровно калибровочная засыпка), у лёгких состав не записан -> вода
            mu_src = MU_O[key] if rho > 1.3 else MU_W[key]
            eps = eps_per_decay(geom, ckey, E, fw, rho, mu_src)
            if not eps:
                continue
            A = rate / eps
            print("%-13s %-8s %9.1f %8.3f %11.4e %9.1f %8.3f"
                  % (geom, nuc, E, rate, eps, A / (mass / 1000.0), A / A0))
            rows.append((geom, nuc, E, rate, eps, A / (mass / 1000.0),
                         aspec, A / A0, dpct, RUNBASE[(geom, ckey)]))
        # Соответствие фона геометрии. У Денты и Петри опорный фон —
        # «пустая защита»: кювета мала, водяной маринелльной защиты нет, и
        # empty_shield — правильный выбор (в имени файла стоит «point5cm»
        # только потому, что тот же сеанс обслуживал точечную геометрию).
        # Настоящая ошибка была бы: маринелльная запись с фоном пустой защиты
        # или наоборот.
        want = "marinelli" if geom == "Marinelli_1L" else "empty_shield"
        bad = (want not in bgref) and ("Фон закр" not in bgref)
        print("   фон: %s%s" % (bgref, "   <- НЕ СООТВЕТСТВУЕТ ГЕОМЕТРИИ" if bad
                                else ""))

    # Числа — файлом, а не только в консоли (пункт 12 протокола проверок).
    if rows:
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "kit_recalc_volume.csv")
        out = os.path.abspath(out)
        with open(out, "w", encoding="utf-8", newline="") as fh:
            fh.write("# Пересчёт объёмных записей комплекта поверки.\n"
                     "# A_изм = R_пика / eps_на_распад; A_пасп приведена на "
                     "дату измерения.\n"
                     "# ratio = A_изм / A_пасп: меньше единицы — модель "
                     "ЗАВЫШАЕТ эффективность.\n"
                     "# run — прогон распада, из которого взята "
                     "eps_на_распад.\n")
            fh.write("geometry,nuclide,E_keV,rate_cps,eps_per_decay,"
                     "A_meas_Bq_kg,A_pass_Bq_kg,ratio,d_pass_pct,run\n")
            for r in rows:
                fh.write("%s,%s,%.3f,%.4f,%.6e,%.1f,%.0f,%.4f,%.1f,%s\n" % r)
        print("\nтаблица: %s (%d строк)" % (out, len(rows)))
        for g in ("Marinelli_1L", "Denta_120mL", "Petri_60mL"):
            rr = [r[7] for r in rows if r[0] == g]
            if rr:
                rr.sort()
                print("   %-13s строк %2d, медиана A/пасп %.3f, "
                      "разброс %.3f..%.3f"
                      % (g, len(rr), rr[len(rr) // 2], rr[0], rr[-1]))

    # Ненайденные записи — вслух. Пустая таблица не должна выглядеть как
    # «всё сошлось»: именно так этот пересчёт и простоял незамеченным.
    if missing:
        print("\nНЕ НАЙДЕНЫ спектры для записей (%d):" % len(missing))
        for m in missing:
            print("   ", m)
        print("Комплект ищется через paths.kit_dir(); проверьте G4MODELS_REF\n"
              "или наличие reference_spectra/reference_kits* в эталонах.")
    if len(missing) == len(VOLUME_RECORDS):
        raise SystemExit("Не найдено НИ ОДНОЙ записи комплекта — считать нечего.")
