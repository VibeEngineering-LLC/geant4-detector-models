"""Мост методик на чистых 4π-моно точечной 5 см (задача 115).

ЧТО ИЗМЕРЯЕТСЯ. Расхождение эффективности модели со штатной кривой прибора
имеет ДВА возможных источника, и до их разделения диагноз невозможен:

  (а) физика — модель поглощает не так, как настоящий кристалл;
  (б) конвенция — «площадь пика» у модели и у аттестации определены
      по-разному. Модель берёт полнопоглощённый пик в узком окне ОСТРОГО
      депозит-спектра; аттестация берёт площадь гаусс-фита по РАЗМЫТОМУ
      спектру, где широкий пик сидит на континууме неполного поглощения и фит
      часть его теряет.

Мост измеряет ровно (б), на одних и тех же событиях:

    B(E) = площадь фита по размытому спектру / истинный пик депозит-спектра

B < 1 означает, что аттестованная конвенция ЗАНИЖАЕТ площадь. Если B(E)
ПЛОСКАЯ, поправка съедает плато расхождения и не трогает разрыв на жёстком
крае; если B(E) ПАДАЕТ с ростом энергии — наоборот, снимает часть разрыва.
Различить эти два ответа и есть смысл задачи; от него зависят задачи 93, 109
и 122.

ПОЧЕМУ ЭТОТ ФАЙЛ, А НЕ `bridge_p5_th228.py`. Тот считался на моно-спектрах
`scat_p5_full_*` от 28.07, снятых ДО правки входного торца (коммит 1a29dee):
прямое A/B на 122,1 кэВ дало +17,8 %. Мост стоял на СМЕСИ ДВУХ ГЕОМЕТРИЙ, и
его числа отозваны (задачи 126, 115). Здесь входы — прогон `bridge_mono_v4.mac`
одним exe с отпечатком провенанса в шапке каждого файла; вердикт входов
пишется в шапку таблицы, так что подмена набора обнаружится сторожем, а не
чтением дат.

НАБОР ЛИНИЙ. Шесть линий, на которых приборная ПШПВ ИЗМЕРЕНА напрямую
(238,632 / 300,087 / 583,187 / 727,330 / 860,557 / 2614,511 — сеанс
30.07.2026), плюс 80,998 / 122,100 / 661,657 для протяжённости диапазона. На
мягком крае результат читать с оговоркой: окно полки конвенции `E−30…E−10`
там отстоит от центра меньше чем на ПШПВ, то есть лежит внутри пика, и в
размытом спектре съедает крыло. Решающими считать 662 и выше.

ШИРИНА РАЗМЫТИЯ — `detector_params.fwhm_measured`, а не `fwhm`: сравнение
идёт с площадями ТОГО ЖЕ сеанса, которым измерены ширины. Оговорка о
десятипроцентном расхождении двух оценок ширины на 662 кэВ — в docstring
`fwhm_measured` и в задаче 138.

ОГОВОРКА О ФОНЕ. Подложка под пиком — ступенька-из-образа плюс полином
степени `BG_DEGREE` = 2, как у прибора (BackgrPower = 2 в слепке lsrm.cnf);
прежняя линейная версия отозвана аудитом 31.07.2026 (P0-1: линейная подложка
не держала комптоновский край внутри зоны и создавала обрыв B сама). Высота
ступеньки не задаётся числом, а подгоняется наравне с остальными
параметрами. Чувствительность к степени полинома и ширине зоны — таблицей
`results/bridge_sensitivity.csv`.

Спектры прогона закоммичены в `results/bridge_spectra/` (девять файлов,
~64 КБ) вместе с макросом `macros/bridge_mono_v4.mac`: мост воспроизводим из
чистого клона без пересчёта Geant4 (аудит 31.07.2026, сверх ТЗ-5). Рабочий
каталог модели, если содержит спектры, имеет приоритет — свежий прогон
перекрывает копию в репозитории.
"""
import glob
import math
import os
import re
import sys

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import detector_params as dp  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
OUT = str(paths.results("Gamma-1S"))
PREFIX = "bridge4pi_E"

# Окно истинного пика в ОСТРОМ депозит-спектре. Шаг гистограммы 1 кэВ, пик
# занимает единицы каналов; ±3 кэВ ловит его целиком вместе с утечкой в
# соседний канал и не захватывает вылет K-рентгена иода (28,6 кэВ).
TRUE_HALF = 3.0
# Полуширина зоны фита в единицах ПШПВ. 2,5 — границы пика ПРИБОРА
# (LeftPeakBound = RightPeakBound = 2.5 в lsrm.cnf, см. peak_to_total.py):
# мост воспроизводит конвенцию аттестации, а не собственную. Прежние 3,2 были
# унаследованы от отозванного bridge_p5_th228.py без вывода, и именно лишние
# 28 % ширины затягивали комптоновский край 2614,5 (2381,8 кэВ) на 114 кэВ
# ВНУТРЬ зоны — линейная подложка с одной ступенькой в центре пика описать
# его не может, фит схлопывал ширину (ws = 0,874) и обрыв B в значительной
# части создавался самим мостом (внутренний аудит 31.07.2026, P0-1/P1-2).
# При 2,5 ПШПВ край на 2614,5 лежит на границе зоны (2343,6), а не внутри.
FIT_ZONE_FWHM = 2.5
# Степень полинома подложки. 2 — конвенция прибора (BackgrPower = 2 в
# слепке); прежняя линейная подложка была вторым отступлением от конвенции.
BG_DEGREE = 2

# Отношение штатная/наша по линиям Th-228 — из results/bridge_attested_ratio.csv
# (задача 151, хвост P1-4): посчитано АНАЛИТИЧЕСКИ по полиномам зон .efa
# (efa_zones), не снято программой и не привязано к дате сеанса СпектраЛайн —
# обновляется автоматически при каждом пересчёте нашей кривой. Раньше здесь
# был словарь трёх значащих цифр, переписанный руками из консоли
# hard_edge_th228.py; см. docstring bridge_attested_ratio.py про замену.
_ATTESTED_PATH = os.path.join(OUT, "bridge_attested_ratio.csv")


def _load_attested():
    if not os.path.exists(_ATTESTED_PATH):
        raise SystemExit(
            "Нет %s — прогоните analysis/bridge_attested_ratio.py"
            " (нужен G1S_LSRM_MASTER_EFA)." % _ATTESTED_PATH)
    out = {}
    for r in csvio.read(_ATTESTED_PATH):
        out[float(r["E_keV"])] = (float(r["ratio_attested_over_ours"]),
                                  float(r["d_ratio_pct"]))
    return out


ATTESTED_RATIO = _load_attested()


def attested_ratio(E, tol=0.05):
    """(отношение штатная/наша; погрешность, %) для линии E; None, если её нет.

    Поиск по БЛИЗОСТИ, а не по точному ключу. Энергия в шапке прогона приходит
    напечатанной с четырьмя знаками (2614,5100), а в таблице аттестации стоит
    2614,511 — точный `dict.get` не находил решающую линию и МОЛЧА выбрасывал
    её из сводки, после чего остаток по трём оставшимся линиям выглядел
    плоским просто потому, что жёсткого края в нём не было.
    """
    for k, v in ATTESTED_RATIO.items():
        if abs(k - E) <= tol:
            return v
    return None


OBS = {
    "quantity": "B(E) — отношение площади гаусс-фита по РАЗМЫТОМУ спектру к"
                " истинному полнопоглощённому пику того же депозит-спектра",
    "area": "числитель — площадь гауссианы фита (ступенька и подложка в"
            " площадь не входят); знаменатель — сумма отсчётов в окне"
            " +-%.1f кэВ депозит-спектра" % TRUE_HALF,
    "window": "фит в зоне +-%.1f ПШПВ (границы пика прибора из lsrm.cnf);"
              " истинный пик в окне +-%.1f кэВ" % (FIT_ZONE_FWHM, TRUE_HALF),
    "shelf": "в числителе подложка — полином степени %d (BackgrPower"
             " прибора) плюс ступенька-из-образа; всё подгоняется;"
             " в знаменателе подложка не вычитается" % BG_DEGREE,
    "blurred": "числитель ДА (приборная ПШПВ по шести измеренным линиям);"
               " знаменатель НЕТ",
    "geometry": "точечный источник на оси z = 91 мм (5 см от наружной"
                " поверхности корпуса); полный 4pi без конуса",
}


def load(path):
    """[(E; отсчёты)] из расчётного спектра; строки шапки пропускаются."""
    out = []
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#") or ln.startswith("E_keV") or not ln.strip():
            continue
        e, c = ln.split(",")
        out.append((float(e), float(c)))
    return out


def true_peak(ec, E0, half=TRUE_HALF):
    return sum(c for e, c in ec if abs(e - E0) < half)


def blur(ec, emax, fwhm_of=None):
    """Свёртка острого спектра с приборным разрешением.

    `emax` приходит аргументом и берётся по САМОЙ ВЕРХНЕЙ энергии набора:
    зашитая константа обрезала бы спектр, а обрезанный спектр молча теряет
    правое крыло пика — тот же класс дефекта, что и потерянные верхние узлы
    сетки (см. `export_curves.curve`). `fwhm_of` — закон ширины; вынесен
    аргументом для таблицы чувствительности (задача 138: две оценки ширины
    расходятся на 10 % и не сведены).
    """
    fwhm_of = dp.fwhm_measured if fwhm_of is None else fwhm_of
    src = np.zeros(emax)
    for e, c in ec:
        i = int(e)
        if 0 <= i < emax:
            src[i] += c
    xs = np.arange(emax, dtype=float)
    out = np.zeros(emax)
    for i in np.nonzero(src)[0]:
        s = fwhm_of(max(xs[i], 10.0)) / 2.3548
        lo, hi = max(0, int(xs[i] - 6 * s)), min(emax, int(xs[i] + 6 * s))
        g = np.exp(-0.5 * ((xs[lo:hi] - xs[i]) / s) ** 2)
        out[lo:hi] += src[i] * g / (s * math.sqrt(2 * math.pi))
    return xs, out


def fit_area(xs, ys, E0, zone=None, bg_degree=None, fwhm_of=None):
    """Площадь гауссианы фита; ступенька и подложка подгоняются вместе с ней.

    Возвращает (площадь; σ площади из ковариации; относительная ширина фита
    к приборной; высота ступеньки). Высота ступеньки — СВОБОДНЫЙ параметр:
    прежде она задавалась тремя числами, подобранными руками по трём линиям,
    и на новых линиях брать было неоткуда.

    `zone`, `bg_degree`, `fwhm_of` вынесены аргументами ради таблицы
    чувствительности: вывод моста обязан переживать выбор зоны фита, степени
    подложки и закона ширины, иначе он принадлежит этим выборам, а не данным
    (внутренний аудит 31.07.2026).
    """
    zone = FIT_ZONE_FWHM if zone is None else zone
    bg_degree = BG_DEGREE if bg_degree is None else bg_degree
    fwhm_of = dp.fwhm_measured if fwhm_of is None else fwhm_of
    fw = fwhm_of(E0)
    s0 = fw / 2.3548
    m = (xs >= E0 - zone * fw) & (xs <= E0 + zone * fw)
    x, y = xs[m], ys[m]

    def model(xx, A, sh, ws, hstep, b0, b1, b2):
        s = ws * s0
        mu = E0 + sh
        amp = A / (s * math.sqrt(2 * math.pi))
        gg = amp * np.exp(-0.5 * ((xx - mu) / s) ** 2)
        st = amp * hstep * 0.5 * erfc((xx - mu) / (s * math.sqrt(2)))
        t = xx - xx.mean()
        bg = b0 + b1 * t + (b2 * t * t if bg_degree >= 2 else 0.0 * t)
        return bg + gg + st

    # При линейной подложке b2 зажимается в ноль границами, а не выкинут из
    # модели: число параметров одно, ковариации сравнимы между вариантами.
    b2lim = 1e-12 if bg_degree < 2 else np.inf
    p, cov = curve_fit(
        model, x, y,
        p0=[y.sum(), 0.0, 1.0, 1e-3, float(np.median(y)), 0.0, 0.0],
        sigma=np.sqrt(np.maximum(y, 1.0)),
        bounds=([0, -0.5 * fw, 0.5, 0.0, -np.inf, -np.inf, -b2lim],
                [np.inf, 0.5 * fw, 2.0, 0.2, np.inf, np.inf, b2lim]),
        maxfev=40000)
    dA = float(math.sqrt(max(cov[0, 0], 0.0)))
    return p[0], dA, p[2], p[3]


def inputs():
    """Файлы прогона в порядке возрастания энергии; энергия — из ИМЕНИ.

    Имя даёт энергию с одним знаком после запятой, а нужна точная — она
    читается из шапки (`E_prim_keV`), имя служит только для сортировки.

    Сначала рабочий каталог модели (свежий прогон), при его пустоте —
    копия спектров, закоммиченная в `results/bridge_spectra/`: мост
    воспроизводим из чистого клона (сверх ТЗ-5).
    """
    src = BUILD
    if not glob.glob(os.path.join(src, PREFIX + "*.csv")):
        src = os.path.join(OUT, "bridge_spectra")
    out = []
    for p in sorted(glob.glob(os.path.join(src, PREFIX + "*.csv"))):
        E = None
        for ln in open(p, encoding="utf-8"):
            if not ln.startswith("#"):
                break
            m = re.match(r"#\s*E_prim_keV\s*=\s*([0-9.]+)", ln)
            if m:
                E = float(m.group(1))
        if E is None:
            raise SystemExit("%s: в шапке нет E_prim_keV" % p)
        out.append((E, p))
    return sorted(out)


def main():
    rows = inputs()
    if not rows:
        raise SystemExit(
            "Не найдены моно-спектры %s*.csv в %s.\n"
            "Прогон: g1s.exe bridge_mono_v4.mac shield" % (PREFIX, BUILD))
    emax = int(max(E for E, _ in rows)) + 400

    print("Мост методик на чистых 4pi-моно точечной 5 см.")
    print("B = площадь фита (размытый спектр) / истинный пик (депозит).")
    print("Зона фита +-%.1f ПШПВ; подложка степени %d — конвенция прибора.\n"
          % (FIT_ZONE_FWHM, BG_DEGREE))
    print("%10s %12s %12s %8s %8s %10s %7s %8s"
          % ("E, кэВ", "истин.пик", "фит", "B", "+-",
             "ПШПВ,кэВ", "ш/ш0", "ступ."))
    table = []
    blurred = {}
    for E0, path in rows:
        ec = load(path)
        npk = true_peak(ec, E0)
        xs, ys = blur(ec, emax)
        blurred[E0] = (ec, xs, ys, npk)
        area, dA, wrel, hstep = fit_area(xs, ys, E0)
        B = area / npk
        # Числитель и знаменатель считаны с ОДНИХ событий и коррелированы,
        # поэтому сумма в квадратурах — ВЕРХНЯЯ оценка σ(B), не точная.
        dB = B * math.hypot(dA / area, 1.0 / math.sqrt(max(npk, 1.0)))
        print("%10.3f %12.0f %12.0f %8.4f %8.4f %10.1f %7.3f %8.4f"
              % (E0, npk, area, B, dB, dp.fwhm_measured(E0), wrel, hstep))
        table.append((E0, npk, area, B, dB, wrel, hstep))

    band = [(E, B) for E, _, _, B, _, _, _ in table if 583.0 <= E <= 900.0]
    edge = [(E, B) for E, _, _, B, _, _, _ in table if E > 900.0]
    soft = [(E, B) for E, _, _, B, _, _, _ in table if E < 583.0]
    print("\nФОРМА B(E).")
    slope = None
    if len(band) >= 2:
        Es = np.array([E for E, _ in band])
        Bs = np.array([B for _, B in band])
        slope, ic = np.polyfit(np.log10(Es), Bs, 1)
        print("  плато 583…861 кэВ (%d линии): B = %.4f…%.4f; наклон"
              " %+.4f на декаду" % (len(band), Bs[0], Bs[-1], slope))
    # Разрыв проверяется ПРОДОЛЖЕНИЕМ ПЛАТО, а не общей прямой по всем точкам:
    # общая прямая по пяти точкам целиком определяется верхней и выдаёт
    # «плавный рост» там, где на самом деле полка и обрыв.
    for E, B in edge:
        if slope is not None:
            pred = ic + slope * math.log10(E)
            print("  %.1f кэВ: B = %.4f; продолжение плато дало бы %.4f —"
                  " отклонение %+.1f п.п." % (E, B, pred, 100 * (B - pred)))
    if soft:
        print("  мягкие линии (%s) — СПРАВОЧНО: окно полки E-30…E-10 там"
              % "; ".join("%.1f" % E for E, _ in soft))
        print("  ближе центра одной ПШПВ; конвенция к ним неприменима"
              " по арифметике (задача 129).")

    # --- Чувствительность к выборам моста (внутренний аудит 31.07, P0-1).
    # Вывод обязан переживать зону фита, степень подложки и закон ширины;
    # проверяется на решающей линии и контрольной середине.
    print("\nЧУВСТВИТЕЛЬНОСТЬ B К ВЫБОРАМ МОСТА (решающая и контрольная):")
    print("%10s %6s %5s %12s %8s %8s"
          % ("E, кэВ", "зона", "полин", "закон ПШПВ", "B", "ш/ш0"))
    sens_rows = []
    variants = [(2.5, 2, "fwhm_measured", dp.fwhm_measured),
                (2.5, 1, "fwhm_measured", dp.fwhm_measured),
                (3.2, 2, "fwhm_measured", dp.fwhm_measured),
                (3.2, 1, "fwhm_measured", dp.fwhm_measured),
                (2.5, 2, "fwhm_sqrt662", dp.fwhm)]
    for E0 in (661.657, 2614.511):
        got = [E for E in blurred if abs(E - E0) <= 0.05]
        if not got:
            continue
        ec, xs, ys, npk = blurred[got[0]]
        for zone, deg, law_name, law in variants:
            # Смена закона ширины меняет и размытие, не только зону фита.
            if law is not dp.fwhm_measured:
                xs2, ys2 = blur(ec, emax, law)
            else:
                xs2, ys2 = xs, ys
            a, dA, w, _h = fit_area(xs2, ys2, got[0], zone, deg, law)
            Bv = a / npk
            print("%10.3f %6.1f %5d %12s %8.4f %8.3f"
                  % (got[0], zone, deg, law_name, Bv, w))
            sens_rows.append(("%.3f" % got[0], "%.1f" % zone, "%d" % deg,
                              law_name, "%.4f" % Bv, "%.3f" % w))

    print("\nОСТАТОК ПОСЛЕ ПРИВЕДЕНИЯ К ОДНОЙ КОНВЕНЦИИ:"
          " (1 + превышение)·B − 1.")
    # Отношение штатная/наша — точное, из bridge_attested_ratio.csv (полином
    # зон .efa, задача 151 P1-4); погрешность d_ratio_rel — своя на линию
    # (квадратура СКО фита обеих зон), не зашитая константа.
    print("%10s %12s %8s %18s" % ("E, кэВ", "превышение", "B", "остаток"))
    res = []
    matched = set()
    for E, _, _, B, dB, _, _ in table:
        rv = attested_ratio(E)
        if rv is None:
            continue
        r, d_ratio_pct = rv
        d_ratio_rel = d_ratio_pct / 100.0
        matched.add(round(E, 1))
        rest = B / r - 1.0
        drest = (B / r) * math.hypot(dB / B, d_ratio_rel)
        res.append((E, rest, drest))
        print("%10.1f %11.1f %% %8.4f %11.1f +- %.1f %%"
              % (E, 100 * (1 / r - 1), B, 100 * rest, 100 * drest))
    missed = sorted(round(E, 1) for E in ATTESTED_RATIO
                    if round(E, 1) not in matched)
    if missed:
        # Отказ, а не предупреждение: без 2614;5 таблица остатка выглядит
        # плоской просто потому, что решающей точки в ней нет. Ровно так и
        # вышло при первом прогоне — ключ 2614;511 против 2614;510 в шапке.
        raise SystemExit(
            "bridge_mono: для линий %s есть измеренное отношение; но нет"
            " моно-прогона.\nВывод о форме остатка без них недействителен —"
            " досчитайте прогон или уберите строку из ATTESTED_RATIO."
            % "; ".join("%.1f" % E for E in missed))
    if len(res) >= 2:
        lo = min(x for _, x, _ in res)
        hi = max(x for _, x, _ in res)
        # Плоскость проверяется против погрешностей, а не зашитым порогом:
        # хи-квадрат согласия остатков с константой (взвешенное среднее).
        w = [1.0 / d ** 2 for _, _, d in res]
        mean = sum(x * ww for (_, x, _), ww in zip(res, w)) / sum(w)
        chi2 = sum(ww * (x - mean) ** 2
                   for (_, x, _), ww in zip(res, w)) / (len(res) - 1)
        print("  разброс остатка: %.1f…%.1f %% (размах %.1f п.п.);"
              " среднее %.1f %%; хи2/ню согласия с константой %.2f"
              % (100 * lo, 100 * hi, 100 * (hi - lo), 100 * mean, chi2))
        if chi2 < 2.0:
            print("  ВЫВОД: остаток СОВМЕСТИМ С ПЛОСКИМ в пределах"
                  " погрешностей. Это не доказательство плоскости:")
            print("  погрешность точки ~1 п.п. от округлённых активностей"
                  " не даёт различить размах меньше ~2 п.п.")
        else:
            print("  ВЫВОД: остаток НЕ плоский — конвенция снимает не весь"
                  " разрыв; остаётся энергетически зависимый вклад.")
    print("  ОГОВОРКА: у СпектраЛайн подложка — ступенька плюс полином"
          " по РЕАЛЬНОМУ спектру с континуумом соседних линий;")
    print("  здесь моно-спектр — континуум только собственный. Перенос"
          " B на сеанс Th-228 остаётся приближением.")

    # Таблица остатка — ФАЙЛОМ, а не только в print. До 01.08.2026 она жила
    # исключительно в консольном выводе, и документы (report.md §5.3, статья)
    # цитировали её оттуда — то есть из отчёта, а не из производящей таблицы.
    # Тот же класс дефекта, что задача 148 у compare_cups.py и задача 109 у
    # perline-таблицы; найден внешним аудитом.
    if res:
        w = [1.0 / d ** 2 for _, _, d in res]
        mean = sum(x * ww for (_, x, _), ww in zip(res, w)) / sum(w)
        chi2 = (sum(ww * (x - mean) ** 2 for (_, x, _), ww in zip(res, w))
                / (len(res) - 1)) if len(res) > 1 else float("nan")
        table_by_E = {E: (n, a, B, dB, wr, h)
                      for E, n, a, B, dB, wr, h in table}
        out_rows = []
        for E, rest, drest in res:
            rv = attested_ratio(E)
            near_E = min(table_by_E, key=lambda k: abs(k - E))
            Bv, dBv = table_by_E[near_E][2], table_by_E[near_E][3]
            out_rows.append(
                ("%.3f" % E, "%+.2f" % (100 * (1 / rv[0] - 1)),
                 "%.4f" % Bv, "%.4f" % dBv,
                 "%+.2f" % (100 * rest), "%.2f" % (100 * drest),
                 "1" if E >= 583.0 else "0"))
        csvio.write(
            os.path.join(OUT, "bridge_residual.csv"),
            ["E_keV", "excess_model_pct", "bridge", "d_bridge",
             "residual_after_bridge_pct", "d_residual_pct", "decisive"],
            out_rows,
            comments=[
                "Остаток после приведения расчёта и аттестации к ОДНОЙ"
                " конвенции съёма площади: (1 + превышение)*B - 1.",
                "excess_model_pct — превышение расчёта над аттестацией"
                " (bridge_attested_ratio.csv); bridge — поправка конвенции"
                " (эта же программа; столбец bridge в bridge_mono.csv).",
                "d_residual_pct — квадратура относительной погрешности"
                " моста и СТРОГОЙ погрешности отношения кривых"
                " (bridge_attested_ratio.csv; с учётом того; что она РАЗНАЯ"
                " на разных линиях).",
                "Согласие остатка с константой: хи2/ню = %.2f при среднем"
                " %.2f %%. Больше ~2 означает; что конвенция снимает не весь"
                " разрыв и остаётся энергетически зависимый вклад."
                % (chi2, 100 * mean),
                "Заведено 01.08.2026: до этого таблица существовала только в"
                " консольном выводе; и документы цитировали её из отчёта; а"
                " не из производящего файла.",
            ],
            stamp=stamp.lines(
                "detectors/Gamma-1S/analysis/bridge_mono.py",
                dict(OBS, quantity="остаток расхождения расчёт/аттестация"
                                   " после приведения к одной конвенции съёма"
                                   " площади; и хи2/ню его согласия с"
                                   " константой"),
                inputs=[p for _, p in rows],
                geometry_dir=str(paths.geometry("Gamma-1S")),
                names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO)))
        print("\nтаблица остатка: %s"
              % os.path.join(OUT, "bridge_residual.csv"))

    csvio.write(
        os.path.join(OUT, "bridge_mono.csv"),
        ["E_keV", "true_peak", "fit_area", "bridge", "d_bridge",
         "fit_width_rel", "step_height", "decisive"],
        [("%.3f" % E, "%.0f" % n, "%.0f" % a, "%.4f" % B, "%.4f" % dB,
          "%.3f" % w, "%.4f" % h, "1" if E >= 583.0 else "0")
         for E, n, a, B, dB, w, h in table],
        comments=[
            "Мост методик: во сколько раз конвенция аттестации (фит по"
            " размытому спектру) занижает площадь против истинного пика.",
            "Зона фита +-%.1f ПШПВ и полином степени %d — конвенция"
            " ПРИБОРА (lsrm.cnf); чувствительность к этим выборам — в"
            " bridge_sensitivity.csv." % (FIT_ZONE_FWHM, BG_DEGREE),
            "decisive=1 — линия внутри области; где конвенция применима;"
            " decisive=0 — мягкий край; читать справочно (задача 129).",
            "fit_width_rel — отношение ширины фита к приборной ПШПВ."
            " Размытый пик гауссов ПО ПОСТРОЕНИЮ blur(); поэтому отход"
            " от 1;0 — отказ модели ПОДЛОЖКИ; не свойство пика.",
            "d_bridge — верхняя оценка: сигма площади из ковариации фита"
            " плюс пуассон знаменателя в квадратурах; числитель и"
            " знаменатель коррелированы (одни события).",
            "step_height — высота ступеньки-из-образа в долях амплитуды;"
            " подгоняется; не задаётся.",
        ],
        stamp=stamp.lines(
            "detectors/Gamma-1S/analysis/bridge_mono.py", OBS,
            inputs=[p for _, p in rows],
            geometry_dir=str(paths.geometry("Gamma-1S")),
            names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO)))
    csvio.write(
        os.path.join(OUT, "bridge_sensitivity.csv"),
        ["E_keV", "fit_zone_fwhm", "bg_degree", "fwhm_law", "bridge",
         "fit_width_rel"],
        sens_rows,
        comments=[
            "Чувствительность B к выборам моста: зона фита (2;5 — прибор;"
            " 3;2 — прежняя реализация); степень подложки; закон ширины.",
            "Введено по внутреннему аудиту 31.07.2026 (P0-1): при зоне 3;2"
            " комптоновский край 2614;5 (2381;8 кэВ) лежит внутри зоны.",
        ],
        stamp=stamp.lines(
            "detectors/Gamma-1S/analysis/bridge_mono.py", OBS,
            inputs=[p for _, p in rows],
            geometry_dir=str(paths.geometry("Gamma-1S")),
            names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO)))
    print("\nтаблицы: %s; %s"
          % (os.path.join(OUT, "bridge_mono.csv"),
             os.path.join(OUT, "bridge_sensitivity.csv")))


if __name__ == "__main__":
    main()
