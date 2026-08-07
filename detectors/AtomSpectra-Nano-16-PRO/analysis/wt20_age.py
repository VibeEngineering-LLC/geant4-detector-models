# -*- coding: utf-8 -*-
"""Возраст ряда тория в пачке WT-20 по отношению площадей опорных пиков.

Определяемая величина — время, прошедшее с ХИМИЧЕСКОЙ ОЧИСТКИ тория (не с
даты изготовления электродов: очистка предшествует изготовлению, и промежуток
между ними этой задачей не определяется).

Наблюдаемая величина
--------------------
Отношение площадей двух пиков, принадлежащих РАЗНЫМ ветвям ряда:

    911,20 кэВ (Ac-228, верхняя ветвь, следует за Ra-228, T½ = 5,75 года)
    238,63 кэВ (Pb-212, нижняя ветвь, следует за Th-228, T½ = 1,91 года)

После очистки Ra-228 нарастает из нуля, а Th-228 остаётся от очистки целиком
(химия отделяет элемент, а не изотоп) и первые годы убывает. Поэтому отношение
ветвей — прямая функция возраста, а не свойство образца.

Почему отношение, а не активности
---------------------------------
Отношение площадей не зависит ни от массы пачки, ни от содержания ThO2, ни от
абсолютной эффективности прибора: все эти множители одинаковы у обоих пиков и
сокращаются. Остаётся зависимость от ОТНОСИТЕЛЬНОГО отклика на 911 и 238 кэВ —
она и составляет основную систематическую погрешность результата.

Как считается
-------------
Возраст ищется не разложением спектра, а прямым сравнением площадей:

    Q(t) = S_мод(окно 911, t) / S_мод(окно 238, t)      — модель при возрасте t
    Q_изм = S_изм(окно 911) / S_изм(окно 238)           — измерение

и решается Q(t) = Q_изм. Модель при этом — полная сумма всех звеньев при
активностях возраста t, поэтому примеси чужих линий в окне (например 240,99 кэВ
Ra-224 внутри окна 238,63) учитываются сами, без отдельных поправок.

Площадь у модели и у измеренного снимается ОДИНАКОВО: сумма в окне ±1 ПШПВ за
вычетом подложки, оценённой по двум внешним подокнам той же ширины. ПШПВ берётся
из кривой прибора, записанной в файле замера, и под спектр не настраивается.

Правки по внешнему адверсариальному аудиту (07.08.2026, вечер; вердикт
«ПРАВКИ ОБЯЗАТЕЛЬНЫ», см. `docs/wt20-methods-appendix.md` §12.6):

  Б1 — подокна подложки исключают бины, лежащие в сигнальном окне ЛЮБОГО
       другого опорного пика (`exclude` в peak_net). Без этого подокно
       911,20 справа на 54 % состояло из дублета Ac-228 965/969 кэВ — сосед
       того же нуклида, а не континуум.
  Б3 — дисперсия суммы по окну считается ПРАВИЛЬНО: не суммой дисперсий по
       каналу (занижает, каналы коррелированы одним МК-розыгрышем, размытым
       ядром на несколько каналов сразу), а по исходным бинам шаблона —
       `window_net_var`.
  Б2 — чувствительность к ширине окна снята прогоном (не предположением) и
       напечатана отдельной диагностикой.
  Б4/Б5 — перенос систематики отклика больше не выдаётся за измеренную
       границу: печатаются ОБЕ оценки (по соседней паре и по паре, кроющей
       весь рабочий интервал), обе явно помечены как экстраполяция.
  Б6 — путь к файлу калибровки печатается, отсутствие — предупреждение, а не
       молчание; применённые коэффициенты идут в выходной CSV.

    python analysis/wt20_age.py <спектр.xml> <шаблоны> [-o <каталог>]
"""
import argparse
import csv
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.environ.get("SPECTRAVIBE_ROOT")
if not _ROOT:
    raise SystemExit("не задана переменная окружения SPECTRAVIBE_ROOT")
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, _HERE)

from wt20_unfold import (E_MAX, E_STEP, TAIL_L, TAIL_R, fwhm_from_file,  # noqa: E402
                         line_shape, read_atomspectra_xml, read_correction,
                         rebin_to_grid)
from wt20_forward import (CHAIN, YEAR_SEC, activities_after_purification,  # noqa: E402
                          branch_tl208, half_lives, load_raw_templates,
                          load_templates, th232_activity, THO2_FRAC_MID)

# Опорные пики. Ветвь указана явно: возраст определяется ТОЛЬКО парами из
# разных ветвей, пары внутри одной ветви служат проверкой отклика.
#
# 238,63 кэВ — Pb-212, выход 43,6 %, самая сильная линия нижней ветви;
# 911,20 и 968,97 кэВ — Ac-228, выходы 25,8 и 15,8 %, верхняя ветвь;
# 583,19 и 2614,51 кэВ — Tl-208, нижняя ветвь, включены как проверка отклика:
#   Tl-208 связан с Pb-212 фиксированным ветвлением Bi-212, поэтому отношение
#   583/238 от возраста почти не зависит и любое его отклонение — свойство
#   отклика, а не ряда.
PEAKS = [
    (238.63,  "Pb-212", "нижняя"),
    (583.19,  "Tl-208", "нижняя"),
    (911.20,  "Ac-228", "верхняя"),
    (968.97,  "Ac-228", "верхняя"),
    (2614.51, "Tl-208", "нижняя"),
]

# Пара, дающая ТОЧЕЧНУЮ оценку возраста: (числитель, знаменатель).
#
# Взята одна — 911,20 / 238,63. Пара 968,97 / 238,63 держится ОТДЕЛЬНО как
# перекрёстная проверка, не как точечная оценка — по причине, изменившейся
# по ходу разработки (записано, чтобы не повторить ошибку):
#
#   ДО фикса Б1 (взаимное исключение подложки): окно 968,97 показывало
#   вклад Tl-208 около минус 85 % площади — отношение мерило не ветви ряда,
#   а точность вычитания подложки под комптоновским краем Tl-208.
#
#   ПОСЛЕ фикса Б1: композиция обоих окон стала чистой (Ac-228 доминирует
#   в обоих), но по КОНСТРУКЦИИ оба пика используют ОДНОСТОРОННЮЮ подложку
#   НА ПРОТИВОПОЛОЖНЫХ сторонах — 911,20 берёт только левое подокно (правое
#   съедено дублетом со стороны 969), 968,97 только правое (левое съедено со
#   стороны 911). Если континуум в этой области имеет наклон, каждая оценка
#   получает свой знак смещения от него, и это правдоподобно объясняет
#   расхождение результатов пары на ~0,7 года при обеих чистых композициях.
#
# Раз причина отбраковки исчезла, а расхождение осталось — пара 968,97/238,63
# участвует в общем диапазоне неопределённости НАРАВНЕ с шириной окна и
# систематикой отклика (см. envelope_lo/envelope_hi), а не отбрасывается
# молча. Первичной остаётся 911,20/238,63 — только по порядку исторического
# выбора, не потому что 968,97/238,63 хуже обоснована.
AGE_PAIRS = [(911.20, 238.63)]
CROSSCHECK_PAIRS = [(968.97, 238.63)]

# Пары-проверки ОТКЛИКА. Общее свойство всех четырёх: обе линии пары связаны
# фиксированным ветвлением/принадлежат одному нуклиду, поэтому отношение от
# возраста практически не зависит, и расхождение с моделью — погрешность
# ОТНОСИТЕЛЬНОГО отклика между двумя энергиями, и ничто иное.
#
#   583,19/238,63  — Tl-208/Pb-212, соседняя пара, ближайшая к рабочему
#                    интервалу 238->911 (доля покрытия по ln E — 67 %);
#   2614,51/583,19 — обе линии Tl-208, самый чистый отклик без завязки на
#                    ветвление;
#   968,97/911,20  — обе линии Ac-228, узкий интервал;
#   2614,51/238,63 — КРОЕТ ВЕСЬ рабочий интервал 238->911 и заходит за него;
#                    добавлена по находке Б5 внешнего аудита — прежняя
#                    редакция сравнивала «устойчивость» с парой на СМЕЖНОМ,
#                    а не на сопоставимом интервале.
CHECK_PAIRS = [(583.19, 238.63), (2614.51, 583.19), (968.97, 911.20),
               (2614.51, 238.63)]
SYS_PAIR = (583.19, 238.63)
CUM_PAIR = (2614.51, 238.63)

# Ширины окна для проверки чувствительности (Б2). 1.0 — принятая; остальные —
# диагностика систематики выбора окна, не альтернативные результаты.
WIDTH_SCALES = [0.70, 0.85, 1.00, 1.15, 1.30]


def peak_windows(fwhm_fn):
    """Сигнальные окна всех опорных пиков — общий список для взаимного
    исключения из подложки (см. peak_net, exclude)."""
    return [(e0, float(fwhm_fn(e0))) for e0, _nuc, _br in PEAKS]


# Минимум бинов на стороне подложки, чтобы её учитывать отдельно; меньше —
# сторона игнорируется (используется только другая, если она надёжна), либо
# окно негодно целиком (обе ниже порога). Поднят с 2 до 5 по находке
# внешнего аудита (07.08.2026, второй проход, Б2): при ЧАСТИЧНОМ взаимном
# исключении (сканы ширины окна между 1,00 и 1,15 ПШПВ) на границе
# исключённой зоны выживает огрызок в 2-3 бина, который при пороге 2
# получал бы тот же вес (0,5), что и надёжная сторона из 25-30 бинов —
# порог 5 не даёт такому огрызку весить наравне с полноценной стороной.
MIN_SIDE_BINS = 5


def _baseline_mode(n_bl, n_br, min_bins=MIN_SIDE_BINS):
    """Как усреднять подложку: 'both' / 'left' / 'right' / None (обе стороны
    ниже порога — окно негодно).

    ОДНА функция на три места (peak_net, net_with_sigma, window_net_var):
    порознь эта логика уже расходилась — МК-дисперсия тайно обнулялась там,
    где взаимное исключение (Б1) оставляло только одну сторону, потому что
    ветвление было продублировано и не перенесено во все три места synchronно
    (находка внешнего аудита, 07.08.2026, второй проход, Б-находка 1).
    """
    ok_l, ok_r = n_bl >= min_bins, n_br >= min_bins
    if ok_l and ok_r:
        return "both"
    if ok_l:
        return "left"
    if ok_r:
        return "right"
    return None


def peak_net(v, centres, e0, fw, exclude=None):
    """Площадь пика: сумма в ±1 ПШПВ минус подложка по внешним подокнам.

    `exclude` — список (E, ПШПВ) чужих пиков: бины подокон подложки,
    попавшие в сигнальное окно (±1 ПШПВ) любого из них, из среднего подложки
    ИСКЛЮЧАЮТСЯ. Без этого подложка близко расположенного сильного пика сама
    оказывается пиком: подокно 911,20 кэВ справа на 54 % состояло из дублета
    Ac-228 965/969 кэВ (разнесены на 58 кэВ при ПШПВ ~50) — вычитание убирало
    70 % нетто самого пика 911,20 (найдено внешним аудитом, 07.08.2026, Б1).

    Подложка усредняется по ДВУМ сторонам, когда обе доступны; если одна
    сторона выбита исключением целиком (911,20 и 968,97 — соседи через
    58 кэВ при ПШПВ ~50: каждый вычищает ВСЮ противоположную сторону
    другого), используется ОДНА уцелевшая сторона, а не отказ — это то же
    решение, которое проверил внешний аудит («только левое подокно»,
    возраст 2,50 против 1,97 при симметричном подсчёте до правки Б1).
    Откaз (None) — только если ОБЕ стороны исключены целиком.

    Возвращает (нетто, брутто, число каналов окна, подложка на канал, маски
    sel/bl/br — bl/br уже с учётом исключения, если оно задано).
    """
    sel = (centres >= e0 - fw) & (centres <= e0 + fw)
    bl = (centres >= e0 - 2.0 * fw) & (centres < e0 - fw)
    br = (centres > e0 + fw) & (centres <= e0 + 2.0 * fw)
    if exclude:
        clean = np.ones(len(centres), dtype=bool)
        for xe, xfw in exclude:
            clean &= ~((centres >= xe - xfw) & (centres <= xe + xfw))
        bl = bl & clean
        br = br & clean
    n_bl, n_br = int(bl.sum()), int(br.sum())
    mode = _baseline_mode(n_bl, n_br)
    if sel.sum() < 3 or mode is None:
        return None
    if mode == "both":
        b = 0.5 * (v[bl].mean() + v[br].mean())
    elif mode == "left":
        b = v[bl].mean()
    else:
        b = v[br].mean()
    net = float(v[sel].sum() - b * sel.sum())
    return net, float(v[sel].sum()), int(sel.sum()), b, sel, bl, br


def net_with_sigma(y_smp, y_bg, centres, e0, fw, exclude=None):
    """Нетто пика по разности «образец − приведённый фон» и его σ.

    Дисперсия считается по исходным отсчётам обоих спектров: у образца это
    сам отсчёт (пуассон), у приведённого фона — отсчёт, умноженный на квадрат
    коэффициента приведения, который уже сидит в переданном массиве. Точное
    значение коэффициента здесь не нужно: приведённый фон много меньше
    образца в окнах опорных пиков, и его вклад в дисперсию оценивается сверху
    как сам приведённый отсчёт.

    Формула ветвится ТЕМ ЖЕ способом, что подложка в peak_net (обе стороны /
    одна уцелевшая) — если бы формулы разошлись, сигма считала бы дисперсию
    не той величины, что реально вычтена из нетто.
    """
    r = peak_net(y_smp - y_bg, centres, e0, fw, exclude=exclude)
    if r is None:
        return None
    net, _gross, n, _b, sel, bl, br = r
    n_bl, n_br = int(bl.sum()), int(br.sum())
    var_in = float((y_smp[sel] + y_bg[sel]).sum())

    def side_sum_var(mask):
        return float(y_smp[mask].sum() + y_bg[mask].sum())

    mode = _baseline_mode(n_bl, n_br)
    if mode == "both":
        # b = 0.5·(ΣL/nL + ΣR/nR); var(n·b) = (n/2)²·(varΣL/nL² + varΣR/nR²)
        var_bg = (0.5 * n / n_bl) ** 2 * side_sum_var(bl) \
            + (0.5 * n / n_br) ** 2 * side_sum_var(br)
    elif mode == "left":
        var_bg = (float(n) / n_bl) ** 2 * side_sum_var(bl)
    else:
        var_bg = (float(n) / n_br) ** 2 * side_sum_var(br)
    var = var_in + var_bg
    return net, math.sqrt(max(var, 0.0))


def window_net_var(e_raw, c_raw, n_prim, centres, sel, bl, br, fwhm_fn,
                    tail_l=TAIL_L, tail_r=TAIL_R):
    """Дисперсия ВЗВЕШЕННОЙ СУММЫ по окну (сигнал минус подложка) от
    статистики Монте-Карло шаблона, а не сумма дисперсий по каналу.

    Один и тот же МК-розыгрыш размывается по многим выходным каналам ядром
    формы линии; суммировать дисперсии независимо по каналам — занижать
    результат (соседние каналы положительно коррелированы одним и тем же
    исходным событием). Правильно: выразить сумму по окну как линейную
    комбинацию исходных (нерасплывшихся) счётчиков c_j и взять дисперсию
    этой комбинации целиком.

    Сумма по окну W (сигнал минус подложка, как в peak_net), ОБЕ стороны:

        S_W = Σ_j (c_j/N) · H_j,   H_j = G_sel,j − (n_sel/2)·(G_bl,j/n_bl + G_br,j/n_br)

    одна уцелевшая сторона (см. `_baseline_mode`, тот же режим, что и в
    peak_net/net_with_sigma — три места ветвили эту логику порознь, и МК-
    дисперсия тайно обнулялась ровно там, где взаимное исключение Б1
    оставляло только одну сторону: находка внешнего аудита, 07.08.2026,
    второй проход):

        H_j = G_sel,j − n_sel·(G_bl,j/n_bl)         [режим 'left']
        H_j = G_sel,j − n_sel·(G_br,j/n_br)         [режим 'right']

    где G_X,j = Σ_{i∈X} K_ij — сумма ядра формы линии от бина j по всем
    каналам множества X (те же ядра, что в `broaden`). Дисперсия по Пуассону
    (Var[c_j] ≈ c_j):

        Var[S_W] = Σ_j c_j · H_j² / N².

    Возвращает Var[S_W] (не сигму), в единицах (отсч./с)² при том же
    масштабе входа, что у `broaden` — умножается на (A·t)² вызывающей
    стороной для получения дисперсии в отсчётах.
    """
    n_sel, n_bl, n_br = int(sel.sum()), int(bl.sum()), int(br.sum())
    mode = _baseline_mode(n_bl, n_br)
    if n_sel == 0 or mode is None:
        return 0.0
    step = centres[1] - centres[0]
    lo0 = centres[0] - 0.5 * step
    reach = 8.0 if (tail_l or tail_r) else 4.0
    var_sum = 0.0
    for e, c in zip(e_raw, c_raw):
        if c <= 0:
            continue
        fw_e = max(float(fwhm_fn(e)), 1.0)
        s = fw_e / 2.3548
        k0 = max(int((e - reach * s - lo0) / step), 0)
        k1 = min(int((e + reach * s - lo0) / step) + 1, len(centres))
        if k1 <= k0:
            continue
        g = line_shape(centres[k0:k1], e, s, tail_l, tail_r)
        ssum = g.sum()
        if ssum <= 0:
            continue
        g = g / ssum
        g_sel = float(g[sel[k0:k1]].sum())
        if mode == "both":
            g_bl = float(g[bl[k0:k1]].sum())
            g_br = float(g[br[k0:k1]].sum())
            h = g_sel - 0.5 * n_sel * (g_bl / n_bl + g_br / n_br)
        elif mode == "left":
            g_bl = float(g[bl[k0:k1]].sum())
            h = g_sel - n_sel * (g_bl / n_bl)
        else:
            g_br = float(g[br[k0:k1]].sum())
            h = g_sel - n_sel * (g_br / n_br)
        var_sum += c * h * h
    return var_sum / (n_prim * n_prim)


def model_at(tmap, acts, centres):
    """Суммарный спектр модели при заданных активностях звеньев, отсч./с."""
    s = np.zeros(len(centres))
    for key, spec_k in tmap.items():
        s += acts[key] * spec_k
    return s


def solve_age(q_target, q_fn, t_lo=0.05, t_hi=8.0, n=4000):
    """Младший возраст на ВОСХОДЯЩЕЙ ветви кривой, при котором Q(t) = q_target.

    Кривая отношения ветвей немонотонна: растёт от нуля, проходит максимум
    около 7,9 года и затем медленно спадает к единице (вековое равновесие).
    Функция ищет корень ТОЛЬКО на восходящем участке [t_lo, t_max]; если
    q_target лежит в диапазоне, достижимом и на падающей ветви (это возможно
    при q_target между значением в максимуме и асимптотикой на больших t),
    второе решение здесь не ищется и не возвращается — вызывающая сторона не
    должна интерпретировать `ok=False` как «возраст равен t_hi»: это означает
    «на восходящей ветви решения нет», а не «решение около 8 лет». Все
    вызовы этой функции в модуле работают со значениями Q, для которых
    восходящая ветвь — единственная физически осмысленная (партия моложе
    десятка лет, см. §12.5 приложения); для проверки больших гипотетических
    возрастов (устойчивость к 20-летнему допущению NUREG/CR-1039) кривая
    Q(t) читается НАПРЯМУЮ (q_fn(20.0)), а не через этот решатель.
    """
    ts = np.linspace(t_lo, t_hi, n)
    qs = np.array([q_fn(t) for t in ts])
    i_max = int(np.argmax(qs))
    ts, qs = ts[:i_max + 1], qs[:i_max + 1]
    if q_target <= qs[0]:
        return float(ts[0]), False
    if q_target >= qs[-1]:
        return float(ts[-1]), False
    j = int(np.searchsorted(qs, q_target))
    t0, t1, q0, q1 = ts[j - 1], ts[j], qs[j - 1], qs[j]
    return float(t0 + (q_target - q0) * (t1 - t0) / (q1 - q0)), True


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("xml", help="файл замера AtomSpectra")
    ap.add_argument("tdir", help="каталог понуклидных МК-шаблонов")
    ap.add_argument("-o", "--out", help="каталог вывода")
    args = ap.parse_args()
    outdir = args.out or os.path.dirname(args.xml)
    os.makedirs(outdir, exist_ok=True)

    spec = read_atomspectra_xml(args.xml)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    t_smp = float(spec.real_time)

    # Калибровка — путём и явным предупреждением, а не молча (Б6 внешнего
    # аудита): прогон с несуществующим каталогом вывода незаметно откатывался
    # на заводскую шкалу и менял Q на 16 %, без единой строки об этом.
    calib_path = os.path.join(outdir, "calibration_fitted.csv")
    corr = read_correction(calib_path)
    if not corr:
        print("ВНИМАНИЕ: поправок калибровки нет (%s не найден) — работаем "
              "по заводской шкале" % calib_path)
    else:
        print("калибровка: %s (образец %s, фон %s)"
              % (calib_path, corr.get("sample"), corr.get("background")))

    edges = np.arange(0.0, E_MAX + E_STEP, E_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])
    y = rebin_to_grid(np.asarray(spec.counts, float), list(spec.energy_cal),
                      corr.get("sample"), edges)
    if bg is not None:
        k_bg = t_smp / float(bg.real_time)
        ybg = rebin_to_grid(np.asarray(bg.counts, float), list(bg.energy_cal),
                            corr.get("background"), edges) * k_bg
    else:
        ybg = np.zeros_like(y)

    fwhm_fn = fwhm_from_file(args.xml)
    if fwhm_fn is None:
        raise SystemExit("в файле замера нет кривой ПШПВ прибора")

    hl = half_lives()
    br_tl = branch_tl208()
    tmap, head0 = load_templates(args.tdir, centres, fwhm_fn)
    raw_tmap = load_raw_templates(args.tdir)   # для дисперсии МК (Б3)
    mass_g = float(head0.get("wt20_mass_g", "0").split()[0])
    if mass_g <= 0:
        raise SystemExit("в шапке шаблона нет массы пачки (wt20_mass_g)")
    # a0 нужна только чтобы прогнать модель по возрасту; отношение Q от неё
    # не зависит (§12.1 приложения, проверено внешним аудитом), поэтому
    # содержание ThO2 здесь берётся средним по марке, без вилки.
    a0 = th232_activity(mass_g, THO2_FRAC_MID)

    windows = peak_windows(fwhm_fn)

    print("замер %.0f с, шаблонов %d, A(Th-232) = %.0f Бк (среднее по марке)"
          % (t_smp, len(tmap), a0))
    print("ПШПВ прибора: 238,63 кэВ -> %.1f кэВ, 911,20 кэВ -> %.1f кэВ"
          % (fwhm_fn(238.63), fwhm_fn(911.20)))

    # Диагностика Б1: сколько бинов подложки реально вырезано взаимным
    # исключением, по сравнению с наивным подсчётом без него.
    print("\n--- взаимное исключение подложки (Б1) ---")
    for e0, nuc, _br in PEAKS:
        fw = float(fwhm_fn(e0))
        r0 = peak_net(y - ybg, centres, e0, fw, exclude=None)
        r1 = peak_net(y - ybg, centres, e0, fw, exclude=windows)
        if r0 is None or r1 is None:
            continue
        n_bl0, n_br0 = int(r0[5].sum()), int(r0[6].sum())
        n_bl1, n_br1 = int(r1[5].sum()), int(r1[6].sum())
        cut = (n_bl0 + n_br0) - (n_bl1 + n_br1)
        flag = "  <-- вырезано" if cut else ""
        print("  %9.2f (%s): подокна %d+%d -> %d+%d бин%s"
              % (e0, nuc, n_bl0, n_br0, n_bl1, n_br1, flag))

    # --- площади опорных пиков в измеренном спектре -------------------------
    print("\n--- площади опорных пиков (образец минус приведённый фон) ---")
    print("  %9s %-8s %-8s %12s %10s %7s"
          % ("E, кэВ", "нуклид", "ветвь", "нетто", "сигма", "отн, %"))
    meas = {}
    for e0, nuc, branch in PEAKS:
        fw = float(fwhm_fn(e0))
        r = net_with_sigma(y, ybg, centres, e0, fw, exclude=windows)
        if r is None:
            print("  %9.2f  окно не помещается в сетку" % e0)
            continue
        net, sig = r
        meas[e0] = (net, sig)
        print("  %9.2f %-8s %-8s %12.0f %10.0f %7.2f"
              % (e0, nuc, branch, net, sig, 100.0 * sig / net if net else 0))

    # --- модель: нетто окна (дёшево — для перебора решателя) ----------------
    def model_net(acts, e0, fw):
        m = model_at(tmap, acts, centres) * t_smp
        r = peak_net(m, centres, e0, fw, exclude=windows)
        return r[0] if r is not None else None

    # --- модель: нетто И его МК-дисперсия (дорого — раз в точке ответа) -----
    # window_net_var перебирает каждую сырую линию каждого нуклида (у Ac-228
    # их 272); solve_age прогоняет q_fn до 4000 раз на вызов — считать
    # дисперсию ВНУТРИ q_fn означало десятки миллионов вычислений ядра формы
    # линии на один прогон программы. Дисперсия нужна только В ТОЧКЕ ОТВЕТА
    # (после того, как t_hat уже найден дешёвым q_fn), поэтому вынесена в
    # отдельную функцию, вызываемую по требованию, не из решателя.
    def model_net_and_var(acts, e0, fw):
        m = model_at(tmap, acts, centres) * t_smp
        r = peak_net(m, centres, e0, fw, exclude=windows)
        if r is None:
            return None
        net, _gross, _n, _b, sel, bl, br = r
        var = 0.0
        for key, (e_raw, c_raw, n_prim) in raw_tmap.items():
            a_k = acts.get(key, 0.0)
            if a_k == 0.0:
                continue
            v = window_net_var(e_raw, c_raw, n_prim, centres, sel, bl, br,
                               fwhm_fn)
            var += (a_k * t_smp) ** 2 * v
        return net, math.sqrt(max(var, 0.0))

    def q_fn(pair):
        num, den = pair

        def q(t_years):
            A = activities_after_purification(t_years * YEAR_SEC, a0, hl, br_tl)
            a = model_net(A, num, float(fwhm_fn(num)))
            b = model_net(A, den, float(fwhm_fn(den)))
            if a is None or b is None or b == 0:
                return float("nan")
            return a / b
        return q

    def rel_mc_at(pair, t_years):
        """Относительная МК-погрешность Q(t) в точке t — статистика шаблонов,
        не счёт замера (Б3)."""
        num, den = pair
        A = activities_after_purification(t_years * YEAR_SEC, a0, hl, br_tl)
        rn = model_net_and_var(A, num, float(fwhm_fn(num)))
        rd = model_net_and_var(A, den, float(fwhm_fn(den)))
        if rn is None or rd is None or rn[0] == 0 or rd[0] == 0:
            return float("nan")
        return math.hypot(rn[1] / rn[0], rd[1] / rd[0])

    # --- из чего сложена площадь окна у модели ------------------------------
    # Отношение площадей меряет ветви ряда лишь настолько, насколько окно
    # принадлежит своему нуклиду. Разбор печатается до всякого вывода о
    # возрасте: он и решает, какая пара годится, а какая нет.
    print("\n--- состав окна у модели (возраст 2 года), доля площади ---")
    A2 = activities_after_purification(2.0 * YEAR_SEC, a0, hl, br_tl)
    parts_by_peak = {}
    for e0, nuc, _branch in PEAKS:
        fw = float(fwhm_fn(e0))
        parts, tot = {}, 0.0
        for key, spec_k in tmap.items():
            r = peak_net(A2[key] * spec_k * t_smp, centres, e0, fw,
                        exclude=windows)
            if r is None:
                continue
            parts[key] = r[0]
            tot += r[0]
        parts_by_peak[e0] = (parts, tot)
        top = sorted(parts.items(), key=lambda p: -abs(p[1]))
        s = ", ".join("%s %+.0f %%" % (k, 100.0 * v / tot)
                      for k, v in top if abs(v) > 0.02 * abs(tot))
        print("  %9.2f (%s): %s" % (e0, nuc, s))

    print("\n--- отношение площадей: измерение и модель по возрастам ---")
    grid_years = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 20.0]
    hdr = "  %-18s %10s" % ("пара", "измерено")
    hdr += "".join("%9s" % ("%g л" % t) for t in grid_years)
    print(hdr)
    rows_out, ages = [], []
    for pair in AGE_PAIRS + CROSSCHECK_PAIRS + CHECK_PAIRS:
        num, den = pair
        if num not in meas or den not in meas:
            continue
        q_meas = meas[num][0] / meas[den][0]
        rel_meas = math.hypot(meas[num][1] / meas[num][0],
                              meas[den][1] / meas[den][0])
        qf = q_fn(pair)
        vals = [qf(t) for t in grid_years]
        print("  %-18s %10.4f%s"
              % ("%.0f/%.0f" % (num, den), q_meas,
                 "".join("%9.4f" % v for v in vals)))
        rows_out.append((num, den, q_meas, rel_meas, vals))
        if pair in AGE_PAIRS or pair in CROSSCHECK_PAIRS:
            t_hat, ok = solve_age(q_meas, qf)
            if not ok:
                print("  [%.0f/%.0f: возраст вне восходящей ветви, "
                      "исключено из сводки]" % (num, den))
                continue
            rel_mc = rel_mc_at(pair, t_hat)
            rel_total = math.hypot(rel_meas, rel_mc)
            # Q(t) растёт с возрастом на восходящей ветви, поэтому БОЛЬШЕЕ Q
            # даёт БОЛЬШИЙ возраст: q·(1+σ) -> t_hi, q·(1−σ) -> t_lo. Прежняя
            # редакция называла их наоборот (t_lo от q·(1+σ)) — на |Δ| это не
            # сказывалось (косметика Ж-отчёта), но имя вводило в заблуждение.
            t_hi_b, ok_hi_b = solve_age(q_meas * (1.0 + rel_total), qf)
            t_lo_b, ok_lo_b = solve_age(q_meas * (1.0 - rel_total), qf)
            ages.append((pair, t_hat, abs(t_hat - t_lo_b), abs(t_hi_b - t_hat),
                         q_meas, rel_meas, rel_mc, rel_total, ok,
                         pair in AGE_PAIRS))

    # --- ширина окна съёма: чувствительность результата (Б2) ---------------
    print("\n--- чувствительность к ширине окна (диагностика, не альтернатива) ---")
    print("  %10s %10s %12s" % ("×ПШПВ", "Q изм", "возраст, лет"))
    num0, den0 = AGE_PAIRS[0]
    fw_num0, fw_den0 = float(fwhm_fn(num0)), float(fwhm_fn(den0))
    width_pairs = []          # (ширина, возраст) — вместе, чтобы не потерять
                               # соответствие при пропуске неудачной ширины
    for wf in WIDTH_SCALES:
        rn = net_with_sigma(y, ybg, centres, num0, wf * fw_num0, exclude=windows)
        rd = net_with_sigma(y, ybg, centres, den0, wf * fw_den0, exclude=windows)
        if rn is None or rd is None or rd[0] == 0:
            print("  %10.2f  окно не помещается в сетку" % wf)
            continue
        q_w = rn[0] / rd[0]

        def q_w_fn(t_years, _wf=wf):
            A = activities_after_purification(t_years * YEAR_SEC, a0, hl, br_tl)
            a = model_net(A, num0, _wf * fw_num0)
            b = model_net(A, den0, _wf * fw_den0)
            if a is None or b is None or b == 0:
                return float("nan")
            return a / b

        t_w, ok_w = solve_age(q_w, q_w_fn)
        mark = "" if wf != 1.00 else "  <- принятая"
        print("  %10.2f %10.4f %12s%s"
              % (wf, q_w, ("%.2f" % t_w) if ok_w else "—", mark))
        if ok_w:
            width_pairs.append((wf, t_w))
    width_lo = min(t for _wf, t in width_pairs) if width_pairs else float("nan")
    width_hi = max(t for _wf, t in width_pairs) if width_pairs else float("nan")

    # --- проверка относительного отклика ------------------------------------
    print("\n--- проверка относительного отклика (от возраста не зависит) ---")
    print("  %-14s %10s %10s %10s  %s"
          % ("пара, кэВ", "измерено", "модель 2 л", "мод/изм", "интервал"))
    resp = {}
    for num, den, q_meas, _rel, vals in rows_out:
        if (num, den) not in CHECK_PAIRS:
            continue
        q_mod = vals[grid_years.index(2.0)]
        resp[(num, den)] = q_mod / q_meas
        print("  %-14s %10.4f %10.4f %10.3f  %.0f-%.0f кэВ"
              % ("%.0f/%.0f" % (num, den), q_meas, q_mod, q_mod / q_meas,
                 min(num, den), max(num, den)))

    # --- возраст ------------------------------------------------------------
    print("\n--- возраст с момента химической очистки ---")
    print("  %-14s %10s %12s %10s %10s  %s"
          % ("пара, кэВ", "Q изм", "возраст, лет", "±года", "МК, %", "статус"))
    for pair, t_hat, dlo, dup, q_meas, rel_meas, rel_mc, rel_total, ok, \
            is_primary in ages:
        note = "первичная" if is_primary else "перекрёстная проверка"
        print("  %-14s %10.4f %12.2f %10.2f %10.2f  %s"
              % ("%.0f/%.0f" % pair, q_meas, t_hat, 0.5 * (dlo + dup),
                 100.0 * rel_mc, note))

    primary_rows = [a for a in ages if a[9]]
    if not primary_rows:
        raise SystemExit("ни одна пара не пригодна для определения возраста")
    pair0, t_hat, dlo, dup, q0, rel_meas0, rel_mc0, rel_total0, _ok, _u = \
        primary_rows[0]
    stat = 0.5 * (dlo + dup)     # полная скобка счёт+МК вместе (rel_total0)
    qf0 = q_fn(pair0)

    # Приближённое разложение общей скобки stat на счётную и МК-часть —
    # линейно по вкладу в rel_total0 (не строгая квадратурная декомпозиция,
    # только для наглядности печати).
    part_meas = stat * rel_meas0 / rel_total0 if rel_total0 else 0.0
    part_mc = stat * rel_mc0 / rel_total0 if rel_total0 else 0.0
    print("\n  на паре %.0f/%.0f: из общей скобки ±%.2f года — счёт даёт "
          "~%.2f, МК-статистика шаблонов ~%.2f (точность отклика МК "
          "%.1f %% отн., разложение приближённое)"
          % (pair0[0], pair0[1], stat, part_meas, part_mc, 100.0 * rel_mc0))

    # --- систематика отклика: ДВЕ оценки, обе — экстраполяция ---------------
    # Прямой меры отклика на паре 238 -> 911 кэВ нет: обе линии Ac-228, а его
    # активность и есть искомая величина. Печатаются ОБЕ доступные оценки —
    # по соседней паре (67 % рабочего интервала) и по паре, кроющей рабочий
    # интервал целиком (238 -> 2614,51, шире 238 -> 911) — вместо того чтобы
    # выдавать одну из них за измеренную границу (находка Б4 внешнего
    # аудита: три альтернативных способа продолжить давали R вне диапазона,
    # который строился по одной лишь соседней паре).
    f_sys = resp.get(SYS_PAIR, float("nan"))
    f_cum = resp.get(CUM_PAIR, float("nan"))
    t_sys, ok_sys = solve_age(q0 * f_sys, qf0) if math.isfinite(f_sys) \
        else (float("nan"), False)
    t_cum, ok_cum = solve_age(q0 * f_cum, qf0) if math.isfinite(f_cum) \
        else (float("nan"), False)
    print("\n  систематика относительного отклика (обе оценки — экстраполяция,")
    print("  не измерение НА рабочей паре):")
    print("    по паре %.0f/%.0f (67 %% интервала): модель/изм = %.2f -> "
          "возраст %s"
          % (SYS_PAIR[0], SYS_PAIR[1], f_sys,
             ("%.1f года" % t_sys) if ok_sys else "вне восходящей ветви"))
    print("    по паре %.0f/%.0f (кроет интервал целиком): модель/изм = "
          "%.2f -> возраст %s"
          % (CUM_PAIR[0], CUM_PAIR[1], f_cum,
             ("%.1f года" % t_cum) if ok_cum else "вне восходящей ветви"))

    sys_candidates = [t for t, ok in ((t_sys, ok_sys), (t_cum, ok_cum)) if ok]
    t_sys_hi = max(sys_candidates) if sys_candidates else t_hat

    # Перекрёстная пара 968,97/238,63 — та же наблюдаемая величина, другая
    # (противоположная по стороне) подложка. Расхождение с первичной парой —
    # правдоподобно, наклон континуума, а не контаминация (см. комментарий
    # у CROSSCHECK_PAIRS) — входит в конверт наравне с шириной окна и
    # систематикой, а не отбрасывается.
    cross_ages = [a[1] for a in ages if not a[9]]

    # --- ИТОГ: точечная оценка + честный конверт неопределённости ----------
    stat_lo, stat_hi = t_hat - stat, t_hat + stat
    lo_candidates = [stat_lo]
    hi_candidates = [stat_hi, t_sys_hi]
    if math.isfinite(width_lo):
        lo_candidates.append(width_lo)
    if math.isfinite(width_hi):
        hi_candidates.append(width_hi)
    lo_candidates.extend(cross_ages)
    hi_candidates.extend(cross_ages)
    envelope_lo = min(lo_candidates)
    envelope_hi = max(hi_candidates)

    print("\n  ТОЧЕЧНАЯ ОЦЕНКА (пара %.0f/%.0f, окно ±1 ПШПВ, подложка без "
          "контаминации дублетом, счёт + статистика шаблонов):"
          % (pair0[0], pair0[1]))
    print("    %.2f ± %.2f года" % (t_hat, stat))
    if cross_ages:
        print("  Перекрёстная пара %.0f/%.0f (та же величина, подложка на "
              "противоположной стороне): %s года"
              % (CROSSCHECK_PAIRS[0][0], CROSSCHECK_PAIRS[0][1],
                 ", ".join("%.2f" % t for t in cross_ages)))
    print("  ОБЩИЙ ДИАПАЗОН (объединяет счёт, статистику шаблонов, выбор")
    print("  ширины окна, перекрёстную пару и обе оценки систематики отклика;")
    print("  не квадратурная сумма — источники не гауссовы и не независимы):")
    print("    %.1f...%.1f года" % (envelope_lo, envelope_hi))
    print("  Точную границу систематика отклика построить не позволяет")
    print("  (см. ниже); качественный вывод — единицы лет, а не десятки —")
    print("  устойчив ко всем перечисленным источникам разброса.")
    print("\n  замер 01.06.2024 -> очистка тория, точечно: %.1f года назад"
          % t_hat)
    print("  [очистка предшествует изготовлению электродов; промежуток между")
    print("   ними этой задачей не определяется]")

    # --- устойчивость к 20-летнему допущению NUREG/CR-1039 -----------------
    # Не решается через solve_age (кривая на больших t уходит за восходящую
    # ветвь — прежняя редакция путала «результат solve_age при скорректиро-
    # ванном Q» с «прямым значением модели при t=20», это разные вещи и
    # давали разные числа; см. Б5 внешнего аудита). Здесь — только прямое
    # сравнение: во сколько раз должен ошибаться отклик, и укладывается ли
    # это в диапазон, который показывают проверочные пары НА СОПОСТАВИМЫХ
    # интервалах, а не по максимуму среди вообще всех пар (среди которых
    # есть заведомо непригодная 968,97/911,20).
    q20 = qf0(20.0)
    f20 = q20 / q0
    print("\n  устойчивость к допущению «20 лет» (NUREG/CR-1039, 1979, срок")
    print("  между разделением и продажей — допущение ИХ дозовой оценки,")
    print("  не измерение):")
    print("    при t=20 лет модель даёт Q = %.4f -> в %.2f раза больше "
          "измеренного" % (q20, f20))
    print("    для сравнения, на проверочных парах отклик расходится в "
          "%.2f (%.0f/%.0f) и %.2f раза (%.0f/%.0f)"
          % (f_sys, SYS_PAIR[0], SYS_PAIR[1], f_cum, CUM_PAIR[0], CUM_PAIR[1]))
    if f20 <= max(f_sys, f_cum):
        print("    %.2f лежит В ПРЕДЕЛАХ уже наблюдаемого разброса отклика —" % f20)
        print("    аргументом «отклик так не ошибается» 20 лет строго не")
        print("    исключить.")
    else:
        print("    %.2f ВЫШЕ обеих наблюдаемых оценок — но их всего две, и")
        print("    строгой границей они не являются.")
    resp_911_969 = resp.get((968.97, 911.20), float("nan"))
    print("    Отдельно: пара 968,97/911,20 (обе линии Ac-228, интервал")
    print("    58 кэВ) даёт мод/изм = %.3f — МЕНЬШЕ, чем у пар с бо́льшим "
          "энергетическим" % resp_911_969)
    print("    охватом (%.2f на 583/239, %.2f на 2615/239). Расхождение "
          "отклика растёт не гладко" % (f_sys, f_cum))
    print("    с энергией: продолжать его на 911 кэВ единой экстраполяцией "
          "формально нечем")
    print("    (эта же пара к тому же сама опирается на две ОДНОСТОРОННИЕ, "
          "взаимно")
    print("    исключающие подложки — см. Б1, — так что как мера отклика "
          "она менее надёжна,")
    print("    чем более широкие проверочные пары). Довод в пользу единиц "
          "лет остаётся")
    print("    прямым решением по измеренному Q, а не аргументом от "
          "систематики отклика.")

    # --- запись -------------------------------------------------------------
    p = os.path.join(outdir, "forward_age.csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["# возраст ряда по отношению площадей опорных пиков;"])
        w.writerow(["# площадь = сумма в +-1 ПШПВ минус подложка по внешним "
                    "подокнам с ВЗАИМНЫМ ИСКЛЮЧЕНИЕМ соседних опорных пиков "
                    "(Б1), одинаково у модели и измеренного"])
        w.writerow(["# калибровка: %s" % calib_path,
                    "найдена" if corr else "НЕ НАЙДЕНА, заводская шкала"])
        if corr:
            w.writerow(["# коэффициенты образца", str(corr.get("sample"))])
            w.writerow(["# коэффициенты фона", str(corr.get("background"))])
        w.writerow(["E_числитель_кэВ", "E_знаменатель_кэВ", "назначение",
                    "Q_измерено", "сигма_Q_счёт_отн", "сигма_Q_МК_отн",
                    "возраст_лет"]
                   + ["Q_модель_%gлет" % t for t in grid_years])
        by_pair = {a[0]: a for a in ages}
        for num, den, q_meas, rel_meas, vals in rows_out:
            a = by_pair.get((num, den))
            role = ("возраст_первичная" if (num, den) in AGE_PAIRS
                    else "возраст_перекрёстная" if (num, den) in CROSSCHECK_PAIRS
                    else "проверка отклика")
            w.writerow(["%.2f" % num, "%.2f" % den, role,
                        "%.4f" % q_meas, "%.4f" % rel_meas,
                        "%.4f" % a[6] if a else "",
                        "%.2f" % a[1] if a else ""]
                       + ["%.4f" % v for v in vals])
        w.writerow([])
        w.writerow(["# чувствительность к ширине окна съёма (Б2)"])
        w.writerow(["×ПШПВ", "возраст_лет"])
        for wf, t_w in width_pairs:
            w.writerow(["%.2f" % wf, "%.2f" % t_w])
        w.writerow([])
        w.writerow(["# состав окна у модели при возрасте 2 года, доля площади"])
        w.writerow(["E_кэВ"] + [k for k, _l, _e, _f, _b in CHAIN
                                if k in tmap])
        for e0, _nuc, _br in PEAKS:
            if e0 not in parts_by_peak:
                continue
            parts, tot = parts_by_peak[e0]
            w.writerow(["%.2f" % e0]
                       + ["%.3f" % (parts.get(k, 0.0) / tot) if tot else ""
                          for k, _l, _e, _f, _b in CHAIN if k in tmap])
        w.writerow([])
        w.writerow(["итог_точка_лет", "%.2f" % t_hat,
                    "итог_счётная_сигма_лет", "%.2f" % stat,
                    "итог_диапазон_от", "%.2f" % envelope_lo,
                    "итог_диапазон_до", "%.2f" % envelope_hi,
                    "систематика_соседняя_пара", "%.2f" % t_sys if ok_sys else "н/д",
                    "систематика_кумулятивная_пара",
                    "%.2f" % t_cum if ok_cum else "н/д",
                    "перекрёстная_пара",
                    ", ".join("%.2f" % t for t in cross_ages) if cross_ages
                    else "н/д"])
    print("\nзаписано: %s" % p)


if __name__ == "__main__":
    main()
