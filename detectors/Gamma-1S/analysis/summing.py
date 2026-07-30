"""Поправка на каскадное суммирование в сосуде Маринелли 1 л.

Идея. В моноэнергетическом прогоне в кристалл приходит один квант, суммироваться
не с чем — это «чистая» эффективность ППП eps_mono(E). В прогоне полного распада
кванты каскада летят одновременно, и если в кристалле сработали два, событие
уходит из своего пика. Отсюда

    eps_decay(E) = A_peak(E) / N_emit(E),      C(E) = eps_mono(E) / eps_decay(E)

где N_emit(E) — число квантов этой энергии, РЕАЛЬНО испущенных за прогон
(файл *_emit.csv). Выход линии p_gamma нигде не вписывается руками: он приходит
из той же базы PhotonEvaporation, что и транспорт. Это принципиально — на
числах из вторых рук в этом проекте уже обжигались.

Проверка метода: Cs-137 и K-40 каскада не имеют, для них C обязан выйти 1,00
(макрос decay_control.mac).

ОГОВОРКА. Geant4 печатает при старте «Enable correlated gamma emission 0»:
угловые корреляции между квантами каскада по умолчанию выключены, кванты
разыгрываются изотропно независимо друг от друга. В маринелльной геометрии
проба охватывает детектор почти со всех сторон, и вклад корреляций мал, но
для точечной геометрии это было бы существеннее.
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
import peakwin  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_lsrm import read_run, WIN, BG0, BG1  # noqa: E402

# Линии, по которым ЛСРМ строил кривую, с указанием прогона распада
LINES = {
    # контроль: каскада нет, C обязан выйти 1,00
    "decay_Cs137": [661.657],
    "decay_K40": [1460.822],
    # каскадные
    "decay_Tl208": [583.187, 2614.511],
    "decay_Bi214": [609.32, 768.36, 1120.294, 1764.491],
}


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


def area(hist, E0, sub=True):
    """Площадь пика с вычетом континуума по левой полке."""
    gross = sum(c for e, c in hist.items() if abs(e - E0) <= WIN)
    side = sum(c for e, c in hist.items() if E0 - BG0 <= e <= E0 - BG1)
    # Центры каналов полуцелые ((i+0,5)·bin в main.cc), поэтому окно шириной
    # 2·WIN содержит 2·WIN каналов, а не 2·WIN+1: считаем по факту, иначе
    # подложка вычитается с лишними 8 %.
    n = math.floor(E0 + WIN - 0.5) - math.ceil(E0 - WIN - 0.5) + 1
    nside = math.floor(E0 - BG1 - 0.5) - math.ceil(E0 - BG0 - 0.5) + 1
    bg = side / nside * n if sub else 0.0
    # D(bg) = (n/nside)^2 * side = (n/nside)*bg; вывод — в export_curves.py
    return gross - bg, math.sqrt(max(gross + (n / nside) * bg, 1.0))


def mono_curve(tag="rho1.60"):
    out = {}
    for p in sorted(glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv"))):
        E0, net, dnet, N = read_run(p)
        if net > 0:
            out[round(E0, 3)] = (net / N, dnet / N)
    return out


def area_var(hist, E0, win, bg0, bg1):
    """Площадь при заданных окне и полке — для проверки устойчивости.

    Единственное место, где параметры окна ЗАКОННО задаются по вызову: смысл
    функции в том, чтобы их менять. Счёт всё равно идёт через `peakwin`, иначе
    сканирование устойчивости мерило бы ещё и разницу двух реализаций окна —
    ровно ту, из-за которой правка окна полки не дошла до объёмного пересчёта.
    """
    if bg0 is None:
        i0, i1, _j0, _j1 = peakwin.channels(E0, win_keV=win)
        return sum(c for e, c in hist.items()
                   if max(0, i0) <= e - 0.5 < i1)
    return peakwin.area(hist, E0, win_keV=win, bg0_keV=bg0, bg1_keV=bg1)


def sensitivity(mono):
    """Насколько C зависит от способа снятия пика.

    В спектре РАСПАДА континуум под линией гораздо богаче, чем в
    моноэнергетическом, и его структура сложнее. Если C сильно поедет при
    смене окна или полки, значит поправка меряет не совпадения, а произвол
    вычитания, и верить ей нельзя.
    """
    VAR = [("окно ±6, полка 30-10", 6.0, 30.0, 10.0),
           ("окно ±6, полка 60-20", 6.0, 60.0, 20.0),
           ("окно ±10, полка 40-15", 10.0, 40.0, 15.0),
           ("окно ±6, БЕЗ вычета", 6.0, None, None)]
    print("\n--- устойчивость C к способу снятия пика ---")
    print("%-8s %9s %s" % ("нуклид", "E, кэВ",
                           "".join("%22s" % v[0] for v in VAR)))
    for tag, lines in LINES.items():
        sp = os.path.join(BUILD, tag + ".csv")
        em = os.path.join(BUILD, tag + "_emit.csv")
        if not (os.path.exists(sp) and os.path.exists(em)):
            continue
        hist, N = load(sp)
        emit, _ = load(em)
        for E in lines:
            nem = sum(c for e, c in emit.items() if abs(e - E) <= 2.0)
            if nem == 0:
                continue
            key = min(mono, key=lambda k: abs(k - E))
            if abs(key - E) > 1.0:
                continue
            row = ""
            for _, w, b0, b1 in VAR:
                # моно-точку снимаем ТЕМ ЖЕ способом — иначе сравнение нечестно
                mp = os.path.join(BUILD, "grid", "rho1.60_E%07.1f.csv" % key)
                if not os.path.exists(mp):
                    row += "%22s" % "-"
                    continue
                mh, mN = load(mp)
                em_ = area_var(mh, key, w, b0, b1) / mN
                ed = area_var(hist, E, w, b0, b1) / nem
                row += "%22.3f" % (em_ / ed if ed > 0 else float("nan"))
            print("%-8s %9.3f %s" % (tag[6:], E, row))


if __name__ == "__main__":
    sub = "--probe" in sys.argv
    root = os.path.join(BUILD, "probe") if sub else BUILD
    mono = mono_curve()

    print("Поправка на каскадное суммирование, Маринелли 1 л, ОИСН-16 1,6\n")
    print("%-8s %9s %9s %10s %10s %7s" %
          ("нуклид", "E, кэВ", "p_gamma", "eps_расп", "eps_моно", "C"))
    res = {}
    for tag, lines in LINES.items():
        sp = os.path.join(root, tag + ".csv")
        em = os.path.join(root, tag + "_emit.csv")
        if not (os.path.exists(sp) and os.path.exists(em)):
            print("  нет файлов для", tag)
            continue
        hist, N = load(sp)
        emit, Ne = load(em)
        assert N == Ne, "разное число распадов в спектре и в испускании"
        for E in lines:
            # испущено квантов этой энергии: линии в испускании острые
            nem = sum(c for e, c in emit.items() if abs(e - E) <= 2.0)
            if nem == 0:
                print("%-8s %9.3f   линия не испускается" % (tag[6:], E))
                continue
            pg = nem / N
            a, da = area(hist, E)
            ed = a / nem
            ded = ed * math.sqrt((da / max(a, 1)) ** 2 + 1.0 / nem)
            key = min(mono, key=lambda k: abs(k - E))
            if abs(key - E) > 1.0:
                print("%-8s %9.3f %9.4f %10.4e   нет моно-точки" % (tag[6:], E, pg, ed))
                continue
            em_, dem = mono[key]
            C = em_ / ed
            dC = C * math.sqrt((ded / ed) ** 2 + (dem / em_) ** 2)
            res[E] = (C, dC)
            mark = ""
            if tag in ("decay_Cs137", "decay_K40"):
                mark = "  <- контроль, ждём 1,00" if abs(C - 1) < 3 * dC \
                       else "  <- КОНТРОЛЬ НЕ ПРОШЁЛ"
            print("%-8s %9.3f %9.4f %10.4e %10.4e %7.3f ± %.3f%s"
                  % (tag[6:], E, pg, ed, em_, C, dC, mark))

    sensitivity(mono)

    print("\nC > 1 — потеря отсчётов из пика на совпадениях;")
    print("C < 1 — приход в пик суммой двух квантов (сумм-пик поверх линии).")
    print("\nЧто это даёт для сверки: измеренная ЛСРМ эффективность каскадных")
    print("линий занижена ровно на этот множитель, значит отношение МК/эксп")
    print("надо делить на C, прежде чем судить о форме кривой.")
