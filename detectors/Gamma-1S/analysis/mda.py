"""МИА установки Гамма-1С против паспорта — по ИЗМЕРЕННОМУ фону.

    МИА = 4*sqrt(2) * sqrt(n_ф) / (sqrt(t) * eta),   eta = eps_ППП * p_gamma
[Методика ЛСРМ «Активность в счётных образцах», 2014, прил. 5, с. 16].
Коэффициент именно 4*sqrt(2) = 5,657 — сверено по оригинальному скану;
в курированной записи RAG-библиотеки стоит ошибочное 4,42.

n_ф — скорость счёта ИЗМЕРЕННОГО фона в окне пика (имп/с). Фон берётся из
того же XML комплекта: там вложен фон той же геометрии (маринелли с водой),
54 000 с. Окно пика — та же процедура, что и для площадей: +-1,25 ПШПВ,
ПШПВ меряется по реальному пику, а не по полиному из заголовка.

eps — РАСЧЁТНАЯ эффективность, посчитанная ПРЯМО В ВОДЕ (сетка grid/water1.00),
а не пересчитанная с ОИСН-16: у ОИСН-16 71 % железа по массе, и перенос
эффективности между такими разными матрицами формулой был бы экстраполяцией.
p_gamma — выход линии НА РАСПАД РОДИТЕЛЯ, из прогонов цепочек того же расчёта.

Паспорт ДЦКИ.412131.001 ПС, п. 2.2: маринелли 1 л с дистиллированной водой,
2 ч: Cs-137 1,5; K-40 25; Ra-226 3; Th-232 3 Бк/кг.

ОГОВОРКА (замечание оператора): МИА Ra-226 считается по линии Bi-214 609,3 и
предполагает равновесие с радоном, то есть верна в меру УДЕРЖАНИЯ радона
пробой. Вода радон удерживает хорошо, сыпучие и аэрируемые пробы — нет.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402
from compare_lsrm import read_run  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)
REF = os.path.join(str(paths.ref("Gamma-1S")), "kit", "Marinelli_1L")
KFACT = 4.0 * math.sqrt(2.0)          # 5,657
TMEAS = 2.0 * 3600.0                  # паспортные 2 часа
MASS_KG = 1.0                         # 1 л дистиллированной воды

# нуклид -> (аналитическая линия, файл прогона для выхода на распад родителя,
#            паспортная МИА Бк/кг)
LINES = [
    ("Cs-137", 661.657, "decay_Cs137_emit.csv", 1.5),
    ("K-40", 1460.822, "decay_K40_emit.csv", 25.0),
    ("Ra-226", 609.32, "chain_Ra226_emit.csv", 3.0),
    ("Th-232", 2614.511, "chain_Th232_emit.csv", 3.0),
    ("Th-232", 583.187, "chain_Th232_emit.csv", 3.0),
]


def yield_per_parent(fn, E):
    p = os.path.join(BUILD, fn)
    if not os.path.exists(p):
        return None
    tot, N = 0, None
    for line in open(p, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            continue
        if line and line[0].isdigit():
            e, c = line.split(",")
            if abs(float(e) - E) <= 2.0:
                tot += int(c)
    return (tot / N) if N else None


def mc_eff(E, tag="water1.00"):
    for p in glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv")):
        E0, net, dnet, N = read_run(p)
        if abs(E0 - E) < 1.0:
            return net / N, dnet / N
    return None


if __name__ == "__main__":
    # фон: вложен в любой файл комплекта Маринелли
    src = paths.find_data("sample_M_cs_легкий_2001-2005.xml")
    if src is None:
        raise SystemExit("не найден спектр комплекта в %s"
                         % paths.ref("Gamma-1S"))
    src = str(src)
    smp, bg = bm.read_checked(src)[:2]
    print("Фон: %s, живое %.0f с, каналов %d\n"
          % ("маринелли с водой (вложен в XML комплекта)", bg.live, len(bg.n)))
    print("%-8s %9s %7s %9s %9s %9s %8s %8s %8s" %
          ("нуклид", "E, кэВ", "ПШПВ%", "n_ф,имп/с", "eps", "p_gamma",
           "МИА", "МДА_К", "паспорт"))
    # ПШПВ: мерить по пику можно только там, где он есть в опорной записи
    # (662 в записи цезия). Для остальных энергий — масштабирование
    # w(E) = w662*sqrt(E/662), стандартный ход разрешения сцинтиллятора.
    w662 = bm.fwhm_at(smp, 661.657)
    for nuc, E, emitf, ref in LINES:
        fw = w662 * math.sqrt(E / 661.657)
        a = bm.peak_area(bg, E, fw)
        if a is None:
            print("%-8s %9.3f   окно вне спектра фона" % (nuc, E))
            continue
        gross, _, _ = a
        # для МИА нужен ПОЛНЫЙ счёт фона в окне, не за вычетом подложки
        w = 1.25 * fw
        graw, _ = bg.counts_between(E - w, E + w)
        nf = graw / bg.live
        eps = mc_eff(E)
        pg = yield_per_parent(emitf, E)
        if eps is None or pg is None:
            print("%-8s %9.3f %7.1f %9.4f   нет %s"
                  % (nuc, E, 100 * fw / E, nf,
                     "eps" if eps is None else "выхода (%s)" % emitf))
            continue
        eta = eps[0] * pg
        mda = KFACT * math.sqrt(nf) / (math.sqrt(TMEAS) * eta) / MASS_KG
        # Вариант Карри/ГОСТ: (2,71 + 4,65*sqrt(B))/(t*eta), B — счёт фона
        B = nf * TMEAS
        mda_c = (2.71 + 4.65 * math.sqrt(B)) / (TMEAS * eta) / MASS_KG
        print("%-8s %9.3f %7.1f %9.4f %9.4e %9.4f %8.2f %8.2f %8.1f"
              % (nuc, E, 100 * fw / E, nf, eps[0], pg, mda, mda_c, ref))
    print("\nt = %.0f с (2 ч), масса пробы %.1f кг" % (TMEAS, MASS_KG))
    print("МИА — формула ЛСРМ, 4*sqrt(2) = 5,657; МДА_К — Карри/ГОСТ, "
          "(2,71+4,65*sqrt(B))/(t*eta).")
    print("Эффективность здесь РАСЧЁТНАЯ (шкала МК). В шкале эксперимента "
          "её надо делить на 1,171 — обе колонки тогда x1,171.")
