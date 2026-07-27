"""Независимая проверка кривой: эффективность прямо по ИЗМЕРЕННЫМ спектрам.

Кривая .efr — уже обработанный ЛСРМ результат. Здесь тот же ответ получается
в обход неё: площадь пика из спектра комплекта, минус измеренный фон той же
геометрии, делённая на активность источника по паспорту с пересчётом на дату
измерения и на выход линии. Если два пути дадут разное, ошибка в разборе .efr;
если одно и то же — сверка МК с экспериментом опирается не на один источник.

Активности и даты — из описи комплекта (README набора) и заголовков XML.
Выходы линий берутся из ТОГО ЖЕ прогона Geant4, что и эффективность
(файлы *_emit.csv), а не из справочника.
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

# Файлы комплекта лежат по подкаталогам нуклидов, поэтому ищем по имени
# (paths.find_data), а не склеиваем путь.
BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)

# Опись комплекта: файл -> (нуклид, уд. активность Бк/кг, погр. %, дата паспорта,
#                           масса г, период полураспада в годах, линии кэВ)
# Периоды: Cs-137 30,08 г; K-40 1,248e9 г; Ra-226 1600 г; Th-232 1,405e10 г.
KIT = {
    "sample_M_cs_легкий_2001-2005.xml":
        ("Cs-137", 1890.0, 5.0, "1997-05-30", 570.0, 30.08, [661.657]),
    "sample_M_k_легкий_2001-2005.xml":
        ("K-40", 2540.0, 10.0, None, 665.0, 1.248e9, [1460.822]),
    "sample_M_ra_легкий_2001-2007.xml":
        ("Ra-226", 1850.0, 10.0, None, 622.0, 1600.0,
         [295.223, 351.932, 609.32, 1120.294, 1764.491]),
    "Th232_420-7-17_Маринелли_0cm.xml":
        ("Th-232", 1940.0, 6.0, "2007-09-17", 1600.0, 1.405e10,
         [238.632, 583.187, 911.204, 2614.511]),
}


def days(a, b):
    from datetime import date
    ya, ma, da = (int(x) for x in a.split("-"))
    yb, mb, db = (int(x) for x in b.split("-"))
    return (date(yb, mb, db) - date(ya, ma, da)).days


def emit_yield(E):
    """Выход линии на распад — из прогонов распада того же расчёта."""
    for f in glob.glob(os.path.join(BUILD, "decay_*_emit.csv")):
        tot, N = 0, None
        for line in open(f, encoding="utf-8"):
            if line.startswith("#"):
                if "N_primaries" in line:
                    N = int(line.split("=")[1])
                continue
            if line and line[0].isdigit():
                e, c = line.split(",")
                if abs(float(e) - E) <= 2.0:
                    tot += int(c)
        if N and tot > 50:
            return tot / N, os.path.basename(f)
    return None, None


def mc_eff(E):
    best = None
    for p in glob.glob(os.path.join(BUILD, "grid", "rho1.60_E*.csv")):
        E0, net, dnet, N = read_run(p)
        if abs(E0 - E) < 1.0:
            best = (net / N, dnet / N)
    return best


if __name__ == "__main__":
    print("Эффективность ППП прямо по измеренным спектрам комплекта\n")
    print("%-8s %9s %8s %10s %11s %11s %8s" %
          ("нуклид", "E, кэВ", "ПШПВ %", "имп/с", "eps_эксп", "eps_МК", "МК/эксп"))
    for fn, (nuc, aspec, dpct, dt0, mass, _t12, lines) in KIT.items():
        path = paths.find_data(fn)
        if path is None:
            print("нет файла:", fn)
            continue
        path = str(path)
        if not os.path.exists(path):
            print("нет файла:", fn)
            continue
        s, b = bm.read(path)
        mdate = None
        for line in open(path, encoding="utf-8"):
            if "<StartTime>" in line:
                mdate = line.split(">")[1][:10]
                break
        # распад к дате измерения
        k = 1.0
        if dt0 and mdate:
            k = 0.5 ** (days(dt0, mdate) / 365.25 / _t12)
        A = aspec * (mass / 1000.0) * k          # Бк в сосуде
        for E in lines:
            fw = bm.fwhm_at(s, E)
            if fw is None or fw <= 0:
                print("%-8s %9.3f   пик не найден" % (nuc, E))
                continue
            r = bm.net_rate(s, b, E, fw)
            if r is None:
                print("%-8s %9.3f   ROI вне спектра" % (nuc, E))
                continue
            rate, drate, rbg, _ = r
            pg, src = emit_yield(E)
            if pg is None:
                print("%-8s %9.3f %8.1f %10.4f   нет выхода линии" %
                      (nuc, E, 100 * fw / E, rate))
                continue
            eps = rate / (A * pg)
            deps = eps * math.sqrt((drate / max(rate, 1e-9)) ** 2 + (dpct / 100) ** 2)
            mc = mc_eff(E)
            if mc is None:
                print("%-8s %9.3f %8.1f %10.4f %11.4e" % (nuc, E, 100 * fw / E, rate, eps))
                continue
            print("%-8s %9.3f %8.1f %10.4f %11.4e %11.4e %8.3f"
                  % (nuc, E, 100 * fw / E, rate, eps, mc[0], mc[0] / eps))
        print("   %s: A = %.0f Бк (паспорт %.0f Бк/кг x %.3f кг x распад %.4f),"
              " живое %.0f с, фон %.0f с"
              % (nuc, A, aspec, mass / 1000, k, s.live, b.live if b else 0))
