"""Трек 1b: приборная калибровка ПШПВ против модельного закона 49,9*sqrt(E).

Три закона ширины на одном приборе:
1. МОДЕЛЬ конвейера: FWHM(E) = 49,9*sqrt(E/661,657) — одно число 49,9,
   снятое по пику 662 кэВ записи цезия комплекта; на нём стоят все окна
   площадей и классификация линий «чистая/бленд».
2. Calibr.cfw (01.11.2024, каталог Data рабочего дерева) — реальная калибровка
   по 25 пикам 11 источников поверки 2024 на 5 см; полином 3-й степени
   ПО z = sqrt(E) (конвенция ЛСРМ, BUG-22 в SpectraVibe) плюс сами
   измеренные точки (E, ПШПВ) — используются здесь как опора.
3. Полином из заголовков .spe поверки (коэффициенты 8,969/-0,898/0,143/
   -0,0017): его происхождение подозрительно — FWHMCALIBRATIONFILE этих
   .spe указывает на файлы ДРУГОГО прибора (рассинхрон серийников в
   конфиге, известная ловушка), и именно ЕГО значение (88,6 кэВ на 1437)
   показывал статус-бар СпектраЛайн в сеансе оператора 29.07.2026.

4. Калибровка ПО САМОМУ спектру Th232_420-7-17_Маринелли_0cm (файл .cfw
   оператора, 29.07.2026): 9 пиков этого спектра, полином той же конвенции.
   Геометрически родная ширина — опора для размытия МК-спектра маринелли
   в Треке 1c. Точки и полином вшиты числами (файл лежал вне SpectraVibe).

Вопросы, на которые отвечает таблица: насколько модельный закон врёт по
диапазону (окна и чистота), и какой из приборных законов активен в ПО.
Ответ по данным: in-situ точки маринелли согласны с Calibr.cfw в ±5-9 %,
«чужой» spe-полином завышает ширину до +21 % в середине (911 кэВ) —
активный в ПО закон не совпадает с реальными пиками этого спектра.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402

CFW = os.path.join("detectors", "Gamma-1S", "raw_lsrm", "Work", "BG",
                   "Gamma-1S", "Data", "Calibr.cfw")
SPE_POLY = (8.969108040517, -0.8981588537219, 0.1433311097464,
            -0.001693803867)      # из заголовков .spe поверки 2024
# Калибровка по пикам самого спектра Th232_420-7-17_Маринелли_0cm
# (.cfw оператора, 29.07.2026): полином по sqrt(E) и 9 измеренных точек.
MAR_POLY = (-19.9307565267979, 3.4810472847871, -0.0516136278516,
            0.0006473958987)
MAR_PTS = [(237.836, 24.563), (337.277, 30.331), (510.769, 41.344),
           (583.188, 45.907), (727.157, 45.905), (911.03, 54.222),
           (964.594, 56.463), (968.797, 56.637), (2616.27, 109.832)]
FWHM662 = 49.9

# опорные энергии: узлы сеток + края
ENERGIES = [59.5, 88.0, 122.1, 165.9, 238.632, 351.932, 583.187, 661.657,
            911.204, 1120.294, 1460.822, 1764.491, 2614.511, 3000.0]


def poly_sqrtE(coefs, E):
    z = math.sqrt(E)
    return sum(c * z ** k for k, c in enumerate(coefs))


def read_cfw(path):
    txt = open(path, "rb").read().decode("cp1251")
    coefs, pts = None, []
    sec = None
    for ln in txt.splitlines():
        ln = ln.strip()
        if ln.startswith("["):
            sec = ln
            continue
        if sec == "[Calibration]" and ln.startswith("Coeff="):
            coefs = tuple(float(x) for x in ln.split("=", 1)[1].split(","))
        if sec == "[CalibrationData]" and "=" in ln:
            e, f = ln.split("=")
            pts.append((float(e), float(f)))
    return coefs, pts


if __name__ == "__main__":
    root = paths.require_spectravibe("чтение Calibr.cfw (приборная ПШПВ)")
    coefs, pts = read_cfw(os.path.join(str(root), CFW))
    print("Calibr.cfw: полином 3-й степени по sqrt(E), %d измеренных точек\n"
          % len(pts))
    print("%9s %9s %9s %9s %9s %11s %11s" %
          ("E, кэВ", "модель", "cfw", "spe-полин", "маринелли",
           "модель/cfw", "spe/cfw"))
    rows = []
    for E in ENERGIES:
        m = FWHM662 * math.sqrt(E / 661.657)
        c = poly_sqrtE(coefs, E)
        s = poly_sqrtE(SPE_POLY, E)
        g = poly_sqrtE(MAR_POLY, E)
        print("%9.1f %9.2f %9.2f %9.2f %9.2f %11.3f %11.3f"
              % (E, m, c, s, g, m / c, s / c))
        rows.append((E, m, c, s, g, m / c, s / c))

    # контроль полинома против измеренных точек самого файла
    worst = max(abs(poly_sqrtE(coefs, e) - f) / f for e, f in pts)
    print("\nконтроль: полином против %d точек файла, худшее отклонение"
          " %.1f%%" % (len(pts), 100 * worst))
    print("in-situ маринелли против законов (точка / cfw / spe):")
    for e, f in MAR_PTS:
        print("   %8.1f  %6.2f  %6.2f (%+5.1f%%)  %6.2f (%+5.1f%%)"
              % (e, f, poly_sqrtE(coefs, e),
                 100 * (poly_sqrtE(coefs, e) / f - 1),
                 poly_sqrtE(SPE_POLY, e),
                 100 * (poly_sqrtE(SPE_POLY, e) / f - 1)))

    out = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results",
        "fwhm_check.csv"))
    csvio.write(
        out,
        ["E_keV", "fwhm_model_keV", "fwhm_cfw_keV", "fwhm_spe_keV",
         "fwhm_marinelli_keV", "model_over_cfw", "spe_over_cfw"],
        [("%.3f" % e, "%.2f" % m, "%.2f" % c, "%.2f" % s, "%.2f" % g,
          "%.4f" % mc, "%.4f" % sc) for e, m, c, s, g, mc, sc in rows],
        comments=[
            "Три закона ПШПВ: модель конвейера 49,9*sqrt(E/661,657);"
            " Calibr.cfw (реальная калибровка 11 источников, 01.11.2024,"
            " полином по sqrt(E));",
            "полином из заголовков .spe поверки (происхождение — файлы"
            " другого прибора, рассинхрон конфига; именно он был активен"
            " в сеансе оператора).",
            "model_over_cfw > 1 — модельное окно шире реального пика.",
        ])
    print("таблица: %s" % out)
