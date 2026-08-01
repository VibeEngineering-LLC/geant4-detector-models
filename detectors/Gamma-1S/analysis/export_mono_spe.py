"""Экспорт моно-МК-спектра в .spe СпектраЛайн для измерения моста (задача 115).

ЗАЧЕМ. Мост методик — во сколько раз аттестационная конвенция съёма (гаусс-
фит со ступенькой и полиномом фона) занижает площадь относительно истинного
полнопоглощённого пика — нельзя надёжно померить своим питон-фитом: он даёт
энергозависимые артефакты (переоценка на низких E, недооценка на 2614 из-за
линейного фона). Правильный ход (предписание аудитора): не воспроизводить фон
СпектраЛайн у себя, а ЭКСПОРТИРОВАТЬ моно-МК-спектр в формат прибора и
прогнать его через СпектраЛайн на параметрах слепка. Тогда B_lsrm(E) =
площадь_фита_СпектраЛайн / N_истинный_пик измеряется в потребительской
конвенции напрямую.

КАК. Депозит-спектр моно-линии полной геометрии (scat_p5_full_E*, 1-кэВ бины)
размывается СОХРАНЁННОЙ ПШПВ-калибровкой прибора (Calibr.cfw, не модельным
законом), пересыпается на канальную сетку прибора по его же энергетической
калибровке и пишется .spe штатным писателем тулкита (write_lsrm_spe).
Заголовок несёт калибровки и метку детектора №GEANT (не заводской номер).

ПРЕДСКАЗАНИЕ ЗАФИКСИРОВАНО ДО ПРОГОНА (задача 115): остаток истинного
завышения модели = B_lsrm(E) * (eps_наша/eps_штат) - 1, где отношение из
hard_edge (2614 -> 1,275; плато -> 1,078). Ожидание: плато 6-7 %, 2614
7-12 %. Порог B_lsrm(2614) ~ 0,84 разделяет «разрыв полностью конвенционный»
и «есть остаточное завышение на жёстком крае».

N_истинный_пик каждой линии печатается — оператор делит площадь фита
СпектраЛайн на него и получает B_lsrm.

Калибровки и канальная сетка берутся из реального измеренного .spe точечной
5 см (через SPECTRAVIBE_ROOT); моно-спектры модели — из рабочего каталога
(G4MODELS_BUILD_GAMMA_1S). Оба вне репозитория.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import detector_params as dp  # noqa: E402

# Моно-линии и их спектры. ИСТОЧНИК СМЕНЁН 01.08.2026: раньше здесь стояли
# `scat_p5_full_E*.csv` — прогон 28.07, ОТОЗВАННЫЙ задачей 126 (снят до
# правки входного торца). Файлы .spe, лежащие в рабочем каталоге прибора с
# 30.07.2026, сделаны из них и потому недействительны: их истинные площади
# (5482 отсчёта на 2614,5) не совпадают с текущими (5946), и прогон таких
# файлов сравнивал бы разные события.
#
# Теперь берутся `results/bridge_spectra/bridge4pi_E*.csv` — те самые
# спектры, по которым считается `bridge_mono.csv`, закоммиченные в
# репозиторий. Значит числитель (площадь, снятая программой) и знаменатель
# (истинный пик) принадлежат одному прогону, а измеренный мост прямо
# сопоставим с модельным из той же таблицы — что и составляет смысл
# проверки.
MONO = {122.100: "bridge4pi_E0122.1.csv",
        165.900: "bridge4pi_E0165.9.csv",
        583.187: "bridge4pi_E0583.2.csv",
        661.657: "bridge4pi_E0661.7.csv",
        727.330: "bridge4pi_E0727.3.csv",
        860.557: "bridge4pi_E0860.6.csv",
        2614.511: "bridge4pi_E2614.5.csv"}
GEOM = "Точечная-5см"
DET = os.environ.get("G1S_LSRM_DETECTOR", "УДС-ГЦ-63х63-USB №GEANT")


def poly(coefs, x):
    return sum(a * x ** i for i, a in enumerate(coefs))


def load_mc(path):
    ec, N = [], None
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"):
            if "N_primaries" in ln:
                N = int(ln.split("=")[1])
            continue
        if ln.strip() and not ln.startswith("E_keV"):
            e, c = ln.split(",")
            ec.append((float(e), float(c)))
    return ec, N


def true_peak(ec, E0, half=3.0):
    return sum(c for e, c in ec if abs(e - E0) < half)


def broaden_to_channels(ec, ecal, fwhm_cal, nchan):
    """Размыть депозит ПШПВ-калибровкой прибора и пересыпать на каналы.

    Для канала c его энергия E(c) = poly(ecal, c); ширина канала —
    dE/dc = poly'(ecal, c). Плотность в энергии сворачивается с гауссианой
    ПШПВ(E) и интегрируется по ширине канала.
    """
    ch = np.arange(nchan, dtype=float)
    E_ch = np.array([poly(ecal, c) for c in ch])
    dEdc = np.gradient(E_ch, ch)                     # кэВ на канал
    out = np.zeros(nchan)
    for e0, n in ec:
        if n <= 0 or e0 <= 0:
            continue
        # ШИРИНА — ТОЙ ЖЕ ФУНКЦИЕЙ, ЧТО В МОСТЕ (исправлено 01.08.2026).
        # Прежде здесь стоял poly(fwhm_cal, e0) — приборный полином ПШПВ,
        # взятый от САМОЙ ЭНЕРГИИ. Он задан от корня из энергии: на 122 кэВ
        # подстановка энергии давала 96 кэВ вместо 15, а выше 500 кэВ —
        # отрицательное значение, которое подменялось единицей. В первом
        # случае линия размазывалась в колокол шириной в сотни кэВ, во
        # втором пик исчезал вовсе. Обе поломки видны глазом в программе
        # обработки (проверено оператором на 122,1 и 583,2 кэВ).
        #
        # Берётся `detector_params.fwhm_measured` — та же ширина, которой
        # размывает `bridge_mono.py`. Это не только чинит ошибку, но и
        # делает сравнение честным: измеренный программой мост и модельный
        # снимаются с ОДИНАКОВО размытого спектра, и различие принадлежит
        # только съёму площади.
        fw = dp.fwhm_measured(e0)
        if fw <= 0:
            fw = 1.0
        s = fw / 2.3548
        # вклад линии e0 в каждый канал: N(канал) = n * g(E_ch; e0,s) * dE
        m = np.abs(E_ch - e0) < 6 * s
        g = np.exp(-0.5 * ((E_ch[m] - e0) / s) ** 2) / (s * math.sqrt(2 * math.pi))
        out[m] += n * g * dEdc[m]
    return np.rint(out).astype(np.int64)


def main():
    root = paths.require_spectravibe("экспорт моно-МК в .spe: калибровки прибора")
    sys.path.insert(0, os.path.join(str(root), "scripts"))
    from gamma.io.lsrm_spe import read_lsrm_spe, write_lsrm_spe
    from gamma.spectrum import Spectrum, StoredFwhmCalibration

    # Калибровки и канальная сетка — из реального измеренного .spe.
    ref = os.environ.get("G1S_LSRM_REF_SPE")
    if not ref or not os.path.exists(ref):
        raise SystemExit(
            "Задайте G1S_LSRM_REF_SPE — путь к реальному .spe точечной 5 см\n"
            "(для калибровок прибора и канальной сетки). В репозиторий не\n"
            "входит: несёт заводской номер и оператора.")
    rsp = read_lsrm_spe(ref)
    ecal = list(rsp.energy_cal)
    fwhm_cal = list(rsp.stored_fwhm_calibration.coefficients)
    nchan = len(rsp.counts)
    tlive, treal = float(rsp.live_time), float(rsp.real_time)
    print("калибровка: %d каналов, ~%.3f кэВ/канал; ПШПВ-калибровка прибора"
          " степени %d\n" % (nchan, (poly(ecal, nchan - 1) - poly(ecal, 0))
                             / (nchan - 1), len(fwhm_cal) - 1))

    build = paths.build("Gamma-1S")
    live = os.environ.get("G1S_LSRM_MONO_DIR")
    outdir = live if (live and os.path.isdir(live)) else str(build)

    print("%9s %11s %11s %s" % ("E, кэВ", "истин.пик", "каналов", "файл .spe"))
    # Спектры моста закоммичены в results/bridge_spectra/; рабочий каталог
    # модели имеет приоритет, если свежий прогон там есть (та же конвенция,
    # что в bridge_mono.inputs()).
    spectra_dir = os.path.join(str(paths.results("Gamma-1S")), "bridge_spectra")
    for E0, fn in sorted(MONO.items()):
        src = os.path.join(str(build), fn)
        if not os.path.exists(src):
            src = os.path.join(spectra_dir, fn)
        if not os.path.exists(src):
            print("%9.1f  нет файла %s" % (E0, fn))
            continue
        ec, N = load_mc(src)
        npeak = true_peak(ec, E0)
        counts = broaden_to_channels(ec, ecal, fwhm_cal, nchan)

        spec = Spectrum(
            counts=counts, live_time=tlive, real_time=treal,
            sample_id="G4MC mono %.0f keV" % E0,
            geometry=GEOM, detector_id=DET,
            comments=("G4MC mono line %.1f keV; true full-absorption peak"
                      " = %.0f counts (bridge B = SL_fit_area / this)"
                      % (E0, npeak)),
            energy_cal=tuple(ecal),
            energy_cal_degree=len(ecal) - 1,
            energy_cal_source="stored",
            stored_fwhm_calibration=StoredFwhmCalibration(
                coefficients=tuple(fwhm_cal)),
        )
        out = os.path.join(outdir, "G4MC_mono_%d.spe" % round(E0))
        write_lsrm_spe(spec, out)
        print("%9.1f %11.0f %11d %s" % (E0, npeak, int((counts > 0).sum()), out))

    print("\nОператору: открыть каждый .spe в СпектраЛайн на параметрах"
          " слепка, снять площадь пика фитом,")
    print("разделить на «истин.пик» -> B_lsrm(E). Предсказание остатка"
          " зафиксировано в задаче 115 ДО прогона.")


if __name__ == "__main__":
    main()
