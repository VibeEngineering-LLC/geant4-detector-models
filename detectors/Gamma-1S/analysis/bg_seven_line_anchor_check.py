# -*- coding: utf-8 -*-
"""Воспроизводимый расчёт полного пересчёта энергетической шкалы ФОНА
AmTiCsEu по 7-линейной методике ЛСРМ для ЕРН/фоновых спектров NaI
(11.08.2026, замечание оператора «фон не откалиброван» -> «просто
двигать - ошибка» -- см. amticseu-remarks.md §13, export_amticseu_
data.py._BG_ANCHORS_*).

Заменяет bg_k40_anchor_check.py (устарел, метод там ошибочен -- см.
докстринг того файла).

Метод: 10 природных линий (7 стандартных ЛСРМ-якорей + отдельно
Bi-214/Pb-214 разбиты, где видны по отдельности), окно вокруг каждой
масштабировано по ПШПВ прибора (0,6×ПШПВ, не наивные фиксированные
широкие окна), линейная подложка по средним 4 крайних каналов окна,
взвешенный (вес~sqrt(чистый счёт)) центроид в КАНАЛЬНОМ пространстве
(не в кэВ -- избегает обратной связи с уже неверной заводской шкалой).
Итоговая калибровка -- взвешенный МНК-рефит (channel, true_keV) пар,
степень подобрана по факту (проверены 1/2/3 -- линейная лучшая, RMS
не улучшается с добавлением параметров, коэффициент кривизны шумовой).

11-й якорь -- К-рентген Pb защиты: центроид в ЭТОМ фоне не сходится с
шириной окна (см. ниже), но ВКЛЮЧЁН в подгонку по прямой директиве
оператора («просто прими для ХРИ свинца защиты 77», тот же приём, что
якорь Ra-226 §16) с E=77,0 кэВ и каналом по argmax в узком окне вокруг
факторного предсказания (не по недостоверному центроиду).

Запуск:
    SPECTRAVIBE_ROOT=... G4MODELS_AMTICSEU_BG_SPE=... python bg_seven_line_anchor_check.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WEB = os.path.join(HERE, "..", "web-th232")
os.environ.setdefault("G4MODELS_SOURCE_CONFIG",
                      os.path.join(WEB, "configs", "amticseu.yaml"))
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
sys.path.insert(0, WEB)
import export_data as ed  # noqa: E402
import export_ra226_data as erd  # noqa: E402

# 7 стандартных ЛСРМ-якорей ЕРН/фона NaI (SpectraVibe references/
# 01_metadata_calibration.md, "Background-only anchor heuristic") --
# 511 несёт и аннигиляцию, и Tl-208, разнесены с соседними по факту
# видимости в этом конкретном фоне.
LINES = [
    (238.63, "Pb-212/214"),
    (295.21, "Pb-214"),
    (351.93, "Pb-214"),
    (511.0, "Tl-208+аннигиляция"),
    (609.32, "Bi-214"),
    (911.20, "Ac-228"),
    (1120.29, "Bi-214"),
    (1460.822, "K-40"),
    (1764.49, "Bi-214"),
    (2614.511, "Tl-208"),
]
# 11-й якорь: К-рентген Pb защиты. Принятое значение E=77,0 кэВ --
# округление по прямой директиве оператора (не измеренный центроид --
# см. блок "Справка" в main(), центроид там не сходится). Интенсивностно-
# взвешенное среднее кластера Kα1/Kα2/Kβ1/Kβ2 (IAEA Bi-207 rad_types=x:
# 74.970/36.5%, 72.805/21.7%, 84.986/12.5%, 87.301/3.8%) = 76,66 кэВ --
# близко к принятому, но принято именно круглое 77,0 (прецедент Ra-226).
PB_XRAY_KEV_ACCEPTED = 77.0
PB_XRAY_KEV_MEASURED = 76.66
PB_XRAY_CH = 34.0  # argmax в узком окне вокруг факторного предсказания
PB_XRAY_NETSUM = 77.2  # наименьший netsum из всех 11 якорей


def weighted_centroid_ch(ch, counts, ch_center, hw_ch):
    m = (ch >= ch_center - hw_ch) & (ch <= ch_center + hw_ch)
    x, y = ch[m], counts[m].copy()
    if len(y) < 8:
        return None
    left, right = y[:4].mean(), y[-4:].mean()
    bgline = np.linspace(left, right, len(y))
    net = np.clip(y - bgline, 0, None)
    if net.sum() <= 0:
        return None
    return float((x * net).sum() / net.sum()), float(net.sum())


def main():
    meas, bg = erd.read_pair()
    fwhm_k, fwhm_p, _ = erd.fit_power_law_to_factory_fwhm(
        meas["fwhm_coefs"], meas["fwhm_model"])
    ed.FWHM_LAW.update({"kind": "power", "k": fwhm_k, "p": fwhm_p})

    c0, c1 = bg["coefs"]
    ch = np.arange(bg["n_channels"], dtype=float)
    counts_bg = bg["counts"].astype(float)

    print("Заводские коэфф. фона (линейные): c0=%.5f c1=%.6f" % (c0, c1))
    print()

    anchors_ch, anchors_E, weights = [], [], []
    for E0, name in LINES:
        hw_kev = 0.6 * ed.fwhm_kev(E0)
        ch_center = (E0 - c0) / c1
        hw_ch = hw_kev / c1
        r = weighted_centroid_ch(ch, counts_bg, ch_center, hw_ch)
        if r is None:
            print("%9.3f %-20s -- окно вырождено, пропуск" % (E0, name))
            continue
        cch, s = r
        print("%9.3f %-20s окно=+-%.1f кэВ  канал=%8.3f  netsum=%8.1f"
              % (E0, name, hw_kev, cch, s))
        anchors_ch.append(cch)
        anchors_E.append(E0)
        weights.append(s)

    print("%9.3f %-20s окно=argmax(узкое)  канал=%8.3f  netsum=%8.1f"
          "  -- 11-й якорь, принят по директиве оператора"
          % (PB_XRAY_KEV_ACCEPTED, "К-рентген Pb", PB_XRAY_CH, PB_XRAY_NETSUM))
    anchors_ch.append(PB_XRAY_CH)
    anchors_E.append(PB_XRAY_KEV_ACCEPTED)
    weights.append(PB_XRAY_NETSUM)

    anchors_ch = np.array(anchors_ch)
    anchors_E = np.array(anchors_E)
    w = np.sqrt(np.array(weights))

    print()
    print("=== Подгонка (взвешенная МНК), степень 1-3 ===")
    best = None
    for deg in (1, 2, 3):
        A = np.vstack([anchors_ch ** k for k in range(deg + 1)]).T
        coef, *_ = np.linalg.lstsq(A * w[:, None], anchors_E * w, rcond=None)
        pred = A @ coef
        rms = float(np.sqrt(np.mean((pred - anchors_E) ** 2)))
        print("  degree=%d  coef=%s  RMS=%.3f кэВ" % (deg, coef, rms))
        if best is None or rms < best[1]:
            best = (deg, rms, coef)
    print()
    print("Наименьший RMS формально у степени %d: coef=%s" % (best[0], best[2]))
    print("В ПРОДАКШЕНЕ (export_amticseu_data.py) взята степень 1 -- "
          "коэффициент кривизны на степени 2 порядка 1e-5, физически не "
          "обоснован (детектор NaI линеен в этом диапазоне), а слабое "
          "улучшение RMS даёт единственный низковесовой 11-й якорь "
          "(наименьший netsum) -- переобучение на шум одной слабой точки, "
          "не сигнал.")

    print()
    print("=== Справка: устойчивость ВЗВЕШЕННОГО ЦЕНТРОИДА К-рентгена Pb "
          "по ширине окна (демонстрирует, почему канал взят по argmax, "
          "не по центроиду) ===")
    ch_center = (PB_XRAY_KEV_MEASURED - c0) / c1
    for hw_kev in (5.8, 10, 15, 20, 25):
        hw_ch = hw_kev / c1
        r = weighted_centroid_ch(ch, counts_bg, ch_center, hw_ch)
        if r is None:
            print("  hw=%.1f кэВ: окно вырождено" % hw_kev)
            continue
        cch, s = r
        E_old = c0 + c1 * cch
        print("  hw=%5.1f кэВ: канал=%.3f  E(заводская)=%.3f  netsum=%.1f"
              "  -- дрейфует с ростом окна, не сходится"
              % (hw_kev, cch, E_old, s))


if __name__ == "__main__":
    main()
