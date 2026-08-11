# -*- coding: utf-8 -*-
"""Воспроизводимый расчёт якоря K-40 для поправки нуля заводской шкалы
ФОНА AmTiCsEu (export_amticseu_data.py, поправка _BG_K40_SHIFT_KEV,
11.08.2026, замечание оператора «фон не откалиброван» -- см.
amticseu-remarks.md §13).

Метод: окно 1350-1560 кэВ на СОБСТВЕННОЙ (некорректированной) шкале
фона, линейная подложка по средним первых/последних 5 каналов окна,
взвешенный центроид остатка. K-40 (1460,822 кэВ) выбран якорем как
универсальная природная линия -- при 13 ч набора фона статистика
хорошая (тысячи чистых отсчётов), в отличие от источника (1 ч).

Запуск:
    SPECTRAVIBE_ROOT=... G4MODELS_AMTICSEU_BG_SPE=... python bg_k40_anchor_check.py
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
import export_ra226_data as erd  # noqa: E402

K40_KEV = 1460.822


def main():
    meas, bg = erd.read_pair()
    e_bg = bg["e_of_ch"]
    counts_bg = bg["counts"].astype(float)

    print("Заводские коэффициенты фона (линейные):", bg["coefs"])
    print("Время набора фона, с:", bg["live_s"])

    lo, hi = 1350.0, 1560.0
    m = (e_bg >= lo) & (e_bg <= hi)
    x = e_bg[m]
    y = counts_bg[m].copy()
    left = y[:5].mean()
    right = y[-5:].mean()
    bgline = np.linspace(left, right, len(y))
    net = np.clip(y - bgline, 0, None)
    centroid = float((x * net).sum() / net.sum())
    shift = K40_KEV - centroid

    print("Окно: %.0f-%.0f кэВ" % (lo, hi))
    print("Взвешенный центроид (линейная подложка по краям окна): %.2f кэВ"
          % centroid)
    print("Сумма чистого счёта в пике: %.1f" % net.sum())
    print("K-40 (табличное): %.3f кэВ" % K40_KEV)
    print("Сдвиг (K40 - центроид): %+.2f кэВ" % shift)


if __name__ == "__main__":
    main()
