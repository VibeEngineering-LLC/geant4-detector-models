# -*- coding: utf-8 -*-
"""Воспроизводимый расчёт якоря Am-241 для локальной поправки шкалы
AmTiCsEu (export_amticseu_data.py, _apply_am241_anchor_correction,
11.08.2026, замечание оператора №8 -- см. amticseu-remarks.md §11).

Найден внешним аудитом (11.08.2026): комментарий в export_amticseu_data.py
заявлял argmax-расчёт канала пика Am-241 (24,0-24,16), но нигде в
репозитории не было сохранённого скрипта, которым он воспроизводится --
голословное утверждение вопреки правилу citation-green. Этот файл
закрывает разрыв: тот же расчёт, что был выполнен в сессии, сохранён
здесь как проверяемый артефакт.

Метод: окно ~0,55×ПШПВ вокруг номинала 59,541 кэВ на фон-вычтенных
отсчётах образца, параболическое уточнение вершины по 3 соседним
каналам вокруг argmax. Наивное окно (фиксированные ±10 кэВ, УЖЕ или
сравнимо с ПШПВ на этой энергии) даёт систематически смещённый
результат на слабых линиях (см. проверку на 964,057 кэВ ниже -- она же
объясняет, почему для АНКЕРОВ ВЫШЕ по энергии в _apply_am241_anchor_
correction канал взят ОБРАТНЫМ пересчётом заводской квадратики, а НЕ
этим же argmax).

Запуск:
    SPECTRAVIBE_ROOT=... G4MODELS_AMTICSEU_BG_SPE=... python am241_anchor_check.py
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


def parab_argmax(e, net, e_center, half_win):
    m = (e >= e_center - half_win) & (e <= e_center + half_win)
    idx = np.where(m)[0]
    imax = idx[np.argmax(net[idx])]
    if imax <= 0 or imax >= len(net) - 1:
        return imax, e[imax]
    y0, y1, y2 = net[imax - 1], net[imax], net[imax + 1]
    denom = y0 - 2 * y1 + y2
    off = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
    off = max(-1.0, min(1.0, off))
    ch_ref = imax + off
    e_ref = np.interp(ch_ref, np.arange(len(e)), e)
    return ch_ref, e_ref


def main():
    meas, bg = erd.read_pair()
    fwhm_k, fwhm_p, _ = erd.fit_power_law_to_factory_fwhm(
        meas["fwhm_coefs"], meas["fwhm_model"])
    ed.FWHM_LAW.update({"kind": "power", "k": fwhm_k, "p": fwhm_p})

    e = meas["e_of_ch"]
    counts = meas["counts"].astype(float)
    bg_on = np.interp(e, bg["e_of_ch"], bg["counts"].astype(float),
                      left=0.0, right=0.0)
    bg_scaled = bg_on * (meas["live_s"] / bg["live_s"])
    net = counts - bg_scaled

    print("Якорь Am-241 (номинал 59,541 кэВ), окно 0,55хПШПВ:")
    hw = 0.55 * ed.fwhm_kev(59.541)
    ch, ee = parab_argmax(e, net, 59.541, hw)
    print("  окно полуширина=%.2f кэВ  канал=%.3f  E(заводская)=%.3f  "
          "дельта к номиналу=%+.3f" % (hw, ch, ee, ee - 59.541))

    print()
    print("Контроль (наивное узкое окно ±10 кэВ -- ЗАВЕДОМО НЕВЕРНЫЙ "
          "метод, приведён для сравнения, почему высокоэнергетичные "
          "анкеры в поправке взяты ОБРАТНЫМ пересчётом заводской "
          "квадратики, а не этим же argmax):")
    for E0 in (344.279, 511.0, 661.657, 964.057):
        ch10, ee10 = parab_argmax(e, net, E0, 10.0)
        chw, eew = parab_argmax(e, net, E0, 0.55 * ed.fwhm_kev(E0))
        print("  E0=%8.3f  окно=10кэВ: канал=%.3f E=%.3f d=%+.3f  |  "
              "окно=0,55хПШПВ: канал=%.3f E=%.3f d=%+.3f"
              % (E0, ch10, ee10, ee10 - E0, chw, eew, eew - E0))


if __name__ == "__main__":
    main()
