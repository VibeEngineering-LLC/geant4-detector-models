# -*- coding: utf-8 -*-
"""Сверка РАСЧЁТНОГО фона внутри домика с ИЗМЕРЕННЫМ — без подгонки.

Расчёт берётся готовым из run_bg_shield.py (b_room_prior.csv): поле ЕРН по
типовым концентрациям бетона, перенос сквозь защиту, отклик прибора. Ни один
параметр здесь не подбирается — скрипт только приводит модель на реальную
канальную сетку прибора и печатает отношение по окнам. Расхождение, если оно
есть, остаётся расхождением.

ПОЧЕМУ НЕЛЬЗЯ ИНТЕРПОЛИРОВАТЬ. У RC-103 1024 канала с КВАДРАТИЧНОЙ
калибровкой, ширина канала 2,4-3,2 кэВ против 1 кэВ у модели. Точечная
np.interp берёт модельную плотность в одной точке на канал и систематически
недосчитывает модель во столько раз, во сколько реальный канал шире
модельного. Правильно — сохраняющий поток ремешок по накопленной сумме; здесь
переиспользуется fit_lines.rebin_model_to_meas(), а не пишется заново.

Запуск:  python compare_bg.py [pb] [имя_файла_измерения]
         python compare_bg.py 50 "Фон домик 23 дня.xml"
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths  # noqa: E402
import read_rcxml  # noqa: E402
import fit_lines as fl  # noqa: E402
import run_bg_shield as bg  # noqa: E402

WINDOWS = bg.WINDOWS


def read_model(pb):
    """-> (E, cps) из b_room_prior.csv, 1 кэВ на канал."""
    path = os.path.join(bg.out_dir(pb), "b_room_prior.csv")
    if not os.path.exists(path):
        raise SystemExit("нет %s — сначала python run_bg_shield.py %.0f"
                         % (path, pb))
    e, c = [], []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line[:1].isdigit():
            continue
        p = line.split(",")
        e.append(float(p[0]))
        c.append(float(p[1]))
    return np.array(e), np.array(c)


def main():
    pb = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
    name = sys.argv[2] if len(sys.argv) > 2 else "Фон домик 23 дня.xml"

    meas_path = str(paths.measured("RadiaCode-103") / name)
    if not os.path.exists(meas_path):
        raise SystemExit(
            "нет файла измерения %s.\nПоложите его в каталог измерений "
            "(переменная G4MODELS_MEASURED)." % meas_path)

    smp = read_rcxml.read(meas_path)[0]
    e_meas = smp.energy
    cps_meas = smp.counts / smp.live

    # P-007 (16.08): ПОСЛЕДНИЙ канал прибора — канал ПЕРЕПОЛНЕНИЯ, в него падает
    # всё, что выше верхней границы шкалы (у RC-103 это 2795,5 кэВ). Замер на
    # «Фон домик 23 дня»: ch1023 = 1,70e-02 cps при соседних ~1,2e-05, то есть
    # в 1400 раз больше — это не физика спектра, а сумма всего хвоста, включая
    # мюонные события с депозитом в единицы МэВ. Сравнивать его с модельным
    # интервалом 2792-2795 кэВ бессмысленно: он давал 70,8 % полосы 1500-3000 и
    # ронял её отношение до 0,258, тогда как все подокна внутри дают 0,73-0,99.
    n_over = cps_meas[-1] * smp.live
    print("канал переполнения ch%d (E=%.1f кэВ) ИСКЛЮЧЁН: %.0f отсчётов, "
          "%.4e cps" % (len(e_meas) - 1, e_meas[-1], n_over, cps_meas[-1]))
    e_meas = e_meas[:-1]
    cps_meas = cps_meas[:-1]

    e_mod, cps_mod = read_model(pb)
    # Модель на 1-кэВ сетке уже свёрнута с разрешением прибора (rcspec.fold в
    # run_bg_shield), поэтому здесь только перенос на каналы прибора.
    cps_mod_ch = fl.rebin_model_to_meas(e_mod, cps_mod, e_meas)

    print("измерение: %s" % meas_path)
    print("живое время: %.0f с (%.2f сут), каналов %d"
          % (smp.live, smp.live / 86400.0, len(e_meas)))
    print("модель:    b_room_prior.csv, Pb %.0f мм, крышка %s, БЕЗ подгонки\n"
          % (pb, "есть" if bg.g.WITH_LID else "НЕТ"))

    # P-006 (16.08): отношение БЕЗ погрешности не позволяет сказать, различаются
    # ли два окна значимо. Пуассоновскую ошибку измерения считаем здесь (счёт в
    # окне известен: cps*live), МК-ошибку модели — из sumw2, если он рядом.
    print("%-24s %13s %13s %9s %9s"
          % ("окно", "измерено, cps", "модель, cps", "модель/изм", "±изм,%"))
    for nm, lo, hi in WINDOWS:
        m = (e_meas >= lo) & (e_meas < hi)
        if not m.any():
            continue
        a = cps_meas[m].sum()
        b = cps_mod_ch[m].sum()
        # Счёт в окне = cps*live; относительная пуассоновская ошибка = 1/sqrt(N).
        n_counts = a * smp.live
        rel_meas = (100.0 / np.sqrt(n_counts)) if n_counts > 0 else float("nan")
        print("%-24s %13.5e %13.5e %9.3f %9.2f"
              % (nm, a, b, (b / a) if a else float("nan"), rel_meas))
    print("\n±изм — только пуассоновская ошибка ИЗМЕРЕНИЯ. Ошибка МОДЕЛИ (при")
    print("importance biasing число отсчётов о точности не говорит) печатается")
    print("драйвером run_bg_shield.py по Sum(w^2) — смотреть там.")

    print("\nОтношение — ОТВЕТ, а не параметр: ни одна величина здесь не")
    print("подбиралась под измерение. Расхождение означает расхождение модели")
    print("с опытом, а не то, что модель надо докрутить.")


if __name__ == "__main__":
    main()