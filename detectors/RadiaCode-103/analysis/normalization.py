# -*- coding: utf-8 -*-
"""Нормировка абсолютной шкалы и проверка светосбора.

Две разные задачи, которые легко перепутать.

1. НОРМИРОВКА. Множитель к публикуемой кривой: измеренное eps_p(662) на реперной
   пробе против того, что даёт та же кривая до нормировки. Кривая и её множитель
   обязаны быть согласованы, иначе расчёт не воспроизводит репер.

2. СВЕТОСБОР. Монте-Карло считает энерговыделение, а в пик полного поглощения
   попадают только события с полностью собранным светом. Значит светосбор должен
   проявиться как РАСХОЖДЕНИЕ множителей для полного счёта и для площади пика.

Во втором пункте важно, что оба множителя берутся из ОДНОГО прогона полного
распада: он даёт и полный счёт, и площадь пика сразу. Если брать полный счёт из
прогона распада, а площадь пика из моноэнергетической сетки, к разнице примешается
статистика двух разных расчётов — и получится дефицит светосбора там, где его нет.
Именно на этом я один раз ошибся: точка сетки на 662 кэВ дала самопоглощение
1,0012, то есть заведомо шумное значение (поглощение не может увеличивать
эффективность), и «светосбор» вышел 0,932 вместо 0,95 +- 0,04.
"""
import os

import numpy as np
from scipy.optimize import curve_fit

import sys
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec
import read_rcxml
from curves import Curve, K_NORM
from fit_peak import BASE, SAMPLE, BG, peak_model, rebin_to

E0, YIELD = 661.657, 0.851
MASS_KG = 0.246
A_TRUE = 3340.0 * MASS_KG      # Бк, по калибровке приложения (сверена с АТОМТЕХ)
NUC = os.path.join("nuclides", "organic_0.49", "nuc_Cs137.csv")
VESSEL = "m500"


def measured():
    """Площадь пика и полный счёт реперной пробы, имп/с, за вычетом фона."""
    smp = read_rcxml.read(SAMPLE)[0]
    bg = read_rcxml.read(BG)[0]
    kbg = smp.live / bg.live
    bgs = rebin_to(smp, bg) * kbg
    net = smp.counts - bgs
    e, dE = smp.energy, np.gradient(smp.energy)

    m = (e > 560) & (e < 790)
    x, y = e[m], net[m] / dE[m]
    # Дисперсия ЧИСТОГО счёта: N_проба + k²·N_фон_сырой, где k —
    # множитель приведения фона по времени. Через уже
    # приведённый фон это N_проба + k·bg_приведённый. Стояло без k,
    # то есть при k < 1 погрешность занижалась, и сверка «модель
    # против аттестации» показывала согласие лучше действительного.
    # Эталон вывода — detectors/Gamma-1S/analysis/export_curves.py.
    err = np.sqrt(np.maximum(smp.counts[m] + kbg * bgs[m], 1)) / dE[m]
    p, cov = curve_fit(peak_model, x, y,
                       p0=[net[m].sum(), E0, 0.084 * 662 / 2.355, y.min(), 0.0],
                       sigma=err, absolute_sigma=True)
    area, d_area = p[0], np.sqrt(cov[0, 0])

    thr = e > 20
    tot = net[thr].sum()
    d_tot = np.sqrt(smp.counts[thr].sum() + kbg * bgs[thr].sum())
    return (area, d_area, tot, d_tot, smp.live)


def main():
    area, d_area, tot, d_tot, live = measured()
    print("реперная проба: %.0f Бк (%.0f Бк/кг), живое время %d с"
          % (A_TRUE, A_TRUE / MASS_KG, live))
    print("  площадь пика 662  %8.0f +- %.0f имп  (%.5f имп/с)"
          % (area, d_area, area / live))
    print("  полный счёт >20   %8.0f +- %.0f имп  (%.4f имп/с)"
          % (tot, d_tot, tot / live))

    # --- 1. множитель к публикуемой кривой -------------------------------
    c = Curve(VESSEL)
    rho = MASS_KG * 1000 / 498.9
    eps_raw = float(c.eps_p(E0, "organic", rho, norm=False))
    eps_meas = area / (A_TRUE * YIELD * live)
    k_curve = eps_meas / eps_raw
    print("\n1. НОРМИРОВКА КРИВОЙ (плотность пробы %.4f г/см³)" % rho)
    print("  eps_p(662) измеренное  %.4e" % eps_meas)
    print("  eps_p(662) расчётное   %.4e  (сглаженная кривая до нормировки)"
          % eps_raw)
    print("  множитель              %.4f     в curves.py стоит %.3f  %s"
          % (k_curve, K_NORM,
             "сходится" if abs(k_curve - K_NORM) < 0.005 else "ОБНОВИТЬ K_NORM"))

    # --- 2. светосбор: всё из одного прогона -----------------------------
    p = rcspec.rdir(*NUC.split(os.sep), v=VESSEL)
    if not os.path.exists(p):
        print("\nнет прогона полного распада: " + p)
        return
    meta, hist = rcspec.read_spec(p)
    n = float(meta["N_primaries"])
    pk = hist[int(E0) - 8:int(E0) + 9].sum()
    tt = hist[20:].sum()
    eps_pk, d_eps_pk = pk / n, pk / n / np.sqrt(pk)
    eps_tt, d_eps_tt = tt / n, tt / n / np.sqrt(tt)

    k_tot = (tot / live) / (A_TRUE * eps_tt)
    d_k_tot = k_tot * np.hypot(d_tot / tot, d_eps_tt / eps_tt)
    k_pk = (area / live) / (A_TRUE * eps_pk)
    d_k_pk = k_pk * np.hypot(d_area / area, d_eps_pk / eps_pk)
    k_pt = k_pk / k_tot
    d_k_pt = k_pt * np.hypot(d_k_pk / k_pk, d_k_tot / k_tot)

    print("\n2. СВЕТОСБОР — оба множителя из одного прогона (%d распадов)" % n)
    print("  %-26s %10s %10s %16s" % ("", "измерено", "расчёт", "отношение"))
    print("  %-26s %10.4f %10.4f %8.3f +- %.3f"
          % ("полный счёт, имп/с", tot / live, A_TRUE * eps_tt, k_tot, d_k_tot))
    print("  %-26s %10.5f %10.5f %8.3f +- %.3f"
          % ("площадь пика, имп/с", area / live, A_TRUE * eps_pk, k_pk, d_k_pk))
    print("  %-26s %10s %10s %8.3f +- %.3f"
          % ("пик / полное = СВЕТОСБОР", "—", "—", k_pt, d_k_pt))
    sig = abs(1.0 - k_pt) / d_k_pt
    print("\n  отличие светосбора от единицы: %.1f сигма — %s"
          % (sig, "значимо" if sig > 2 else "НЕ ЗНАЧИМО, дефицита в данных нет"))
    print("  верхняя граница дефицита по этой сверке: %.0f %%"
          % (100 * max(0.0, 1.0 - (k_pt - 2 * d_k_pt))))
    print("\n  вывод: систематика сводится к одному множителю %.3f +- %.3f,"
          % (0.5 * (k_tot + k_pk), 0.5 * abs(k_tot - k_pk) + d_k_tot))
    print("  одинаковому для пика и полного счёта; отвечает эффективному объёму")
    print("  кристалла %.2f см³ вместо номинального 1,00." % k_tot)


if __name__ == "__main__":
    main()
