# -*- coding: utf-8 -*-
"""Сверка расчётных кривых с независимым расчётом LSRM (из состава BecqMoni).

Зачем: до сих пор абсолютная шкала опирална одно измерение аттестованной пробы.
LSRM считает эффективность своим кодом по своей модели того же прибора в тех же
сосудах — это ВТОРОЕ независимое мнение, причём его расхождения с моим расчётом
объясняются конкретными различиями моделей, а не «поправкой из ниоткуда».

ЧТО У LSRM ИНАЧЕ (из .in-файлов, reference/lsrm):
  * кристалл — ЦИЛИНДР Ø10 x 10 мм, объём 0.785 см³ вместо куба 1.000 см³;
  * колодец осесимметричный Ø20 мм, тогда как реальная полость 36.7 x 19.8 мм
    (m500) и 37.2 x 20.7 мм (m200) — у LSRM проба подходит к кристаллу ближе;
  * корпуса прибора нет: только обёртка кристалла (1 мм отражателя + 1 мм
    оболочки + 1 мм крепления);
  * проба — ВОДА 1.0 г/см³, поэтому сверять надо с моей конфигурацией water_1.00;
  * описание обёртки у самого LSRM не согласовано между файлами: в авторской
    маринелли это полиэтилен 0.93 + TiO2 4.26, в классической — Al 2.7 + MgO 2.25.

Запуск:  python compare_lsrm.py m500
"""
import csv
import math
import os

import numpy as np

import sys
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec

HERE = os.path.dirname(os.path.abspath(__file__))
# Выгрузки LSRM лежат прямо в reference/. Стоял путь reference/lsrm — каталога
# с таким именем здесь нет вовсе (он есть у Гамма-1С), поэтому сверка тихо
# ничего не находила, а curves.py так же тихо пропускал блок lsrm: перегенерация
# curves.json убрала бы график сверки из статьи, и никто бы не заметил.
REF = os.path.abspath(os.path.join(HERE, "..", "reference"))

# сосуд -> файл выгруженной кривой LSRM
LSRM_FILE = {
    "m200": "RadiaCode - author marinelli 0.2.txt",
    "m500": "RadiaCode - author marinelli 0.5.txt",
}
CFG = "full_water_1.00"      # у LSRM проба — вода 1.0 г/см³

V_CYL = np.pi * 0.5 ** 2 * 1.0   # см³, кристалл LSRM
V_CUBE = 1.0                     # см³, кристалл по паспорту (мой)


def read_lsrm(path):
    """-> E, eps, отн. погрешность (доли)."""
    E, eps, unc = [], [], []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.replace(",", " ").split()
            if len(p) < 3:
                continue
            try:
                e, y, u = float(p[0]), float(p[1]), float(p[2])
            except ValueError:
                continue          # строка заголовка
            E.append(e)
            eps.append(y)
            unc.append(u / 100.0)
    o = np.argsort(E)
    return np.array(E)[o], np.array(eps)[o], np.array(unc)[o]


def read_mine(vessel, cfg=CFG):
    E, ep, dep = [], [], []
    with open(rcspec.rdir("efficiency.csv", v=vessel), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["config"] == cfg:
                E.append(float(r["E_keV"]))
                ep.append(float(r["eps_p"]))
                dep.append(float(r["d_eps_p"]))
    o = np.argsort(E)
    return np.array(E)[o], np.array(ep)[o], np.array(dep)[o]


def logint(x, xs, ys):
    return np.exp(np.interp(np.log(np.clip(x, xs[0], xs[-1])),
                            np.log(xs), np.log(ys)))


def measured_eps(E0=661.657):
    """eps_p из аттестованного измерения: площадь пика / (A*p*t).

    Считается из спектров, а не берётся числом, чтобы сверка опиралась на данные.
    Если спектров нет под рукой — возвращает None, сверка тогда пропускается.
    """
    try:
        import read_rcxml
        from fit_peak import BASE, SAMPLE, BG, YIELD_CS, APP_BQ_PER_KG, \
            peak_model, rebin_to
        from scipy.optimize import curve_fit
    except Exception:
        return None
    if not os.path.exists(SAMPLE) or not os.path.exists(BG):
        return None

    smp = read_rcxml.read(SAMPLE)[0]
    bg = read_rcxml.read(BG)[0]
    bgs = rebin_to(smp, bg) * (smp.live / bg.live)
    net = smp.counts - bgs
    e, dE = smp.energy, np.gradient(smp.energy)
    m = (e > 560) & (e < 790)
    x, y = e[m], net[m] / dE[m]
    err = np.sqrt(np.maximum(smp.counts[m] + bgs[m], 1)) / dE[m]
    p, _ = curve_fit(peak_model, x, y,
                     p0=[net[m].sum(), E0, 0.084 * 662 / 2.355, y.min(), 0.0],
                     sigma=err, absolute_sigma=True)
    A = APP_BQ_PER_KG * smp.weight          # Бк, аттестованная привязка
    return p[0] / (A * YIELD_CS * smp.live)


def plot(v, Em, epm, depm, El, epl, uncl):
    """Две панели: сами кривые и их отношение с полосой погрешности."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as ex:
        print("график не построен:", ex)
        return

    name = {"m200": "200 мл", "m500": "500 мл"}[v]
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))

    ax[0].errorbar(Em, epm, yerr=depm, fmt="o-", ms=4, lw=1.2, color="#1f5fa8",
                   label="этот расчёт (Geant4, куб 1 см³)")
    ax[0].plot(El, epl, "-", lw=1.4, color="#c2521a",
               label="LSRM / BecqMoni (цилиндр 0,785 см³)")
    ax[0].fill_between(El, epl * (1 - uncl), epl * (1 + uncl), color="#c2521a",
                       alpha=0.18, lw=0)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("энергия, кэВ")
    ax[0].set_ylabel(r"$\varepsilon_p$, отсчёт в пике на фотон")
    ax[0].set_title("Маринелли %s, проба вода 1,00 г/см³" % name, fontsize=10)
    ax[0].legend(fontsize=8, loc="lower left")
    style(ax[0])

    yl = np.array([logint(e, El, epl) for e in Em])
    ul = np.array([logint(e, El, uncl) for e in Em])
    r = epm / yl
    dr = r * np.hypot(depm / epm, ul)
    ax[1].axhline(1.0, color="0.4", lw=1.0, ls="--")
    ax[1].errorbar(Em, r, yerr=dr, fmt="o-", ms=4, lw=1.2, color="#1f5fa8")
    ax[1].axhline(V_CUBE / V_CYL, color="#2e7d32", lw=1.0, ls=":",
                  label="ожидание по размеру кристалла (1,27)")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("энергия, кэВ")
    ax[1].set_ylabel("этот расчёт / LSRM")
    ax[1].set_title("Отношение (полоса — обе погрешности)", fontsize=10)
    ax[1].legend(fontsize=8, loc="lower left")
    style(ax[1])

    fig.tight_layout()
    d = rcspec.rdir("figures", v=v)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "lsrm_compare.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print("график:", p)


def style(ax):
    ax.grid(True, which="both", alpha=0.25, lw=0.6)
    ax.tick_params(labelsize=9)


def main():
    v = rcspec.vessel()
    path = os.path.join(REF, LSRM_FILE[v])
    if not os.path.exists(path):
        raise SystemExit("нет файла LSRM: " + path)

    El, epl, uncl = read_lsrm(path)
    Em, epm, depm = read_mine(v)

    # LSRM даёт первые точки с погрешностью в сотни процентов — они бессмысленны
    good = uncl < 0.30
    print("сосуд %s, конфигурация %s (у LSRM проба — вода 1.0)" % (v, CFG))
    print("кривая LSRM: %s" % LSRM_FILE[v])
    print("точек у LSRM %d, из них с погрешностью < 30 %%: %d (%.0f..%.0f кэВ)"
          % (len(El), good.sum(), El[good][0], El[good][-1]))

    lo, hi = El[good][0], El[good][-1]
    sel = (Em >= lo) & (Em <= hi)

    print("\n%8s %12s %7s %12s %7s %8s %8s" %
          ("E, кэВ", "мой eps_p", "мой ±%", "LSRM eps_p", "LSRM ±%",
           "мой/LSRM", "±%"))
    rat = []
    for e, y, dy in zip(Em[sel], epm[sel], depm[sel]):
        yl = logint(e, El[good], epl[good])
        ul = logint(e, El[good], uncl[good])
        r = y / yl
        um = dy / y if y > 0 else float("nan")
        rat.append((e, r))
        print("%8.1f %12.4e %7.1f %12.4e %7.1f %8.3f %8.1f"
              % (e, y, 100 * um, yl, 100 * ul, r,
                 100 * math.hypot(um, ul)))

    rat = np.array(rat)
    band = (rat[:, 0] >= 100) & (rat[:, 0] <= 1500)
    rmid = rat[band, 1]
    print("\nотношение мой/LSRM в 100..1500 кэВ: среднее %.3f, разброс %.3f..%.3f"
          % (rmid.mean(), rmid.min(), rmid.max()))

    # ФОРМА кривой — то, что не зависит от объёма кристалла и нормировки
    r662 = logint(661.7, rat[:, 0], rat[:, 1])
    print("\nФОРМА (отношение, приведённое к 1 на 662 кэВ) — проверка хода кривой,")
    print("свободная от объёма кристалла и от абсолютной нормировки:")
    print("%8s %10s" % ("E, кэВ", "форма"))
    for e, r in rat:
        print("%8.1f %10.3f" % (e, r / r662))
    sh = rat[band, 1] / r662
    print("отклонение формы в 100..1500 кэВ: %.3f..%.3f (идеал 1.000)"
          % (sh.min(), sh.max()))

    # Сравнение с аттестованным измерением. Кривая LSRM посчитана для воды
    # 1.0 г/см³, а репер — ягоды 0.49 г/см³, поэтому LSRM надо пересчитать на
    # плотность пробы. Отношение беру из своего же расчёта: оно определяется
    # самопоглощением, то есть тем, что у обеих моделей общее.
    if v == "m500":
        Ew, epw, _ = read_mine(v, "full_water_1.00")
        Eo, epo, _ = read_mine(v, "full_organic_0.50")
        E0 = 661.657
        kden = logint(E0, Eo, epo) / logint(E0, Ew, epw)
        lsrm_water = logint(E0, El[good], epl[good])
        lsrm_berry = lsrm_water * kden
        mine_berry = logint(E0, Eo, epo)
        meas = measured_eps(E0)
        print("\nсверка на 662 кэВ с аттестованной пробой (ягоды 0.49 г/см³):")
        print("  поправка на плотность 1.00 -> 0.49 (из моего расчёта): x%.3f"
              % kden)
        print("  %-34s %.4e" % ("LSRM, вода 1.00", lsrm_water))
        print("  %-34s %.4e" % ("LSRM, пересчёт на 0.49", lsrm_berry))
        print("  %-34s %.4e" % ("мой расчёт, 0.49", mine_berry))
        if meas:
            print("  %-34s %.4e" % ("ИЗМЕРЕНИЕ (аттестовано)", meas))
            print("  завышение: LSRM в %.3f раза, мой расчёт в %.3f раза"
                  % (lsrm_berry / meas, mine_berry / meas))
            print("  то есть НЕЗАВИСИМЫЙ расчёт тоже выше измерения; общая")
            print("  причина не в объёме кристалла (у LSRM он на 21 % меньше)")

    print("\nучёт разной формы кристалла:")
    print("  LSRM: цилиндр диам.10x10 = %.4f см³; мой: куб = %.4f см³, отн. %.3f"
          % (V_CYL, V_CUBE, V_CUBE / V_CYL))
    print("  если привести LSRM к 1 см³ по объёму, отношение мой/LSRM станет"
          " %.3f" % (rmid.mean() * V_CYL / V_CUBE))
    # Множитель берётся ИЗ curves.py, а не переписывается сюда числом. Здесь
    # стояло 0,80 — значение, ОБЪЯВЛЕННОЕ ОШИБОЧНЫМ в самом curves.py
    # (множитель к точке сетки, у которой самопоглощение на 662 кэВ вышло
    # 1,0012, то есть заведомо шумное). Одно число в двух видах внутри
    # одного отчёта. Тот же дефект и то же лечение, что на Гамма-1С.
    from curves import K_NORM
    print("  нормировка модели по аттестованной пробе даёт множитель %.3f;"
          % K_NORM)
    print("  с ним отношение к LSRM: %.3f" % (K_NORM * rmid.mean()))

    plot(v, Em[sel], epm[sel], depm[sel], El[good], epl[good], uncl[good])

    out = rcspec.rdir("lsrm_compare.csv", v=v)
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E_keV", "eps_p_mine", "eps_p_lsrm", "ratio", "shape"])
        for e, r in rat:
            w.writerow(["%.1f" % e, "%.6e" % logint(e, Em, epm),
                        "%.6e" % logint(e, El[good], epl[good]),
                        "%.4f" % r, "%.4f" % (r / r662)])
    print("\nтаблица:", out)


if __name__ == "__main__":
    main()
