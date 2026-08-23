# -*- coding: utf-8 -*-
"""Картинка: РАСЧЁТНЫЙ фон против ИЗМЕРЕННОГО, без подгонки.

Две панели. Сверху — спектры на реальной канальной сетке прибора: измерение,
полная модель, и раздельно её слагаемые (гамма ЕРН помещения и космические
мюоны). Снизу — отношение модель/измерение по каналам, сглаженное скользящим
окном, с горизонталью на единице.

ЧТО ЗДЕСЬ НЕЗАВИСИМО, А ЧТО НЕТ (исправлено 16.08.2026, P-006 — прежний текст
утверждал «ни один параметр не подбирался» и типовые UNSCEAR 400/40/30, что
было неверно в обеих частях).

Концентрации в бетоне (K-40 315,92, Ra-226 7,65, Th-232 9,97 Бк/кг) не типовые,
а ИЗМЕРЕННЫЕ по линиям ОТКРЫТОГО фона — и измеренные ЭТОЙ ЖЕ моделью: площадь
линии делилась на модельный отклик той же линии. Поэтому любой равномерный
множитель, общий для отклика при извлечении и при предсказании, сокращается
тождественно, и вертикальное положение кривой по линиям сойтись ОБЯЗАНО —
это не результат проверки, а следствие построения.

Независимы здесь: ФОРМА спектра (кривая эффективности между линиями), перенос
сквозь свинец, и мюоны (априорный поток PDG 0,0167 см⁻²с⁻¹ на горизонтальную
поверхность, с активностями бетона никак не связанный).

Модель приводится на канальную сетку сохраняющим поток ремешком
(fit_lines.rebin_model_to_meas): у прибора 1024 канала с квадратичной
калибровкой, ширина канала 2,4-3,2 кэВ против 1 кэВ у модели, и точечная
интерполяция систематически недосчитывала бы модель во столько же раз.

Запуск:  python plot_bg_compare.py [pb] [out.png]"""
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths  # noqa: E402
import rcspec  # noqa: E402
import read_rcxml  # noqa: E402
import fit_lines as fl  # noqa: E402

J_PDG = 0.0167          # мюон/(см² с), горизонтальная поверхность, уровень моря
R_DISK_MM = 1000.0      # насыщенный радиус диска (проверено развёрткой)
MU_DIR_NAME = "musat_box"
# P-006: словарь типовых UNSCEAR 400/40/30 здесь БЫЛ, но нигде не использовался
# (единственное вхождение в файле), а подпись на картинке ссылалась на него как
# на источник кривой. Удалён: реальные активности живут в
# run_bg_shield.PRIOR_BQKG и попадают сюда уже вшитыми в b_room_prior.csv.


def read_two_col(path, col=1):
    e, c = [], []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line[:1].isdigit():
            continue
        p = line.split(",")
        e.append(float(p[0]))
        c.append(float(p[col]))
    return np.array(e), np.array(c)


def model_gamma(pb):
    # Каталог спрашиваем у САМОГО драйвера, а не собираем копией шаблона имени:
    # в имя вошла посадка прибора (P-012/P-013), и локальная копия «pb50_nolid»
    # молча читала бы прежний вертикальный прогон, пока драйвер пишет новый.
    import run_bg_shield as bg
    path = os.path.join(bg.out_dir(pb), "b_room_prior.csv")
    if not os.path.exists(path):
        raise SystemExit("нет %s — сначала run_bg_shield.py" % path)
    # P-014: колонки файла — E,cps_total,cps_K,cps_Ra,cps_Th,cps_mu. Брать
    # колонку 1 (cps_total) НЕЛЬЗЯ: с 16.08 в неё входят и мюоны (#SHIELD-13),
    # а ниже model_muon() добавляет их ВТОРОЙ раз. Читаем гамма-компоненты
    # поимённо — тогда двойной учёт невозможен по построению.
    out = np.zeros(rcspec.NBINS)
    for col in (2, 3, 4):
        e, c = read_two_col(path, col)
        idx = e.astype(int)
        ok = (idx >= 0) & (idx < rcspec.NBINS)
        out[idx[ok]] += c[ok]
    return out


def model_muon():
    """-> (cps на канал, сырые отсчёты) — мюоны при априорном потоке PDG."""
    build = str(paths.build("RadiaCode-103"))
    d = os.path.join(build, MU_DIR_NAME)
    hist = np.zeros(rcspec.NBINS)
    n_tot = 0.0
    for name in sorted(os.listdir(d)):
        if not name.startswith("mu3_") or not name.endswith(".csv"):
            continue
        if name.endswith(".sumw2.csv"):
            continue
        meta, h = rcspec.read_spec(os.path.join(d, name))
        hist += h
        n_tot += float(meta["N_primaries"])
    if n_tot <= 0:
        raise SystemExit("нет мюонных прогонов в " + d)
    area = math.pi * (R_DISK_MM / 10.0) ** 2
    return rcspec.fold(hist, "103") * (J_PDG * area / n_tot), hist


def model_pb210():
    """-> cps на 1 Бк/кг Pb-210 в свинце домика (5 ячеек, prov: bg_budget.py).

    Активность свинца ЭТОГО домика не измерена; на графике компонент рисуется
    ВИЛКОЙ по литературным значениям (67 OPERA / 91 GeMSE Бк/кг), не линией.
    """
    cells = {"sh0_Pb_bot": 35.5, "sh0_Pb_xhi": 54.6, "sh0_Pb_xlo": 54.6,
             "sh0_Pb_yhi": 32.8, "sh0_Pb_ylo": 32.8}
    build = str(paths.build("RadiaCode-103"))
    out = np.zeros(rcspec.NBINS)
    for cell, mass in cells.items():
        meta, h = rcspec.read_spec(os.path.join(build, "pb210_%s.csv" % cell))
        out += rcspec.fold(h, "103") * mass / float(meta["N_primaries"])
    return out


PB210_LO, PB210_HI = 67.0, 91.0   # Бк/кг, литературная вилка


def smooth(y, k=9):
    if k < 3:
        return y
    ker = np.ones(k) / k
    return np.convolve(y, ker, mode="same")


def main():
    pb = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
    out = sys.argv[2] if len(sys.argv) > 2 else "bg_compare.png"
    # Диапазон по энергии — аргументами: та же картинка годится и как обзор
    # 20-3000, и как увеличенный низ, где сидят линии самого свинца.
    emin = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
    emax = float(sys.argv[4]) if len(sys.argv) > 4 else 3000.0

    meas_path = str(paths.measured("RadiaCode-103") / "Фон домик 23 дня.xml")
    smp = read_rcxml.read(meas_path)[0]
    # P-007: последний канал — переполнение, не энергетический интервал.
    # P-017/оператор 18.08 («там калибровка по Pb210 и по суммарной хри не совпадает.
    # Прими эти пики как реперные»): реперы — линия Pb-210 46,539 кэВ (ch 17,96,
    # фит гауссианой на нетто после SNIP) И суммарный XRF-комплекс свинца (Ka2 72,80 +
    # Ka1 74,97 + Kb1,3 84,9 + Kb2 87,3, ЭФФЕКТИВНАЯ энергия = взвешенный по
    # интенсивностям центроид 76,68 кэВ, ch 30,81) — прежде якорь стоял на голой
    # Ka1 74,97, хотя прибор комплекс не разрешает и видит только суммарный центроид.
    # 4 якоря (оба реперных пика + K-40 1460,82 + Tl-208 2614,51) на 3 параметра
    # квадратичной -> rms 0,68 кэВ, 1 степень свободы, проверяемо.
    CAL_BG = [0.888094236, 2.484743967, 0.000221908]
    _ch = np.arange(len(smp.counts))
    e_meas = sum(c * _ch ** i for i, c in enumerate(CAL_BG))[:-1]
    cps_meas = (smp.counts / smp.live)[:-1]
    # Пуассоновская ошибка измерения — по СЫРЫМ отсчётам канала.
    sd_meas = (np.sqrt(np.maximum(smp.counts, 0.0)) / smp.live)[:-1]

    e_mod = np.arange(rcspec.NBINS) + 0.5
    gam = model_gamma(pb)
    mu, mu_raw = model_muon()
    pb210 = model_pb210()

    gam_ch = fl.rebin_model_to_meas(e_mod, gam, e_meas)
    mu_ch = fl.rebin_model_to_meas(e_mod, mu, e_meas)
    pb_ch = fl.rebin_model_to_meas(e_mod, pb210, e_meas)
    base_ch = gam_ch + mu_ch
    tot_lo = base_ch + PB210_LO * pb_ch
    tot_hi = base_ch + PB210_HI * pb_ch
    tot_ch = 0.5 * (tot_lo + tot_hi)   # середина вилки — для нижней панели

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(11.5, 8.6), dpi=160, sharex=True,
        gridspec_kw=dict(height_ratios=[3, 1], hspace=0.06))

    ax.step(e_meas, cps_meas, where="mid", lw=1.1, color="#1b1b1b",
            label="измерено: домик, 23,0 сут")
    ax.fill_between(e_meas, np.maximum(cps_meas - sd_meas, 1e-12),
                    cps_meas + sd_meas, step="mid", color="#1b1b1b", alpha=0.18,
                    lw=0)
    ax.fill_between(e_meas, np.maximum(tot_lo, 1e-12), np.maximum(tot_hi, 1e-12),
                    step="mid", color="#c1121f", alpha=0.30, lw=0,
                    label="модель ПОЛНАЯ: ЕРН + мюоны + Pb-210 (67…91 Бк/кг)")
    ax.step(e_meas, tot_hi, where="mid", lw=1.1, color="#c1121f")
    ax.step(e_meas, gam_ch, where="mid", lw=1.0, color="#0a6ebd", alpha=0.9,
            label="модель: гамма ЕРН помещения")
    ax.step(e_meas, mu_ch, where="mid", lw=1.0, color="#2a9d8f", alpha=0.9,
            label="модель: космические мюоны")
    ax.step(e_meas, PB210_HI * pb_ch, where="mid", lw=1.0, color="#e07b00",
            alpha=0.9, label="модель: Pb-210 свинца домика (91 Бк/кг)")

    ax.set_yscale("log")
    ax.set_ylabel("скорость счёта, имп/с на канал")
    ax.set_xlim(emin, emax)
    vis = (e_meas >= emin) & (e_meas <= emax)
    lo = max(1e-8, np.min(cps_meas[vis][cps_meas[vis] > 0]) * 0.5)
    ax.set_ylim(lo, max(cps_meas[vis].max(), tot_ch[vis].max()) * 3)
    ax.grid(alpha=0.22, lw=0.4)
    ax.legend(fontsize=8.4, loc="upper right", framealpha=0.95)
    ax.set_title(
        "Фон внутри свинцового домика: расчёт против измерения, БЕЗ ПОДГОНКИ\n"
        "Pb %.0f мм, полость 150×150×385 мм, верх открыт, маринелли m200 пустая"
        % pb, fontsize=11)

    LINES = ((46.5, "Pb-210\n(свинец домика)"), (72.8, "Pb Kα2"), (75.0, "Pb Kα1"),
             (84.9, "Pb Kβ1"), (87.3, "Pb Kβ2"), (238.6, "Pb-212"),
             (295.2, "Pb-214"), (351.9, "Pb-214"), (609.3, "Bi-214"),
             (1120.3, "Bi-214"), (1460.8, "K-40"), (1764.5, "Bi-214"),
             (2614.5, "Tl-208"))
    for x, nm in LINES:
        if not (emin <= x <= emax):
            continue
        pb_own = x < 100.0        # линии самого свинца — выделяем цветом
        ax.axvline(x, color=("#c1121f" if pb_own else "#8d99ae"), lw=0.6,
                   ls=":", alpha=0.85)
        ax.annotate(nm, xy=(x, ax.get_ylim()[1]),
                    xytext=(x, ax.get_ylim()[1] * 0.30),
                    fontsize=6.6, rotation=90,
                    color=("#c1121f" if pb_own else "#495057"),
                    ha="center", va="top")

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(cps_meas > 0, tot_ch / cps_meas, np.nan)
    axr.step(e_meas, smooth(np.nan_to_num(ratio), 9), where="mid", lw=1.0,
             color="#c1121f")
    axr.axhline(1.0, color="#1b1b1b", lw=0.8)
    axr.axhline(0.5, color="#8d99ae", lw=0.6, ls="--")
    axr.set_xlim(emin, emax)
    axr.set_ylim(0, 2.0)
    axr.set_ylabel("модель / изм.")
    axr.set_xlabel("энергия, кэВ")
    axr.grid(alpha=0.22, lw=0.4)

    # P-006 (16.08): прежняя подпись была ЛОЖНОЙ — заявляла типовые UNSCEAR
    # 400/40/30 и «свободных параметров нет», тогда как файл читает
    # b_room_prior.csv, построенный на ПОДОБРАННЫХ 315,92/7,65/9,97, причём
    # подобранных по измерению открытого фона этой же моделью. Совпадение по
    # жёстким линиям поэтому не является независимым подтверждением: активности
    # определены ровно так, чтобы эти линии сошлись.
    fig.text(0.5, 0.012,
             "Концентрации в бетоне ИЗМЕРЕНЫ по линиям открытого фона той же "
             "моделью (K-40 316, Ra-226 7,7, Th-232 10,0 Бк/кг); совпадение "
             "линий K-40 1461 и Tl-208 2615 ожидаемо по построению.\n"
             "Pb-210 в свинце домика — ЛИТЕРАТУРНАЯ вилка 67 (OPERA) … 91 (GeMSE) "
             "Бк/кг, активность этого свинца не измерялась и не подбиралась. "
             "Мюоны — поток PDG 0,0167 см⁻²·с⁻¹. Подогнанных параметров нет.",
             ha="center", fontsize=7.2, color="#343a40")
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(out)
    print("записано", out)


if __name__ == "__main__":
    main()