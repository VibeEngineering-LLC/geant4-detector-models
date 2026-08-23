# -*- coding: utf-8 -*-
"""ПОЛНЫЙ бюджет фона в домике против «Фон домик 23 дня.xml».

Слагаемые: ЕРН бетона (K/Ra/Th, измеренные активности) + мюоны (поток PDG)
+ Pb-210 в свинце самого домика (литературные активности, НЕ подгонка).
ЕРН и мюоны берутся ГОТОВЫМИ из b_room_prior.csv того каталога, который
называет сам драйвер run_bg_shield (посадка прибора зашита в имя — P-012);
Pb-210 — из прогонов pbself по пяти ячейкам свинца."""
import math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths, rcspec, read_rcxml
import fit_lines as fl
import run_bg_shield as bg

# --- измерение --------------------------------------------------------------
meas = str(paths.measured("RadiaCode-103") / "Фон домик 23 дня.xml")
smp = read_rcxml.read(meas)[0]
# P-007: последний канал — переполнение (всё выше шкалы), сравнивать его как
# энергетический интервал нельзя.
# P-018: НЕ smp.energy — штатная шкала прибора врёт (rms 6,1 кэВ на якорях).
# Та же калибровка, что в plot_bg_compare.CAL_BG: реперы Pb-210 46,539 и суммарный
# XRF-комплекс свинца (эффективный центроид 76,68) + K-40 1460,82 + Tl-208 2614,51,
# квадратичная по 4 якорям, rms 0,68 кэВ.
CAL_BG = [0.888094236, 2.484743967, 0.000221908]
_ch = np.arange(len(smp.counts))
e_meas = sum(c * _ch ** i for i, c in enumerate(CAL_BG))[:-1]
cps_meas = (smp.counts / smp.live)[:-1]

# --- модель, ЕРН + мюоны — ГОТОВЫМИ из b_room_prior.csv ---------------------
# Формат: E_keV,cps_total,cps_K,cps_Ra,cps_Th,cps_mu; мюоны уже в total
# (#SHIELD-13). Складывать их вторично из musat_box было бы двойным учётом.
g = os.path.join(bg.out_dir(50.0), "b_room_prior.csv")
if not os.path.exists(g):
    raise SystemExit("нет %s — сначала run_bg_shield.py" % g)
gamma_full = np.zeros(rcspec.NBINS)
mu_full = np.zeros(rcspec.NBINS)
for line in open(g, encoding="utf-8"):
    if line.startswith("#") or not line[:1].isdigit():
        continue
    p = line.split(",")
    i = int(float(p[0]))
    if 0 <= i < rcspec.NBINS:
        gamma_full[i] = float(p[2]) + float(p[3]) + float(p[4])
        mu_full[i] = float(p[5])

# --- Pb-210 в свинце домика -------------------------------------------------
# Прогоны pbself: 1 распад Pb-210 (цепочка до Bi-210, P-010) на историю,
# равномерно по объёму ячейки. При активности a Бк/кг ячейки массой m идёт
# a*m распадов/с; N историй = T = N/(a*m) секунд, cps = fold(hist)*a*m/N.
# Массы ячеек — из геометрии слоя, сумма 210,3 кг сверена против
# massPb_kg=210.259 в шапке прогона (0,02 %).
PB_CELLS = {"sh0_Pb_bot": 35.5, "sh0_Pb_xhi": 54.6, "sh0_Pb_xlo": 54.6,
            "sh0_Pb_yhi": 32.8, "sh0_Pb_ylo": 32.8}
pb_per_bqkg = np.zeros(rcspec.NBINS)   # cps на 1 Бк/кг (вся защита)
for cell, mass in PB_CELLS.items():
    p = os.path.join(str(paths.build("RadiaCode-103")), "pb210_%s.csv" % cell)
    meta, h = rcspec.read_spec(p)
    pb_per_bqkg += rcspec.fold(h, "103") * mass / float(meta["N_primaries"])

# Литературные активности Pb-210 в коммерческом свинце (НЕ подгонка):
# 67 Бк/кг — OPERA (низкофоновый Boliden); 91 Бк/кг — GeMSE (обычный
# коммерческий). Активность свинца ЭТОГО домика не измерена и не подбирается.
PB210_BQKG = (67.0, 91.0)

# --- на канальную сетку прибора ---------------------------------------------
e_mod = np.arange(rcspec.NBINS) + 0.5
gam_ch = fl.rebin_model_to_meas(e_mod, gamma_full, e_meas)
mu_ch = fl.rebin_model_to_meas(e_mod, mu_full, e_meas)
pb_ch = fl.rebin_model_to_meas(e_mod, pb_per_bqkg, e_meas)

WIN = (("полный 20-3000", 20, 3000),
       ("20-60", 20, 60), ("60-100", 60, 100), ("100-300", 100, 300),
       ("300-700", 300, 700), ("700-1500", 700, 1500),
       ("1500-3000", 1500, 3000),
       ("K-40 1400-1520", 1400, 1520), ("Tl-208 2560-2670", 2560, 2670),
       ("2700-2790", 2700, 2790))
print("измерение: %s, живое %.2f сут (канал переполнения исключён, P-007)"
      % (os.path.basename(meas), smp.live / 86400))
print("модель: %s" % g)
print("ЕРН: " + ", ".join("%s=%.2f Бк/кг" % (k, v)
                          for k, v in bg.PRIOR_BQKG.items())
      + " (измеренные по открытому фону)")
print("Pb-210: %s Бк/кг — литературные (OPERA/GeMSE), НЕ подгонка\n"
      % " и ".join("%.0f" % a for a in PB210_BQKG))
hdr = ("%-17s %10s %10s %8s" % ("окно, кэВ", "измерено", "ЕРН+мю", "м/и")
       + "".join(" %10s %7s" % ("+Pb %g" % a, "м/и") for a in PB210_BQKG))
print(hdr)
for nm, lo, hi in WIN:
    m = (e_meas >= lo) & (e_meas < hi)
    a = cps_meas[m].sum()
    base = gam_ch[m].sum() + mu_ch[m].sum()
    row = "%-17s %10.4e %10.4e %8.3f" % (nm, a, base,
                                         base / a if a else float("nan"))
    for act in PB210_BQKG:
        tot = base + act * pb_ch[m].sum()
        row += " %10.4e %7.3f" % (tot, tot / a if a else float("nan"))
    print(row)