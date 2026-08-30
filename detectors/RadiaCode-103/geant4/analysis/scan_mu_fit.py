# -*- coding: utf-8 -*-
"""Сканирование фиксированной амплитуды мюонов: сколько космики убрать,
чтобы форма ОСТАВШЕГОСЯ спектра в домике совпала с известными активностями
нуклидных цепочек. Не подгонка целиком (mu не варьируется), а диагностика
вопроса "не мюонный ли шаблон завышен" (по указанию оператора 30.08.2026).

Метод: амплитуда мюонов ФИКСИРУЕТСЯ на доле scale от известной (из открытого
фона), её вклад вычитается из измерения в домике, оставшиеся 3 колонки
(K40, Ra226_chain, Th232_chain) подгоняются свободно. Если при каком-то
scale отношения подогнанных амплитуд к известным сходятся к 1 одновременно
для всех трёх цепочек — шаблон мюонов даёт систематику именно такого масштаба.
Если не сходятся ни при каком scale — дело не в мюонах.

Запуск: python scan_mu_fit.py
"""
import os, sys, io, contextlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc
import predict_shield as ps

cnt_o, e_o, live_o, _ = ps.read_meas(ftc.MEAS_NAME, ftc.CAL_ROOM)
cps_o, var_o, hm_o, _ = ps.load_cps(ps.FMT_OPEN, ftc.MUON_CSV)
names, A_o, V_o = ps.columns_on_grid(cps_o, var_o, e_o, hm_o)
sel = (e_o >= ftc.E_LO) & (e_o < ftc.E_HI)
y_o = cnt_o[sel]
w = 1.0 / np.sqrt(np.maximum(y_o, 1.0))
with contextlib.redirect_stdout(io.StringIO()):
    amp, _sd, _p, _c, _s = ftc.fit(A_o[sel] * live_o, y_o, w, names, "", "",
                                   var_counts=V_o[sel] * live_o * live_o)
amp_map = dict(zip(names, amp))

cnt_s, e_s, live_s, _ = ps.read_meas(ps.MEAS_SHIELD, None)
cps_s, var_s, hm_s, posture_s = ps.load_cps(ps.FMT_SHIELD, ps.MUON_SHIELD_CSV)
if not hm_s:
    sys.exit("нет мюонного шаблона с домиком - сканировать нечего")
names_s, A_s, V_s = ps.columns_on_grid(cps_s, var_s, e_s, hm_s)
sel_s = (e_s >= ftc.E_LO) & (e_s < ftc.E_HI)
y_s = cnt_s[sel_s]
A_s_sel = A_s[sel_s, :]
V_s_sel = V_s[sel_s, :]
w_s = 1.0 / np.sqrt(np.maximum(y_s, 1.0))

i_mu = names_s.index("mu")
other = [n for n in names_s if n != "mu"]
oi = [i for i, n in enumerate(names_s) if n != "mu"]
mu_known = amp_map["mu"]

print(f"Постановка: {posture_s}. mu_known (из открытого фона) = {mu_known:.3f} с^-1\n")
hdr = f"{'scale':>6} {'mu_fix':>9} " + " ".join(f"{n:>13}" for n in other) + f" {'chi2/ndf':>10}"
print(hdr)
print("-" * len(hdr))
for scale in [1.5, 1.25, 1.0, 0.9, 0.8, 0.77, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]:
    mu_fix = mu_known * scale
    y_sub = y_s - A_s_sel[:, i_mu] * mu_fix * live_s
    with contextlib.redirect_stdout(io.StringIO()):
        amp_f, _sd, pred_f, _c, _s = ftc.fit(A_s_sel[:, oi] * live_s, y_sub, w_s, other, "", "",
                                             var_counts=V_s_sel[:, oi] * live_s * live_s)
    pred_total = pred_f + A_s_sel[:, i_mu] * mu_fix * live_s
    chi2ndf, _shape = ftc.metrics(pred_total, y_s, len(other))
    ratios = " ".join(f"{(amp_f[i] / amp_map[n] if amp_map[n] else float('nan')):>13.3f}"
                      for i, n in enumerate(other))
    print(f"{scale:>6.2f} {mu_fix:>9.2f} {ratios} {chi2ndf:>10.3f}")