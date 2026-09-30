# -*- coding: utf-8 -*-
"""#GS-37, #SA-10: шкала образца тория по ФОРМЕ (cal_shape.json) против центроидов таблицы пиков СпектраЛайн
(независимый путь: программа прибора, README-референсы.md) и точки горба суммирования (GS_CAL_SUM).
Печатает E_shape(канал) − E_lib для каждой точки. Запуск: GS_CAL_SHAPE=1 python gs37_shape_vs_sl.py"""
import os, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "cal"))
os.environ["GS_CAL_SHAPE"] = "1"
import gs2020_calib as C
SL_PEAKS = [(543.253, 238.632), (1297.795, 583.187), (1986.459, 911.204), (2098.076, 964.766),
            (2106.831, 968.971), (3129.695, 1460.822), (5537.873, 2614.511), (6785.0, 3187.0)]   # копия fit_gs2020_th232_m1.py:73-74 + GS_CAL_SUM
n = len(C.read("sample")["counts"])
e_shape = C.energy_axis("sample", n)
e_file = C.s.channel_to_energy(np.arange(n, dtype=float), C.read("sample")["coeffs"])
print("канал      E_lib   E_файл−E_lib  E_форма−E_lib")
for ch, lib in SL_PEAKS:
    ef = np.interp(ch, np.arange(n), e_file); es = np.interp(ch, np.arange(n), e_shape)
    print("%8.1f %9.3f %+12.2f %+13.2f" % (ch, lib, ef - lib, es - lib))
