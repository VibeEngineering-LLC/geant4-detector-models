# -*- coding: utf-8 -*-
"""Пропускание фона сосудом T по фоновым линиям (их нет в Th-232): T = площадь(образец − модель М1) / площадь(k·фон).
#CAL-1: пик ищется В КАЖДОМ спектре отдельно (свой центроид и ПШПВ), окно фона на образец не переносится.
Запуск: python bg_transmission_k40.py [путь к fit_m1.json]"""
import sys, json, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
import sa10_th232 as s
FIT = sys.argv[1] if len(sys.argv) > 1 else r"C:\g4work\gs2020\run_marinelli\out_v5\fit_m1.json"
r = s.read_atomspectra_xml(s.SAMPLE_XML); b = s.read_atomspectra_xml(s.BKG_XML)
k = r["live_time"] / b["live_time"]
e = s.channel_to_energy(np.arange(len(r["counts"])), r["coeffs"])   # общая координата каналов; центроиды — свои у каждого
M = np.array(json.load(open(FIT, encoding="utf-8"))["chain"]["model"])
LINES = [("Pb-214", 351.93), ("Bi-214", 609.31), ("Bi-214", 1120.29), ("K-40", 1460.82), ("Bi-214", 1764.49), ("Bi-214", 2204.1)]
for nuc, E in LINES:
    hw = 3 * np.sqrt(max(-518.78 + 4.48609 * E, 100))
    fb, fs = s.fit_peak(e, k * b["counts"], E, hw), s.fit_peak(e, r["counts"] - M, E, hw)
    if fb is None or fs is None: print(nuc, E, "фит не сошёлся:", "фон" if fb is None else "образец"); continue
    a_b = s.net_peak_area(e, k * b["counts"], fb["mu"], fb["fwhm"]); a_s = s.net_peak_area(e, r["counts"] - M, fs["mu"], fs["fwhm"])
    T = a_s["net"] / a_b["net"]; dT = T * np.hypot(a_s["sigma"] / a_s["net"], a_b["sigma"] / a_b["net"])
    print("%-6s %8.2f  центроид фон %8.2f / образец %8.2f  k·фон %8.0f ± %4.0f  образец−М1 %8.0f ± %4.0f  T = %.3f ± %.3f"
          % (nuc, E, fb["mu"], fs["mu"], a_b["net"], a_b["sigma"], a_s["net"], a_s["sigma"], T, dT))
