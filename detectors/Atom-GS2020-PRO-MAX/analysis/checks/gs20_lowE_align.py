# -*- coding: utf-8 -*-
"""Мягкая зона GS2020: положение порогового фронта (полувысота слева) и горба (максимум 40-150 кэВ) в КАНАЛАХ и кэВ по
своей шкале каждого спектра (cal_own.json). Сдвиг шкал ниже 238 кэВ даёт провал «измерение − фон». Запуск: python gs20_lowE_align.py"""
import os, sys; import numpy as np
sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cal"))
import gs2020_calib as C
for tag in ("kcl", "bgw", "sample", "bg"):
    d = C.read(tag); c = np.convolve(np.asarray(d["counts"], float), np.ones(9) / 9, mode="same"); e = C.energy_axis(tag, len(c))
    win = (e > 40) & (e < 150); im = np.flatnonzero(win)[np.argmax(c[win])]
    il = im
    while il > 0 and c[il] > 0.5 * c[im]: il -= 1
    fr = il + (0.5 * c[im] - c[il]) / (c[il + 1] - c[il])   # доля канала до полувысоты
    print("%-6s живое %8.0f с | горб: канал %5d, %6.1f кэВ | фронт 1/2: канал %7.1f, %6.1f кэВ | E(канал 200) %6.1f кэВ"
          % (tag, d["live_time"], im, e[im], fr, np.interp(fr, np.arange(len(e)), e), e[200]))
