# -*- coding: utf-8 -*-
"""Выгрузка сетки run_shield_grid.py в data.js для интерактивной страницы.

Источники (не коммитятся, GRID_DIR — рабочий каталог задачи №4):
  b_room_gamma_pb<N>.csv  — K+Ra+Th сквозь защиту, cps/канал (fit_lines.py)
  b_room_mu_pb<N>.csv     — мюон сквозь защиту, cps/канал (a_mu из fit_lines.py)
  s_K40_per_bqkg_pb<N>.csv, s_Cs137_per_bqkg_pb<N>.csv — собственная
                             активность пробы, cps/канал НА 1 Бк/кг

Спектры перебиваются в 10 кэВ/канал (было 1 кэВ) — для веса страницы и
потому что 1-кэВ детализация всё равно смыта разрешением прибора 8,4%.
Пересчитывается СУММОЙ по исходным каналам внутри нового бина (сохраняет
поток) — не подвыборкой по одной точке.

Запуск: python export_data.py
"""
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths  # noqa: E402

BUILD = str(paths.build("RadiaCode-103"))
GRID_DIR = os.path.join(BUILD, "shield_grid")
OUT_JS = os.path.join(HERE, "data.js")

PB_NODES = [5, 10, 15, 20, 25, 30, 40, 50, 100]
RCAV, HZCAV = 50.0, 90.0            # ShieldGeom defaults, PbShield.hh
PB_DENSITY_G_CM3 = 11.35            # G4_Pb NIST, сверено с логом shieldrun geom
                                     # (pb=20 -> 22.3926 кг = 1972.92 см³ × 11.35)
K_BASELINE_BQKG = 250.0             # чистая ягода, постановка задачи
FWHM_662_REL = 0.084                # RC-103, R662 в rcspec.py
NBIN = 10.0                         # кэВ, экспортный бин
NCH_OUT = 320                       # 0..3200 кэВ


def rebin(e, c, step=NBIN, nch=NCH_OUT):
    """Сумма исходных (1 кэВ) каналов внутри каждого нового бина ширины step."""
    out = [0.0] * nch
    for ei, ci in zip(e, c):
        k = int(ei // step)
        if 0 <= k < nch:
            out[k] += ci
    return out


def read_csv(path):
    e, c = [], []
    if not os.path.exists(path):
        return e, c
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line[:1].isdigit():
            continue
        a, b = line.split(",")
        e.append(float(a))
        c.append(float(b))
    return e, c


def pb_mass_kg(pb):
    """Масса свинца, Cd=Cu=0 — та же формула, что PbShield.cc BuildShield
    (аналитическая разность объёмов двух цилиндров), сверено с логом
    shieldrun geom pb=20 -> 22.3926 кг."""
    r_out, hz_out = (RCAV + pb) / 10.0, (HZCAV + pb) / 10.0   # см
    r_cav, hz_cav = RCAV / 10.0, HZCAV / 10.0
    v_cm3 = 2 * math.pi * (r_out ** 2 * hz_out - r_cav ** 2 * hz_cav)
    return v_cm3 * PB_DENSITY_G_CM3 / 1000.0


def fwhm_kev(e_kev):
    return FWHM_662_REL * 661.657 * math.sqrt(max(e_kev, 1.0) / 661.657)


def git_stamp():
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.join(HERE, ".."), stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        sha = "не закоммичено"
    return sha


def main():
    nodes = []
    for pb in PB_NODES:
        e_g, c_g = read_csv(os.path.join(GRID_DIR, "b_room_gamma_pb%d.csv" % pb))
        e_m, c_m = read_csv(os.path.join(GRID_DIR, "b_room_mu_pb%d.csv" % pb))
        e_k, c_k = read_csv(os.path.join(GRID_DIR, "s_K40_per_bqkg_pb%d.csv" % pb))
        e_c, c_c = read_csv(os.path.join(GRID_DIR, "s_Cs137_per_bqkg_pb%d.csv" % pb))
        if not e_g:
            raise SystemExit("нет данных для pb=%d в %s — сначала run_shield_grid.py"
                             % (pb, GRID_DIR))
        nodes.append({
            "pb": pb,
            "mass_kg": round(pb_mass_kg(pb), 3),
            "gamma": [round(x, 8) for x in rebin(e_g, c_g)],
            "mu": [round(x, 8) for x in rebin(e_m, c_m)] if e_m else [0.0] * NCH_OUT,
            "k40_per_bqkg": [round(x, 10) for x in rebin(e_k, c_k)] if e_k else [0.0] * NCH_OUT,
            "cs137_per_bqkg": [round(x, 10) for x in rebin(e_c, c_c)] if e_c else [0.0] * NCH_OUT,
        })

    data = {
        "generated": "2026-08-12",
        "git": git_stamp(),
        "nch": NCH_OUT,
        "bin_kev": NBIN,
        "k_baseline_bqkg": K_BASELINE_BQKG,
        "eps_p_662": 6.635e-4,
        "i_cs137": 0.851,
        "mass_kg_sample": 0.1001,
        "fwhm_662_rel": FWHM_662_REL,
        "n_min": 20.0,   # эвристический порог достоверности N (mda_shield.py)
        "cs_window": [632.2, 691.2],   # 662±1.25 sigma, mda_shield.py/fit_lines.py
        "nodes": nodes,
    }

    with open(OUT_JS, "w", encoding="utf-8") as f:
        f.write("// Автосгенерировано export_data.py — не редактировать руками.\n")
        f.write("const SHIELD_DATA = ")
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"))
        f.write(";\n")

    print("записано", OUT_JS, "%.0f КБ" % (os.path.getsize(OUT_JS) / 1024.0))
    for n in nodes:
        print("  pb=%3d мм  масса=%6.2f кг  gamma_sum=%.4e  mu_sum=%.4e"
              % (n["pb"], n["mass_kg"], sum(n["gamma"]), sum(n["mu"])))


if __name__ == "__main__":
    main()
