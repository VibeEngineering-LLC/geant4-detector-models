# -*- coding: utf-8 -*-
"""Предельная оценка: может ли распределение тория по сечению объяснить 583/2614.

Измеренное отношение площадей линий 583,19 и 2614,51 кэВ равно 7,845, модель на
тех же окнах даёт 3,493 — расхождение 2,25 раза. Одно из объяснений, которое
приходится проверять: торий распределён по стержню не равномерно, а собран у
поверхности, и тогда мягкая линия теряет в самопоглощении меньше, чем считает
модель.

Проверка не требует перебора профилей, потому что у эффекта есть ПОТОЛОК.
Самопоглощение всегда сильнее для мягкой линии, поэтому отношение выходов
f(583)/f(2614) не может превысить единицу ни при каком распределении: в пределе,
когда весь торий лежит на самой поверхности, обе линии выходят одинаково и
отношение стремится к 1. Значит максимум, который распределение способно дать, —
это 1 / (f(583)/f(2614))_равномерное.

Здесь этот потолок считается прямо: тем же прослеживанием, но с розыгрышем в
тонком поверхностном слое вместо равномерного по объёму.

    python analysis/wt20_surface_limit.py [каталог вывода]
"""
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import wt20_source_scatter as S      # noqa: E402

LINES = [238.63, 583.19, 911.20, 2614.51]
XRAY = 77.30                          # Ka висмута — K-серия дочерних
DEPTHS_MM = [1.60, 0.50, 0.10, 0.02]  # 1,60 = равномерно по всему радиусу


def trace_shell(E0_kev, depth_cm, n=150000):
    """То же прослеживание, но розыгрыш в слое толщиной depth от поверхности."""
    R = S.ROD_R
    r_in = max(R - depth_cm, 0.0)
    k = S.RNG.integers(0, S.N_RODS, n)
    # равномерно по площади кольца r_in..R
    u = S.RNG.random(n)
    r = np.sqrt(r_in ** 2 + u * (R ** 2 - r_in ** 2))
    ph = 2 * math.pi * S.RNG.random(n)
    x = S.CENTRES[k] + r * np.cos(ph)
    y = r * np.sin(ph)
    cos_t = 2 * S.RNG.random(n) - 1.0
    sin_t = np.sqrt(np.clip(1 - cos_t ** 2, 1e-12, None))
    psi = 2 * math.pi * S.RNG.random(n)
    ux, uy, uz = sin_t * np.cos(psi), sin_t * np.sin(psi), cos_t
    E = np.full(n, E0_kev / 1000.0)
    alive = np.ones(n, bool)
    out_E = np.zeros(n)
    out_uy = np.zeros(n)
    escaped = np.zeros(n, bool)

    for _ in range(60):
        idx = np.where(alive)[0]
        if idx.size == 0:
            break
        Ei = E[idx]
        mu_c = S.mu_of(Ei, S.COH); mu_i = S.mu_of(Ei, S.INC); mu_p = S.mu_of(Ei, S.PE)
        mu_t = mu_c + mu_i + mu_p
        s = -np.log(S.RNG.random(idx.size)) / mu_t
        esc, t_hit = S.metal_path_to_interaction(x[idx], y[idx], ux[idx], uy[idx], s)
        gone = idx[esc]
        escaped[gone] = True
        out_E[gone] = E[gone]
        out_uy[gone] = uy[gone]
        alive[gone] = False
        hit = idx[~esc]
        if hit.size == 0:
            continue
        th = t_hit[~esc]
        # t — трёхмерная длина, смещение прямо на компоненты направления
        x[hit] += ux[hit] * th
        y[hit] += uy[hit] * th
        uu = S.RNG.random(hit.size)
        mc = S.mu_of(E[hit], S.COH); mi = S.mu_of(E[hit], S.INC); mp = S.mu_of(E[hit], S.PE)
        mt = mc + mi + mp
        is_pe = uu < mp / mt
        is_coh = (~is_pe) & (uu < (mp + mc) / mt)
        is_inc = ~(is_pe | is_coh)
        alive[hit[is_pe]] = False
        j = hit[is_coh]
        if j.size:
            c = 2 * S.RNG.random(j.size) - 1.0
            ux[j], uy[j], uz[j] = S.rotate(ux[j], uy[j], uz[j], c)
        j = hit[is_inc]
        if j.size:
            c = S.sample_klein_nishina(E[j], j.size)
            E[j] = E[j] / (1.0 + (E[j] / 0.510998950) * (1.0 - c))
            ux[j], uy[j], uz[j] = S.rotate(ux[j], uy[j], uz[j], c)
            alive[j[E[j] < S.E_CUT]] = False

    fwhm = math.sqrt(S.FWHM_F0 + S.FWHM_F1 * E0_kev)
    inwin = escaped & (np.abs(out_E * 1000.0 - E0_kev) <= 0.5 * fwhm) & (out_uy > 0)
    return float(inwin.sum()) / n


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(outdir, exist_ok=True)
    print("Выход в окне линии, вверх, при розыгрыше в поверхностном слое")
    print("(1,60 мм = весь радиус, то есть равномерно по объёму)")
    print()
    hdr = "  слой, мм " + "".join("%10.2f" % e for e in LINES) + "%10.2f" % XRAY
    print(hdr)
    rows = []
    for d_mm in DEPTHS_MM:
        d = d_mm / 10.0
        vals = [trace_shell(e, d) for e in LINES]
        vx = trace_shell(XRAY, d)
        rows.append((d_mm, vals, vx))
        print("  %8.2f " % d_mm + "".join("%10.4f" % v for v in vals)
              + "%10.4f" % vx)

    print()
    print("  Отношения, которые проверяются:")
    print("  слой, мм   f(583)/f(2614)   потолок эффекта   f(77)/f(238)")
    base = None
    for d_mm, vals, vx in rows:
        r = vals[1] / vals[3]
        if base is None:
            base = r
        print("  %8.2f       %8.4f          x %.3f         %8.4f"
              % (d_mm, r, r / base, vx / vals[0]))

    need = 2.25
    ceil = 1.0 / base
    print()
    print("  Расхождение с измерением, которое надо объяснить: x%.2f" % need)
    print("  Потолок, доступный распределению тория:           x%.2f" % ceil)
    print("  Вывод: распределение по сечению %s"
          % ("объясняет" if ceil >= need else "НЕ объясняет расхождение"))

    p = os.path.join(outdir, "wt20_surface_limit.csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        f.write("# выход в окне линии вверх при розыгрыше в поверхностном слое\n")
        f.write("слой_мм;" + ";".join("E%.2f" % e for e in LINES)
                + ";E%.2f\n" % XRAY)
        for d_mm, vals, vx in rows:
            f.write("%.2f;" % d_mm + ";".join("%.5f" % v for v in vals)
                    + ";%.5f\n" % vx)
    print()
    print("записано:", p)


if __name__ == "__main__":
    main()
