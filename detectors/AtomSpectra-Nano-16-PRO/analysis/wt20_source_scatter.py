# -*- coding: utf-8 -*-
"""Перенос излучения ВНУТРИ источника: сколько линия теряет и что возвращается.

`wt20_source_forward.py` считает выход по узкому пучку — exp(-mu*L) с ПОЛНЫМ
коэффициентом ослабления. Там любое взаимодействие считается потерей, и это
верно ровно для одной величины: доли квантов, вышедших БЕЗ единого
взаимодействия. Для комптоновского рассеяния такой счёт неполон: квант не
исчезает, он уходит вниз по энергии и меняет направление. Часть таких квантов
всё равно выходит из пачки, и часть из них попадает обратно в окно своей линии,
потому что при малом угле потеря энергии меньше разрешения прибора.

Здесь это считается прослеживанием: фотоэффект — поглощение, некогерентное
рассеяние — Клейн–Нишина с изменением энергии и направления, когерентное —
изменение направления без потери энергии (форм-фактор не вводится, угловое
распределение берётся томсоновским: при доле когерентного 3–8 % огрубление
углов сдвигает результат меньше статистики).

Геометрия — та же пачка: 10 цилиндров ⌀3,20 мм с шагом 4,85 мм, воздух между
ними не ослабляет. Цилиндры считаются бесконечными по длине: при длине 175 мм
и пробеге порядка сантиметра вклад торцов ниже процента.

Считаются три величины на каждую линию:
  * **выход без взаимодействий** — то же, что даёт узкопучковая формула;
  * **выход в окне линии** — вышло с энергией в пределах ±ПШПВ/2 от исходной,
    то есть попадёт в тот же пик; именно эта величина отвечает площади ППП;
  * **выход всего** — вышло с любой энергией; разность с предыдущим уходит в
    непрерывную часть спектра.

    python analysis/wt20_source_scatter.py [каталог вывода]
"""
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DET = os.path.dirname(_HERE)
_XCOM = os.path.join(_DET, "reference", "xcom_wt20_alloy.csv")

N_RODS = 10
ROD_R = 0.160          # см
PITCH = 0.485          # см
RHO = 18.925           # г/см3
E_CUT = 0.015          # МэВ, ниже — считаем поглощённым
N_HIST = 200000

# разрешение прибора, ПШПВ² = f0 + f1·E — из подгонки разложения
FWHM_F0, FWHM_F1 = 200.0, 2.20

LINES = [238.63, 300.09, 583.19, 727.33, 860.56, 911.20, 968.97, 2614.51]

RNG = np.random.default_rng(20260807)


def load_xcom():
    """E (МэВ) и линейные коэффициенты по процессам, 1/см."""
    E, coh, inc, pe = [], [], [], []
    with io.open(_XCOM, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") or ln.startswith("E_"):
                continue
            v = [float(x) for x in ln.strip().split(";")]
            E.append(v[0]); coh.append(v[1]); inc.append(v[2]); pe.append(v[3])
    E = np.array(E)
    o = np.argsort(E)
    return (E[o], np.array(coh)[o] * RHO, np.array(inc)[o] * RHO,
            np.array(pe)[o] * RHO)


EG, COH, INC, PE = load_xcom()


def mu_of(E, table):
    """Лог-лог интерполяция. Скачки на K-краях сетка XCOM содержит узлами."""
    return np.exp(np.interp(np.log(E), np.log(EG), np.log(np.maximum(table, 1e-30))))


CENTRES = (np.arange(N_RODS) - (N_RODS - 1) / 2.0) * PITCH


def metal_path_to_interaction(x, y, ux, uy, sin_t, s_need):
    """Идём по лучу через 10 цилиндров, пока не набрано s_need металла.

    Возвращает (escaped, t_hit) — вышел ли квант и параметр точки
    взаимодействия вдоль луча (в трёхмерной длине).
    """
    n = len(x)
    ur2 = np.maximum(ux ** 2 + uy ** 2, 1e-12)
    t_in = np.full((n, N_RODS), np.inf)
    t_out = np.full((n, N_RODS), np.inf)
    for j, xc in enumerate(CENTRES):
        dx = x - xc
        b = dx * ux + y * uy
        c = dx ** 2 + y ** 2 - ROD_R ** 2
        disc = b ** 2 - ur2 * c
        ok = disc > 0
        sq = np.sqrt(np.clip(disc, 0, None))
        t1 = np.where(ok, (-b - sq) / ur2, np.inf)
        t2 = np.where(ok, (-b + sq) / ur2, np.inf)
        t1 = np.maximum(t1, 0.0)                 # назад по лучу не идём
        t2 = np.where(t2 > 0, t2, np.inf)
        bad = ~ok | (t2 <= t1)
        t_in[:, j] = np.where(bad, np.inf, t1)
        t_out[:, j] = np.where(bad, np.inf, t2)
    order = np.argsort(t_in, axis=1)
    t_in = np.take_along_axis(t_in, order, axis=1)
    t_out = np.take_along_axis(t_out, order, axis=1)
    # inf - inf в np.where считался бы всё равно (обе ветви вычисляются),
    # поэтому разность берётся только там, где отрезок существует
    fin = np.isfinite(t_in) & np.isfinite(t_out)
    seg = np.zeros_like(t_in)                                 # длина в плоскости
    np.subtract(t_out, t_in, out=seg, where=fin)
    seg3 = seg / sin_t[:, None]                                # в трёхмерной длине
    cum = np.cumsum(seg3, axis=1)
    total = cum[:, -1]
    escaped = s_need >= total
    idx = np.argmax(cum >= s_need[:, None], axis=1)
    prev = np.where(idx > 0, np.take_along_axis(cum, np.maximum(idx - 1, 0)[:, None],
                                                axis=1)[:, 0], 0.0)
    t_seg_in = np.take_along_axis(t_in, idx[:, None], axis=1)[:, 0]
    t_hit = t_seg_in + (s_need - prev) * sin_t                 # обратно в параметр луча
    return escaped, t_hit


def sample_klein_nishina(E_mev, n):
    """Косинус угла рассеяния по Клейну–Нишине, отбраковкой."""
    a = E_mev / 0.510998950
    cos_t = np.empty(n)
    todo = np.arange(n)
    while todo.size:
        m = todo.size
        c = 2 * RNG.random(m) - 1.0
        eps = 1.0 / (1.0 + a[todo] * (1.0 - c))
        # dsigma/dOmega ~ eps^2 (eps + 1/eps - 1 + c^2), максимум при c=1 равен 2
        f = eps ** 2 * (eps + 1.0 / eps - 1.0 + c ** 2) / 2.0
        acc = RNG.random(m) < f
        cos_t[todo[acc]] = c[acc]
        todo = todo[~acc]
    return cos_t


def rotate(ux, uy, uz, cos_t):
    """Поворот направления на угол с косинусом cos_t и случайным азимутом."""
    sin_t = np.sqrt(np.clip(1 - cos_t ** 2, 0, None))
    phi = 2 * math.pi * RNG.random(len(cos_t))
    # ортонормированный базис вокруг (ux,uy,uz)
    ax = np.where(np.abs(uz) < 0.9, 0.0, 1.0)
    vx = uy * ax - uz * 0.0
    vy = uz * 0.0 - ux * ax
    vz = ux * 0.0 - uy * 0.0
    # проще: строим перпендикуляр через векторное произведение с ортом
    hx = np.where(np.abs(uz) < 0.9, 0.0, 1.0)
    hy = np.zeros_like(ux)
    hz = np.where(np.abs(uz) < 0.9, 1.0, 0.0)
    vx = uy * hz - uz * hy
    vy = uz * hx - ux * hz
    vz = ux * hy - uy * hx
    nv = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    vx, vy, vz = vx / nv, vy / nv, vz / nv
    wx = uy * vz - uz * vy
    wy = uz * vx - ux * vz
    wz = ux * vy - uy * vx
    cx = np.cos(phi) * sin_t
    cy = np.sin(phi) * sin_t
    nx = cos_t * ux + cx * vx + cy * wx
    ny = cos_t * uy + cx * vy + cy * wy
    nz = cos_t * uz + cx * vz + cy * wz
    nn = np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
    return nx / nn, ny / nn, nz / nn


def trace_line(E0_kev, n=N_HIST):
    E0 = E0_kev / 1000.0
    k = RNG.integers(0, N_RODS, n)
    r = ROD_R * np.sqrt(RNG.random(n))
    ph = 2 * math.pi * RNG.random(n)
    x = CENTRES[k] + r * np.cos(ph)
    y = r * np.sin(ph)
    cos_t = 2 * RNG.random(n) - 1.0
    sin_t = np.sqrt(np.clip(1 - cos_t ** 2, 1e-12, None))
    psi = 2 * math.pi * RNG.random(n)
    ux, uy, uz = sin_t * np.cos(psi), sin_t * np.sin(psi), cos_t
    E = np.full(n, E0)
    virgin = np.ones(n, bool)
    alive = np.ones(n, bool)
    out_E = np.zeros(n)
    out_uy = np.zeros(n)
    out_virgin = np.zeros(n, bool)
    escaped_any = np.zeros(n, bool)

    for _ in range(60):
        idx = np.where(alive)[0]
        if idx.size == 0:
            break
        Ei = E[idx]
        mu_c = mu_of(Ei, COH); mu_i = mu_of(Ei, INC); mu_p = mu_of(Ei, PE)
        mu_t = mu_c + mu_i + mu_p
        s = -np.log(RNG.random(idx.size)) / mu_t
        st = np.sqrt(np.clip(1 - uz[idx] ** 2, 1e-12, None))
        esc, t_hit = metal_path_to_interaction(x[idx], y[idx], ux[idx], uy[idx],
                                               st, s)
        gone = idx[esc]
        escaped_any[gone] = True
        out_E[gone] = E[gone]
        out_uy[gone] = uy[gone]
        out_virgin[gone] = virgin[gone]
        alive[gone] = False

        hit = idx[~esc]
        if hit.size == 0:
            continue
        th = t_hit[~esc]
        x[hit] += ux[hit] * th * np.sqrt(np.clip(1 - uz[hit] ** 2, 1e-12, None)) / \
            np.maximum(np.sqrt(ux[hit] ** 2 + uy[hit] ** 2), 1e-12)
        y[hit] += uy[hit] * th * np.sqrt(np.clip(1 - uz[hit] ** 2, 1e-12, None)) / \
            np.maximum(np.sqrt(ux[hit] ** 2 + uy[hit] ** 2), 1e-12)

        u = RNG.random(hit.size)
        mc = mu_of(E[hit], COH); mi = mu_of(E[hit], INC); mp = mu_of(E[hit], PE)
        mt = mc + mi + mp
        p_pe = mp / mt
        p_coh = (mp + mc) / mt
        is_pe = u < p_pe
        is_coh = (~is_pe) & (u < p_coh)
        is_inc = ~(is_pe | is_coh)

        alive[hit[is_pe]] = False                       # фотоэффект — поглощение

        j = hit[is_coh]
        if j.size:
            c = 2 * RNG.random(j.size) - 1.0            # томсоновское огрубление
            nx, ny, nz = rotate(ux[j], uy[j], uz[j], c)
            ux[j], uy[j], uz[j] = nx, ny, nz
            virgin[j] = False

        j = hit[is_inc]
        if j.size:
            c = sample_klein_nishina(E[j], j.size)
            E[j] = E[j] / (1.0 + (E[j] / 0.510998950) * (1.0 - c))
            nx, ny, nz = rotate(ux[j], uy[j], uz[j], c)
            ux[j], uy[j], uz[j] = nx, ny, nz
            virgin[j] = False
            dead = j[E[j] < E_CUT]
            alive[dead] = False

    fwhm = math.sqrt(FWHM_F0 + FWHM_F1 * E0_kev)
    inwin = escaped_any & (np.abs(out_E * 1000.0 - E0_kev) <= 0.5 * fwhm)
    up = out_uy > 0
    return dict(
        f_virgin=float(out_virgin.sum()) / n,
        f_virgin_up=float((out_virgin & up).sum()) / n,
        f_peak=float(inwin.sum()) / n,
        f_peak_up=float((inwin & up).sum()) / n,
        f_any=float(escaped_any.sum()) / n,
        f_any_up=float((escaped_any & up).sum()) / n,
        fwhm=fwhm,
    )


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_DET, "results")
    os.makedirs(outdir, exist_ok=True)
    print("Перенос в источнике, %d историй на линию" % N_HIST)
    print("сплав W 0,980000 / Th 0,017576 / O 0,002424, XCOM;"
          " пачка 10 x ⌀3,20 мм, шаг 4,85 мм")
    print()
    print(" E, кэВ  ПШПВ  без вз-вий  в окне линии  всего вышло"
          "   вклад рассеяния в окно")
    rows = []
    for e in LINES:
        r = trace_line(e)
        gain = (r["f_peak"] - r["f_virgin"]) / max(r["f_virgin"], 1e-12)
        rows.append((e, r))
        print(" %7.2f %5.1f   %8.4f     %8.4f     %8.4f        %+6.1f %%"
              % (e, r["fwhm"], r["f_virgin"], r["f_peak"], r["f_any"],
                 100 * gain))
    print()
    print(" то же, только в верхнюю полусферу (сторона прибора):")
    print(" E, кэВ  без вз-вий  в окне линии  всего вышло")
    for e, r in rows:
        print(" %7.2f   %8.4f     %8.4f     %8.4f"
              % (e, r["f_virgin_up"], r["f_peak_up"], r["f_any_up"]))

    p = os.path.join(outdir, "wt20_source_scatter.csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        f.write("# перенос внутри источника, %d историй на линию\n" % N_HIST)
        f.write("# f_без — вышло без единого взаимодействия (узкий пучок)\n")
        f.write("# f_окно — вышло с энергией в пределах ±ПШПВ/2 от линии\n")
        f.write("# f_всего — вышло с любой энергией\n")
        f.write("E_кэВ;ПШПВ_кэВ;f_без;f_окно;f_всего;"
                "f_без_вверх;f_окно_вверх;f_всего_вверх\n")
        for e, r in rows:
            f.write("%.2f;%.1f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f\n"
                    % (e, r["fwhm"], r["f_virgin"], r["f_peak"], r["f_any"],
                       r["f_virgin_up"], r["f_peak_up"], r["f_any_up"]))
    print()
    print("записано:", p)


if __name__ == "__main__":
    main()
