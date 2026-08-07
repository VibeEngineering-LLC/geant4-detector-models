# -*- coding: utf-8 -*-
"""Активность ПО ОТДЕЛЬНЫМ ЛИНИЯМ: оконный съём измеренного и модельного пика.

Зачем при наличии шаблонной подгонки. Подгонка по всему спектру опирается на
форму континуума, а континуум в контактной геометрии зависит от того, что лежит
ВОКРУГ (стол, стена, рука) — то есть от вещей, которых в модели нет. Оконный
съём этого не касается: и в измерении, и в модели берётся ОДНО И ТО ЖЕ окно
вокруг линии, из обоих вычитается подложка по одному правилу, и отношение даёт
активность. Если два способа расходятся, расхождение — это результат, а не шум.

Величина определена только вместе с окном. Каталог линий и границы окон берутся
из веб-конструктора ROI (`reference/roi/wizard_lines_iaea.xml`,
https://vibeengineering-llc.github.io/becqmoni-roi-wizard/). Оператор
подобрал их по форме измеренного спектра: где пик хорошо разделяется — узкое
окно, где мультиплет — окно шире (в имени линии в скобках). Если окна в имени
нет, берётся ±1 ПШПВ. Подложка — среднее в двух окнах шириной 0,5 ПШПВ,
отставленных на 1,5 ПШПВ от края главного окна.

    python analysis/wt20_lines.py <спектр.xml> <каталог шаблонов> [каталог вывода]
"""
import csv
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DET = os.path.dirname(_HERE)
_ROOT = os.environ.get("SPECTRAVIBE_ROOT")
if not _ROOT:
    raise SystemExit("не задана переменная окружения SPECTRAVIBE_ROOT")
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, _HERE)

import wt20_unfold as U                                            # noqa: E402
import roi_lines as R                                              # noqa: E402


# Границы пригодности линии для оценки АКТИВНОСТИ. Ниже E_MIN_ACT модель
# заведомо занижена в разы (K-серия дочерних из внутренней конверсии в шаблонах
# Geant4 не выходит), и деление измеренного на такую модель даёт величины, не
# имеющие смысла: линия 84,37 кэВ Th-228 давала 1,7·10⁷ Бк при номинале
# 1,9·10⁴. ERR_MAX_PCT отсекает случаи, где счётная погрешность сравнима с
# самой величиной.
E_MIN_ACT = 150.0
ERR_MAX_PCT = 50.0
# Доля СВОЕГО нуклида в окне, ниже которой линия непригодна. Считается по самим
# шаблонам: в окно ±1 ПШПВ попадает всё, что туда попадает, и если больше
# четверти нетто даёт чужой нуклид, отношение «измеренное к модели своего
# нуклида» активностью уже не является. Пример: 240,99 кэВ (Ra-224) отстоит от
# 238,63 (Pb-212) на 2,35 кэВ при ПШПВ около 24 — линии не разделяются вовсе, и
# оконный съём приписывал Ra-224 в 30 раз завышенную активность.
PURITY_MIN = 0.75


def fwhm_keV(e, f0=200.0, f1=2.0):
    """ПШПВ² = f0 + f1·E — та же модель, что в wt20_unfold.broaden."""
    return math.sqrt(max(f0 + f1 * e, 1.0))


def window_area(counts_e, centres, e0, lo=None, hi=None, half_fwhm=1.0):
    """Площадь пика в окне с вычетом линейной подложки.

    Если lo/hi заданы (окно из каталога) — используется явное окно; иначе
    берётся ±half_fwhm·ПШПВ. Подложка строится по двум внешним подокнам
    шириной 0,5 ПШПВ, отставленным на 1,5 ПШПВ от границ главного окна.
    """
    fw = fwhm_keV(e0)
    if lo is None or hi is None:
        lo = e0 - half_fwhm * fw
        hi = e0 + half_fwhm * fw
    sel = (centres >= lo) & (centres <= hi)
    if sel.sum() < 3:
        return None
    bl = (centres >= lo - 0.5 * fw) & (centres < lo)
    br = (centres > hi) & (centres <= hi + 0.5 * fw)
    if bl.sum() < 2 or br.sum() < 2:
        return None
    b = 0.5 * (counts_e[bl].mean() + counts_e[br].mean())
    gross = float(counts_e[sel].sum())
    back = b * sel.sum()
    return dict(gross=gross, back=back, net=gross - back, n=int(sel.sum()),
                fwhm=fw, lo=lo, hi=hi)


def load_spectrum(src, outdir):
    """Спектр образца с вычетом фона на энергетической сетке (та же, что unfold)."""
    from gamma.io.atomspectra_xml import read_atomspectra_xml
    spec = read_atomspectra_xml(src)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    t = float(spec.real_time)
    corr = U.read_correction(os.path.join(outdir, "calibration_fitted.csv"))
    edges = np.arange(0.0, U.E_MAX + U.E_STEP, U.E_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])
    y = U.rebin_to_grid(np.asarray(spec.counts, float),
                        list(spec.energy_cal), corr.get("sample"), edges)
    if bg is not None:
        yb = U.rebin_to_grid(np.asarray(bg.counts, float),
                             list(bg.energy_cal), corr.get("background"),
                             edges) * (t / float(bg.real_time))
        y = y - yb
    return y, centres, t


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, tdir = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(src)
    os.makedirs(outdir, exist_ok=True)

    y, centres, t = load_spectrum(src, outdir)
    print("время образца %.0f с, вычтенных отсчётов в сетке %.0f"
          % (t, y.sum()))

    lines = R.parse_xml(R.DEFAULT_XML)
    print("линий из конструктора: %d" % len(lines))

    # Кэш модельных шаблонов на сетке разложения. broaden применяется здесь
    # один раз на нуклид, а окно берётся из результата.
    # Хранится и уширенный шаблон «на распад», и число разыгранных историй:
    # статистическая погрешность модели определяется числом СОБЫТИЙ в окне, а
    # не величиной, отнесённой к одному распаду. Прежняя редакция брала
    # 1/√(нетто на распад) — при нетто порядка 10⁻⁴ это давало погрешность в
    # тысячи процентов, и все линии отсеивались как непригодные.
    tmap = {}
    nprim = {}
    for key in set(r["key"] for r in lines):
        if not key:
            continue
        p = os.path.join(tdir, "%s.csv" % key)
        if not os.path.exists(p):
            continue
        head, te, tc = U.read_template(p)
        n = float(head["N_primaries"])
        tmap[key] = U.broaden(te, tc / n, centres, 200.0, 2.0)
        nprim[key] = n

    print("\n  %-6s %10s  %-11s  %9s  %6s  %9s  %10s  %8s  %-8s"
          % ("нуклид", "E, кэВ", "окно, кэВ", "нетто", "выход%",
             "мод/расп", "A, Бк", "±, %", "критерий"))
    rows = []
    for r in lines:
        e0 = r["E"]
        win = r["window"]
        lo, hi = (win if win else (None, None))
        wm_key = r["key"] if r["key"] in tmap else None
        wm = window_area(tmap[wm_key], centres, e0, lo, hi) if wm_key else None
        wy = window_area(y, centres, e0, lo, hi)
        if not wy:
            continue
        # чистота окна: доля нетто, которую даёт СОБСТВЕННЫЙ нуклид линии
        purity = float("nan")
        if wm and wm["net"] > 0:
            tot = 0.0
            for k2, sp in tmap.items():
                if k2 in ("XI", "XW", "XD1", "XD2"):
                    continue
                w2 = window_area(sp, centres, e0, lo, hi)
                if w2 and w2["net"] > 0:
                    tot += w2["net"]
            purity = wm["net"] / tot if tot > 0 else float("nan")
        wlbl = ("%.1f-%.1f" % (win[0], win[1])) if win else "±ПШПВ"
        y_pct = r["yield_pct"]
        # критерий Currie: линия видна, если нетто >= 3σ подложки
        # σ_нетто ≈ √(гросс + подложка), т.к. подложка тоже пуассонова
        sigma = math.sqrt(max(wy["gross"] + wy["back"], 1.0))
        detected = wy["net"] >= 3.0 * sigma
        crit = "видна" if detected else "ниже LD"

        a = float("nan"); err = float("nan"); mstr = "—"; astr = "—"; estr = "—"
        if wm and wm["net"] > 0 and r["key"] not in ("XI", "XW"):
            mstr = "%.3e" % wm["net"]
            if detected:
                a = wy["net"] / (wm["net"] * t)
                du = sigma / max(wy["net"], 1.0)
                # число событий Монте-Карло в окне = нетто на распад × число
                # разыгранных распадов
                n_mc = wm["net"] * nprim.get(r["key"], 1.0)
                dm = 1.0 / math.sqrt(max(n_mc, 1.0))
                err = 100.0 * math.hypot(du, dm)
                astr = "%.0f" % a
                estr = "%.1f" % err
        # Пригодность линии для оценки АКТИВНОСТИ — отдельно от того, видна ли
        # она. Три отсева, каждый по известной причине:
        #   E < E_MIN_ACT — участок, где модель заведомо занижена в разы
        #     (K-серия дочерних из внутренней конверсии в шаблонах не выходит);
        #     линия 84,37 кэВ Th-228 давала оттуда 1,7·10⁷ Бк;
        #   XI, XW — характеристический рентген, не член цепочки: поле «выход»
        #     у них означает нормированную интенсивность СЕРИИ, а не выход на
        #     распад, и делить на него нельзя;
        #   err > ERR_MAX_PCT — счётная погрешность больше самой величины.
        usable = bool(detected and math.isfinite(a) and a > 0
                      and e0 >= E_MIN_ACT
                      and r["key"] not in ("XI", "XW")
                      and math.isfinite(err) and err <= ERR_MAX_PCT
                      and math.isfinite(purity) and purity >= PURITY_MIN)
        why = ("годна" if usable else
               "ниже предела" if not detected else
               "ХРИ, не член цепочки" if r["key"] in ("XI", "XW") else
               "модель занижена ниже %.0f кэВ" % E_MIN_ACT if e0 < E_MIN_ACT else
               "погрешность выше %.0f %%" % ERR_MAX_PCT
               if not (math.isfinite(err) and err <= ERR_MAX_PCT) else
               "в окне чужие линии (своих %.0f %%)" % (100 * purity))
        print("  %-6s %10.3f  %-11s  %9.0f  %6.3g  %9s  %10s  %8s  %s"
              % (r["key"] or "-", e0, wlbl, wy["net"], y_pct, mstr,
                 astr, estr, why))
        rows.append(dict(nuc=r["key"] or "-", E=e0, name=r["name"],
                         win=wlbl, net=wy["net"], y_pct=y_pct,
                         model=wm["net"] if wm else float("nan"),
                         A=a, err_pct=err, detected=detected, purity=purity,
                         usable=usable, why=why))

    p = os.path.join(outdir, "line_activities.csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        f.write("# оконный съём по каталогу линий конструктора ROI\n")
        f.write("# https://vibeengineering-llc.github.io/becqmoni-roi-wizard/\n")
        f.write("# колонка «пригодна» — годится ли линия для оценки активности;\n")
        f.write("# отсев: ниже %.0f кэВ модель занижена в разы, ХРИ не член\n"
                % E_MIN_ACT)
        f.write("# цепочки, погрешность выше %.0f %% бессмысленна\n"
                % ERR_MAX_PCT)
        f.write("нуклид;E_кэВ;окно;выход_%;нетто_изм;чистота_окна_%;A_Бк;"
                "неопр_%;пригодна;причина\n")
        for r in rows:
            f.write("%s;%.3f;%s;%.4g;%.0f;%s;%s;%s;%s;%s\n" %
                    (r["nuc"], r["E"], r["win"], r["y_pct"], r["net"],
                     ("%.0f" % (100 * r["purity"]))
                     if math.isfinite(r["purity"]) else "",
                     ("%.4g" % r["A"]) if r["usable"] else "",
                     ("%.3g" % r["err_pct"]) if r["usable"] else "",
                     "да" if r["usable"] else "нет", r["why"]))
    ok = [r for r in rows if r["usable"]]
    print("\nгодных для оценки активности линий: %d из %d"
          % (len(ok), len(rows)))
    if ok:
        # Активность ЗВЕНА и активность РЯДА — разные величины. Tl-208
        # образуется лишь в 35,94 % распадов Bi-212, поэтому его активность
        # приводится к ряду делением на ветвление; остальные звенья в вековом
        # равновесии равны активности ряда напрямую.
        br = U.branch_to_tl208()
        chain = [(r, r["A"] / br if r["nuc"] == "Tl208" else r["A"])
                 for r in ok]
        print("  ветвление Bi-212 -> Tl-208: %.2f %% (МАГАТЭ)" % (100 * br))
        print("  %-6s %10s %12s %12s" % ("нуклид", "E, кэВ", "A звена, Бк",
                                         "A ряда, Бк"))
        for r, ac in sorted(chain, key=lambda z: z[0]["E"]):
            print("  %-6s %10.2f %12.0f %12.0f" % (r["nuc"], r["E"],
                                                   r["A"], ac))
        aa = sorted(a for _, a in chain)
        print("  разброс по ряду: %.0f … %.0f Бк (в %.2f раза)"
              % (aa[0], aa[-1], aa[-1] / aa[0]))
    print("записано:", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
