"""Пересчёт активности точечных записей комплекта по ШТАТНОЙ кривой ЛСРМ.

ЗАЧЕМ. У point_recalc.py эффективность — наша Geant4-модель. Здесь та же
скорость счёта (rate) и та же паспортная активность (A0), но эффективность
берётся из измеренной кривой .efr прибора (то, чем реально считает сам
ЛСРМ). Три числа на одну линию: A/пасп по НАШЕЙ модели, A/пасп по ШТАТНОЙ
кривой, и — по построению — паспорт. Если штатная кривая тоже даёт разброс
по энергиям, это потолок точности, с которым не совестно сравнивать нашу
модель: ЛСРМ откалиброван по этому же комплекту не идеально.

ЕДИНИЦЫ. `.efr` даёт эффективность НА ИСПУЩЕННЫЙ КВАНТ (см. fetch_efr.py и
формат секции: `E=eff,dpct,нуклид,площадь,dплощадь`) — ту же величину, что
и наша `eps_mono_point`, НЕ эффективность на распад. Чтобы получить
сравнимую с `rate/ed` величину, эффективность из .efr домножается на выход
линии `yield_5cm(tag, E)` — ТОТ ЖЕ выход, что использует наша модель:
сравнение эффективностей регистрации, а не библиотек нуклидов.

ИНТЕРПОЛЯЦИЯ — та же зонная (`curvefit.local_quad`), что в
compare_effcalcmc.py и compare_point.py: густая кривая ЛСРМ (24 узла)
интерполируется в энергии наших аналитических линий.

Переиспользует rate/purity/паспорт/выход из point_recalc.py — не
дублирует их. Различается только источник эффективности и печать/свод.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

sys.path.insert(0, str(paths.tools()))
from fetch_efr import parse_efr  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402
from contam import dirty_shelves  # noqa: E402
from curvefit import local_quad  # noqa: E402
from kit_recalc import lsrm_average as kr_average  # noqa: E402
from point_recalc import (
    NUC, OBS, PASSPORT, TAU_SHAPE, USED, decay_factor, eps_decay_5cm,
    eps_mono_point, fwhm, purity, record, yield_5cm,
)  # noqa: E402


_LSRM_EV = {}


MAX_GAP = 60.0   # кэВ: дальше НИ С ОДНОЙ стороны — интерполяция ненадёжна
NEAR_NODE = 5.0  # кэВ: ближе — считается «на узле», дальний край не важен

# Наблюдаемая — та же, что у point_recalc.py (эффективность модели считается
# ЕГО функциями), но НА СУЖЕННОМ НАБОРЕ ЛИНИЙ: остаются только те, где кривая
# .efr попадает в надёжную область интерполяции. Столбец ratio_ours поэтому НЕ
# равен сводке point_recalc.py и подменять её собой не может.
#
# Сужение выписано в штамп не для порядка. Расхождение сводок 1,19 против
# 1,0116 три дня числилось за «третьей цепочкой расчёта» — а на деле сводка
# point_recalc.py в репозитории просто не была перегенерирована после правок
# отбора линий и паспортов. Единственная настоящая разница между таблицами —
# вот это сужение набора, и на 5 см она не сдвигает ничего (1,012 против
# 1,012), на 25 см даёт 1,027 против 1,040 при 4 рядах вместо 6.
OBS_LSRM = dict(
    OBS,
    quantity="A_изм/A_пасп ПАРОЙ: по нашей модели и по аттестованной кривой"
             " .efr; rate и паспорт общие",
    lineset="СУЖЕН: только линии внутри надёжной области интерполяции .efr"
            " (узел ближе %.0f кэВ либо оба анкера ближе %.0f кэВ)"
            % (NEAR_NODE, MAX_GAP),
    warning="ratio_ours НЕ сводка point_recalc.py — набор линий другой;"
            " для сводки брать kit_activity_point.csv",
)


def eps_lsrm_point(geom, E):
    """Эффективность НА РАСПАД по .efr — интерполяция * выход линии.

    ПРОВЕРКА РАЗРЫВА: у .efr узлы стоят там, где у ЛСРМ есть калибровочные
    источники, а не равномерно по энергии. Большинство наших линий (NUC)
    исторически взяты С ТЕХ ЖЕ узлов — разрыв до БЛИЖНЕГО узла там доли кэВ
    (238,63 при узле 238,632), и то, что ДАЛЬНИЙ сосед за сто кэВ, роли не
    играет: `local_quad` в такой точке фактически воспроизводит сам узел.
    Первая версия проверки требовала близости С ОБЕИХ сторон и по ошибке
    выкидывала такие «точно на узле» линии (Th-228 238,63/583,19/2614,51)
    вместе с реально ненадёжными.

    Настоящая дыра — как у Ce-139 (165,9 кэВ): ни с одной стороны узла
    ближе NEAR_NODE, интерполяция через 116 кэВ дала eps в 3,4 раза ниже
    нашей модели (активность 3,388 против паспорта при модельных 0,999).
    Линии типа Bi-207 (569,7/1063,66) — не узлы ЛСРМ вовсе, ближайший сосед
    за 13–52 кэВ с одной стороны и за 100+ с другой: не интерполяция, а
    экстраполяция от одной точки. Оба случая — вне надёжной области,
    возвращаем None, линия не идёт в сравнение.
    """
    if geom not in _LSRM_EV:
        efr_name = "Точечная-5см" if geom == "Point_5cm" else "Точечная-25см"
        path = paths.efficiency_curve(efr_name)
        if path is None:
            _LSRM_EV[geom] = None
        else:
            pts = sorted(p for s in parse_efr(paths.read_text(path))
                        for p in s["points"])
            Ec = [p[0] for p in pts]
            yc = [p[1] for p in pts]
            _LSRM_EV[geom] = (Ec, local_quad(Ec, yc))
    entry = _LSRM_EV[geom]
    if entry is None:
        return None
    Ec, ev = entry
    lo, hi = Ec[0], Ec[-1]
    if E < lo or E > hi:
        return None
    below = [e for e in Ec if e <= E]
    above = [e for e in Ec if e >= E]
    gap_lo = E - below[-1] if below else 1e9
    gap_hi = above[0] - E if above else 1e9
    near_node = min(gap_lo, gap_hi) <= NEAR_NODE
    if not near_node and (gap_lo > MAX_GAP or gap_hi > MAX_GAP):
        # не рядом с реальным узлом, и хотя бы один анкер интерполяции
        # дальше MAX_GAP — надёжной пары анкеров нет
        return None
    return ev(E)


if __name__ == "__main__":
    print("Пересчёт активности по ШТАТНОЙ кривой ЛСРМ (.efr), сравнение с "
          "нашей моделью и паспортом.\n")
    print("%-11s %-8s %9s %11s %11s %8s %8s" %
          ("геометрия", "нуклид", "E, кэВ", "A наша/пасп", "A ЛСРМ/пасп",
           "наша", "ЛСРМ"))
    rows_model, rows_lsrm = [], []
    for geom in ("Point_5cm", "Point_25cm"):
        for nuc, (tag, lines) in NUC.items():
            key = (geom, nuc)
            if key not in PASSPORT:
                continue
            p = record(geom, nuc)
            if not p:
                continue
            s, b = bm.read_checked(p)[:2]
            txt = open(p, encoding="utf-8", errors="replace").read()
            md = re.search(r"<StartTime>(\d{4}-\d{2}-\d{2})", txt)
            md = md.group(1) if md else None
            A0raw, dpct, d0 = PASSPORT[key]
            A0 = A0raw * decay_factor(nuc, d0, md)
            R = float(s.n.sum()) / s.live
            pile = math.exp(2 * TAU_SHAPE * R)
            for E in lines:
                fw = fwhm(E)
                frac, dirt = purity(tag, E, fw)
                usable = frac is not None and frac >= 0.95
                if not usable:
                    continue
                if geom == "Point_25cm":
                    bad = dirty_shelves(nuc, E, fw)
                    if bad:
                        continue
                r = bm.net_rate(s, b, E, fw, roi=1.0, side=1.0)
                if r is None or r[0] <= 0:
                    continue
                rate = r[0] * pile

                # НАША МОДЕЛЬ — тем же способом, что point_recalc.py: для
                # 5 см эффективность на распад прямо из прогона распада, для
                # 25 см — через сетку моноэнергий и выход.
                pg = yield_5cm(tag, E)
                if geom == "Point_5cm":
                    ed_model, _ = eps_decay_5cm(tag, E)
                else:
                    em25 = eps_mono_point("p25cm", E)
                    ed_model = em25 * pg if (em25 and pg) else None

                eps_l = eps_lsrm_point(geom, E)
                ed_lsrm = eps_l * pg if (eps_l and pg) else None

                if not ed_model or not ed_lsrm:
                    continue
                A_model = rate / ed_model
                A_lsrm = rate / ed_lsrm
                dA = rate * (r[1] / r[0])
                rows_model.append((geom, nuc, E, A_model / A0,
                                   (dA / ed_model) / A0, dpct))
                rows_lsrm.append((geom, nuc, E, A_lsrm / A0,
                                  (dA / ed_lsrm) / A0, dpct))
                print("%-11s %-8s %9.1f %11.3f %11.3f %8.3f %8.3f"
                      % (geom, nuc, E, A_model / A0, A_lsrm / A0,
                         A_model / A0, A_lsrm / A0))

    def summarize(rows, label):
        print("\n=== %s ===" % label)
        summary = []
        for g in ("Point_5cm", "Point_25cm"):
            v = [(x[3], x[4]) for x in rows if x[0] == g]
            if not v:
                continue
            per_nuc = []
            for nuc in sorted({x[1] for x in rows if x[0] == g}):
                sel = [x for x in rows if x[0] == g and x[1] == nuc]
                an = kr_average([(x[3], x[4]) for x in sel])
                if not an:
                    continue
                dp = sel[0][5] / 100.0
                dtot = math.hypot(an[1], dp * an[0])
                per_nuc.append((an[0], dtot))
                print("   %-11s %-8s %d лин.  %.3f ± %.3f"
                      % (g, nuc, len(sel), an[0], dtot))
            av = kr_average(per_nuc)
            if av:
                print("   %-11s ИТОГО  A/пасп = %.3f ± %.3f (%s, %d рядов)"
                      % (g, av[0], av[1], av[2], av[3]))
                summary.append((g, av[0], av[1]))
        return summary

    sm_model = summarize(rows_model, "НАША МОДЕЛЬ")
    sm_lsrm = summarize(rows_lsrm, "ШТАТНАЯ КРИВАЯ ЛСРМ")

    print("\n=== СВОДКА: наша модель vs штатная кривая ЛСРМ, обе против "
          "паспорта ===")
    print("%-11s %11s %11s %9s" % ("геометрия", "наша/пасп", "ЛСРМ/пасп",
                                   "наша/ЛСРМ"))
    dm = {g: (v, dv) for g, v, dv in sm_model}
    dl = {g: (v, dv) for g, v, dv in sm_lsrm}
    out_rows = []
    for g in ("Point_5cm", "Point_25cm"):
        if g in dm and g in dl:
            vm, dvm = dm[g]
            vl, dvl = dl[g]
            print("%-11s %7.3f±%.3f %7.3f±%.3f %9.3f"
                  % (g, vm, dvm, vl, dvl, vm / vl))
            out_rows.append((g, "%.4f" % vm, "%.4f" % dvm,
                             "%.4f" % vl, "%.4f" % dvl, "%.4f" % (vm / vl)))

    if out_rows:
        op = os.path.abspath(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "results",
            "kit_activity_point_vs_lsrm.csv"))
        csvio.write(
            op,
            ["geometry", "ratio_ours", "d_ratio_ours", "ratio_lsrm",
             "d_ratio_lsrm", "ours_over_lsrm"],
            out_rows,
            comments=[
                "Активность по нашей Geant4-модели и по штатной кривой ЛСРМ"
                " (.efr), обе против паспорта источника.",
                "Один и тот же rate и паспорт A0 для обеих колонок —"
                " различается только источник эффективности.",
                "ours_over_lsrm показывает, сколько нашей систематики"
                " остаётся ПОСЛЕ вычета того, что не сходится",
                "  даже у штатной, аттестованной кривой прибора.",
                "ratio_ours посчитан на СУЖЕННОМ наборе линий (см. obs.lineset"
                " в штампе) и НЕ равен",
                "  сводке kit_activity_point.csv. Сводка модели — там; здесь"
                " пара на общем наборе.",
            ],
            stamp=stamp.lines(
                "detectors/Gamma-1S/analysis/point_recalc_lsrm.py", OBS_LSRM,
                inputs=sorted(USED),
                geometry_dir=str(paths.geometry("Gamma-1S")),
                names=stamp.SRC_LISTS["Gamma-1S"],
                repo_dir=str(paths.REPO)))
        print("\nсводка: %s" % op)
