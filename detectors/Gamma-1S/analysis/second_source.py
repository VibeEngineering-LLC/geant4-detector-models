"""Второй независимый источник Th-232 (партия 420-17031) в трёх геометриях.

Задача #107, Фаза А1 плана стабилизации. Партия 420-17031 (паспорт
860 Бк/кг ±6 %, аттестация 05.06.2017, матрица ОИСН-06 ро=0,64) измерена в
поверке 2024 во всех трёх объёмных геометриях. Это другая партия и другая
удельная активность, чем у 420-7-17 (1940 Бк/кг, ОИСН-16 ро=1,6), которой
калибрована штатная кривая ЛСРМ, — поэтому запись 420-17031 свободна от
цикличности «кривая подгонялась к этому же паспорту».

Расчёт — ТЕМ ЖЕ конвейером, что kit_recalc.py для комплекта: A_изм = R/eps,
eps на распад из прогона полного распада своей геометрии (выход, каскадное
суммирование, бленды и континуум внутри), пересчёт самопоглощения на
плотность источника через f(mu*ро*d) с подогнанной d_eff. Отбор линий тот же:
чистота по спектру испускания, порог 0,95 — у Th-232 годна одна линия 2614,5.

Чего здесь ждать. Если модель завышает эффективность на 2614,5 одинаково по
обеим партиям — систематика в модели (или в общем для партий описании
геометрии). Если 420-17031 согласуется с паспортом, а 420-7-17 нет —
подозрение переносится на конкретную засыпку/паспорт 420-7-17 (задача #83).

Спектры и фон. Три .spe поверки 2024 читаются штатным читателем SpectraVibe
(SPECTRAVIBE_ROOT), встроенного фона в .spe нет — берётся фон из kit-XML
записи Th-232 той же геометрии того же периода (тот же экран, та же
обстановка). Матрица ОИСН-06 в мю-таблицах не посчитана; как и в kit_recalc,
лёгкую среду представляет вода (ро=0,64 < 1,3). На линиях 583-2614 кэВ
сечение комптоновское, состав входит слабо.

Допущение (как и во всём разборе комплекта): цепочка Th-232 в вековом
равновесии. Для природного тория партии 2017 года к 2024 это выполняется;
распад самого Th-232 (T12 = 1,4e10 лет) пренебрежим.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402
import kit_recalc as kr  # noqa: E402


def spe_reader():
    root = paths.require_spectravibe("чтение .spe поверки 2024 (партия 420-17031)")
    sys.path.insert(0, os.path.join(str(root), "scripts"))
    from gamma.io.lsrm_spe import read_lsrm_spe
    return read_lsrm_spe, root


# (геометрия kit_recalc, путь от корня SpectraVibe, масса г, объём мл)
# Массы — из паспортов «Эталон_*_Аспект2017_.src», блок 420-17031; совпадают
# с SAMPLEMASS в заголовках самих .spe (640,0/38,4/76,8).
SPE_DIR = os.path.join("detectors", "Gamma-1S", "raw_lsrm", "Work", "BG",
                       "Gamma-1S", "Spe - поверки", "Поверка 2024")
RECORDS = [
    ("Marinelli_1L", os.path.join("Маринелли",
     "Th-232_420-17031_Маринелли_0cm.spe"), 640.0, 1000.0),
    ("Petri_60mL", os.path.join("Петри-60мл",
     "Th-232_420-17031_Петри-60мл_0cm.spe"), 38.4, 60.0),
    ("Denta_120mL", os.path.join("Дента120мл",
     "Th-232_420-17031_Дента-120мл_0cm.spe"), 76.8, 120.0),
]
ASPEC, DPCT = 860.0, 6.0     # Бк/кг, ±% (паспорт 05.06.2017)
NUC, CKEY = "Th-232", "Th232chain"
LINES = [583.187, 911.204, 2614.511]     # как в kit_recalc.VLINES


def kit_background(geom):
    """Фон из kit-XML записи Th-232 той же геометрии (тот же экран/период)."""
    kd = paths.kit_dir(geom)
    files = sorted(str(p) for p in kd.rglob("*Th232*")) if kd else []
    if not files:
        return None, None
    s, b, _cal = bm.read_checked(files[0])
    return b, os.path.basename(files[0])


if __name__ == "__main__":
    read_spe, sv_root = spe_reader()
    print("Второй источник Th-232 (420-17031, 860 Бк/кг, ОИСН-06 ро=0,64):\n"
          "тот же конвейер, что kit_recalc, годная линия — 2614,5.\n")
    print("%-13s %9s %8s %11s %9s %8s" %
          ("геометрия", "E, кэВ", "имп/с", "eps/распад", "A, Бк/кг", "A/пасп"))
    rows = []
    for geom, rel, mass, vol in RECORDS:
        p = os.path.join(str(sv_root), SPE_DIR, rel)
        if not os.path.exists(p):
            print("НЕТ ФАЙЛА: %s" % p)
            continue
        sp = read_spe(p)
        s = bm.Spectrum(sp.counts, list(sp.energy_cal),
                        float(sp.live_time), float(sp.real_time), geom)
        b, bgname = kit_background(geom)
        rho = mass / vol
        A0 = ASPEC * mass / 1000.0            # распад Th-232 пренебрежим
        R = float(s.n.sum()) / s.live
        pile = math.exp(2 * kr.TAU_SHAPE * R)
        for E in LINES:
            fw = kr.FWHM662 * math.sqrt(E / 661.657)
            base = kr.RUNBASE[(geom, CKEY)]
            frac, dirt = kr.purity(base, E, fw)
            usable = frac is not None and frac >= kr.CLEAN_FRAC
            r = bm.net_rate(s, b, E, fw, roi=1.0, side=1.0)
            if r is None or r[0] <= 0:
                print("%-13s %9.1f    нет пика" % (geom, E))
                continue
            rate = r[0] * pile
            key = min(kr.MU_O, key=lambda k: abs(k - E))
            mu_src = kr.MU_O[key] if rho > 1.3 else kr.MU_W[key]
            eps = kr.eps_per_decay(geom, CKEY, E, fw, rho, mu_src)
            if not eps:
                print("%-13s %9.1f    нет eps (прогон распада?)" % (geom, E))
                continue
            A = rate / eps
            dA = A * math.hypot(r[1] / r[0], DPCT / 100.0)
            tag = "  чистота %.2f%s" % (
                frac if frac is not None else 0.0,
                "" if usable else " — В АКТИВНОСТЬ НЕ ИДЁТ")
            print("%-13s %9.1f %8.3f %11.4e %9.1f %8.3f%s"
                  % (geom, E, rate, eps, A / (mass / 1000.0), A / A0, tag))
            rows.append((geom, E, rate, eps, A / (mass / 1000.0), ASPEC,
                         A / A0, frac, usable, A, dA, A0, bgname or "?"))

    if not rows:
        raise SystemExit("не посчитано ни одной строки")

    out = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results",
        "second_source_th232.csv"))
    csvio.write(
        out,
        ["geometry", "E_keV", "rate_cps", "eps_per_decay", "A_meas_Bq_kg",
         "A_pass_Bq_kg", "ratio", "purity", "usable", "bg_ref"],
        [(g, "%.3f" % E, "%.4f" % rate, "%.6e" % eps, "%.1f" % am,
          "%.0f" % ap, "%.4f" % rt, "%.3f" % (fr if fr is not None else 0.0),
          "%d" % (1 if us else 0), bg)
         for g, E, rate, eps, am, ap, rt, fr, us, _A, _dA, _A0, bg in rows],
        comments=[
            "Второй независимый источник Th-232: партия 420-17031"
            " (860 Бк/кг ±6%, ОИСН-06 ро=0,64), поверка 2024.",
            "A_изм = R/eps_на_распад тем же конвейером, что kit_recalc"
            " (прогон распада, f(mu*ро*d), вода как лёгкая среда).",
            "ratio = A_изм/A_пасп; usable=0 — линия плохо разделена, в"
            " активность не идёт (порог чистоты %.2f)." % kr.CLEAN_FRAC,
            "Фон — из kit-XML записи Th-232 той же геометрии (bg_ref).",
        ])
    print("\nтаблица: %s (%d строк)" % (out, len(rows)))

    # Сопоставление двух партий по годным линиям (2614,5): согласие между
    # геометриями и между партиями — прямой судья систематики модели.
    #
    # Компаратор — ТОЛЬКО Th-232 старой партии (kit_activity_volume, строки
    # nuclide=Th-232 с их d_ratio), а не сводка сосуда по всем нуклидам:
    # сравнивать надо одноимённые величины. Замечание аудитора 29.07.2026 —
    # первая версия сводила новую партию со средним сосуда и получала
    # «совпадение 1,5 %» с неверным компаратором; честная мера — разность
    # в единицах суммарной сигмы двух рядов (паспорта партий независимы).
    print("\nГодные линии (2614,5), обе партии (Th-232 против Th-232):")
    print("   %-13s %18s %18s %8s" % ("геометрия", "420-17031 A/пасп",
                                      "420-7-17 A/пасп", "разн."))
    old = {}
    kitcsv = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results",
        "kit_activity_volume.csv"))
    if os.path.exists(kitcsv):
        for line in open(kitcsv, encoding="utf-8"):
            if line.startswith("#") or line.startswith("geometry"):
                continue
            fld = line.strip().split(",")
            if len(fld) > 6 and fld[1] == "Th-232":
                old[fld[0]] = (float(fld[5]), float(fld[6]))
    vals = []
    for g, E, rate, eps, am, ap, rt, fr, us, A, dA, A0, bg in rows:
        if not us:
            continue
        vals.append((rt, dA / A0))
        if g in old:
            o, do = old[g]
            sig = abs(rt - o) / math.hypot(dA / A0, do)
            cmp = "%11.3f±%.3f %7.2fσ" % (o, do, sig)
        else:
            cmp = "%18s" % "—"
        print("   %-13s %11.3f±%.3f %s" % (g, rt, dA / A0, cmp))
    if len(vals) > 1:
        av = kr.lsrm_average([(a, d) for a, d in vals])
        c2 = (sum((av[0] - a) ** 2 / d ** 2 for a, d in vals)
              / (len(vals) - 1))
        print("   сводно по геометриям: %.3f±%.3f (%s), хи2/ню = %.2f%s"
              % (av[0], av[1], av[2], c2,
                 "" if c2 <= 2.0 else "  <- ГЕОМЕТРИИ НЕСОГЛАСОВАНЫ"))
