# -*- coding: utf-8 -*-
"""Сколько защита добавляет в пик полного поглощения.

ЗАЧЕМ. Два разных вопроса упираются в одно число.

1. В точечных прогонах источник разыгрывается КОНУСОМ на детектор, а результат
   делится на долю телесного угла. Кванты вне конуса не рождаются, значит путь
   «вылетел в сторону — рассеялся на свинце — вернулся в кристалл» не
   воспроизводится. Насколько это занижает ППП?
2. Опорные расчёты сторонних кодов защиту могут не моделировать вовсе. Тогда
   сравнивать наш расчёт с их расчётом можно только зная этот вклад.

КАК. Один и тот же точечный источник, полный 4π (не конус — конус вырезал бы
именно те боковые пути, ради которых всё измеряется), два прогона: со сборкой
защиты и без неё (bare). Геометрия детектора, расстояние и статистика
одинаковы, отличается только защита. Разница площадей пика — целиком её
вклад.

ДВЕ ГЕОМЕТРИИ, РАЗНЫЕ ПУТИ (замечание аудитора 28.07.2026). На 5 см источник
внутри полости защиты, крышка ЗАКРЫТА (mode=shield): путь «вбок → стенка
полости → назад в кристалл» короткий и воспроизводится хорошо. На 25 см
источник НАД защитой, крышка ОТКРЫТА (mode=open): геометрия рассеяния другая,
и вклад защиты там не обязан быть тем же числом. Меряются обе.

    python detectors/Gamma-1S/analysis/shield_role.py 5     # точечная 5 см
    python detectors/Gamma-1S/analysis/shield_role.py 25    # точечная 25 см

Прогоны делаются макросами scat_test.mac/scat_bare.mac (5 см) или
scat_test_25.mac/scat_bare_25.mac (25 см) — см. REPORT, раздел про роль
защиты.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import peakwin  # noqa: E402
import stamp  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
CASES = [(45.3, "E45.3"), (56.1, "E56.1"), (59.5, "E59.5"), (88.0, "E88.0"),
        (122.1, "E122.1"), (165.9, "E165.9"), (661.657, "E661.7"),
        (2614.511, "E2614.5")]
# Мягкие энергии добавлены 28.07.2026 (замечание оператора): ЛСРМ/EffCalcMC
# защиту не моделирует вовсе, а именно на 45–90 кэВ у модели самый большой
# дефицит против измерения. Раньше проверялось только 122 кэВ и выше.
# 165,9 — КОНТРОЛЬ ЯМЫ (замечание аудитора): по остатку от уровня это дно
# провала между избытком мягкого края и избытком жёсткого; если вклад
# защиты монотонно гаснет с энергией, здесь он должен быть уже мал.


def area(path, E, symmetric=True):
    """(N, net, d, tot). symmetric=True — подложка СРЕДНЕЙ плотностью с
    ДВУХ полок (`peakwin.area(side="both")`); False — одностороннее окно
    слева, штатное `peakwin.area`.

    ОДНОСТОРОННЕЕ ОКНО СИСТЕМАТИЧНО ВРЁТ НА ПАДАЮЩЕМ КОНТИНУУМЕ (замечание
    аудитора): слева от пика континуум ВЫШЕ, чем справа (падает с ростом E),
    и левое окно завышает подложку, занижая чистую площадь. Хуже того — для
    сравнения shield/bare эта ошибка НЕ сокращается в отношении: защита
    меняет форму континуума именно СЛЕВА (обратное рассеяние встаёт туда,
    где рождается энергия ниже E), значит смещение разное в двух прогонах.

    Прежде здесь жила СВОЯ копия правила (окно в кэВ, полка E-30..E-10 с
    вылетом иода внутри) — тот самый дефект «одно правило в N местах»,
    против которого заведён peakwin; мигрировано задачей 151.
    """
    N, hist = None, {}
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"):
            if "N_primaries" in ln:
                N = int(ln.split("=")[1])
            continue
        if not ln[:1].isdigit():
            continue
        e, c = ln.split(",")
        hist[float(e)] = int(c)
    tot = sum(hist.values())
    det = {}
    net = peakwin.area(hist, E, detail=det,
                       side="both" if symmetric else "left")
    # статистическая оценка: дисперсия gross плюс дисперсия вычтенной полки
    if symmetric:
        var_bg = (0.5 * det["n_peak"]) ** 2 \
            * (det["side"] / det["n_side"] ** 2
               + det["side_r"] / det["n_side_r"] ** 2)
    else:
        var_bg = (det["n_peak"] / det["n_side"]) ** 2 * det["side"]
    d = math.sqrt(max(det["gross"] + var_bg, 1.0))
    return N, net, d, tot


def _obs():
    """Объявление наблюдаемой: отношение чистых площадей ППП shield/bare,
    основная оценка — симметричная полка (peakwin.declare(side='both'))."""
    return dict(
        {
            "quantity":
                "отношение чистых площадей ППП с защитой и без (shield/bare);"
                " точечный источник 5 см; полный 4π",
            "area":
                "чистая площадь пика; одно правило на обе стороны отношения"
                " (common/py/peakwin.area)",
        },
        **peakwin.declare(side="both"))


def main():
    print("Вклад защиты в пик полного поглощения; точечный источник 5 см, "
          "полный 4π\n")
    print("%9s %10s %10s %9s %9s" %
          ("E, кэВ", "симметр.", "слева", "разница", "полный счёт"))
    bad = 0
    rows, used = [], []
    for E, tag in CASES:
        ps = os.path.join(BUILD, "scat_p5_full_%s.csv" % tag)
        pb = os.path.join(BUILD, "scat_p5_bare_%s.csv" % tag)
        if not (os.path.exists(ps) and os.path.exists(pb)):
            print("%9.1f  нет пары прогонов (%s)" % (E, tag))
            bad += 1
            continue
        used += [ps, pb]
        # Симметричная подложка — основная оценка (замечание аудитора: левое
        # окно на падающем континууме смещает shield и bare РАЗНО, а не
        # сокращается в отношении). Одностороннее — рядом, чтобы видеть,
        # чувствителен ли вывод к способу.
        Ns, ns_, dns, ts = area(ps, E, symmetric=True)
        Nb, nb_, dnb, tb = area(pb, E, symmetric=True)
        Ns1, ns1, _, _ = area(ps, E, symmetric=False)
        Nb1, nb1, _, _ = area(pb, E, symmetric=False)
        if Ns != Nb:
            print("   ВНИМАНИЕ: статистика разная, %d против %d" % (Ns, Nb))
        r = ns_ / nb_ if nb_ else float("nan")
        dr = r * math.hypot(dns / max(ns_, 1), dnb / max(nb_, 1))
        r1 = ns1 / nb1 if nb1 else float("nan")
        print("%9.1f %10.4f %10.4f %+9.4f %9s"
              % (E, r, r1, r - r1, "%.4f/%.4f" % (ts / Ns, tb / Nb)))
        rows.append(("%.1f" % E, "%.4f" % r, "%.4f" % dr, "%.4f" % r1,
                     "%.5f" % (ts / Ns), "%.5f" % (tb / Nb)))
    if bad:
        # Сводка пишется и при неполном наборе: мягкие пары (45-88 кэВ) не
        # прогнаны ни разу, и ранний выход навсегда оставлял бы числа только
        # в консоли (дыра задачи 148). Код возврата 1 сохраняется как сигнал.
        print("\nПрогоны: g1s.exe scat_test.mac shield и "
              "g1s.exe scat_bare.mac bare")
    print("\n«симметр.» — %s;" % peakwin.declare(side="both")["shelf"])
    print("«слева» — %s (может врать на падающем континууме)."
          % peakwin.declare()["shelf"])
    print("Отношение больше единицы — защита ДОБАВЛЯЕТ в пик обратно\n"
          "рассеянные кванты; меньше — отнимает больше, чем добавляет.")
    if rows:
        op = os.path.abspath(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "results",
            "shield_role.csv"))
        csvio.write(
            op,
            ["E_keV", "ratio_sym", "d_ratio", "ratio_left",
             "total_per_prim_shield", "total_per_prim_bare"],
            rows,
            comments=[
                "Вклад защиты в ППП: отношение чистых площадей shield/bare;"
                " точечный источник 5 см; полный 4pi.",
                "ratio_sym - основная оценка (симметричная полка);"
                " ratio_left - контроль чувствительности к способу.",
                "Больше единицы - защита добавляет в пик обратнорассеянные"
                " кванты.",
            ],
            stamp=stamp.lines(
                "detectors/Gamma-1S/analysis/shield_role.py", _obs(),
                inputs=used,
                geometry_dir=str(paths.geometry("Gamma-1S")),
                names=stamp.SRC_LISTS["Gamma-1S"],
                repo_dir=str(paths.REPO)))
        print("\nсводка: %s" % op)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
