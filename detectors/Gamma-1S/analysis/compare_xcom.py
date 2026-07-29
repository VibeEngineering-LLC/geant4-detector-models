"""Сверка сечений самого Geant4 с NIST XCOM на материалах входного торца.

ЗАЧЕМ. Разбор мягкого края кривой (толщина и плотность слоёв торца) молча
опирается на то, что коэффициенты ослабления, которыми считает тулкит, верны.
Это предположение до сих пор не проверялось: расхождение модели с измерением
искали в геометрии, ни разу не спросив, не смещены ли сами сечения. Здесь
вопрос закрывается прямо — mu/ro из G4EmCalculator против первоисточника.

ЧТО С ЧЕМ СРАВНИВАЕТСЯ. Основной расчёт ослабления в mucalc.cc идёт БЕЗ
когерентного рассеяния (compt+phot+conv, см. шапку mucalc.cc), поэтому колонка
no_coh сверяется с XCOM «without coherent», а with_coh — с «with coherent».
Смешивать нельзя: на 59,5 кэВ эти две величины расходятся на 12 %, и сверка
не с той колонкой дала бы ложную тревогу такого же порядка.

ИСТОЧНИК ЭТАЛОНА. NIST XCOM, https://physics.nist.gov/PhysRefData/Xcom/html/
xcom1.html — расчёт на точно запрошенных энергиях (база интерполирует лог-лог
по своей внутренней сетке, то есть это расчёт, а не считывание фиксированного
узла; для сверки двух расчётов на одной энергии это ровно то, что нужно).
Составы: MgO и NaI — по формуле; резина — состав ICRU-37 (NIST STAR, материал
243: H 0,118371 / C 0,881629), готовой таблицы mu/ro для резины у NIST нет.
Независимый контроль для Al: таблица Хаббелла–Зельцера (NISTIR 5632, Table 3,
z13) даёт 2,778E-01 при 60 кэВ — сходится с колонкой WITH coherent, что и
подтверждает, какая это величина.

ВХОД: mu_xcom_check.csv из каталога расчётов (пишется mucalc.exe).
"""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402

# NIST XCOM, см²/г: (без когерентного, с когерентным).
# Энергии — реперные линии комплекта: Am-241, Cd-109, Co-57, Cs-137.
XCOM = {
    ("MgO", 59.5): (2.075e-1, 2.326e-1),
    ("MgO", 88.0): (1.621e-1, 1.743e-1),
    ("MgO", 122.1): (1.430e-1, 1.496e-1),
    ("MgO", 661.7): (7.655e-2, 7.679e-2),
    ("Al", 59.5): (2.466e-1, 2.810e-1),
    ("Al", 88.0): (1.697e-1, 1.865e-1),
    ("Al", 122.1): (1.429e-1, 1.520e-1),
    ("Al", 661.7): (7.435e-2, 7.467e-2),
    ("NaI", 59.5): (6.370e0, 6.595e0),
    ("NaI", 88.0): (2.216e0, 2.331e0),
    ("NaI", 122.1): (9.387e-1, 1.004e0),
    ("NaI", 661.7): (7.397e-2, 7.664e-2),
    ("rubber", 59.5): (1.848e-1, 1.936e-1),
    ("rubber", 88.0): (1.699e-1, 1.741e-1),
    ("rubber", 122.1): (1.571e-1, 1.593e-1),
    ("rubber", 661.7): (8.609e-2, 8.617e-2),
}

# Порог, выше которого расхождение перестаёт быть «согласием»: сечения такого
# уровня зрелости обязаны сходиться заметно лучше процента, и всё, что выше,
# требует разбора, а не списания на округление.
TOL_PCT = 1.0


def main():
    src = paths.build("Gamma-1S") / "mu_xcom_check.csv"
    if not os.path.exists(src):
        raise SystemExit(
            "Нет %s.\nОн пишется mucalc.exe — соберите цель mucalc и запустите"
            " её в каталоге расчётов." % src)

    with open(src, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    print("Сечения Geant4 против NIST XCOM, mu/ro в см²/г.\n")
    print("%-8s %8s %11s %11s %8s   %11s %11s %8s"
          % ("материал", "E, кэВ", "G4 без ког", "XCOM", "откл, %",
             "G4 с ког", "XCOM", "откл, %"))

    out, worst, worst_lbl = [], 0.0, ""
    for r in rows:
        key = (r["material"], round(float(r["E_keV"]), 1))
        if key not in XCOM:
            continue
        g_no, g_wi = float(r["no_coh"]), float(r["with_coh"])
        x_no, x_wi = XCOM[key]
        d_no = 100.0 * (g_no - x_no) / x_no
        d_wi = 100.0 * (g_wi - x_wi) / x_wi
        print("%-8s %8.1f %11.4e %11.4e %+8.2f   %11.4e %11.4e %+8.2f"
              % (key[0], key[1], g_no, x_no, d_no, g_wi, x_wi, d_wi))
        out.append((key[0], "%.1f" % key[1], "%.5e" % g_no, "%.5e" % x_no,
                    "%+.3f" % d_no, "%.5e" % g_wi, "%.5e" % x_wi,
                    "%+.3f" % d_wi))
        for d, lab in ((d_no, "без ког."), (d_wi, "с ког.")):
            if abs(d) > abs(worst):
                worst, worst_lbl = d, "%s %.1f кэВ, %s" % (key[0], key[1], lab)

    if not out:
        raise SystemExit("В mu_xcom_check.csv нет ни одной сверяемой строки.")

    print("\nмаксимальное отклонение: %+.2f %% (%s)" % (worst, worst_lbl))
    if abs(worst) <= TOL_PCT:
        print("Сечения тулкита согласуются с XCOM: расхождение модели с"
              " измерением\nобъяснить ими не удаётся, искать надо в другом"
              " месте.")
    else:
        print("ВНИМАНИЕ: расхождение выше %.1f %% — разобрать до того, как"
              " списывать\nнесогласие модели на геометрию." % TOL_PCT)

    op = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results",
        "mu_geant4_vs_xcom.csv"))
    csvio.write(
        op,
        ["material", "E_keV", "g4_no_coh", "xcom_no_coh", "dev_no_coh_pct",
         "g4_with_coh", "xcom_with_coh", "dev_with_coh_pct"],
        out,
        comments=[
            "Массовый коэффициент ослабления, см²/г: Geant4 11.2.1"
            " (EmStandardPhysics_option4) против NIST XCOM.",
            "no_coh = compt+phot+conv, сверяется с XCOM without coherent;"
            " with_coh добавляет Rayleigh.",
            "Эталон: https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html",
            "Резина — состав ICRU-37 (NIST STAR, материал 243), готовой"
            " таблицы mu/ro для неё у NIST нет.",
        ])
    print("\nсводка: %s" % op)


if __name__ == "__main__":
    main()
