"""Пайлап в аттестации точечных кривых: скорости счёта записей поверки 2024.

Задача 112 (предписание аудитора). Случайные наложения выносят счёт из пика
в сумм-континуум; живое время это НЕ компенсирует (оно чинит только мёртвое
время тракта). Если аттестация площадей не вводила поправку на наложения,
аттестованная eps занижена на ~2*тау*R, причём на 5 см скорость счёта в
десятки раз выше, чем на 25 см, — кандидат на ту часть асимметрии отношений
модель/ЛСРМ (1,223 против 1,111), которую не объясняет TCC.

Скорости берутся из заголовков .spe самой поверки (CPS, TLIVE, TREAL) —
ни одного прогона не требуется. Заодно из (TREAL-TLIVE)/TREAL и R
восстанавливается фактическая постоянная мёртвого времени тракта — контроль
принятого в kit_recalc тау_форм = 3 мкс.

ВНИМАНИЕ к формату .efr: третье число в строке линии — это интенсивность
I%, а не время измерения; времена есть только в .spe. (Уточнение к
формулировке предписания.)

Номера источников из SHIFR в выход не пишутся — в таблице только нуклид и
геометрия (правило обезличивания репозитория).
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

TAU_SHAPE = 3.0e-6      # мкс, как в kit_recalc
POVERKA = os.path.join("detectors", "Gamma-1S", "raw_lsrm", "Work", "BG",
                       "Gamma-1S", "Spe - поверки", "Поверка 2024")
DIRS = {"p5cm": "Точечная-5см", "p25cm": "Точечная-25см"}


# Объявление наблюдаемой — что именно за число лежит в таблице. Без него
# таблицу нельзя сравнивать ни с какой другой: за один вечер 30.07.2026
# подмена определения стоила вывода четыре раза (method-rules §5).
OBS = {
    "quantity":
        "скорость счёта поверочной записи и оценка доли пика; потерянной на наложения; как 2*тау*R",
    "area":
        "не применимо — площади не снимаются; берутся полный счёт и живое время из шапки записи",
    "window":
        "не применимо — величина относится к записи целиком; не к линии",
    "shelf":
        "не применимо",
    "blurred":
        "не применимо — измеренные записи как есть",
}


def _stamp(inputs=None):
    return stamp.lines("detectors/Gamma-1S/analysis/attestation_pileup.py", OBS,
                       inputs=inputs,
                       geometry_dir=str(paths.geometry("Gamma-1S")),
                       names=stamp.SRC_LISTS["Gamma-1S"],
                       repo_dir=str(paths.REPO))


def header(path):
    raw = open(path, "rb").read()
    i = raw.find(b"SPECTR=")
    txt = raw[:i].decode("cp1251", errors="replace")
    g = {}
    for ln in txt.splitlines():
        if "=" in ln:
            k, v = ln.split("=", 1)
            g[k] = v
    return g


if __name__ == "__main__":
    root = paths.require_spectravibe("скорости счёта записей поверки 2024")
    rows = []
    print("Пайлап аттестации: 2*тау*R при тау=3 мкс; тау_мёртв — из самих"
          " записей.\n")
    print("%-6s %-8s %9s %8s %9s %9s" %
          ("геом.", "нуклид", "CPS", "мёртв,%", "потеря,%", "тау_м,мкс"))
    for tag, d in DIRS.items():
        dp = os.path.join(str(root), POVERKA, d)
        for fn in sorted(os.listdir(dp)):
            if not fn.lower().endswith(".spe"):
                continue
            g = header(os.path.join(dp, fn))
            m = re.match(r"([A-Z][a-z]?-\d+)", fn)
            nuc = m.group(1) if m else "?"
            cps = float(g.get("CPS", "0"))
            tl, tr = float(g.get("TLIVE", "0")), float(g.get("TREAL", "0"))
            dead = (tr - tl) / tr * 100 if tr else 0.0
            loss = 2 * TAU_SHAPE * cps * 100
            tau_d = (dead / 100) / cps * 1e6 if cps else 0.0
            print("%-6s %-8s %9.1f %8.2f %9.2f %9.2f"
                  % (tag, nuc, cps, dead, loss, tau_d))
            rows.append((tag, nuc, cps, tl, dead, loss, tau_d))
    out = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results",
        "attestation_pileup.csv"))
    csvio.write(
        out,
        ["geometry", "nuclide", "cps", "tlive_s", "dead_pct",
         "pileup_loss_pct", "tau_dead_us"],
        [(t, n, "%.1f" % c, "%.1f" % tl, "%.3f" % dd, "%.3f" % lo,
          "%.2f" % td) for t, n, c, tl, dd, lo, td in rows],
        comments=[
            "Скорости счёта записей поверки 2024 (точечные геометрии) и"
            " оценка потерь пика на наложения 2*тау*R (тау=3 мкс).",
            "Живое время компенсирует мёртвое время, но НЕ потери пика в"
            " сумм-континуум: если аттестация не вводила поправку,",
            "аттестованная eps занижена на pileup_loss_pct, сильнее всего"
            " на 5 см. tau_dead_us = dead/(R) - контроль постоянной тракта.",
        ],
        stamp=_stamp())
    print("\nтаблица: %s (%d строк)" % (out, len(rows)))
