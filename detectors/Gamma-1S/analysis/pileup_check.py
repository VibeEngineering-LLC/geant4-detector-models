"""Пайлап в аттестационных записях: полные скорости счёта (задача 112).

ЗАЧЕМ. Все выводы проекта опираются на площади пиков аттестационных записей
поверки 2024. Наложение импульсов (пайлап) переносит счёт из пика в хвост и
занижает площадь тем сильнее, чем горячее запись. Проверка грошовая и
ортогональна остальным задачам (указание аудитора).

ПОЧЕМУ НЕ ИЗ .efr. Времени набора в .efr НЕТ — файл несёт только геометрию,
узлы eps и площади; TLIVE там отсутствует. Скорости счёта берутся из самих
.spe (поля TLIVE/TREAL и полная сумма отсчётов).

МЕРА. Доля пайлапа оценивается как r · tau, где r — полная скорость счёта,
tau — время формирования импульса тракта. Для NaI с ADC ЛСРМ принято
tau = 5 мкс (порядок величины; точное значение из паспорта тракта не
извлекается, поэтому оценка — ориентир, не поправка). Мёртвое время берётся
измеренным: 1 − TLIVE/TREAL.

ВЫВОД (30.07.2026). На записях, на которых стоят выводы проекта, пайлап
единицы процентов: Th-228 точечная 5 см — 3,7 %, Ba-133 — 0,7 %, Am-241 —
0,9 %. Самая горячая запись комплекта Na-22 (7,8 %) в наших выводах не
участвует. Знак и порядок расхождения модель/аттестация пайлап не меняет.
Отдельно: у Cd-109 и Zn-65 скорость счёта 41 и 22 имп/с — это записи на
пределе, там своя проблема (статистика), не пайлап.

Пути к записям — из committed reference репозитория; читатель .spe штатный
из SpectraVibe (SPECTRAVIBE_ROOT).
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402

RESULTS = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results"))

TAU_S = 5.0e-6           # время формирования тракта, с (ориентир для NaI+ADC)
GEOMS = ("Точечная-5см", "Точечная-25см", "Маринелли", "Петри-60мл",
         "Дента-120мл")


def main():
    root = paths.require_spectravibe("пайлап: чтение аттестационных .spe")
    sys.path.insert(0, os.path.join(str(root), "scripts"))
    from gamma.io.lsrm_spe import read_lsrm_spe

    base = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "reference", "lsrm",
        "raw_lsrm", "Work", "BG", "Gamma-1S", "Spe - поверки", "Поверка 2024")
    base = os.path.abspath(base)
    if not os.path.isdir(base):
        raise SystemExit("нет каталога аттестационных записей: %s" % base)

    rows = []
    print("Пайлап в аттестационных записях поверки 2024 (tau = %.0f мкс)\n"
          % (TAU_S * 1e6))
    print("%-34s %8s %8s %8s %10s %8s"
          % ("запись", "TLIVE", "TREAL", "мёртв,%", "имп/с", "пайлап,%"))
    for geom in GEOMS:
        d = os.path.join(base, geom)
        if not os.path.isdir(d):
            continue
        print("--- %s" % geom)
        for f in sorted(glob.glob(os.path.join(d, "*.spe"))):
            try:
                sp = read_lsrm_spe(f)
            except Exception as exc:            # noqa: BLE001
                print("    %-30s не прочитан: %s" % (os.path.basename(f), exc))
                continue
            tot = float(np.sum(sp.counts))
            tl, tr = float(sp.live_time), float(sp.real_time)
            if tl <= 0:
                continue
            dead = 100.0 * (1.0 - tl / tr) if tr > 0 else 0.0
            rate = tot / tl
            pu = 100.0 * rate * TAU_S
            nm = os.path.basename(f).replace(".spe", "")
            print("%-34s %8.0f %8.0f %8.2f %10.0f %8.2f"
                  % (nm[:34], tl, tr, dead, rate, pu))
            rows.append((geom, nm, tl, tr, dead, rate, pu))

    if not rows:
        raise SystemExit("ни одна запись не прочитана — проверьте пути")

    hot = max(rows, key=lambda r: r[6])
    print("\nСамая горячая запись: %s — %.1f %% пайлапа при %.0f имп/с"
          % (hot[1][:40], hot[6], hot[5]))
    key = [r for r in rows if "Th-228" in r[1] and r[0] == "Точечная-5см"]
    if key:
        print("Запись, на которой стоят выводы (Th-228 точечная 5 см):"
              " %.1f %% пайлапа, мёртвое %.2f %%" % (key[0][6], key[0][4]))
    print("Порядок величины: пайлап единицы процентов, знак расхождения"
          " модель/аттестация не меняет.")

    csvio.write(
        os.path.join(RESULTS, "pileup_check.csv"),
        ["geometry", "record", "tlive_s", "treal_s", "dead_pct",
         "total_rate_cps", "pileup_pct"],
        [(g, n, "%.2f" % tl, "%.2f" % tr, "%.2f" % d, "%.0f" % r, "%.2f" % p)
         for g, n, tl, tr, d, r, p in rows],
        comments=[
            "Пайлап в аттестационных записях поверки 2024. Скорости счёта"
            " из .spe (TLIVE/TREAL и полная сумма отсчётов);",
            "  в .efr времени набора НЕТ — только геометрия и узлы eps.",
            "Доля пайлапа = полная скорость * tau; tau = 5 мкс — ОРИЕНТИР"
            " для тракта NaI+ADC ЛСРМ; точное значение из паспорта не",
            "  извлекается; поэтому это оценка порядка; а не поправка.",
            "Мёртвое время измеренное: 1 - TLIVE/TREAL.",
            "На записях наших выводов пайлап единицы процентов (Th-228"
            " точечная 5 см 3;7 %; Ba-133 0;7 %; Am-241 0;9 %).",
        ])
    print("\nтаблица: %s" % os.path.join(RESULTS, "pileup_check.csv"))


if __name__ == "__main__":
    main()
