"""Поправки на загрузку тракта по всем записям комплекта.

Два разных эффекта, которые нельзя смешивать [Будыка §8.1; ГОСТ 26874-86]:

1. МЁРТВОЕ ВРЕМЯ АЦП. Уже учтено делением на живое время, но стоит проверить
   его честность: эффективное тау на импульс
       тау_АЦП = (T_real - T_live) / N_импульсов
   обязано быть примерно ОДИНАКОВЫМ у всех записей одного тракта. Если у
   какой-то записи оно выпадает — живому времени этой записи веры меньше.

2. СЛУЧАЙНЫЕ НАЛОЖЕНИЯ (pile-up). Два импульса в пределах времени формирования
   сливаются в один с неверной амплитудой; событие уходит из пика. Живое время
   этого НЕ учитывает. Поправка к скорости счёта пика:
       A_ист = A * exp(2 * тау_форм * R),
   R — полная скорость счёта [ЛСРМ; Будыка §8.1]. Для NaI тау_форм ~ 3 мкс.

Обе величины печатаются по каждой записи; помечаются те, где поправка на
наложения превышает 0,1 %.
"""
import glob
import os
import re
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402

KIT = str(paths.ref("Gamma-1S"))
TAU_SHAPE = 3.0e-6     # с; ДОПУЩЕНИЕ для NaI, см. шапку

if __name__ == "__main__":
    rows = []
    # Ищем рекурсивно: в скачанном наборе XML лежат на два уровня ниже
    # корня, в закоммиченном — на три (геометрия/нуклид/файл).
    KIT_GEOM = ("Marinelli_1L", "Denta_120mL", "Petri_60mL",
                "Point_25cm", "Point_5cm")
    for p in sorted(str(q) for q in paths.ref("Gamma-1S").rglob("*.xml")):
        s, b = bm.read(p)
        parts = p.replace("\\", "/").split("/")
        geom = next((x for x in reversed(parts) if x in KIT_GEOM),
                    os.path.basename(os.path.dirname(p)))
        name = os.path.basename(p)[7:41]
        N = float(s.n.sum())
        R = N / s.live
        dt = s.real - s.live
        tau_adc = dt / N if N else 0.0
        pile = 2 * TAU_SHAPE * R          # ~ exp(x)-1 при малых x
        rows.append((geom, name, R, 100 * dt / s.real, 1e6 * tau_adc,
                     100 * pile))

    print("Загрузка тракта: %d записей. тау_форм = %.1f мкс (допущение)\n"
          % (len(rows), 1e6 * TAU_SHAPE))
    print("%-13s %-34s %8s %7s %10s %10s" %
          ("геометрия", "запись", "R, имп/с", "мёрт,%", "тауАЦП,мкс",
           "налож., %"))
    for g, n, R, d, t, pi in rows:
        mark = "  <-" if pi > 0.1 else ""
        print("%-13s %-34s %8.1f %7.2f %10.1f %10.3f%s"
              % (g, n, R, d, t, pi, mark))

    taus = sorted(t for _, _, _, _, t, _ in rows if t > 0)
    if not taus:
        raise SystemExit(
            "не найдено ни одного спектра комплекта (*.xml) в %s.\n"
            "Комплект в формате BecqMoni лежит в "
            "reference/lsrm/reference_spectra/reference_kits_becqmoni;\n"
            "переменная G4MODELS_REF, если задана, должна указывать на "
            "каталог с ним." % paths.ref("Gamma-1S"))
    med = taus[len(taus) // 2]
    print("\nмедианное тау_АЦП = %.0f мкс; разброс %0.0f..%0.0f"
          % (med, taus[0], taus[-1]))
    print("Постоянство тау по записям подтверждает честность живого времени;")
    print("выпадающие записи см. по столбцу.")
