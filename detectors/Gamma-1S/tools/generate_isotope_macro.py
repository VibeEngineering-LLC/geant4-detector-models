#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Генератор макроса Geant4 для прогонов «один нуклид на прогон» ветви
распада, заданной в configs/<источник>.yaml (поле geant4.isotope_runs).

Зачем. Язык макрокоманд Geant4 не знает переменных и циклов -- раньше
блок (nucleusLimits/gps ion/outFile/beamOn) на каждый нуклид набивался
руками, GPS-регион пробы дублировался буквально в семи файлах ветви
Th-232, а добор статистики слабых звеньев (R66) жил ВТОРЫМ файлом,
молча перезаписывающим вывод первого для части нуклидов -- ловушка
именно того класса, что уже дважды случалась в этом дереве макросов
(decay.mac/decay_all.mac и decay_control.mac/decay_all.mac, см. их
комментарии). Здесь одно число n_events на нуклид в конфиге, один
генерируемый файл, коллизия невозможна структурно.

Инвентаризация (задача #176, агент A) отдельно отметила: файл, который
пропускает обязательную строку `thresholdForVeryLongDecayTime`
(decay.mac), даёт молчаливый нулевой выход на долгоживущих ядрах.
Генератор пишет эту строку ВСЕГДА -- пропустить её невозможно так же,
как нельзя разойтись в тексте GPS-региона.

Область действия (сознательно, этап 3 задачи #175): только
decay_th232_isotopes.mac -- он используется ИСКЛЮЧИТЕЛЬНО этим
конвейером (export_data.py читает iso_<Nuc>.csv и больше ничего).
Полноцепочечный прогон chain_Th232.csv лежит в decay_all.mac вместе с
Cs-137/K-40 -- общая инфраструктура с другими потребителями, генератор
её не трогает; это отдельная задача, не молчаливое упущение.

Использование:
    python tools/generate_isotope_macro.py [--config PATH] [--out PATH]
Без аргументов: configs/th232.yaml -> macros/decay_th232_isotopes.mac
(тот же путь, что раньше писался руками).
"""
import argparse
import os

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG = os.path.join(HERE, "..", "web-th232", "configs", "th232.yaml")
DEFAULT_OUT = os.path.join(HERE, "..", "macros", "decay_th232_isotopes.mac")


def render(cfg):
    g4 = cfg["geant4"]
    reg = g4["gps_region"]
    cx, cy, cz = reg["centre_mm"]
    src_id = cfg["source"]["id"]

    lines = []
    lines.append("# Отдельные прогоны каждого γ-испускающего звена ветви распада"
                  " источника %s." % src_id)
    lines.append("# СГЕНЕРИРОВАНО: tools/generate_isotope_macro.py из"
                  " web-th232/configs/%s.yaml (geant4.isotope_runs)." % src_id)
    lines.append("# Руками не редактировать -- правки потеряются при следующей"
                  " генерации; менять n_events/состав в конфиге.")
    lines.append("#")
    lines.append("# Ключ разделения: nucleusLimits в Geant4 ограничивает диапазон"
                  " (A_min,A_max,")
    lines.append("# Z_min,Z_max) распадающихся ядер. Задавая один нуклид на прогон,"
                  " получаем")
    lines.append("# ТОЛЬКО его собственный вклад -- цепочка ниже не пойдёт, ничего"
                  " к нему не")
    lines.append("# примешается. n_events выбран по каждому звену отдельно (не"
                  " «побольше")
    lines.append("# всем») -- слабая гамма-ветвь требует на порядки больше"
                  " распадов на тот же")
    lines.append("# отсчёт в кристалле (см. комментарий в конфиге, R66).")
    lines.append("/run/initialize")
    lines.append("/control/verbose 0")
    lines.append("/run/verbose 0")
    lines.append("/process/had/rdm/verbose 0")
    lines.append("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns")
    lines.append("")
    lines.append("/gps/particle ion")
    lines.append("/gps/energy 0 keV")
    lines.append("/gps/pos/type Volume")
    lines.append("/gps/pos/shape Cylinder")
    lines.append("/gps/pos/centre %g %g %g mm" % (cx, cy, cz))
    lines.append("/gps/pos/radius %g mm" % reg["radius_mm"])
    lines.append("/gps/pos/halfz %g mm" % reg["halfz_mm"])
    lines.append("/gps/pos/confine %s" % reg["confine_volume"])
    lines.append("/gps/ang/type iso")

    for r in g4["isotope_runs"]:
        lines.append("")
        lines.append("# --- %s (Z=%d, A=%d) ---" % (r["key"], r["z"], r["a"]))
        lines.append("/process/had/rdm/nucleusLimits %d %d %d %d"
                      % (r["a"], r["a"], r["z"], r["z"]))
        lines.append("/gps/ion %d %d 0 0" % (r["z"], r["a"]))
        lines.append("/g1s/outFile %s" % r["out_file"])
        lines.append("/run/beamOn %d" % r["n_events"])
    lines.append("")
    return "\r\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()

    with open(args.config, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    text = render(cfg)
    with open(args.out, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    print("написано: %s (%d байт, %d нуклидов)"
          % (args.out, len(text.encode("utf-8")), len(cfg["geant4"]["isotope_runs"])))


if __name__ == "__main__":
    main()
