# -*- coding: utf-8 -*-
r"""Библиотека метода 2 для Th-232 с порогом 0,5 % (оператор 25.09, вариант «а»): конфиг донора th232.yaml,
library.lines заменён γ-линиями ENSDF (data\ensdf_th232_chain_lines.csv донора) с I >= 0,5 % в окне 150–2900 кэВ.
Линии прежней библиотеки берутся из конфига как есть. Сумм-пики — без изменений.
#XR-1 (оператор 25.09 «рентген K L учитываем всегда во всех методах»): вторым проходом добавлены K/L-линии
рентгена всех 8 звеньев БЕЗ порога по интенсивности (`add_xray_lines`); их энергии должны совпасть с узлами
сетки #XR-1 (`C:\g4work\gs2020\run_marinelli\out_v5\jobs_xray.txt`, 42 энергии) до 0,01 кэВ.
Выход: scripts\configs\th232_gs2020_lib05.yaml. Спека: scripts\specs\SPEC-make_th232_lib05.md"""

import csv
import os
import sys
import yaml

sys.stdout.reconfigure(encoding="utf-8")

DONOR = r"D:\repos-folder\repos\geant4-detector-models\detectors\Gamma-1S\web-th232"
SRC_CFG = os.path.join(DONOR, "configs", "th232.yaml")
ENSDF = os.path.join(DONOR, "data", "ensdf_th232_chain_lines.csv")
KEYS = {"Ac228", "Ra224", "Pb212", "Bi212", "Tl208"}
# #M2FULL-ALL (оператор 27.09 «делай сейчас», таблица #CFG-1 согласована): THR=0 — «метод 2 полный», без порога по
# интенсивности, аналог режима донора Гамма-1С. Второй параметр argv — суффикс имени выходного файла.
THR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
OUT_SUF = sys.argv[2] if len(sys.argv) > 2 else "lib05"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs", "th232_gs2020_%s.yaml" % OUT_SUF)
LO, HI = 150.0, 2900.0
# #XR-1 (оператор 25.09 «рентген K L учитываем всегда во всех методах»): K- и L-рентген ВСЕХ звеньев цепочки,
# без порога по интенсивности (значения ENSDF местами < 0,01 %, но правило — «всегда», не «выше порога»).
# Энергии должны совпасть с узлами сетки run_marinelli grid_mar_E*.csv до 0,01 кэВ (без интерполяции, см.
# export_amticseu_data.grid_response докстринг) — список узлов #XR-1 (jobs_xray.txt) построен из ЭТОГО же CSV.
XRAY_KEYS = {"Ac228", "Ra224", "Pb212", "Bi212", "Tl208", "Th228", "Rn220", "Th232"}


def add_xray_lines(new_lines, seen_keys):
    """#XR-1: K/L-рентген всех звеньев, без порога THR — энергии совпадают с узлами #XR-1 (jobs_xray.txt)."""
    n = 0
    with open(ENSDF, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["line_type"] != "xray" or row["nuclide"] not in XRAY_KEYS:
                continue
            try:
                e_val, i_val = float(row["E_keV"]), float(row["I_percent"])
            except ValueError:
                continue
            key = round(e_val, 3)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            new_lines.append({"e_kev": e_val, "i_pct": i_val, "nuclide": row["nuclide"],
                               "note": "ENSDF, K/L-рентген, #XR-1"})
            n += 1
    return n


def main():
    with open(SRC_CFG, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    old = {}
    for l in cfg["library"]["lines"]:
        key = round(float(l["e_kev"]), 3)
        old[key] = l

    new_lines = []
    seen_keys = set()

    with open(ENSDF, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["line_type"] != "gamma":
                continue
            if row["nuclide"] not in KEYS:
                continue
            i_str = row["I_percent"].strip()
            if not i_str:
                continue
            try:
                i_val = float(i_str)
            except ValueError:
                continue
            if i_val < THR:
                continue
            try:
                e_val = float(row["E_keV"])
            except ValueError:
                continue
            if not (LO <= e_val <= HI):
                continue

            key = round(e_val, 3)
            if key in seen_keys:
                continue
            seen_keys.add(key)

            if key in old:
                new_lines.append(old[key])
            else:
                new_lines.append({
                    "e_kev": e_val,
                    "i_pct": i_val,
                    "nuclide": row["nuclide"],
                    "note": "ENSDF, порог 0,5 %"
                })

    missing = []
    for k in old:
        found = False
        for nl in new_lines:
            if abs(float(nl["e_kev"]) - k) < 0.001:
                found = True
                break
        if not found:
            missing.append(k)

    # Линии прежней библиотеки вне окна 150–2900 (129,065 кэВ) сохраняются как есть; внутри окна пропажа — отказ.
    lost = [k for k in missing if LO <= k <= HI]
    if lost:
        raise SystemExit("ОТКАЗ: линии прежней библиотеки в окне не попали в выборку: " + str(lost))
    for k in missing:
        new_lines.append(old[k])
        print("сохранена линия прежней библиотеки вне окна: %.3f кэВ" % k)

    n_xray = add_xray_lines(new_lines, seen_keys)
    print(f"#XR-1: добавлено рентгеновских линий {n_xray}")

    new_lines.sort(key=lambda x: float(x["e_kev"]))
    cfg["library"]["lines"] = new_lines
    cfg["library"]["intensity_threshold_pct"] = THR
    if THR <= 0:
        cfg["library"]["note_threshold"] = "M2FULL-ALL: без порога по интенсивности (все линии ENSDF окна LO-HI)"

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)

    print(f"Старое количество линий: {len(old)}")
    print(f"Новое количество линий: {len(new_lines)}")

    added = []
    for nl in new_lines:
        k = round(float(nl["e_kev"]), 3)
        if k not in old:
            added.append(nl)

    for a in added:
        print(f"Добавлено: {a['nuclide']} {float(a['e_kev']):.3f} кэВ, интенсивность {float(a['i_pct'])}%")

    print(f"Выходной файл: {OUT}")


if __name__ == "__main__":
    main()
