#!/usr/bin/env bash
# Слияние отрезков + 6 подгонок GS2020 Th-232 (М1, М2, М2-2 × фон воды / лаборатории) СТРОГО ПОСЛЕДОВАТЕЛЬНО
# (гонка на fwhm_points_gs2020.csv). Параметры — как этапы 1–2 (SESSION-STATE). Запуск: run_gs2020_fits.sh <папка out_*>
S="<WORKDIR>/GEANT4/scripts"; N="$1"; L=/c/g4work/gs2020/run_marinelli/$N
[ -d "$L/chunks" ] || { echo "ОТКАЗ: нет $L/chunks"; exit 2; }
export GS_OUT="C:\\g4work\\gs2020\\run_marinelli\\$N"
. "$S/gs2020_fit_env.sh"   # W-155: параметры подгонки — из одного файла с KCl
python "$S/merge_templates_gs2020.py" "$L/chunks" "$L" > "$L/merge.log" 2>&1; echo "merge rc=$?"
FULL="$S/configs/th232_gs2020_full_noThresh.yaml"
run() {  # $1 лог, $2 скрипт, $3 GS_BG_WATER (1|пусто), $4 GS_M2_CONFIG (путь|пусто)
  ( [ -n "$3" ] && export GS_BG_WATER=1; [ -n "$4" ] && export GS_M2_CONFIG="$4"
    python "$S/$2" > "$L/$1.log" 2>&1 ); echo "$1 rc=$?"
}
run m1_bgw fit_gs2020_th232_m1.py 1 "" && run m1_lines fit_gs2020_th232_m1.py "" "" &&
run m2_bgw fit_gs2020_th232_m2.py 1 "" && run m2_lines fit_gs2020_th232_m2.py "" "" &&
run m2f_bgw fit_gs2020_th232_m2.py 1 "$FULL" && run m2f_lines fit_gs2020_th232_m2.py "" "$FULL"
for f in m1_bgw m1_lines m2_bgw m2_lines m2f_bgw m2f_lines; do echo "== $f"; grep -E "ФОРМА ПИКОВ|ЦЕПОЧКА В РАВНОВЕСИИ" "$L/$f.log"; done
