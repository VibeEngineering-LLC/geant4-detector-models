#!/usr/bin/env bash
# Этап 2 (W-153): М2 и М2-2 на исправленной геометрии, фон воды и лаборатории — 4 прогона СТРОГО ПОСЛЕДОВАТЕЛЬНО
# (гонка на fwhm_points_gs2020.csv). Параметры — как у М1 этапа 1 (SESSION-STATE, хендофф CMP-2).
S="<WORKDIR>/GEANT4/scripts"
export GS_OUT='C:\g4work\gs2020\run_marinelli\out_v5_lvl1v2' PYTHONIOENCODING=utf-8
export GS_TAIL=0 GS_CAL_SL=4 GS_CAL_SUM="6785:3187"
export GS_FWHM_SCALE="238.632:1.05,583.187:1.05,911.204:1.05,964.766:1.05,968.971:1.05,2614.511:1.05"
L=/c/g4work/gs2020/run_marinelli/out_v5_lvl1v2
FULL="$S/configs/th232_gs2020_full_noThresh.yaml"
run() {  # $1 имя лога, $2 GS_BG_WATER (1|пусто), $3 GS_M2_CONFIG (путь|пусто)
  ( [ -n "$2" ] && export GS_BG_WATER=1; [ -n "$3" ] && export GS_M2_CONFIG="$3"
    python "$S/fit_gs2020_th232_m2.py" > "$L/$1.log" 2>&1 ); echo "$1 rc=$?"
}
run m2_bgw 1 "" && run m2_lines "" "" && run m2f_bgw 1 "$FULL" && run m2f_lines "" "$FULL"
for f in m2_bgw m2_lines m2f_bgw m2f_lines; do echo "== $f"; grep -E "ФОРМА ПИКОВ|ЦЕПОЧКА В РАВНОВЕСИИ" "$L/$f.log"; done
