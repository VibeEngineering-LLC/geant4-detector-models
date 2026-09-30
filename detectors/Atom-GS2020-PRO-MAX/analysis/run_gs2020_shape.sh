#!/usr/bin/env bash
# #GS-37: пересчёт GS2020 на шкале по ФОРМЕ (GS_CAL_SHAPE=1, cal_shape.json). Строго последовательно (#М1/М2).
# 1) копия текущих результатов (опубликованный гибрид) в base_hybrid_2026-09-28\ — ничего не затирается без копии;
# 2) калибровка по форме (синтетика — отдельно, см. SESSION-STATE); 3) 6 подгонок Th; 4) 3 подгонки KCl.
S="<WORKDIR>/GEANT4/scripts"; L=/c/g4work/gs2020/run_marinelli/out_v5_oisn10; K=/c/g4work/gs2020/kcl
for D in "$L" "$K"; do mkdir -p "$D/base_hybrid_2026-09-28"; cp -n "$D"/*.json "$D"/*.log "$D/base_hybrid_2026-09-28/" 2>/dev/null; done
echo "копия: $(ls "$L/base_hybrid_2026-09-28" | wc -l) + $(ls "$K/base_hybrid_2026-09-28" | wc -l) файлов"
export GS_OUT='C:\g4work\gs2020\run_marinelli\out_v5_oisn10'; . "$S/gs2020_fit_env.sh"
( cd "$S/cal" && python gs2020_calib_shape.py bg bgw kcl sample > "$L/cal_shape.log" 2>&1 ) || { echo "ОТКАЗ: калибровка rc=$?"; exit 1; }
export GS_CAL_SHAPE=1
FULL="$S/configs/th232_gs2020_full_noThresh.yaml"
run() {  # $1 лог, $2 скрипт, $3 GS_BG_WATER (1|пусто), $4 GS_M2_CONFIG (путь|пусто)
  ( [ -n "$3" ] && export GS_BG_WATER=1; [ -n "$4" ] && export GS_M2_CONFIG="$4"
    python "$S/$2" > "$L/$1.log" 2>&1 ); r=$?; echo "$1 rc=$r"; return $r
}
run m1_bgw fit_gs2020_th232_m1.py 1 "" && run m1_lines fit_gs2020_th232_m1.py "" "" &&
run m2_bgw fit_gs2020_th232_m2.py 1 "" && run m2_lines fit_gs2020_th232_m2.py "" "" &&
run m2f_bgw fit_gs2020_th232_m2.py 1 "$FULL" && run m2f_lines fit_gs2020_th232_m2.py "" "$FULL" || exit 1
for f in m1_bgw m1_lines m2_bgw m2_lines m2f_bgw m2f_lines; do echo "== $f"; grep -E "ШКАЛА|#GS-37|ФОРМА ПИКОВ|ЦЕПОЧКА В РАВНОВЕСИИ" "$L/$f.log"; done
bash "$S/run_gs2020_kcl.sh"; echo "kcl rc=$?"
