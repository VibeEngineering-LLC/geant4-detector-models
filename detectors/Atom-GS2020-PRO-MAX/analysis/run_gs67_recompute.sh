#!/usr/bin/env bash
# #GS-67 (оператор 10.10 «Гладкая шкала + новый множитель»): все подгонки GS2020 на шкале по форме БЕЗ изломов
# (cal_shape.json = cal_shape_poly3_2026-10-09.json, GS_SHAPE_KNOTS="" DEG 3) → выгрузки → сборка. Схема — run_gs63_recompute.sh.
# Прежние JSON/логи/данные страницы — <папка>/keep_pre_gs67/. Множитель ширины 1,05 подтверждён сканом #SHAPE-1 на новой шкале.
S="<WORKDIR>/GEANT4/scripts"; X='C:\g4work\gs2020\ext'; P="../web/gs2020-th232-page"
K=/c/g4work/gs2020/kcl_1l_v4w85_83; T=/c/g4work/gs2020/run_marinelli/out_v5_oisn10
cmp -s /c/g4work/gs2020/run_marinelli/out_v5/cal_shape.json /c/g4work/gs2020/run_marinelli/out_v5/cal_shape_poly3_2026-10-09.json || { echo "ОТКАЗ: рабочая шкала не poly3"; exit 1; }
for D in "$K" "$T"; do mkdir -p "$D/keep_pre_gs67"; cp -p "$D"/fit_*.json "$D"/*.log "$D/keep_pre_gs67/" 2>/dev/null; cp -rp "$D/gs60_bgr" "$D/keep_pre_gs67/" 2>/dev/null; done
cd "$S" || exit 1; cp -p $P/gs2020_th232_data.json $P/gs2020_k40_data.json $K/keep_pre_gs67/
nice bash run_gs60_fits_bgr.sh > "$T/gs67_fits.out" 2>&1; echo "fits rc=$?"
. ./gs2020_fit_env.sh; export GS_CAL_SHAPE=1 GS_EXTRA=1 PYTHONIOENCODING=utf-8; unset GS_FWHM_OWN
(GS_OUT='C:\g4work\gs2020\run_marinelli\out_v5_oisn10' GS_BG_WATER=1 GS_BG_R="$X\r_ext_4pi_th.csv" python export_page_gs2020.py > $T/gs67_export_th.log 2>&1; echo "th rc=$?")
(export GS_OUT='C:\g4work\gs2020\kcl_1l_v4w85_83' GS_BG_WATER=1 GS_BG_R="$X\r_ext_4pi_kcl1l.csv"; python export_page_gs2020_k40.py > $K/gs67_export_k40.log 2>&1; echo "k40 rc=$?"; python export_sources_gs2020.py > $K/gs67_export_sources.log 2>&1; echo "src rc=$?")
python build_page_gs2020.py > $K/gs67_build.log 2>&1; echo "build rc=$?"
for f in "$K/fit_kcl_main.log" "$K/fit_kcl_m2.log" "$T/gs60_bgr/m1_bgw.log" "$T/gs60_bgr/m2_bgw.log" "$T/gs60_bgr/m2f_bgw.log"; do echo "== $f"; grep -E "ФОРМА ПИКОВ|K-40 \(метод|χ²/ν|ЦЕПОЧКА В РАВНОВЕСИИ" "$f" | head -6; done
