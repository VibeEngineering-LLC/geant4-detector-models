#!/usr/bin/env bash
# #GS-77 (оператор 10.10 «делай»: ядро разброса усиления 0,14/0,028 из gs2020_fit_env.sh во всех подгонках; схема — run_gs77_recompute.sh): все подгонки GS2020 на шкале по форме БЕЗ изломов
# (cal_shape.json = cal_shape_poly3_2026-10-09.json, GS_SHAPE_KNOTS="" DEG 3) → выгрузки → сборка. Схема — run_gs63_recompute.sh.
# Прежние JSON/логи/данные страницы — <папка>/keep_pre_gs77b/. Множитель ширины 1,05 подтверждён сканом #SHAPE-1 на новой шкале.
S="<WORKDIR>/GEANT4/scripts"; X='C:\g4work\gs2020\ext'; P="../web/gs2020-th232-page"
K=/c/g4work/gs2020/kcl_1l_v4w85_83; T=/c/g4work/gs2020/run_marinelli/out_v5_oisn10
# W-175: с #GS-70 в cal_shape.json есть ключ berry (GS_SHAPE_MERGE=1) — сверка с poly3 по ключам sample/bg/bgw/kcl, не файлом целиком
PYTHONIOENCODING=utf-8 python -c "import json,sys;d='C:/g4work/gs2020/run_marinelli/out_v5/';a,b=(json.load(open(d+f,encoding='utf-8')) for f in ('cal_shape.json','cal_shape_poly3_2026-10-09.json'));sys.exit(any(a.get(k)!=b.get(k) for k in ('sample','bg','bgw','kcl')))" || { echo "ОТКАЗ: рабочая шкала (sample/bg/bgw/kcl) не poly3"; exit 1; }
for D in "$K" "$T"; do mkdir -p "$D/keep_pre_gs77b"; cp -p "$D"/fit_*.json "$D"/*.log "$D/keep_pre_gs77b/" 2>/dev/null; cp -rp "$D/gs60_bgr" "$D/keep_pre_gs77b/" 2>/dev/null; done
cd "$S" || exit 1; cp -p $P/gs2020_th232_data.json $P/gs2020_k40_data.json $K/keep_pre_gs77b/
nice bash run_gs60_fits_bgr.sh > "$T/gs77_fits.out" 2>&1; echo "fits rc=$?"
. ./gs2020_fit_env.sh; export GS_CAL_SHAPE=1 GS_EXTRA=1 PYTHONIOENCODING=utf-8; unset GS_FWHM_OWN
(GS_OUT='C:\g4work\gs2020\run_marinelli\out_v5_oisn10' GS_BG_WATER=1 GS_BG_R="$X\r_ext_4pi_th.csv" python export_page_gs2020.py > $T/gs77_export_th.log 2>&1; echo "th rc=$?")
(export GS_OUT='C:\g4work\gs2020\kcl_1l_v4w85_83' GS_BG_WATER=1 GS_BG_R="$X\r_ext_4pi_kcl1l.csv"; python export_page_gs2020_k40.py > $K/gs77_export_k40.log 2>&1; echo "k40 rc=$?"; python export_sources_gs2020.py > $K/gs77_export_sources.log 2>&1; echo "src rc=$?")
python build_page_gs2020.py > $K/gs77_build.log 2>&1; echo "build rc=$?"
for f in "$K/fit_kcl_main.log" "$K/fit_kcl_m2.log" "$T/gs60_bgr/m1_bgw.log" "$T/gs60_bgr/m2_bgw.log" "$T/gs60_bgr/m2f_bgw.log"; do echo "== $f"; grep -E "ФОРМА ПИКОВ|K-40 \(метод|χ²/ν|ЦЕПОЧКА В РАВНОВЕСИИ" "$f" | head -6; done
