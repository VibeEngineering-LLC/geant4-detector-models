#!/usr/bin/env bash
# #GS-62: полный пересчёт KCl на выровненной сумме 88,3 ч: калибровка по форме (все теги) → своя ширина → М1/М2 → выгрузка → сборка страницы.
S="<WORKDIR>/GEANT4/scripts"; O=/c/g4work/gs2020/run_marinelli/out_v5; K=/c/g4work/gs2020/kcl_1l_v4w85_83; P="../web/gs2020-th232-page"; L=/c/g4work/gs2020/kcl_parts/aligned
mkdir -p $L; cd "$S" || exit 1; cp -p ../results/gs61_own_fwhm/kcl.csv ../results/gs61_own_fwhm/kcl_88h_raw.csv
( . ./gs2020_fit_env.sh && unset GS_CAL_SHAPE && nice python cal/gs2020_calib_shape.py > $L/cal_shape.log 2>&1 ); echo "shape rc=$?"; cp -p $O/cal_shape.json $L/
export GS_CAL_SHAPE=1 PYTHONIOENCODING=utf-8 SPECTRAVIBE_ROOT="<DOSIM>/ИИ/1 Скилы/0_Work/gamma-spectrum-analysis"
( . ./gs2020_fit_env.sh; python gs61_own_fwhm.py kcl > $L/own_fwhm.log 2>&1 ); echo "own rc=$?"
( export GS_EXTRA=1 GS_BG_R='C:\g4work\gs2020\ext\r_ext_4pi_kcl1l.csv'; nice bash run_gs2020_kcl.sh > $L/run_fits.out 2>&1 ); echo "fits rc=$?"
( . ./gs2020_fit_env.sh && export GS_OUT='C:\g4work\gs2020\kcl_1l_v4w85_83' GS_BG_WATER=1 GS_EXTRA=1 GS_BG_R='C:\g4work\gs2020\ext\r_ext_4pi_kcl1l.csv' GS_FWHM_OWN='<WORKDIR>\GEANT4\results\gs61_own_fwhm\kcl.csv' && python export_page_gs2020_k40.py > $L/export_k40.log 2>&1; echo "exp rc=$?"; python export_sources_gs2020.py > $L/export_sources.log 2>&1; echo "src rc=$?"; unset GS_OUT GS_BG_WATER GS_EXTRA GS_BG_R GS_FWHM_OWN; python build_page_gs2020.py > $L/build.log 2>&1; echo "build rc=$?" )
cp -p "$P/gs2020_k40_data.json" $L/
