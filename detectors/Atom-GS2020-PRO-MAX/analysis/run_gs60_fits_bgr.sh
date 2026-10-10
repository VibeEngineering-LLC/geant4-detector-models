#!/usr/bin/env bash
# #GS-60 (08.10, оператор «Встроить и опубликовать»): все подгонки GS2020 как run_gs2020_extra_accept.sh (GS_EXTRA=1), но фон × r(E):
# KCl — r_ext_4pi_kcl1l.csv, торий — r_ext_4pi_th.csv (gs44_ext_ratio.py). Опубликованные JSON перезаписываются (копии keep_pre_gs60_2026-10-08).
# Копии новых JSON и логов — <папка>/gs60_bgr/. Запуск: bash run_gs60_fits_bgr.sh
set -u
S="<WORKDIR>/GEANT4/scripts"; X='C:\g4work\gs2020\ext'
K=/c/g4work/gs2020/kcl_1l_v4w85_83; T=/c/g4work/gs2020/run_marinelli/out_v5_oisn10
cd "$S" || exit 1
. ./gs2020_fit_env.sh; export GS_CAL_SHAPE=1 GS_EXTRA=1
DK="$K/gs60_bgr"; DT="$T/gs60_bgr"; mkdir -p "$DK" "$DT"
GS_BG_R="$X\\r_ext_4pi_kcl1l.csv" bash run_gs2020_kcl.sh > "$DK/run_kcl.out" 2>&1 || { echo "ОТКАЗ KCl rc=$?"; exit 1; }
cp -p "$K"/fit_kcl_bgw.json "$K"/fit_kcl_bgw_tail.json "$K"/fit_kcl_m2_bgw.json "$K"/fit_kcl_main.log "$K"/fit_kcl_tail.log "$K"/fit_kcl_m2.log "$DK"/
export GS_OUT='C:\g4work\gs2020\run_marinelli\out_v5_oisn10' GS_BG_WATER=1 GS_BG_R="$X\\r_ext_4pi_th.csv"
python fit_gs2020_th232_m1.py > "$DT/m1_bgw.log" 2>&1 || { echo "ОТКАЗ М1 тория rc=$?"; exit 1; }
python fit_gs2020_th232_m2.py > "$DT/m2_bgw.log" 2>&1 || { echo "ОТКАЗ М2 тория rc=$?"; exit 1; }
GS_M2_CONFIG='configs\th232_gs2020_full_noThresh.yaml' python fit_gs2020_th232_m2.py > "$DT/m2f_bgw.log" 2>&1 || { echo "ОТКАЗ М2 full rc=$?"; exit 1; }
cp -p "$T"/fit_m1_bgw_tail0_fwscale_calshape_calsum.json "$T"/fit_m2_bgw_tail0.json "$T"/fit_m2_full_noThresh_bgw_tail0.json "$DT"/
grep -h "ФОН × r" "$K"/fit_kcl_main.log "$DT"/m1_bgw.log | head -2
echo "ГОТОВО → $DK, $DT"
