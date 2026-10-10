#!/usr/bin/env bash
# #GS-42 приёмка: все подгонки GS2020 СТРОГО последовательно при заданном GS_EXTRA (KCl: run_gs2020_kcl.sh = tail → main → m2;
# торий, фон воды: М1, М2, М2 full_noThresh). Выходы JSON и логи копируются в <папка>/gs42_extra<N>_<метка>/.
# Запуск: bash run_gs2020_extra_accept.sh 0|1 <метка>
set -u
X="$1"; TAG="$2"; S="<WORKDIR>/GEANT4/scripts"
K=/c/g4work/gs2020/kcl_1l_v4w85_83; T=/c/g4work/gs2020/run_marinelli/out_v5_oisn10
cd "$S" || exit 1
. ./gs2020_fit_env.sh; export GS_CAL_SHAPE=1 GS_EXTRA="$X"
DK="$K/gs42_extra${X}_$TAG"; DT="$T/gs42_extra${X}_$TAG"; mkdir -p "$DK" "$DT"
bash run_gs2020_kcl.sh > "$DK/run_kcl.out" 2>&1 || { echo "ОТКАЗ KCl rc=$?"; exit 1; }
cp -p "$K"/fit_kcl_bgw.json "$K"/fit_kcl_bgw_tail.json "$K"/fit_kcl_m2_bgw.json "$K"/fit_kcl_main.log "$K"/fit_kcl_tail.log "$K"/fit_kcl_m2.log "$DK"/
export GS_OUT='C:\g4work\gs2020\run_marinelli\out_v5_oisn10' GS_BG_WATER=1
python fit_gs2020_th232_m1.py > "$DT/m1_bgw.log" 2>&1 || { echo "ОТКАЗ М1 тория rc=$?"; exit 1; }
python fit_gs2020_th232_m2.py > "$DT/m2_bgw.log" 2>&1 || { echo "ОТКАЗ М2 тория rc=$?"; exit 1; }
GS_M2_CONFIG='configs\th232_gs2020_full_noThresh.yaml' python fit_gs2020_th232_m2.py > "$DT/m2f_bgw.log" 2>&1 || { echo "ОТКАЗ М2 full rc=$?"; exit 1; }
cp -p "$T"/fit_m1_bgw_tail0_fwscale_calshape_calsum.json "$T"/fit_m2_bgw_tail0.json "$T"/fit_m2_full_noThresh_bgw_tail0.json "$DT"/
echo "ГОТОВО GS_EXTRA=$X → $DK, $DT"
