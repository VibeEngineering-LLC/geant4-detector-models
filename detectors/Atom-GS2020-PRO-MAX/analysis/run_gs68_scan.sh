#!/usr/bin/env bash
# #GS-68 перебор вариантов BecqMoni на тории М1 (оператор 10.10 «в бекмони несколько вариантов кривой калибровки ПШПВ и формы пика.
# Нужно все проверять»): кривая ПШПВ × форма пика × множитель ширины. Итог — из логов (ФОРМА ПИКОВ ±1,5 и ±3). ≤20 потоков.
S="<WORKDIR>/GEANT4/scripts"; T=/c/g4work/gs2020/run_marinelli/out_v5_oisn10; O=/c/g4work/gs2020/npsm68/scan; mkdir -p $O/keep
cp -p $T/fit_m1_bgw_tail0_fwscale_calshape_calsum.json $O/keep/   # рабочий JSON тория — восстановить в конце
J=$O/jobs.txt; NP=${NP:-20}
if [ "${RETRY:-0}" = 1 ]; then   # дозапуск: только задания без итога в логе (10.10: 314 упали по памяти при -P 20)
  while read -r C K P a b; do t="${C}_k${K}_${P}_${a:-0}_${b:-0}"; grep -q "ЦЕПОЧКА В РАВНОВЕСИИ" "$O/$t.log" 2>/dev/null || echo "$C $K $P $a $b"; done < $J > $O/jobs_retry.txt; J=$O/jobs_retry.txt
else : > $J
for C in interp sqrt sqrtpoly power sqrtpoly+k40 power+k40 cfw; do for K in 1.00 1.025 1.05 1.075; do [ $C = cfw ] && [ $K != 1.00 ] && continue
  for P in "gauss" "ege 1.0 1.0" "ege 1.6 1.6" "ege 2.4 2.4" "ege 1.0 2.4" "ege 2.4 1.0" "ege 1.6 2.4" "ege 2.4 1.6" "ege 1.0 1.6" "ege 1.6 1.0" \
           "pvoigt 0.05" "pvoigt 0.1" "pvoigt 0.2" "gain 0.10 0.020" "gain 0.18 0.019" "gain 0.14 0.028" "gain 0.25 0.015"; do
    echo "$C $K $P" >> $J; done; done; done
fi
export S O; xargs -P $NP -L 1 bash -c '
  C=$0; K=$1; P=$2; a=${3:-0}; b=${4:-0}; tag="${C}_k${K}_${P}_${a}_${b}"
  cd "$S" && . ./gs2020_fit_env.sh && export GS_CAL_SHAPE=1 GS_EXTRA=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  export GS_OUT="C:\g4work\gs2020\run_marinelli\out_v5_oisn10" GS_BG_WATER=1 GS_BG_R="C:\g4work\gs2020\ext\r_ext_4pi_th.csv" GS_FWHM_CURVE=${C%+k40} GS_FWHM_CSV="C:/g4work/gs2020/npsm68/scan/fw_$tag.csv"
  [ "${C#*+}" = k40 ] && export GS_FWHM_K40=1; [ $C = cfw ] && export GS_FWHM_CFW=1 GS_FWHM_CURVE=interp
  export GS_FWHM_SCALE="238.632:$K,583.187:$K,911.204:$K,964.766:$K,968.971:$K,1460.822:$K,2614.511:$K"
  case $P in ege) export GS_PEAK=ege GS_EGE_KL=$a GS_EGE_KR=$b;; pvoigt) export GS_PEAK=pvoigt GS_PV_ETA=$a;; gain) export GS_GAIN_W=$a GS_GAIN_S=$b;; esac
  python fit_gs2020_th232_m1.py > "$O/$tag.log" 2>&1; echo "$tag rc=$?"' < $J > $O/run.log 2>&1
cp -p $O/keep/fit_m1_bgw_tail0_fwscale_calshape_calsum.json $T/ && echo "ГОТОВО (рабочий JSON тория восстановлен)" >> $O/run.log
