#!/usr/bin/env bash
# #GS-70 шаблоны метода 1 для черники (#CFG-1 «Да, по таблице» 10.10): Cs-137 (до Ba-137m) и K-40 по 2·10^7, Sr-90 и Y-90 по 10^8
# кусками 5·10^6; NP=12 (память Geant4 мала, рядом идёт перебор #GS-68 на 8 потоках — итого ≤20). Слияние — merge_templates_gs2020.py.
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
G=C:/g4work/gs2020/berry/GS2020_marinelli_berry_dry_1000ml_v4_w85_83.gdml; X=C:/g4work/build/gs2020-marinelli-npsm/gs2020_marinelli.exe
O=/c/g4work/gs2020/berry; mkdir -p $O/chunks; J=$O/jobs_templates.txt; : > $J; s=83000
for j in "Cs137 55:137 56 4" "K40 19:40 19 4" "Sr90 38:90 38 20" "Y90 39:90 39 20"; do set -- $j
  for i in $(seq 1 $4); do s=$((s+1)); echo "$1 $2 $3 $s" >> $J; done; done
export G X O; xargs -P ${NP:-12} -L 1 bash -c 'GS2020_ZMAX=$2 GS2020_GDML=$G $X ion:$1 5000000 $O/chunks/mix_$0_s$3.csv $3 > $O/chunks/mix_$0_s$3.log 2>&1; echo "$0 s$3 rc=$?"' < $J > $O/run_templates.log 2>&1
echo "ГОТОВО" >> $O/run_templates.log
