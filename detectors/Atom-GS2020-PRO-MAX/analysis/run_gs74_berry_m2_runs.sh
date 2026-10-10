#!/usr/bin/env bash
# #GS-74 метод 2 черники: сетка монолиний (γ) и β-компоненты Cs-137/K-40 в геометрии черники. Таблица: CFG1-2026-10-10-berry-m2.md.
# Куски 5·10^6 → chunks_m2\, слияние merge_templates_gs2020.py → C:\g4work\gs2020\berry\ (grid_mar_E*.csv, beta_merged\beta_*.csv).
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
G=C:/g4work/gs2020/berry/GS2020_marinelli_berry_dry_1000ml_v4_w85_83.gdml; O=/c/g4work/gs2020/berry; C=$O/chunks_m2; mkdir -p $C
XG=C:/g4work/build/gs2020-marinelli-npsm/gs2020_marinelli.exe; XB=C:/g4work/build/gs2020-marinelli-beta/gs2020_marinelli.exe
J=$O/jobs_m2.txt; : > $J; s=84000
for E in 661.657 31.817 32.194 36.304 36.378 37.255 4.466 4.828 5.156 1460.822 0.265 2.956 2.958 3.19; do
  for i in 1 2 3 4; do s=$((s+1)); echo "g grid_mar_E$E $E 0 $s" >> $J; done; done
for n in "Cs137 55:137 56" "K40 19:40 19"; do set -- $n
  for i in $(seq 1 10); do s=$((s+1)); echo "b beta_$1 ion:$2 $3 $s" >> $J; done; done
export G C XG XB
xargs -P ${NP:-12} -L 1 bash -c 'if [ $0 = g ]; then GS2020_GDML=$G $XG $2 5000000 $C/$1_s$4.csv $4 > $C/$1_s$4.log 2>&1;
  else GS2020_BETA_ONLY=1 GS2020_ZMAX=$3 GS2020_GDML=$G $XB $2 5000000 $C/$1_s$4.csv $4 > $C/$1_s$4.log 2>&1; fi; echo "$1 s$4 rc=$?"' < $J > $O/run_m2.log 2>&1
grep -c "rc=0" $O/run_m2.log; grep -v "rc=0" $O/run_m2.log
PYTHONIOENCODING=utf-8 python "<WORKDIR>/GEANT4/scripts/merge_templates_gs2020.py" "$(cygpath -w $C)" "$(cygpath -w $O/m2_merged)" > $O/merge_m2.log 2>&1; echo "merge rc=$?"
echo ГОТОВО >> $O/run_m2.log
