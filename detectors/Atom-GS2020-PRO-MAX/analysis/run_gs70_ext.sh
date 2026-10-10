#!/usr/bin/env bash
# #GS-70 внешнее поле для r(E) черники (#CFG-1 «Да, по таблице» 10.10): 21 энергия × up/down × 5·10^7, как #GS-60 (вода — готовые #GS-44).
# Ждёт конца шаблонов (run_gs70_templates.sh), затем NP=12. Выход: C:\g4work\gs2020\ext\chunks\ext_berry_E<E>_<h>.csv
S="<WORKDIR>/GEANT4/scripts"; O=/c/g4work/gs2020/ext; G=C:/g4work/gs2020/berry/GS2020_marinelli_berry_dry_1000ml_v4_w85_83.gdml
J=$O/jobs_ext_gs70.txt; : > $J; s=84001
for e in 30 40 50 60 80 100 120 150 200 250 300 400 500 609.312 800 1000 1250 1460.822 1764.494 2000 2614.511; do
  for h in up down; do echo "$G $e 50000000 C:/g4work/gs2020/ext/chunks/ext_berry_E${e}_${h}.csv $s $h" >> $J; s=$((s+1)); done; done
until grep -q "ГОТОВО" /c/g4work/gs2020/berry/run_templates.log 2>/dev/null; do sleep 60; done
cd "$S" && xargs -P ${NP:-12} -L 1 bash run_gs20_ext_job.sh < $J > $O/run_gs70_ext.log 2>&1; echo "ГОТОВО" >> $O/run_gs70_ext.log
