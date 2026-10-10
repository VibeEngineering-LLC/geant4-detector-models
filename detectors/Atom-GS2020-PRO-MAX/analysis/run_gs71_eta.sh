#!/usr/bin/env bash
# #GS-71 этап 2 (#CFG-1 «Да, по таблице» 10.10): скан η NPSM (прочие — Payne) на монолиниях 661,657 и 1460,822 кэВ, 2·10^7, геометрия mx_oisn10;
# 5 заданий параллельно (как стенд #GS-68). Разбор: python gs71_npsm_edge.py <csv> <E0> ...
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
G=C:/g4work/gs2020/gs20_matrix/mx_oisn10.gdml; X=C:/g4work/build/gs2020-marinelli-npsm/gs2020_marinelli.exe; O=/c/g4work/gs2020/npsm71; mkdir -p $O
s=86000; : > $O/jobs.txt
for eta in 0.6 0.7 0.8 0.9 1.0; do for E in 661.657 1460.822; do s=$((s+1)); echo "$eta $E $O/line_E${E}_eta${eta}.csv $s" >> $O/jobs.txt; done; done
export G X; xargs -L 1 -P 5 bash -c 'GS2020_NPSM=1 GS2020_NPSM_ETA=$0 GS2020_GDML=$G $X $1 20000000 $2 $3 > ${2%.csv}.log 2>&1; echo "$2 rc=$?"' < $O/jobs.txt > $O/run.log 2>&1
echo "ГОТОВО" >> $O/run.log
