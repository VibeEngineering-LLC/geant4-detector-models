#!/usr/bin/env bash
# #GS-70 дым шаблонов черники (#CFG-1 «Да, по таблице» 10.10): по 10^5 распадов Cs-137 (до Ba-137m), K-40, Sr-90, Y-90; время и пики.
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
G=C:/g4work/gs2020/berry/GS2020_marinelli_berry_dry_1000ml_v4_w85_83.gdml; X=C:/g4work/build/gs2020-marinelli-npsm/gs2020_marinelli.exe
O=/c/g4work/gs2020/berry/smoke; mkdir -p $O; cd $O
for j in "Cs137 55:137 56" "K40 19:40 19" "Sr90 38:90 38" "Y90 39:90 39"; do set -- $j
  s=$(date +%s); GS2020_ZMAX=$3 GS2020_GDML=$G $X ion:$2 100000 smoke_$1.csv 82000 > smoke_$1.log 2>&1; r=$?
  echo "$1: rc=$r, $(( $(date +%s) - s )) с, $(grep CHAIN_CUT smoke_$1.log)"
  awk -F, -v n=$1 'f{t+=$2; if($1>=640&&$1<=680)a+=$2; if($1>=28&&$1<=38)x+=$2; if($1>=1440&&$1<=1480)k+=$2; if($1>1000)h+=$2} /^bin_keV/{f=1}
    END{printf "  %s: всего %d; 640-680 %d; 28-38 %d; 1440-1480 %d; >1000 %d\n", n, t, a, x, k, h}' smoke_$1.csv
done
