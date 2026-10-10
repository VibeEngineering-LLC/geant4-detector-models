#!/usr/bin/env bash
# #GS-68 стенд NPSM на монолиниях GS2020 (#CFG-1 «Да, по таблице» 10.10): регрессия → парный замер 10^5 → 5 линий × 2·10^7.
# Payne по умолчанию в exe; геометрия референса mx_oisn10; ≤20 потоков (5 заданий), родитель BelowNormal.
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
G=C:/g4work/gs2020/gs20_matrix/mx_oisn10.gdml; X=C:/g4work/build/gs2020-marinelli-npsm/gs2020_marinelli.exe
OLD=C:/g4work/build/gs2020-marinelli-cfg1/gs2020_marinelli.exe; O=/c/g4work/gs2020/npsm68; mkdir -p $O; cd $O
# 1) регрессия: новый exe без модели == старый (побитово; старый cfg1 собран до 4 строк шапки beta_only…table_drawn — их исключаем); с моделью столбец edep == без модели
GS2020_GDML=$G $OLD 1460.8 20000 reg_old.csv 81900 > reg_old.log 2>&1
GS2020_GDML=$G $X 1460.8 20000 reg_off.csv 81900 > reg_off.log 2>&1
GS2020_NPSM=1 GS2020_GDML=$G $X 1460.8 20000 reg_on.csv 81900 > reg_on.log 2>&1
diff --strip-trailing-cr <(cat reg_old.csv) <(grep -v -E "^(beta_only|killed_decay_photons|primary_table|table_drawn)," reg_off.csv) > /dev/null && echo "РЕГРЕССИЯ off==old: ДА" || { echo "РЕГРЕССИЯ off==old: НЕТ"; exit 1; }
d=$(diff --strip-trailing-cr <(grep -A9999 '^bin_keV' reg_off.csv | cut -d, -f1,2) <(grep -A9999 '^bin_keV' reg_on.csv | cut -d, -f1,2) | wc -l)
echo "edep on==off: строк расхождения $d"; [ "$d" = 0 ] || exit 1
# 2) парный замер стоимости, одно зерно
for m in 0 1; do s=$(date +%s); GS2020_NPSM=$m GS2020_GDML=$G $X 2614.5 100000 time_npsm$m.csv 81901 > time_npsm$m.log 2>&1; r=$?; echo "время npsm=$m: $(( $(date +%s) - s )) с rc=$r"; done
# 3) 5 линий × 2·10^7
n=0; for E in 238.6 583.2 911.2 1460.8 2614.5; do n=$((n+1)); echo "$G $E 20000000 $O/line_E${E}_npsm.csv 8100$n"; done > jobs.txt
xargs -L 1 -P 5 bash -c 'GS2020_NPSM=1 GS2020_GDML=$0 '"$X"' $1 $2 $3 $4 > ${3%.csv}.log 2>&1; echo "$3 rc=$?"' < jobs.txt
echo "ГОТОВО"
