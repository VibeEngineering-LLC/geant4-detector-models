#!/usr/bin/env bash
# #GS-32/A6: один прогон фотонов из табличного спектра (IB по KUB), GS2020_PRIMARY_TABLE. Поля строки jobs: GDML TABLE N OUT SEED.
G="$1"; TAB="$2"; N="$3"; OUT="$4"; SD="$5"
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
GS2020_PRIMARY_TABLE="$TAB" GS2020_GDML="$G" "C:/g4work/build/gs2020-marinelli-beta/gs2020_marinelli.exe" 0 "$N" "$OUT" "$SD" > "${OUT%.csv}.log" 2>&1
echo "$OUT rc=$?"
