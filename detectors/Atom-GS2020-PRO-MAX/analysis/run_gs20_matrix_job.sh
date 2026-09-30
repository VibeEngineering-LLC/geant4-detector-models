#!/usr/bin/env bash
# #GS-20: один прогон варианта матрицы пробы. Поля строки jobs: GDML E N OUT SEED (xargs -L 1).
G="$1"; E="$2"; N="$3"; OUT="$4"; SD="$5"
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
GS2020_GDML="$G" "C:/g4work/build/gs2020-marinelli-cfg1/gs2020_marinelli.exe" "$E" "$N" "$OUT" "$SD" > "${OUT%.csv}.log" 2>&1
echo "$OUT rc=$?"
