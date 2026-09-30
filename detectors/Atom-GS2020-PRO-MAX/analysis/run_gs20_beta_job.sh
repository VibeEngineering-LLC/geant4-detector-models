#!/usr/bin/env bash
# #GS-32: один прогон β/e⁻-компоненты (GS2020_BETA_ONLY=1, сборка gs2020-marinelli-beta). Поля строки jobs: GDML E N OUT SEED.
G="$1"; E="$2"; N="$3"; OUT="$4"; SD="$5"
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
GS2020_BETA_ONLY=1 GS2020_GDML="$G" "C:/g4work/build/gs2020-marinelli-beta/gs2020_marinelli.exe" "$E" "$N" "$OUT" "$SD" > "${OUT%.csv}.log" 2>&1
echo "$OUT rc=$?"
