#!/usr/bin/env bash
# #GS-44 (29.09): one external-field run (ext_source.hh). Job fields: GDML E N OUT SEED HEMI (xargs -L 1).
# Sphere radius/centre: env GS_EXT_R_MM (default 240), GS_EXT_CZ_MM (default 0). Exe: separate ext build (cfg1 untouched).
G="$1"; E="$2"; N="$3"; OUT="$4"; SD="$5"; H="$6"
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
GS2020_GDML="$G" GS2020_EXT_R_MM="${GS_EXT_R_MM:-240}" GS2020_EXT_CZ_MM="${GS_EXT_CZ_MM:-0}" GS2020_EXT_HEMI="$H" \
  "C:/g4work/build/gs2020-marinelli-ext/gs2020_marinelli.exe" "$E" "$N" "$OUT" "$SD" > "${OUT%.csv}.log" 2>&1
echo "$OUT rc=$?"
