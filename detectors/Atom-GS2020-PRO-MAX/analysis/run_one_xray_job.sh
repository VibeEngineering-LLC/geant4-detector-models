#!/usr/bin/env bash
# Один прогон сетки для #XR-1: xargs -L 1 передаёт 4 поля строки jobs_xray.txt как $1..$4.
E="$1"; N="$2"; OUT="$3"; SD="$4"
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"   # 07.10: без него rc=127 (xerces-c_3_3.dll), как в run_gs20_matrix_job.sh
"C:/g4work/build/gs2020-marinelli/gs2020_marinelli.exe" "$E" "$N" "$OUT" "$SD" > "${OUT%.csv}.log" 2>&1
echo "$OUT rc=$?"
