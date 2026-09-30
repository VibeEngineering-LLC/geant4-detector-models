#!/usr/bin/env bash
# Один прогон сетки для #XR-1: xargs -L 1 передаёт 4 поля строки jobs_xray.txt как $1..$4.
E="$1"; N="$2"; OUT="$3"; SD="$4"
"C:/g4work/build/gs2020-marinelli/gs2020_marinelli.exe" "$E" "$N" "$OUT" "$SD" > "${OUT%.csv}.log" 2>&1
echo "$OUT rc=$?"
