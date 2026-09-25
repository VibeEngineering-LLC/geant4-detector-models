#!/usr/bin/env bash
# Один прогон сетки для #XR-1: xargs -L 1 передаёт 4 поля строки jobs_xray.txt как $1..$4.
E="$1"; N="$2"; OUT="$3"; SD="$4"
EXE="${GS2020_EXE:?Переменная окружения GS2020_EXE не установлена (см. README.md)}"
"$EXE" "$E" "$N" "$OUT" "$SD" > "${OUT%.csv}.log" 2>&1
echo "$OUT rc=$?"
