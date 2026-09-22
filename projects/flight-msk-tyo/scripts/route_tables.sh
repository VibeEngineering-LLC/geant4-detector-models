#!/bin/bash
# 77 tables PARMA: 11 route points x 7 components (ip 0,1,29,30,31,32,33), W=75, d=231.4654, g=10. Run from project root.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EXE="$ROOT/build/parma/parma_tables.exe"
OUT="$ROOT/tables/route"
mkdir -p "$OUT"
RC=(2.195 1.710 1.413 1.376 1.561 2.075 2.929 4.294 5.976 9.217 11.235)
fail=0
cd "$ROOT/parma" || exit 3
for i in "${!RC[@]}"; do
  p=$(printf "%02d" $((i+1)))
  for ip in 0 1 29 30 31 32 33; do
    "$EXE" $ip 75 ${RC[$i]} 231.4654 10 "$OUT/w75_p${p}_ip${ip}.tab" > "$OUT/w75_p${p}_ip${ip}.log" 2> "$OUT/w75_p${p}_ip${ip}.err"
    rc=$?
    if [ $rc -ne 0 ]; then echo "FAIL p=$p ip=$ip rc=$rc"; fail=$((fail+1)); fi
  done
done
echo "route tables done, failures=$fail, files=$(ls "$OUT"/*.tab 2>/dev/null | wc -l)"
exit $fail
