#!/bin/bash
# 7 tables PARMA for DAD-SGN cruise: W=76 (Oulu 20.09.2026: 75.8), Rc=17.284 GV (midpoint; DAD 17.073, SGN 17.447), d=278.719 (11 km), g=10.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; OUT="$ROOT/tables/dadsgn"; mkdir -p "$OUT"; cd "$ROOT/parma" || exit 3
fail=0
for ip in 0 1 29 30 31 32 33; do
  "$ROOT/build/parma/parma_tables.exe" $ip 76 17.284 278.719 10 "$OUT/w76_h98_ip$ip.tab" > "$OUT/w76_h98_ip$ip.log" 2> "$OUT/w76_h98_ip$ip.err" || { echo "FAIL ip=$ip"; fail=$((fail+1)); }
  sed -n 2p "$OUT/w76_h98_ip$ip.tab" | awk -v ip=$ip '{print "ip="ip" total_flux="$6}'
done
exit $fail
