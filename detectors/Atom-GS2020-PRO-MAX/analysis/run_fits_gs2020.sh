#!/usr/bin/env bash
# Шесть подгонок GS2020 Th-232: {фон как снят, фон по реальному замеру Маринелли+вода} x {М1, М2 библ. 0,5 %, М2 библ. 2 %}.
# Логи — C:\g4work\gs2020\run_marinelli\out_v5\fit_*.log; печатает код возврата каждой и итоговые строки цепочки.
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
OUT="${GS2020_OUT:?Переменная окружения GS2020_OUT не установлена (см. README.md)}"
BGT_ATT="${BGT_ATT:-0.6,6.02}"
for T in "" "$BGT_ATT"; do
  sfx=""; [ -n "$T" ] && sfx="_bgT${T/,/_}"
  GS_BG_T="$T" python fit_gs2020_th232_m1.py > "$OUT/fit_m1$sfx.log" 2>&1; echo "m1$sfx $?"
  GS_BG_T="$T" GS_M2_CONFIG=configs/th232_gs2020_lib05.yaml python fit_gs2020_th232_m2.py > "$OUT/fit_m2_lib05$sfx.log" 2>&1; echo "m2l05$sfx $?"
  GS_BG_T="$T" python fit_gs2020_th232_m2.py > "$OUT/fit_m2$sfx.log" 2>&1; echo "m2$sfx $?"
done
grep -h "#CAL-1 шкала фона\|ЦЕПОЧКА В РАВНОВЕСИИ" "$OUT"/fit_m1.log "$OUT"/fit_m2_lib05.log "$OUT"/fit_m2.log \
  "$OUT"/fit_m1_bgT*"${BGT_ATT/,/_}".log "$OUT"/fit_m2_lib05_bgT*"${BGT_ATT/,/_}".log "$OUT"/fit_m2_bgT*"${BGT_ATT/,/_}".log
