#!/usr/bin/env bash
# #GS-24 KCl: подгонка K-40 в конфигурации тория (gs2020_fit_env.sh) + вариант с донорским хвостом ядра (чувствительность, W-155).
# Выход: <GS_KDIR, по умолчанию C:\g4work\gs2020\kcl_1l_v4w85_83>\fit_kcl_bgw.json (основной) и fit_kcl_bgw_tail.json (вариант). Запуск: bash run_gs2020_kcl.sh
S="<WORKDIR>/GEANT4/scripts"; K="${GS_KDIR:-/c/g4work/gs2020/kcl_1l_v4w85_83}"   # GS_KDIR: папка шаблонов (v4 — kcl_1l_v4w85_83, конус колодца)
export GS_OUT="$(cygpath -w "$K")" GS_BG_WATER=1; . "$S/gs2020_fit_env.sh"
export GS_CAL_SHAPE="${GS_CAL_SHAPE:-1}"
# #GS-61 (оператор 09.10 «давай попробуем»): своя ширина пика 1460 из спектра KCl (scripts/gs61_own_fwhm.py kcl); GS_FWHM_OWN= — прежнее
# #GS-63 (оператор 09.10 «калибровку из тория надо брать. там все пики сходятся»): по умолчанию свёртка KCl — ширина тория;
# GS_FWHM_OWN=<csv> — своя ширина (#GS-61, только явно)
export GS_FWHM_OWN="${GS_FWHM_OWN-}"   # 30.09: опубликованные подгонки шли через run_gs2020_extra_accept.sh:9 с GS_CAL_SHAPE=1; отдельный запуск молча брал cal_own
[ "${GS_DRY:-0}" = 1 ] && { echo "K=$K GS_OUT=$GS_OUT"; exit 0; }
( unset GS_TAIL GS_FWHM_SCALE; python "$S/fit_gs2020_kcl.py" > "$K/fit_kcl_tail.log" 2>&1 ) || { echo "ОТКАЗ: вариант rc=$?"; exit 1; }
cp "$K/fit_kcl_bgw.json" "$K/fit_kcl_bgw_tail.json"
python "$S/fit_gs2020_kcl.py" > "$K/fit_kcl_main.log" 2>&1 || { echo "ОТКАЗ: основной rc=$?"; exit 1; }
python "$S/fit_gs2020_kcl_m2.py" > "$K/fit_kcl_m2.log" 2>&1 || { echo "ОТКАЗ: метод 2 rc=$?"; exit 1; }   # #GS-24: метод 2, та же конфигурация
for f in main tail m2; do echo "== $f"; grep -E "ФОРМА|K-40|Чувствительность|узлов" "$K/fit_kcl_$f.log"; done
