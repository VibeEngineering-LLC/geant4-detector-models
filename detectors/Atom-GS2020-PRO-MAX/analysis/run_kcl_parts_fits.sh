#!/usr/bin/env bash
# Дрейф KCl 88,3 ч: своя ширина + М1/М2 + выгрузка K-40 по частям (p1 = 41,5 ч, p2 = 88,3 − 41,5 ч) на их калибровках по форме.
# Выход: /c/g4work/gs2020/kcl_parts/<p>/. В конце восстанавливает рабочее состояние 88,3 ч (cal_shape, выходы подгонок, данные страницы).
S="<WORKDIR>/GEANT4/scripts"; O=/c/g4work/gs2020/run_marinelli/out_v5; D=/c/g4work/gs2020/kcl_parts; K=/c/g4work/gs2020/kcl_1l_v4w85_83
P="<WORKDIR>/GEANT4/web/gs2020-th232-page"; R="<DOSIM>\Спектры\Atom GS2020 PRO MAX\Референсы"
declare -A X=([p1]="$R\KCl ч 1 л 1085 г (41,5 ч, 05-07.10).xml" [p2]="<WORKDIR>\GEANT4\results\kcl_parts\KCl ч 1 л 1085 г часть 2 (88,3 − 41,5 ч).xml")
cd "$S" && . ./gs2020_fit_env.sh
export GS_CAL_SHAPE=1 PYTHONIOENCODING=utf-8 GS_EXTRA=1 GS_BG_R='C:\g4work\gs2020\ext\r_ext_4pi_kcl1l.csv' SPECTRAVIBE_ROOT="<DOSIM>/ИИ/1 Скилы/0_Work/gamma-spectrum-analysis"
for p in p1 p2; do
  mkdir -p "$D/$p/x" && cp "$D/cal_shape_$p.json" "$O/cal_shape.json" && export GS2020_KCL_XML="${X[$p]}"
  ( cd "$D/$p/x" && python "$S/gs61_own_fwhm.py" kcl > "$D/$p/own_fwhm.log" 2>&1 ); echo "$p own rc=$?"
  export GS_FWHM_OWN="$(cygpath -w "$D/$p/results/gs61_own_fwhm/kcl.csv")"
  nice bash "$S/run_gs2020_kcl.sh" > "$D/$p/run_fits.out" 2>&1; echo "$p fits rc=$?"; cp -p $K/fit_kcl_*.json $K/fit_kcl_*.log $K/fwhm_points_gs2020.csv "$D/$p/"
  ( export GS_OUT="$(cygpath -w $K)" GS_BG_WATER=1; python "$S/export_page_gs2020_k40.py" > "$D/$p/export_k40.log" 2>&1 ); echo "$p exp rc=$?"; cp -p "$P/gs2020_k40_data.json" "$D/$p/"
done
cp "$O/cal_shape_kcl88h_2026-10-09.json" "$O/cal_shape.json"; cp -p $K/keep_kcl88h_final/* $K/; cp -p "$D/gs2020_k40_data_88h.json" "$P/gs2020_k40_data.json"; echo "восстановлено 88,3 ч"
