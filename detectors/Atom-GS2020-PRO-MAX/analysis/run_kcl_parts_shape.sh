#!/usr/bin/env bash
# Дрейф KCl 88,3 ч: калибровка по форме всеми четырьмя тегами отдельно для части 1 (41,5 ч) и части 2 (88,3 − 41,5 ч).
# Рабочий cal_shape.json (88,3 ч) сохранён в cal_shape_kcl88h_2026-10-09.json и восстанавливается в конце.
S="<WORKDIR>/GEANT4/scripts"; O=/c/g4work/gs2020/run_marinelli/out_v5; D=/c/g4work/gs2020/kcl_parts
R="<DOSIM>\Спектры\Atom GS2020 PRO MAX\Референсы"
declare -A X=([p1]="$R\KCl ч 1 л 1085 г (41,5 ч, 05-07.10).xml" [p2]="<WORKDIR>\GEANT4\results\kcl_parts\KCl ч 1 л 1085 г часть 2 (88,3 − 41,5 ч).xml")
cd "$S" && . ./gs2020_fit_env.sh && unset GS_CAL_SHAPE
for p in p1 p2; do
  GS2020_KCL_XML="${X[$p]}" PYTHONIOENCODING=utf-8 nice python cal/gs2020_calib_shape.py > "$D/cal_shape_$p.log" 2>&1; r=$?
  echo "$p rc=$r"; [ $r = 0 ] && cp "$O/cal_shape.json" "$D/cal_shape_$p.json"
done
cp "$O/cal_shape_kcl88h_2026-10-09.json" "$O/cal_shape.json" && echo "рабочий cal_shape (88,3 ч) восстановлен"
