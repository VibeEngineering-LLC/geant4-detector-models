#!/bin/bash
# Независимый пересчёт T_sim = N / (flux * pi * R^2) из таблиц PARMA и сверка с .meta этапа I. Код возврата 1, если расхождение > 1e-4 (отн.).
# Использование: audit_tsim.sh <каталог stage1> <префикс таблиц> ; R = 480 см.
D="${1:-results/stage1/prod1}"; P="${2:-tables/dadsgn/w76_h98_ip}"; bad=0
declare -A IP=([neutron]=0 [proton]=1 [mup]=29 [mum]=30 [em]=31 [ep]=32 [gamma]=33)
for m in "$D"/*.psp.meta; do
  comp=$(basename "$m" | sed 's/_[0-9]*\.psp\.meta//')
  if [ "$comp" = k40 ]; then   # K-40: T_sim = N / активность (потока PARMA нет)
    n=$(grep '^N=' "$m" | cut -d= -f2); t=$(grep '^T_sim_s=' "$m" | cut -d= -f2); a=$(grep '^k40_activity_Bq=' "$m" | cut -d= -f2)
    res=$(awk -v n="$n" -v a="$a" -v t="$t" 'BEGIN{e=n/a; d=(e-t)/e; if(d<0)d=-d; printf "%.6e %.6e rel=%.2e %s", e, t, d, (d>1e-4?"FAIL":"OK")}')
    echo "$(basename "$m") act=$a N=$n expected/meta: $res"; case "$res" in *FAIL) bad=1;; esac; continue
  fi
  flux=$(sed -n 2p "${P}${IP[$comp]}.tab" | awk '{print $6}')
  n=$(grep '^N=' "$m" | cut -d= -f2); t=$(grep '^T_sim_s=' "$m" | cut -d= -f2)
  res=$(awk -v n="$n" -v f="$flux" -v t="$t" 'BEGIN{e=n/(f*3.14159265358979*480*480); d=(e-t)/e; if(d<0)d=-d; printf "%.6e %.6e rel=%.2e %s", e, t, d, (d>1e-4?"FAIL":"OK")}')
  echo "$(basename "$m") flux=$flux N=$n expected/meta: $res"; case "$res" in *FAIL) bad=1;; esac
done
exit $bad
