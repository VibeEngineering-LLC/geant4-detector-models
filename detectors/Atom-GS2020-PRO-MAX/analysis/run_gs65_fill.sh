#!/usr/bin/env bash
# #GS-65: дозаполнение кривой эффективности 16–68 кэВ (22 узла, #CFG-1 «да» 09.10). Ждёт конца сетки #GS-64, ≤20 потоков; родитель — BelowNormal.
cd "$(dirname "$0")"
until grep -q 'ГОТОВО' /c/g4work/gs2020/density_corr/run.log 2>/dev/null; do sleep 60; done
xargs -L 1 -P 20 bash run_gs20_matrix_job.sh < /c/g4work/gs2020/density_corr/jobs_fill.txt > /c/g4work/gs2020/density_corr/run_fill.log 2>&1
echo "ГОТОВО rc=$?" >> /c/g4work/gs2020/density_corr/run_fill.log
