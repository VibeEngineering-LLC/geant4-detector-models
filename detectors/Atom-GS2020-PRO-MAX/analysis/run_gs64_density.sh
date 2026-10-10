#!/usr/bin/env bash
# #GS-64: сетка поправки на плотность/уровень (#CFG-1 «да» 09.10). 324 задания, ≤20 потоков; запускать из процесса с BelowNormal.
cd "$(dirname "$0")"
sleep 5   # время выставить приоритет родителю до первого дочернего процесса
xargs -L 1 -P 20 bash run_gs20_matrix_job.sh < /c/g4work/gs2020/density_corr/jobs.txt > /c/g4work/gs2020/density_corr/run.log 2>&1
echo "ГОТОВО rc=$?" >> /c/g4work/gs2020/density_corr/run.log
