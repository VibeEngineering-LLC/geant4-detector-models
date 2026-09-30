#!/usr/bin/env bash
# #GS-49 (30.09): учёт состава распада (DECAY_TALLY) по звеньям цепочки Th-232. Запуск: bash run_gs49_tally.sh [N=1000000]
# Выход: C:\g4work\gs2020\gs49\tally_<нуклид>.txt + лог. Exe: отдельная сборка gs2020-marinelli-tally (прежние не тронуты).
N="${1:-1000000}"; D="/c/g4work/gs2020/gs49"; mkdir -p "$D"
export PATH="/c/geant4/bin:/c/g4work/thirdparty/xerces-install/bin:$PATH"
i=0
for job in Ra228:88:228 Ac228:89:228 Th228:90:228 Ra224:88:224 Rn220:86:220 Pb212:82:212 Bi212:83:212 Tl208:81:208; do
  nm="${job%%:*}"; za="${job#*:}"; i=$((i+1)); SD=$((5000+i))
  GS2020_GDML="C:/g4work/gs2020/GS2020_marinelli_th232.gdml" GS2020_DECAY_TALLY="C:/g4work/gs2020/gs49/tally_${nm}.txt" \
    "C:/g4work/build/gs2020-marinelli-tally/gs2020_marinelli.exe" "ion:${za}" "$N" "C:/g4work/gs2020/gs49/${nm}.csv" "$SD" > "$D/${nm}.log" 2>&1
  echo "$nm N=$N seed=$SD rc=$? $(grep -h 'DECAY_TALLY:' "$D/${nm}.log" | tail -1)"
done
