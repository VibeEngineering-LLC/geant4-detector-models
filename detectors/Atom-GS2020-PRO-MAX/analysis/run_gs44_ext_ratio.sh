#!/usr/bin/env bash
# #GS-44 (29.09): r(E) внешнего фона KCl/вода (gs44_ext_ratio.py) в окружении подгонки KCl (gs2020_fit_env.sh, как run_gs2020_kcl.sh).
# GS_EXT_SELFTEST=1 — KCl подменён водой (r обязано быть 1). Выход: $GS_EXT_DIR (по умолчанию C:\g4work\gs2020\ext).
S="<WORKDIR>/GEANT4/scripts"; K="${GS_KDIR:-/c/g4work/gs2020/kcl_1l_v4w85_83}"
export GS_OUT="$(cygpath -w "$K")" GS_BG_WATER=1; . "$S/gs2020_fit_env.sh"
export GS_CAL_SHAPE=1 GS_EXTRA=1   # как в run_gs2020_extra_accept.sh:9 (опубликованные подгонки KCl): без GS_CAL_SHAPE шкала cal_own ≠ cal_shape → сетка e не совпадёт с fit_kcl_bgw.json
PYTHONIOENCODING=utf-8 python "$S/gs44_ext_ratio.py"
