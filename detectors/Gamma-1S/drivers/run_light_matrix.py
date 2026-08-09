"""Прямой прогон распада Th-232 в маринелли с лёгкой матрицей ОИСН-06.

Задача 107, проверка метода. Активность партии 420-17031 (ро = 0,64,
ОИСН-06) считалась через eps штатного прогона (ро = 1,6, ОИСН-16) с
пересчётом самопоглощения f(мю*ро*d_eff), где лёгкую среду представляла
вода. Пересчёт — приближение сразу в трёх местах: экспоненциальная формула
вместо переноса, эффективная толщина d_eff вместо реальной геометрии, вода
вместо состава ОИСН-06 (15 % железа). Здесь то же самое считается прямым
транспортом: тот же макрос цепочки, но засыпка ро = 0,64 из ОИСН-06.

Сравнение eps_прямая/eps_пересчитанная по линиям — цена f-приближения.
Если она мала, батч-разность маринелли (1,32 сигма, сама по себе незначимая)
поправкой не объясняется и вопрос закрыт; если велика — весь пересчёт между
плотностями в kit_recalc/second_source требует замены на прямые прогоны.

Выход: chain_Th232_oisn06.csv (+ _emit) в BUILD. Статистика 600 000 распадов
даёт ~3 % по площади 2614,5 — достаточно против ожидаемых процентов эффекта.
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gps_region  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
EXE = os.path.join(BUILD, "g1s.exe")
if not os.path.exists(EXE):
    raise SystemExit("нет g1s.exe в %s" % BUILD)

ARGS = ["vessel", "0.64", "OISN06", "1000"]

# Тело розыгрыша подставляется из выгрузки построенной геометрии: записанное
# константами (73/45/16) оно перестало покрывать пробу, когда сосуд перешёл на
# чертёж изготовителя, и теряло наружное кольцо и верх пробы (R68, R75).
MAC = """/run/initialize
/control/verbose 0
/run/verbose 0
/process/had/rdm/verbose 0
/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns
/gps/particle ion
/gps/energy 0 keV
%s
/gps/ang/type iso
/process/had/rdm/nucleusLimits 208 232 81 90
/gps/ion 90 232 0 0
/g1s/outFile %s
/run/beamOn 600000
"""

if __name__ == "__main__":
    mp = os.path.join(BUILD, "light_matrix.mac")
    open(mp, "w", encoding="utf-8").write(
        MAC % ("\n".join(gps_region.gps_lines(BUILD, ARGS[0], ARGS[1:])),
               os.path.join(BUILD, "chain_Th232_oisn06.csv")))
    print("=== Th-232 цепочка, маринелли, ОИСН-06 ро=0,64, 600k ===",
          flush=True)
    r = subprocess.run([EXE, mp, "vessel", "0.64", "OISN06", "1000"],
                       cwd=BUILD, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    for ln in (r.stdout or "").splitlines():
        if ln.startswith(("RESULT", "EMIT", "SETUP")) or "проба" in ln:
            print("  ", ln.strip(), flush=True)
    if r.returncode != 0:
        print("!! код возврата", r.returncode)
        print((r.stderr or "")[-1500:])
        sys.exit(1)
    print("готово", flush=True)
