"""Влияние угловых корреляций каскада: прогон с флагом против прогона без него.

ЗАЧЕМ. `G4RadioactiveDecay` по умолчанию разыгрывает каждый квант каскада
независимо: `G4DeexPrecoParameters` задаёт `fCorrelatedGamma = false`.
Включается это `/process/had/deex/correlatedGamma true`. Пока флаг выключен,
угловые корреляции в расчёте не участвуют вообще, и списать на них расхождение
модели с измерением нельзя ни в какую сторону — списывать не на что.

ПОЧЕМУ НЕ НА РАБОЧЕЙ ГЕОМЕТРИИ 5 см. Оценка сделана ДО прогона, а не после.
Собственный TCC модели на 5 см равен 0,970 ± 0,014, то есть суммирование
составляет 3,0 %. Если корреляции усиливают совпадения на 7,7 % (величина по
W(θ) для каскада Co-60 при характерном угле 25°), суммирование станет 3,23 %,
а TCC сместится на 0,23 процентного пункта. Чтобы различить такой сдвиг на
трёх сигмах, нужна статистика в 330 раз больше нынешней — около 130 млн
распадов; по сумм-пику 2505 кэВ, где сейчас 18 отсчётов на 400 тыс., — около
34 млн. Это десятки часов счёта ради одного числа.

Поэтому флаг проверяется там, где эффект велик: источник вплотную к торцу.
Телесный угол вырастает примерно в шесть раз, суммирование — вместе с ним,
и эффект корреляций становится различим на обычной статистике. Проверяется
при этом ровно то, что нужно: работает ли флаг и какого порядка вносимая им
поправка. Перенос на 5 см — через известную зависимость суммирования от
телесного угла, а не прямым сравнением.

ЛОВУШКА, РАДИ КОТОРОЙ ЗАВЕДЕНА ПЕРЕМЕННАЯ ОКРУЖЕНИЯ. Из макроса флаг включить
нельзя: это параметр деэксцитации, принимается только до инициализации, а
макрос выполняется после неё. Команда отвергается с «Illegal application
state», НО КОД ВОЗВРАТА ОСТАЁТСЯ НУЛЕВЫМ. Прогон выглядит успешным, корреляции
молча не включаются, и сравнение «с флагом против без» показывает отсутствие
эффекта по причине, не имеющей отношения к физике. Поэтому main.cc читает
G1S_CORRELATED_GAMMA до `Initialize()`, а этот драйвер проверяет, что строка
подтверждения включения в выводе действительно есть, и отказывается считать,
если её нет.
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
EXE = os.path.join(BUILD, "g1s.exe")

if not os.path.exists(EXE):
    raise SystemExit("Не найдена собранная модель %s" % EXE)

ZFACE = 41.0
ZSRC = ZFACE + 1.0      # мм: вплотную к крышке — максимум телесного угла
NDECAY = int(os.environ.get("G1S_CORR_NDECAY", "400000"))
MARK = "SETUP correlatedGamma = true"


def macro(out):
    t = ["/control/verbose 0", "/run/verbose 0", "/process/had/rdm/verbose 0",
         "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns",
         "/process/had/rdm/nucleusLimits 60 60 27 28",
         "/gps/particle ion", "/gps/energy 0 keV",
         "/gps/pos/type Point", "/gps/pos/centre 0 0 %.1f mm" % ZSRC,
         "/gps/ang/type iso", "/gps/ion 27 60 0 0",
         "/g1s/outFile %s" % os.path.join(BUILD, out),
         "/run/beamOn %d" % NDECAY]
    return "\n".join(t) + "\n"


def run(out, correlated):
    mp = os.path.join(BUILD, "corr_%s.mac" % ("on" if correlated else "off"))
    open(mp, "w", encoding="utf-8").write(macro(out))
    env = dict(os.environ)
    if correlated:
        env["G1S_CORRELATED_GAMMA"] = "1"
    else:
        env.pop("G1S_CORRELATED_GAMMA", None)

    print("=== корреляции %s, %d распадов Co-60 на %.0f мм от торца ==="
          % ("ВКЛ" if correlated else "выкл", NDECAY, ZSRC - ZFACE), flush=True)
    r = subprocess.run([EXE, mp, "shield"], cwd=BUILD, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", env=env)
    out_txt = r.stdout or ""
    for ln in out_txt.splitlines():
        s = ln.strip()
        if s.startswith(("RESULT", "SETUP")) or "Illegal" in s \
                or "WARNING" in s or "ERROR" in s:
            print("  ", s, flush=True)
    if r.returncode != 0:
        print("!! код возврата", r.returncode)
        print((r.stderr or "")[-1500:])
        sys.exit(1)

    # Главная защита: без этой строки прогон с «включёнными» корреляциями на
    # деле идёт без них, и сравнение молча теряет смысл.
    if correlated and MARK not in out_txt:
        raise SystemExit(
            "Корреляции НЕ включились: в выводе нет строки «%s».\n"
            "Проверьте, что main.cc читает G1S_CORRELATED_GAMMA ДО"
            " Initialize() и что запущен свежий бинарь." % MARK)
    if not correlated and MARK in out_txt:
        raise SystemExit(
            "Корреляции включились в контрольном прогоне, где не должны:"
            " переменная окружения протекла.")


if __name__ == "__main__":
    run("corr_off.csv", correlated=False)
    run("corr_on.csv", correlated=True)
    print("готово: corr_off.csv и corr_on.csv в %s" % BUILD)
