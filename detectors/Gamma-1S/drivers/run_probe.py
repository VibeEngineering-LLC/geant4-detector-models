"""Два целевых прогона по итогам сверки с ЛСРМ.

ПРОГОН 1 — ГЛУБИНА КОЛОДЦА МАРИНЕЛЛИ. Точечная геометрия дала МК/эксп = 0,971,
то есть модель детектора верна на 3 %, а маринелли — 1,165. Значит расхождение
в геометрии объёмного источника. Из документов у сосуда известны только Ø150 и
H = 110; стенка 2 мм, колодец Ø80 и глубина 74 мм — допущения. Мельче колодец —
проба дальше от кристалла — эффективность ниже. Вопрос: существует ли
правдоподобная глубина, при которой отношение садится на единицу.

ПРОГОН 2 — ПЛОТНОСТЬ ОТРАЖАТЕЛЯ MgO. Мягкий край точечной кривой занижен
(МК/эксп 0,78 на 59,5 кэВ при ~0,95 в середине). Принятая насыпная плотность
2,0 г/см³ — допущение из диапазона 1,5–2,4. Проверяем 1,3 / 1,5 / 2,0 на
точечной геометрии 5 см. КРИТЕРИЙ: мягкий край должен выправиться, НЕ ломая
середину и жёсткий край — иначе это подгонка, а не уточнение.

Оба прогона короткие: по нескольким энергиям, а не всей сеткой.
"""
import math
import os
import subprocess
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402


BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)

EXE = "g1s.exe"

if not os.path.exists(os.path.join(BUILD, EXE)):
    raise SystemExit(
        "Не найдена собранная модель %s.\n"
        "Соберите её (см. common/cmake и README детектора) или укажите\n"
        "G4MODELS_BUILD_GAMMA_1S на каталог, где она уже лежит."
        % os.path.join(BUILD, EXE))
OUT = os.path.join(BUILD, "probe2")
os.makedirs(OUT, exist_ok=True)

N = 400000
WELLS = [74.0, 65.0, 55.0, 45.0]          # мм
MGOS = [1.30, 1.50, 2.00]                 # г/см³
E_WELL = [661.657, 1460.822]              # где сверка маринелли надёжнее всего
E_MGO = [59.5, 88.0, 122.1, 661.657, 2614.511]
THETA = 60.0                              # конус точечной 5 см
FRAC = (1 - math.cos(math.radians(THETA))) / 2


# --- Область розыгрыша объёмного источника ----------------------------------
# /gps/pos/confine отбрасывает точки вне пробы, но НЕ добавляет точки, которых
# в области розыгрыша не было. Значит цилиндр обязан объемлеть пробу целиком,
# иначе часть её просто не облучается — молча, без единого предупреждения.
#
# Именно это здесь и происходило: цилиндр был записан константами (r=73,
# halfz=45, центр z=16) под ШТАТНУЮ глубину колодца 74 мм. Но весь смысл этого
# прогона — менять глубину. Чем мельче колодец, тем выше поднимается уровень
# засыпки, и тем большая часть пробы оказывается за верхней границей области:
# 2,7 % при 65 мм, 13,9 % при 55, 25,1 % при 45.
#
# Хуже того, теряется ВЕРХ пробы — самый дальний от кристалла и потому наименее
# эффективный. Выбрасывая его, мы завышаем эффективность мелкого колодца, то
# есть смещаем результат ПРОТИВ проверяемой гипотезы (мельче колодец — ниже
# эффективность). Скан нашёл бы оптимум там, где его нет.
#
# Поэтому габарит считается из тех же формул, что и G1SDetector::BuildVessel,
# и печатается в лог: расхождение с «объём пробы» из вывода модели сразу видно.
WALL, Z_FACE, R_IN, R_WELL_OUT, SAMPLE_CM3 = 2.0, 41.0, 73.0, 42.0, 1000.0


def sample_span(well):
    """(центр, полувысота) цилиндра, объемлющего пробу, мм. См. BuildVessel."""
    z_lo = Z_FACE - well + 2 * WALL              # низ пробы (над дном сосуда)
    ring = math.pi * (R_IN ** 2 - R_WELL_OUT ** 2) * (well - WALL) / 1000.0
    top_h = (SAMPLE_CM3 - ring) * 1000.0 / (math.pi * R_IN ** 2)
    z_hi = Z_FACE + WALL + top_h                 # уровень засыпки
    return 0.5 * (z_lo + z_hi), 0.5 * (z_hi - z_lo)


def mac_vol(tag, lines, well):
    zc, hz = sample_span(well)
    hz += 0.5                                    # запас на округления
    print("    колодец %.0f мм: розыгрыш z = %.2f +- %.2f мм"
          % (well, zc, hz), flush=True)
    t = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
         "/gps/particle gamma", "/gps/pos/type Volume",
         "/gps/pos/shape Cylinder", "/gps/pos/centre 0 0 %.2f mm" % zc,
         "/gps/pos/radius %.1f mm" % R_IN, "/gps/pos/halfz %.2f mm" % hz,
         "/gps/pos/confine Sample", "/gps/ang/type iso"]
    for e in lines:
        t += ["/gps/energy %.3f keV" % e,
              "/g1s/outFile %s" % os.path.join(OUT, "%s_E%07.1f.csv" % (tag, e)),
              "/run/beamOn %d" % N]
    return "\n".join(t) + "\n"


def mac_point(tag, lines):
    t = ["/run/initialize", "/control/verbose 0", "/run/verbose 0",
         "/gps/particle gamma", "/gps/pos/type Point",
         # 50 мм от наружной плоскости торца (Z_FACE); прежнее «93» держало
         # старую плоскость z=43 (протектор), убранную правкой торца.
         "/gps/pos/centre 0 0 %.1f mm" % (Z_FACE + 50.0), "/gps/ang/type iso",
         "/gps/ang/maxtheta %.1f deg" % THETA]
    for e in lines:
        t += ["/gps/energy %.3f keV" % e,
              "/g1s/outFile %s" % os.path.join(OUT, "%s_E%07.1f.csv" % (tag, e)),
              "/run/beamOn %d" % N]
    return "\n".join(t) + "\n"


def run(mac, args, label):
    p = os.path.join(BUILD, "probe2.mac")
    open(p, "w", encoding="utf-8").write(mac)
    print("=== %s ===" % label, flush=True)
    r = subprocess.run([os.path.join(BUILD, "g1s.exe"), p] + args, cwd=BUILD,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace")
    for ln in (r.stdout or "").splitlines():
        if ln.startswith("RESULT") or "MgO" in ln or "проба" in ln:
            print("   ", ln.strip(), flush=True)
    if r.returncode != 0:
        print("!! код", r.returncode, (r.stderr or "")[-800:])
        sys.exit(1)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "well"):
        for w in WELLS:
            tag = "well%.0f" % w
            # маринелли, ОИСН-16 ро=1,6 — как в сверке с .efr
            run(mac_vol(tag, E_WELL, w),
                ["vessel:marinelli", "1.6", "OISN16", "1000", "2.0", str(w)],
                "колодец %.0f мм" % w)
    if what in ("all", "mgo"):
        for m in MGOS:
            tag = "mgo%.2f" % m
            run(mac_point(tag, E_MGO),
                ["shield", "1.6", "OISN16", "1000", str(m), "74"],
                "MgO %.2f г/см³, точечный 5 см" % m)
    print("доля телесного угла точечной: %.5f" % FRAC)
    print("готово")
