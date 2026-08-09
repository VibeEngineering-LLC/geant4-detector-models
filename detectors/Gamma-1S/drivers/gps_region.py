"""Область розыгрыша объёмного источника — по ВЫГРУЗКЕ построенной геометрии.

ЗАЧЕМ ОТДЕЛЬНЫМ МОДУЛЕМ. `/gps/pos/confine` отбрасывает точки вне тома, но НЕ
добавляет точек, которых в теле розыгрыша не было. Значит тело обязано
объемлеть пробу целиком, иначе часть её не облучается — молча, без единого
предупреждения, и теряется как раз дальняя от кристалла часть, отчего
эффективность выходит завышенной.

Числа этого тела жили в каждом драйвере своими константами (r=73, halfz=45,
центр z=16 — под сосуд, габариты которого брались из таблицы ЛСРМ). Когда
сосуд Маринелли переставили на чертёж изготовителя (R68: наружный Ø154,
колодец Ø97 глубиной 65), проба выросла до r = 75 и до z = +65,2, а константы
драйверов остались прежними: терялось наружное кольцо 2 мм и верхние 4 мм
(R75). Дублирование размеров модели в скриптах и есть причина — здесь оно
снимается: тело считается из того, что Geant4 ПОСТРОИЛ.

Как: exe запускается с пустым макросом и переменной окружения G1S_DUMP_GEOM,
которая заставляет его выгрузить построенное дерево (имя, материал, r и z
каждого тела). Из выгрузки берутся строки тома `Sample`. Прогон дешёвый:
инициализация без единого события.

Сторож той же проверки стоит и В МОДЕЛИ (main.cc, CheckConfinement) и роняет
прогон, если тело розыгрыша не покрывает том. Здесь — чтобы драйвер строил
правильно с первого раза, там — чтобы неправильное не посчиталось.

    from gps_region import sample_region
    r, zc, hz = sample_region(BUILD, "vessel:marinelli", ["1.6", "OISN16"])
"""
import os
import subprocess
import tempfile

# Запас к габариту пробы, мм. Нужен на округления при печати команд макроса
# (%.1f) и на то, что уровень засыпки считается из объёма и «дышит» с
# плотностью и матрицей.
MARGIN = 3.0

_CACHE = {}

NOP = ("# Пустой макрос: задаёт режим позиционным аргументом.\n"
       "# Выгрузка геометрии происходит после Initialize(), до макроса.\n"
       "/control/verbose 0\n/run/verbose 0\n")


def sample_bounds(build, mode, args=(), exe="g1s.exe"):
    """(rmax, zlo, zhi) тома пробы в мм для данного режима, по выгрузке.

    Возвращает None, если тома `Sample` в этом режиме нет (bare/open/shield).
    """
    key = (build, mode, tuple(args))
    if key in _CACHE:
        return _CACHE[key]

    tmp = tempfile.mkdtemp(prefix="g1sgeom")
    mac = os.path.join(tmp, "nop.mac")
    dump = os.path.join(tmp, "geom.csv")
    with open(mac, "w", encoding="utf-8") as fh:
        fh.write(NOP)
    env = dict(os.environ, G1S_DUMP_GEOM=dump)
    r = subprocess.run([os.path.join(build, exe), mac, mode] + list(args),
                       cwd=build, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=env)
    if r.returncode != 0 or not os.path.exists(dump):
        raise SystemExit(
            "не удалось снять габарит пробы: %s %s -> код %d\n%s"
            % (exe, mode, r.returncode, (r.stderr or r.stdout or "")[-800:]))

    rmax = zlo = zhi = None
    with open(dump, encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("Sample,"):
                continue
            p = ln.rstrip("\n").split(",")
            if "?" in p[2:6]:
                raise SystemExit("габарит пробы не снят: тело неизвестного "
                                 "класса в выгрузке (%s)" % ln.strip())
            ro, z0, z1 = float(p[3]), float(p[4]), float(p[5])
            rmax = ro if rmax is None else max(rmax, ro)
            zlo = z0 if zlo is None else min(zlo, z0)
            zhi = z1 if zhi is None else max(zhi, z1)

    out = None if rmax is None else (rmax, zlo, zhi)
    _CACHE[key] = out
    return out


def sample_region(build, mode, args=(), exe="g1s.exe", margin=MARGIN):
    """(radius, centre_z, halfz) цилиндра розыгрыша, мм, с запасом."""
    b = sample_bounds(build, mode, args, exe)
    if b is None:
        raise SystemExit(
            "в режиме %s тома `Sample` нет — /gps/pos/confine Sample в этом "
            "режиме будет молча снят (разбор R61)." % mode)
    rmax, zlo, zhi = b
    return (rmax + margin, 0.5 * (zlo + zhi), 0.5 * (zhi - zlo) + margin)


def gps_lines(build, mode, args=(), exe="g1s.exe", margin=MARGIN):
    """Готовые команды макроса: тело розыгрыша и ограничение по пробе."""
    r, zc, hz = sample_region(build, mode, args, exe, margin)
    return ["/gps/pos/type Volume", "/gps/pos/shape Cylinder",
            "/gps/pos/centre 0 0 %.2f mm" % zc,
            "/gps/pos/radius %.2f mm" % r,
            "/gps/pos/halfz %.2f mm" % hz,
            "/gps/pos/confine Sample"]
