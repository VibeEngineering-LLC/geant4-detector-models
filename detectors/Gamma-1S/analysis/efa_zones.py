"""Полином зон .efa СпектраЛайн: разбор, пересчёт под свои eps, сборка.

ЗАЧЕМ. Программа считает активность НЕ по точкам файла эффективности, а по
полиному зон. Подмена одних точек ничего не меняет: три прогона на «нашей»
кривой (Th-228, точечная 5 см, 30.07.2026) дали фактически применённую
эффективность 3,152E-02 на 238,6 кэВ — узел ЛСРМ, а не наш 3,472E-02.
Пока секция Zones остаётся мастерской, кривая остаётся штатной, какой бы
детектор ни стоял в заголовке.

КОНВЕНЦИЯ (восстановлена по рабочему .efa прибора, не документирована).
Для каждой зоны:

    Zone_k   = степень, log10(E_min), log10(E_max), СКО фита
    Curve_k_j= базисный полином номер j (j = 1..степень+1),
               коэффициенты по УБЫВАЮЩИМ степеням x = log10(E)
    Curve_k  = коэффициенты c_j при базисных полиномах

    eps(E) = 10 ** sum_j c_j * P_j(log10 E)

Базис ортогонален не на отрезке, а НА НАБОРЕ УЗЛОВ с весами 1/σ²: корень
P_2 первой зоны (x = 2,5920) совпадает со средневзвешенным log10(E) её
девяти узлов (2,5897), а не с серединой отрезка (2,2325).

ПРОВЕРКА РЕЦЕПТА. Взвешенный МНК по узлам мастера в его же базисе
воспроизводит записанные коэффициенты до последнего знака (зона 1:
-1.706245 -0.155679 -0.045806) — значит, обратная задача решена верно и
пересчёт под свои eps законен.

ПОБОЧНО. Сам штатный полином отходит от своих калибровочных точек до
3,4 % в первой зоне и до 10,7 % во второй: «кривая ЛСРМ» в расчёте и
«узлы ЛСРМ» в файле — разные вещи.

Базис при пересчёте СОХРАНЯЕТСЯ: он определяется энергиями узлов и их
весами, а мы меняем только значения eps. Пересчитываются лишь Curve_k.
"""
import math

# Перевод «разы <-> десятичный логарифм». Множитель 230 = 100*ln10 в формуле
# погрешности ЛСРМ (см. fit_zone) — это он же, вместе с процентами.
LN10 = math.log(10.0)


def horner(cs, x):
    """Полином, коэффициенты по убывающим степеням."""
    r = 0.0
    for c in cs:
        r = r * x + c
    return r


def parse_block(lines):
    """Разбор блока .efa: узлы, зоны, базис, коэффициенты."""
    nodes, zones, basis, coef = [], {}, {}, {}
    for ln in lines:
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        if k[:1].isdigit():
            f = v.split(",")
            nodes.append((float(k), float(f[0]), float(f[1])))
        elif k.startswith("Zone_"):
            zones[int(k[5:])] = [float(x) for x in v.split(",")]
        elif k.startswith("Curve_") and "_" in k[6:]:
            z, j = k[6:].split("_")
            basis[(int(z), int(j))] = [float(x) for x in v.split(",")]
        elif k.startswith("Curve_"):
            coef[int(k[6:])] = [float(x) for x in v.split(",")]
    return nodes, zones, basis, coef


def _solve(A, b):
    """Гаусс с выбором главного элемента; матрица малая (3x3)."""
    n = len(b)
    A = [row[:] for row in A]
    b = b[:]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]
        b[i], b[p] = b[p], b[i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
            b[r] -= f * b[i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]
    return x


def fit_zone(zone, basis_rows, nodes, eps_of):
    """Коэффициенты Curve_k под свою кривую eps_of(E) на узлах зоны.

    nodes — [(E, eps, d_eps_%)]; веса берутся из d_eps штатного файла:
    они задают, насколько узел затягивает фит, и должны остаться теми же,
    иначе изменится и базис.

    ШКАЛА ВЕСА И РАЗБРОСА — ДЕСЯТИЧНЫЙ ЛОГАРИФМ (исправлено 01.08.2026 по
    внешнему аудиту). Отклик фита есть y = lg(eps), поэтому и вес узла, и
    разброс обязаны быть в той же шкале: заявленная погрешность узла d %
    переводится в lg как d/100/ln10. Прежняя реализация брала вес
    1/(d/100)^2 и разброс curve/eps-1 — то есть относительную шкалу при
    логарифмическом отклике. На коэффициенты это не влияло (общий
    множитель ln10^2 в весах сокращается), но 4-е поле Zone_k, куда
    записывается разброс, конвенцией ЛСРМ определено в единицах lg: по
    формуле производителя погрешность кривой есть 230*D*||Q(X)||, где
    множитель 230 = 100*ln10 переводит D из lg в проценты. Запись
    линейного разброса в это поле завышала объявленную погрешность в
    ln10 ~ 2,3 раза.

    Возвращает (коэффициенты, СКО в единицах lg, макс. отклонение в
    единицах lg, число узлов).
    """
    deg, xmin, xmax = int(zone[0]), zone[1], zone[2]
    inz = [(E, d) for E, _e, d in nodes
           if xmin - 1e-9 <= math.log10(E) <= xmax + 1e-9]
    n = deg + 1
    A = [[0.0] * n for _ in range(n)]
    b = [0.0] * n
    for E, d in inz:
        x = math.log10(E)
        w = 1.0 / (d / 100.0 / LN10) ** 2
        ph = [horner(basis_rows[j], x) for j in range(n)]
        y = math.log10(eps_of(E))
        for i in range(n):
            b[i] += w * ph[i] * y
            for j in range(n):
                A[i][j] += w * ph[i] * ph[j]
    c = _solve(A, b)

    def curve(E):
        x = math.log10(E)
        return 10.0 ** sum(c[j] * horner(basis_rows[j], x) for j in range(n))

    dev = [math.log10(curve(E) / eps_of(E)) for E, _d in inz]
    rms = math.sqrt(sum(v * v for v in dev) / len(dev))
    return c, rms, max(abs(v) for v in dev), len(inz)


def fmt_coef(v):
    """Формат коэффициента как в файле: .123456789012E+1"""
    s = "%.12E" % abs(v)
    m, e = s.split("E")
    m = m.replace(".", "")
    return "%s.%sE%s%d" % ("-" if v < 0 else "", m,
                           "+" if int(e) + 1 >= 0 else "-", abs(int(e) + 1))


def rewrite_block(lines, eps_of):
    """Блок .efa со своими eps в точках И пересчитанным полиномом зон.

    Возвращает (строки, отчёт по зонам) — отчёт печатается вызывающим,
    чтобы качество фита не проходило молча. В отчёте разброс приведён к
    ПРОЦЕНТАМ (умножением на 100*ln10, конвенция ЛСРМ), а в поле Zone_k
    записывается исходное значение в единицах lg — см. fit_zone.
    """
    nodes, zones, basis, _coef = parse_block(lines)
    report = []
    new_coef = {}
    rms_lg = {}
    for z, zf in sorted(zones.items()):
        rows = [basis[(z, j + 1)] for j in range(int(zf[0]) + 1)]
        c, rms, dmax, k = fit_zone(zf, rows, nodes, eps_of)
        new_coef[z] = c
        rms_lg[z] = rms
        report.append((z, 10 ** zf[1], 10 ** zf[2], k,
                       100 * LN10 * rms, 100 * LN10 * dmax))

    out = []
    for ln in lines:
        if "=" in ln:
            k, rest = ln.split("=", 1)
            if k[:1].isdigit():
                parts = rest.split(",")
                parts[0] = "%.6E" % eps_of(float(k))
                ln = k + "=" + ",".join(parts)
            elif k.startswith("Curve_") and "_" not in k[6:]:
                z = int(k[6:])
                if z in new_coef:
                    ln = k + "=" + ",".join(fmt_coef(v) for v in new_coef[z])
            elif k.startswith("Zone_"):
                z = int(k[5:])
                if z in new_coef:
                    f = ln.split("=", 1)[1].split(",")
                    # 4-е поле Zone_k — параметр D конвенции ЛСРМ, в
                    # единицах lg (не в разах и не в процентах): по формуле
                    # производителя погрешность есть 230*D*||Q(X)|| %.
                    f[3] = "%.11f" % rms_lg[z]
                    ln = k + "=" + ",".join(f)
        out.append(ln)
    return out, report
