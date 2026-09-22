import sys
sys.stdout.reconfigure(encoding="utf-8")

def tol_of(E, tol_a, tol_b):
    """Допуск совпадения в кэВ."""
    return tol_a + tol_b * E

def nearest_level(levels, Lt):
    """(idx, dist): индекс уровня, ближайшего по энергии к Lt, и расстояние |E_level - Lt|. levels не пуст."""
    idx, level = min(levels.items(), key=lambda item: abs(item[1]["E"] - Lt))
    return (idx, abs(level["E"] - Lt))  # первая версия возвращала словарь уровня вместо расстояния

def is_secondary_candidate(levels, E, tol):
    """True, если у КАКОГО-ЛИБО уровня есть переход с |Eg - E| <= tol."""
    for level in levels.values():
        for _, Eg, _, _ in level["trans"]:
            if abs(Eg - E) <= tol:
                return True
    return False

def classify_lines(levels, Sn, M, lines, tol_a, tol_b):
    """Список записей (E, y, is_primary, level_idx_or_None) в порядке lines."""
    result = []
    for E, y in lines:
        if y <= 0:
            result.append((E, y, False, None))
            continue
        E0 = E
        Lt_val = Sn - E0 - E0 * E0 / (2 * M)  # энергия уровня, который питала бы первичная линия E0 (один шаг, без итераций)
        j, dist = nearest_level(levels, Lt_val)
        tol = tol_of(E0, tol_a, tol_b)
        primary_candidate = dist <= tol
        secondary = is_secondary_candidate(levels, E0, tol)
        if E0 >= 0.5 * Sn:
            is_primary = primary_candidate
        else:
            is_primary = primary_candidate and not secondary
        result.append((E0, y, is_primary, j if is_primary else None))
    return result

def main():
    levels = {0: {"E": 0.0, "trans": []}, 1: {"E": 300.0, "trans": [(0, 300.0, 1.0, 0.0)]}, 2: {"E": 1000.0, "trans": [(1, 700.0, 0.75, 0.0), (0, 1000.0, 0.25, 0.0)]},
              3: {"E": 1003.0, "trans": []},   # соседний уровень: отличим от уровня 2 только с учётом отдачи (Lt = 999,7 против 1002,1 без неё)
              4: {"E": 5000.0, "trans": [(0, 5000.0, 1.0, 0.0)]},
              5: {"E": 3000.0, "trans": [(0, 3000.0, 1.0, 0.0)]}}   # переход 3000 кэВ делает линию 3000 вторичным кандидатом
    Sn = 8000.0
    M = 11 * 931494.0
    tol_a = 1.5
    tol_b = 3e-4
    lines = [(6997.9, 0.5), (7697.8, 0.2), (700.0, 0.3), (50.0, 0.1), (3000.0, 0.3)]
    expected = [
        (6997.9, 0.5, True, 2),
        (7697.8, 0.2, True, 1),
        (700.0, 0.3, False, None),
        (50.0, 0.1, False, None),
        (3000.0, 0.3, False, None)   # кандидат в первичные (Lt = 4999,6 ~ уровень 4), но и вторичный, а E < 0,5*Sn
    ]
    result = classify_lines(levels, Sn, M, lines, tol_a, tol_b)
    if result == expected:
        print("SELFTEST PASS")
        return 0
    else:
        print("SELFTEST FAIL", result)
        return 1

if __name__ == "__main__":
    exit(main())
