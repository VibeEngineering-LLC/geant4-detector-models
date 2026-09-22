import sys
import os
from collections import deque

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from capture_data import merge_lines  # после sys.path.insert: иначе модуль не импортируется из другого каталога

def reachable_levels(levels, primary):
    """Множество индексов уровней, достижимых из ключей primary по переходам (включая сами ключи)."""
    reachable = set()
    queue = deque(j for j in primary if j in levels)
    while queue:
        j = queue.popleft()
        if j in reachable: continue
        reachable.add(j)
        for d, _, _, _ in levels[j]["trans"]:
            if d in levels and d not in reachable:
                queue.append(d)
    return reachable

def propagate_yields(levels, primary, Sn, M):
    """Список (E_keV, yield_per_capture) всех гамма модели: первичных и каскадных, склеенный merge_lines."""
    v = {j: 0.0 for j in levels}
    for j, w in primary.items():
        if j in levels:
            v[j] += w
    reachable = reachable_levels(levels, primary)
    yields = []
    for j in sorted(reachable, key=lambda x: levels[x]["E"], reverse=True):
        for d, Eg, prob, alpha in levels[j]["trans"]:
            if prob <= 0.0: continue
            flow = v[j] * prob
            yield_gamma = flow / (1.0 + alpha)
            yields.append((Eg, yield_gamma))
            if d in v:
                v[d] += flow
    for j, w in primary.items():
        if j not in levels: continue
        E0 = Sn - levels[j]["E"]
        E0 = E0 - E0 * E0 / (2.0 * M)
        if E0 > 0:
            yields.append((E0, w))
    return merge_lines(yields)

def selftest():
    levels = {
        0: {"E": 0.0, "trans": []},
        1: {"E": 300.0, "trans": [(0, 300.0, 1.0, 1.0)]},
        2: {"E": 1000.0, "trans": [(1, 700.0, 0.75, 0.0), (0, 1000.0, 0.25, 0.0)]}
    }
    primary = {2: 0.6, 1: 0.2}
    Sn = 8000.0
    M = 11 * 931494.0
    result = propagate_yields(levels, primary, Sn, M)
    expected_lines = [
        (700.0, 0.45),
        (1000.0, 0.15),
        (300.0, 0.325),
        (6997.9, 0.6),
        (7697.5, 0.2)
    ]
    reachable = reachable_levels(levels, {2: 0.6})
    if reachable != {0, 1, 2}:
        print("SELFTEST FAIL reachable_levels")
        return
    for E, y in expected_lines:
        found = False
        for E0, y0 in result:
            if abs(E - E0) < 2.0 and abs(y - y0) < 1e-9:
                found = True
                break
        if not found:
            print(f"SELFTEST FAIL expected line {E} keV with yield {y}")
            return
    print("SELFTEST PASS")

if __name__ == "__main__":
    selftest()
