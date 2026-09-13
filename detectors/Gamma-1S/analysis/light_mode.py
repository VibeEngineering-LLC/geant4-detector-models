import sys
import json
import argparse

sys.stdout.reconfigure(encoding="utf-8")

def load_scale(path):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    try:
        a = data["a"]
        b = data["b"]
        if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
            raise SystemExit(f"Ошибка в файле {path}: поле 'a' или 'b' не является числом")
        return (a, b)
    except KeyError as e:
        raise SystemExit(f"Ошибка в файле {path}: отсутствует поле {e}")

def light_to_energy(hist, a, b, e_of_ch, n_channels, l_min=1.0):
    result = {}
    stats = {
        "n_in": 0,
        "n_used": 0,
        "n_below_lmin": 0,
        "n_out_of_range": 0,
        "sum_in": 0,
        "sum_used": 0
    }
    
    for light, cnt in hist.items():
        stats["n_in"] += 1
        stats["sum_in"] += cnt
        
        if light < l_min:
            stats["n_below_lmin"] += 1
            continue
            
        ch = a + b * light
        if ch < 0 or ch > n_channels - 1:
            stats["n_out_of_range"] += 1
            continue
            
        E = e_of_ch(ch)
        if E in result:
            result[E] += cnt
        else:
            result[E] = cnt
        stats["n_used"] += 1
        stats["sum_used"] += cnt
        
    return (result, stats)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    
    if not args.selftest:
        raise SystemExit("Не реализовано: только самопроверка")
        
    # Проверка 1: точность перевода
    a, b = 1.5, 0.6
    e_of_ch = lambda c: 2.5 * c + 10.0
    n_channels = 100
    
    hist = {100.0: 7.0}
    result, stats = light_to_energy(hist, a, b, e_of_ch, n_channels)
    
    if len(result) != 1:
        raise SystemExit("SELFTEST FAIL: перевод: ожидается один ключ")
        
    E = list(result.keys())[0]
    expected_E = 163.75
    if abs(E - expected_E) > 1e-9:
        raise SystemExit(f"SELFTEST FAIL: перевод: ожидается {expected_E}, получено {E}")
        
    if result[E] != 7.0:
        raise SystemExit(f"SELFTEST FAIL: перевод: значение не совпадает")
    
    # Проверка 2: отбрасывание
    hist = {0.5: 3.0, 100.0: 7.0, 1000.0: 5.0}
    result, stats = light_to_energy(hist, a, b, e_of_ch, n_channels)
    
    if (stats["n_used"] != 1 or 
        stats["n_below_lmin"] != 1 or 
        stats["n_out_of_range"] != 1 or
        stats["sum_in"] != 15.0 or
        stats["sum_used"] != 7.0):
        raise SystemExit("SELFTEST FAIL: отбрасывание")
    
    # Проверка 3: слияние ключей. Приёмка 12.09: было L2 = L1 + 1e-13 — такие ключи РАЗЛИЧАЮТСЯ
    # (шаг float у 100.0 около 1,4e-14), в результате оказывалось два ключа, и ветка слияния не
    # исполнялась вовсе. Слияние возможно только при шкале, которая квантует канал, — так и берём.
    e_step = lambda c: 2.5 * round(c) + 10.0          # 100,0 и 100,4 дают каналы 61,5 и 61,74 → один ключ
    result, stats = light_to_energy({100.0: 3.0, 100.4: 4.0}, a, b, e_step, n_channels)
    if len(result) != 1 or abs(sum(result.values()) - 7.0) > 1e-9 or stats["n_used"] != 2:
        raise SystemExit(f"SELFTEST FAIL: слияние: ключей {len(result)}, "
                         f"сумма {sum(result.values())}, n_used {stats['n_used']}")

    # Проверка 4: монотонность — ПО ВЫХОДУ функции. Приёмка 12.09: прежняя версия считала канал и
    # энергию сама и проверяла собственную арифметику, то есть не краснела ни на одной мутации.
    lights = [10.0, 20.0, 50.0, 120.0]
    result, stats = light_to_energy({L: 1.0 for L in lights}, a, b, e_of_ch, n_channels)
    energies = sorted(result)
    if len(energies) != len(lights):
        raise SystemExit(f"SELFTEST FAIL: монотонность: ключей {len(energies)} вместо {len(lights)}")
    if any(energies[i] <= energies[i - 1] for i in range(1, len(energies))):
        raise SystemExit(f"SELFTEST FAIL: монотонность: {energies}")
            
    print("SELFTEST OK")

if __name__ == "__main__":
    main()
