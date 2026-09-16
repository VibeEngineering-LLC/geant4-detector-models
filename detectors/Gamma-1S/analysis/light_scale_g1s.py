import sys, os, re, json, argparse, numpy as np
A = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.dirname(os.path.dirname(A))
sys.path.insert(0, A + "/analysis"); sys.path.insert(0, A + "/geant4/run_g1s_npsm")
import mix_unfold_g1s as g1s
from analyze_stage5 import find_peak_position

sys.stdout.reconfigure(encoding="utf-8")

def fit_line(x, y):
    # np.polyfit возвращает [наклон, свободный член] — генератор ставил их наоборот (приёмка 12.09)
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    return a, b, resid

def read_grid_file(path):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    header_line = None
    for i, line in enumerate(lines):
        if line.strip().startswith('bin_keV,count_edep,count_light'):
            header_line = i
            break
    
    if header_line is None:
        return None
    
    data = []
    for line in lines[header_line+1:]:
        parts = line.strip().split(',')
        if len(parts) >= 3:
            # битая строка таблицы — громкий ValueError, не молчаливый пропуск (спека: без глотания ошибок)
            bin_keV, count_edep, count_light = map(float, parts[:3])
            data.append((bin_keV, count_light))
    
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spe', default=A + "/reference/lsrm/raw_lsrm/Work/BG/Gamma-1S/Spe - поверки/Поверка 2016/Маринелли/Смесь_AmTiCsEu_Маринелли.spe")
    parser.add_argument('--grid', default=_REPO_ROOT + "/build/Gamma-1S-npsm-1142/out")
    parser.add_argument('--ch-offset', type=float, default=1.0)
    parser.add_argument('--out', default=A + "/web-th232/data/light_scale_amticseu.json")
    parser.add_argument('--selftest', action='store_true')
    
    args = parser.parse_args()
    
    if args.selftest:
        # Тест подгонки
        a_true, b_true = -3.0, 0.55
        L = [30, 70, 150, 200, 300, 400, 500, 650, 800]
        c = [a_true + b_true * l for l in L]
        
        a_rest, b_rest, resid = fit_line(L, c)
        if abs(a_rest - a_true) > 1e-9 or abs(b_rest - b_true) > 1e-9:
            raise SystemExit("SELFTEST FAIL: подгонка не совпала")
        
        # Имитация ошибки
        c[0] += 1.0
        _, _, resid = fit_line(L, c)
        if max(abs(r) for r in resid) <= 0.5:
            raise SystemExit("SELFTEST FAIL: максимальная невязка не превышает 0.5")
        
        print("SELFTEST OK")
        return
    
    # Считываем спектр
    spec = g1s.read_lsrm_spe(args.spe)
    e_of, n = g1s.recalibrate_energy(spec, verbose=False, ch_offset=args.ch_offset)
    pairs = e_of.pairs
    
    # Собираем кривую света
    light_curve = {}
    requires_interpretation = []
    
    for file in os.listdir(args.grid):
        if not file.startswith('gridon_E') or not file.endswith('.csv'):
            continue
        
        match = re.search(r'gridon_E([0-9.]+)\.csv', file)
        if not match:
            continue
            
        E = float(match.group(1))
        
        path = os.path.join(args.grid, file)
        data = read_grid_file(path)
        if data is None:
            requires_interpretation.append(f"Файл {file} не содержит данных")
            continue
        
        bins = [row[0] for row in data]
        count_light = [row[1] for row in data]
        
        peak_pos = find_peak_position(bins, count_light, E)
        if peak_pos is None:
            requires_interpretation.append(f"Файл {file} не содержит пика")
            continue
            
        light_curve[E] = peak_pos
    
    # Интерполяция
    sorted_E = sorted(light_curve.keys())
    if len(sorted_E) < 2:
        raise SystemExit("Недостаточно точек для интерполяции")
    
    def L_at(E):
        if E in light_curve:
            return light_curve[E]
        
        requires_interpretation.append(f"Интерполирована точка {E} кэВ")
        return np.interp(E, sorted_E, [light_curve[e] for e in sorted_E])
    
    # Реперы
    refs = []
    for c_i, E_i in pairs:
        L_i = L_at(E_i)
        refs.append({'E': E_i, 'c': c_i, 'L': L_i})
    
    if len(refs) < 3:
        raise SystemExit("Недостаточно реперных точек")
    
    # Подгонка прямых
    L_vals = [r['L'] for r in refs]
    c_vals = [r['c'] for r in refs]
    E_vals = [r['E'] for r in refs]
    
    a, b, resid_light = fit_line(L_vals, c_vals)
    a2, b2, resid_energy = fit_line(E_vals, c_vals)
    
    # Вывод результатов
    print("Реперные точки:")
    print("E (кэВ), c (каналы), предсказанный канал, невязка (каналы)")
    for r in refs:
        pred_c_light = a + b * r['L']
        pred_c_energy = a2 + b2 * r['E']
        resid_light_val = r['c'] - pred_c_light
        resid_energy_val = r['c'] - pred_c_energy
        print(f"{r['E']:.3f}, {r['c']:.3f}, {pred_c_light:.3f}, {resid_light_val:.3f}")
    
    std_light = np.std(resid_light)
    max_resid_light = max(abs(r) for r in resid_light)
    print(f"СКО невязок по свету: {std_light:.4f}")
    print(f"Максимальная невязка по свету: {max_resid_light:.4f}")
    
    std_energy = np.std(resid_energy)
    max_resid_energy = max(abs(r) for r in resid_energy)
    print(f"СКО невязок по энергии: {std_energy:.4f}")
    print(f"Максимальная невязка по энергии: {max_resid_energy:.4f}")
    
    # Проверка ниже реперов
    checks = []
    
    # Cs K-комплекс
    cs_lines = [(30.625, 32.6), (30.973, 60.2), (35.089, 17.6), (35.818, 4.29)]
    E_c = sum(e * i for e, i in cs_lines) / sum(i for e, i in cs_lines)
    L_c = sum(L_at(e) * i for e, i in cs_lines) / sum(i for e, i in cs_lines)
    c_meas_cs = 12.30
    spread_cs = 0.35
    
    c_light_cs = a + b * L_c
    c_energy_cs = a2 + b2 * E_c
    
    diff_light_cs = c_light_cs - c_meas_cs
    diff_energy_cs = c_energy_cs - c_meas_cs
    
    checks.append({
        'name': 'Cs K',
        'E_c': E_c,
        'L_c': L_c,
        'c_meas': c_meas_cs,
        'spread': spread_cs,
        'c_light': c_light_cs,
        'c_energy': c_energy_cs
    })
    
    # Sm K-комплекс
    sm_lines = [(39.522, 20.87), (40.117, 37.8), (45.523, 11.81), (46.575, 3.05)]
    E_c = sum(e * i for e, i in sm_lines) / sum(i for e, i in sm_lines)
    L_c = sum(L_at(e) * i for e, i in sm_lines) / sum(i for e, i in sm_lines)
    c_meas_sm = 14.62
    spread_sm = 0.29
    
    c_light_sm = a + b * L_c
    c_energy_sm = a2 + b2 * E_c
    
    diff_light_sm = c_light_sm - c_meas_sm
    diff_energy_sm = c_energy_sm - c_meas_sm
    
    checks.append({
        'name': 'Sm K',
        'E_c': E_c,
        'L_c': L_c,
        'c_meas': c_meas_sm,
        'spread': spread_sm,
        'c_light': c_light_sm,
        'c_energy': c_energy_sm
    })
    
    print("\nПроверка ниже реперов:")
    for check in checks:
        diff_light = check['c_light'] - check['c_meas']
        diff_energy = check['c_energy'] - check['c_meas']
        print(f"{check['name']}:")
        print(f"  Центр комплекса: {check['E_c']:.3f} кэВ, {check['L_c']:.1f} света")
        print(f"  Измеренный канал: {check['c_meas']:.2f}")
        print(f"  Предсказано по свету: {check['c_light']:.2f}, разность: {diff_light:.3f} ({diff_light/check['spread']:.2f}σ)")
        print(f"  Предсказано по энергии: {check['c_energy']:.2f}, разность: {diff_energy:.3f}")
    
    # Условия спеки шага 7, которых генератор не написал (приёмка 12.09); невязки — в JSON (шаг 6)
    for r, rl, r_e in zip(refs, resid_light, resid_energy):
        r['resid_light'], r['resid_energy'] = float(rl), float(r_e)
        if abs(rl) > 0.5:
            requires_interpretation.append(f"Репер {r['E']:.3f} кэВ: невязка по свету {rl:+.3f} кан > 0,5")
    for chk in checks:
        d = chk['c_light'] - chk['c_meas']
        if abs(d) > chk['spread']:
            requires_interpretation.append(f"{chk['name']}: разность по свету {d:+.3f} кан > 1 разброса ({chk['spread']})")
    if std_light > std_energy:
        requires_interpretation.append(f"СКО по свету {std_light:.4f} > СКО по энергии {std_energy:.4f}")

    # Требует толкования
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    if not requires_interpretation:
        print("пусто")
    else:
        for item in requires_interpretation:
            print(item)
    
    # Запись в JSON
    result = {
        'curve': sorted([[e, light_curve[e]] for e in light_curve.keys()]),
        'a': a,
        'b': b,
        'a2': a2,
        'b2': b2,
        'refs': refs,
        'checks': checks,
        'ch_offset': args.ch_offset,
        'spe': args.spe
    }
    
    with open(args.out, 'w', encoding='utf-8') as f:
        json.dump(result, f, ensure_ascii=False, indent=1)

if __name__ == '__main__':
    main()
