import sys
sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

def read_table(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f.readlines()]
    
    # Строка 1: комментарий
    line1 = lines[0]
    if not line1.startswith("#"):
        raise ValueError("Первая строка должна быть комментарием")
    ip = None
    for token in line1[1:].split():
        if token.startswith("ip="):
            ip = token[3:]
            break
    if ip is None:
        raise ValueError("Не найден параметр ip в комментарии")
    
    # Строка 2: данные
    parts = lines[1].split()
    if len(parts) != 6:
        raise ValueError("Неверное количество полей во второй строке")
    nebin, nabin, ie511 = map(int, parts[:3])
    flux511_cont, flux511_line, total_flux = map(float, parts[3:])
    
    # Строка 3: границы энергетических интервалов
    ehigh = list(map(float, lines[2].split()))
    if len(ehigh) != nebin + 1:
        raise ValueError("Неверное количество границ энергий")
    
    # Далее строки с данными D[k][ia]
    D = []
    for i in range(nebin):
        if i + 3 >= len(lines):
            raise ValueError("Недостаточно строк данных")
        row = list(map(float, lines[i + 3].split()))
        if len(row) != nabin:
            raise ValueError("Неверное количество столбцов в строке данных")
        D.append(row)
    
    return {
        "nebin": nebin,
        "nabin": nabin,
        "ie511": ie511,
        "flux511_cont": flux511_cont,
        "flux511_line": flux511_line,
        "total_flux": total_flux,
        "ehigh": ehigh,
        "D": D,
        "ip": ip
    }

def main():
    if len(sys.argv) < 3:
        print("Недостаточно аргументов", file=sys.stderr)
        sys.exit(2)
    
    out_file = sys.argv[1]
    args = sys.argv[2:]
    
    tables = []
    weights = []
    
    for arg in args:
        colon_idx = arg.rfind(":")
        if colon_idx == -1:
            print(f"Неверный формат аргумента: {arg}", file=sys.stderr)
            sys.exit(2)
        path = arg[:colon_idx]
        try:
            weight = float(arg[colon_idx+1:])
        except ValueError:
            print(f"Неверный вес в аргументе: {arg}", file=sys.stderr)
            sys.exit(2)
        try:
            table = read_table(path)
        except Exception as e:
            print(f"Ошибка чтения файла {path}: {e}", file=sys.stderr)
            sys.exit(2)
        tables.append(table)
        weights.append(weight)
    
    # Проверка согласованности
    first = tables[0]
    for i, t in enumerate(tables):
        if t["nebin"] != first["nebin"]:
            print("nebin не совпадают", file=sys.stderr)
            sys.exit(2)
        if t["nabin"] != first["nabin"]:
            print("nabin не совпадают", file=sys.stderr)
            sys.exit(2)
        if t["ie511"] != first["ie511"]:
            print("ie511 не совпадают", file=sys.stderr)
            sys.exit(2)
        if t["ip"] != first["ip"]:
            print("ip не совпадают", file=sys.stderr)
            sys.exit(2)
        if any(abs(a - b) > 1e-8 * max(abs(a), abs(b), 1e-300) for a, b in zip(first["ehigh"], t["ehigh"])):
            print("энергетические границы не совпадают", file=sys.stderr)
            sys.exit(2)
    
    nebin = first["nebin"]
    nabin = first["nabin"]
    ie511 = first["ie511"]
    ip = first["ip"]
    
    # Вычисление смеси
    D_mix = [[0.0] * nabin for _ in range(nebin)]
    flux511_cont_mix = 0.0
    flux511_line_mix = 0.0
    total_flux_mix = 0.0
    
    for i, (t, w) in enumerate(zip(tables, weights)):
        flux511_cont_mix += w * t["flux511_cont"]
        flux511_line_mix += w * t["flux511_line"]
        total_flux_mix += w * t["total_flux"]
        for k in range(nebin):
            for ia in range(nabin):
                D_mix[k][ia] += w * t["D"][k][ia]
    
    # Вывод
    with open(out_file, "w", encoding="utf-8", newline="\n") as f:
        # Строка 1
        weights_str = ",".join(f"{w:.9f}" for w in weights)
        sum_w = sum(weights)
        f.write(f"# mix ip={ip} n_inputs={len(tables)} weights={weights_str} sum_w={sum_w:.9f}\n")
        
        # Строка 2
        f.write(f"{nebin} {nabin} {ie511} {flux511_cont_mix:.9e} {flux511_line_mix:.9e} {total_flux_mix:.9e}\n")
        
        # Строка 3: границы энергий
        f.write(" ".join(f"{x:.9e}" for x in first["ehigh"]) + "\n")
        
        # Данные D
        for row in D_mix:
            f.write(" ".join(f"{x:.9e}" for x in row) + "\n")
    
    # Вывод в stdout
    print(f"total_flux_mix={total_flux_mix:.9e}")
    
    # Вычисление континуума
    continuum_check = 0.0
    ehigh = first["ehigh"]
    for k in range(nebin):
        de = ehigh[k+1] - ehigh[k]
        for ia in range(nabin):
            continuum_check += D_mix[k][ia] * de
    
    print(f"continuum_check={continuum_check:.9e}")

if __name__ == "__main__":
    main()
