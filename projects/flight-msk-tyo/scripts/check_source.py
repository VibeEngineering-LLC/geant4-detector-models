import sys
import math
sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

def read_table(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            lines = [line.strip() for line in f.readlines()]
    except Exception as e:
        print(f"Ошибка чтения файла таблицы: {e}", file=sys.stderr)
        sys.exit(2)

    if len(lines) < 3:
        print("Недостаточно строк в файле таблицы", file=sys.stderr)
        sys.exit(2)

    # Пропускаем комментарий
    line2 = lines[1].split()
    if len(line2) != 6:
        print("Неверный формат строки 2 файла таблицы (ожидалось 6 полей)", file=sys.stderr)
        sys.exit(2)
    nebin, nabin, ie511, flux511_cont, flux511_line, total_flux = map(float, line2)
    nebin = int(nebin)
    nabin = int(nabin)
    ie511 = int(ie511)

    # Читаем границы энергетических интервалов
    edges = list(map(float, lines[2].split()))
    if len(edges) != nebin + 1:
        print("Неверное количество границ энергии", file=sys.stderr)
        sys.exit(2)

    # Читаем матрицу D[k][ia]
    D = []
    for i in range(nebin):
        if i + 3 >= len(lines):
            print("Недостаточно строк в файле таблицы для данных D", file=sys.stderr)
            sys.exit(2)
        row = list(map(float, lines[i+3].split()))
        if len(row) != nabin:
            print(f"Неверное количество столбцов в строке {i+4} файла таблицы", file=sys.stderr)
            sys.exit(2)
        D.append(row)

    # Вычисляем массы и вероятности
    mass = []
    total_mass = 0.0
    for k in range(nebin):
        m_k = sum(D[k]) * (edges[k+1] - edges[k])
        if k == ie511 - 1 and ie511 > 0:
            m_k += flux511_line
        mass.append(m_k)
        total_mass += m_k

    p = [m / total_mass for m in mass]

    # Вычисляем условные вероятности q[k][ia]
    q = []
    for k in range(nebin):
        row_sum = sum(D[k])
        if row_sum == 0:
            q.append([0.0] * nabin)
        else:
            q.append([D[k][ia] / row_sum for ia in range(nabin)])

    # Вычисляем косинусные распределения c_ia
    c = [0.0] * nabin
    for ia in range(nabin):
        for k in range(nebin):
            c[ia] += p[k] * q[k][ia]

    return nebin, nabin, edges, D, p, c

def read_csv(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            lines = [line.strip() for line in f.readlines()]
    except Exception as e:
        print(f"Ошибка чтения CSV файла: {e}", file=sys.stderr)
        sys.exit(2)

    if len(lines) < 3:
        print("Недостаточно строк в CSV файле", file=sys.stderr)
        sys.exit(2)

    # Парсим заголовок
    header = lines[0]
    tokens = {}
    for part in header[2:].split():
        if '=' in part:
            key, value = part.split('=', 1)
            tokens[key] = value

    try:
        N = int(tokens['N'])
        hits = int(tokens['hits'])
        underflow = int(tokens['underflow'])
        overflow = int(tokens['overflow'])
        R_cm = float(tokens['R_cm'])
        probe_r_cm = float(tokens['probe_r_cm'])
    except Exception as e:
        print(f"Ошибка парсинга заголовка CSV: {e}", file=sys.stderr)
        sys.exit(2)

    # Читаем E_hist
    if not lines[1].startswith("E_hist,"):
        print("Неверный формат строки 2 CSV файла", file=sys.stderr)
        sys.exit(2)
    E_hist = list(map(int, lines[1][7:].split(',')))
    if len(E_hist) < 20:
        print("E_hist слишком короткая", file=sys.stderr)
        sys.exit(2)

    # Читаем C_hist
    if not lines[2].startswith("C_hist,"):
        print("Неверный формат строки 3 CSV файла", file=sys.stderr)
        sys.exit(2)
    C_hist = list(map(int, lines[2][7:].split(',')))

    return N, hits, underflow, overflow, R_cm, probe_r_cm, E_hist, C_hist

def main():
    if len(sys.argv) != 3:
        print("Использование: python check_source.py <table_file> <src_check_csv>", file=sys.stderr)
        sys.exit(2)

    table_path = sys.argv[1]
    csv_path = sys.argv[2]

    nebin, nabin, edges, D, p, c = read_table(table_path)
    N, hits, underflow, overflow, R_cm, probe_r_cm, E_hist, C_hist = read_csv(csv_path)

    # Проверка H (hit fraction)
    f = (probe_r_cm / R_cm) ** 2
    exp = N * f
    if N * f * (1 - f) == 0:
        z = float('inf')
    else:
        z = (hits - exp) / math.sqrt(N * f * (1 - f))
    h_pass = abs(z) < 4.5
    print(f"CHECK H hits={hits} exp={exp:.2f} z={z:.2f} {'PASS' if h_pass else 'FAIL'}")

    # Проверка S (sums)
    sum_E = sum(E_hist) + underflow + overflow
    sum_C = sum(C_hist)
    s_pass = (sum_E == hits) and (sum_C == hits)
    print(f"CHECK S sum_E={sum_E} sum_C={sum_C} hits={hits} {'PASS' if s_pass else 'FAIL'}")

    # Проверка E (energy groups)
    group_size = 20  # по спецификации: группы по 20 соседних тонких бинов
    ngroups = nebin // group_size
    if ngroups * group_size != nebin:
        ngroups += 1

    chi2_E = 0.0
    max_z_E = 0.0
    ndf_E = 0
    for g in range(ngroups):
        start_bin = g * group_size
        end_bin = min((g + 1) * group_size, nebin)
        exp_g = hits * sum(p[k] for k in range(start_bin, end_bin))  # p[k] уже вероятность ячейки
        obs_g = sum(E_hist[k] for k in range(start_bin, end_bin))
        if exp_g >= 25:
            ndf_E += 1
            if exp_g > 0:
                z = (obs_g - exp_g) / math.sqrt(exp_g)
                max_z_E = max(max_z_E, abs(z))
                chi2_E += z * z
    if ndf_E == 0:
        chi2_ndf_E = 0.0
    else:
        chi2_ndf_E = chi2_E / ndf_E
    e_pass = (max_z_E < 4.5) and (chi2_ndf_E < 1.6)
    print(f"CHECK E groups={ndf_E} max_abs_z={max_z_E:.2f} chi2_ndf={chi2_ndf_E:.2f} {'PASS' if e_pass else 'FAIL'}")

    # Проверка C (cosine bins)
    chi2_C = 0.0
    max_z_C = 0.0
    ndf_C = 0
    for ia in range(nabin):
        exp_ia = hits * c[ia]
        obs_ia = C_hist[ia]
        if exp_ia >= 25:
            ndf_C += 1
            if exp_ia > 0:
                z = (obs_ia - exp_ia) / math.sqrt(exp_ia)
                max_z_C = max(max_z_C, abs(z))
                chi2_C += z * z
    if ndf_C == 0:
        chi2_ndf_C = 0.0
    else:
        chi2_ndf_C = chi2_C / ndf_C
    c_pass = (max_z_C < 4.5) and (chi2_ndf_C < 1.6)
    print(f"CHECK C bins={ndf_C} max_abs_z={max_z_C:.2f} chi2_ndf={chi2_ndf_C:.2f} {'PASS' if c_pass else 'FAIL'}")

    all_pass = h_pass and s_pass and e_pass and c_pass
    print("SOURCE_CHECK " + ("PASS" if all_pass else "FAIL"))
    sys.exit(0 if all_pass else 1)

if __name__ == "__main__":
    main()
