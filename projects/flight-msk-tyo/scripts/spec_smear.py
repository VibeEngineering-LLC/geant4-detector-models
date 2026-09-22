import sys, argparse
import numpy as np
import math

sys.stdout.reconfigure(encoding="utf-8")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input CSV file")
    parser.add_argument("output", help="Output CSV file")
    parser.add_argument("--fwhm662", type=float, default=41.6, help="FWHM at 661.657 keV in keV")
    parser.add_argument("--emin", type=int, default=20, help="Minimum energy in keV")
    parser.add_argument("--emax", type=int, default=10000, help="Maximum energy in keV")
    parser.add_argument("--light", action="store_true", help="Use light output instead of deposit")

    args = parser.parse_args()

    try:
        data = np.loadtxt(args.input, delimiter=",", skiprows=2)
    except Exception as e:
        print(f"Ошибка чтения файла: {e}", file=sys.stderr)
        sys.exit(2)

    # Разделение данных
    bin_keV = data[:, 0]
    sumw = data[:, 1]
    sumw2 = data[:, 2]
    if args.light:
        sumw = data[:, 3]
        sumw2 = data[:, 4]

    # Удалить переполненный бин
    overflow_rate_per_h = sumw[-1] * 3600
    bin_keV = bin_keV[:-1]
    sumw = sumw[:-1]
    sumw2 = sumw2[:-1]

    # Фильтрация нулевых значений
    nonzero = sumw > 0
    bin_keV = bin_keV[nonzero]
    sumw = sumw[nonzero]
    sumw2 = sumw2[nonzero]

    # Выходная сетка: бины [j, j+1) кэВ, j = emin..emax-1, центры j+0.5 (та же сетка, что во входном файле)
    lo, hi = int(args.emin), int(args.emax)
    if lo >= hi:
        print("Нет пересечения диапазонов", file=sys.stderr)
        sys.exit(2)
    nout = hi - lo
    edges = np.arange(lo, hi + 1, dtype=float)
    verf = np.vectorize(math.erf)
    smeared_rates = np.zeros(nout)
    smeared_vars = np.zeros(nout)
    for E_i, r_i, v_i in zip(bin_keV + 0.5, sumw, sumw2):
        sig = args.fwhm662 * math.sqrt(E_i / 661.657) / 2.35482
        j0 = max(0, int(math.floor(E_i - 5 * sig)) - lo)
        j1 = min(nout, int(math.ceil(E_i + 5 * sig)) - lo + 1)
        if j1 <= j0:
            continue
        cdf = 0.5 * (1.0 + verf((edges[j0:j1 + 1] - E_i) / (sig * math.sqrt(2.0))))
        w = np.diff(cdf)
        smeared_rates[j0:j1] += w * r_i
        smeared_vars[j0:j1] += w * w * v_i
    output_centers = edges[:-1] + 0.5
    counts_per_h = smeared_rates * 3600
    sigma_counts_per_h = np.sqrt(smeared_vars) * 3600

    # Запись в файл
    with open(args.output, "w") as f:
        f.write("E_keV,counts_per_h_per_keV,sigma_counts_per_h_per_keV\n")
        for E, c, s in zip(output_centers, counts_per_h, sigma_counts_per_h):
            f.write(f"{E:.1f},{c:.6f},{s:.6f}\n")

    # Подсчет суммарных значений
    total_smeared = np.sum(counts_per_h)
    total_unsmeared = np.sum(sumw[(bin_keV >= lo) & (bin_keV < hi)] * 3600)

    ratio = total_smeared / total_unsmeared if total_unsmeared > 0 else 1.0
    uncertainty_total = np.sqrt(np.sum(smeared_vars)) * 3600

    print(f"Суммарные счета (смазанные): {total_smeared:.2f}")
    print(f"Суммарные счета (не смазанные): {total_unsmeared:.2f}")
    print(f"Отношение: {ratio:.4f}")
    print(f"Стат. неопределенность: {uncertainty_total:.2f}")
    print(f"Переполнение: {overflow_rate_per_h:.2f}")

if __name__ == "__main__":
    main()
