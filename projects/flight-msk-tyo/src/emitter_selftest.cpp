#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <iomanip>
#include <cmath>
#include "capture_emitter.hh"

int main(int argc, char* argv[]) {
    if (argc < 6 || argc > 7) {
        std::cerr << "Usage: emitter_selftest <db_file> <Z> <A_target> <N> <seed> [dump_prefix]\n";
        return 2;
    }

    const std::string db_file = argv[1];
    const int Z = std::stoi(argv[2]);
    const int A_target = std::stoi(argv[3]);
    const long long N = std::stoll(argv[4]);
    const unsigned seed = std::stoul(argv[5]);
    const bool dump = argc > 6;
    const std::string dump_prefix = dump ? argv[6] : "";

    CaptureEmitter emitter;
    std::string err;
    if (!emitter.Load(db_file, err)) {
        std::cerr << err << "\n";
        return 2;
    }

    if (!emitter.Has(Z, A_target)) {
        std::cerr << "No data for Z=" << Z << ", A=" << A_target << "\n";
        return 2;
    }

    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    auto rand01 = [&]() { return u01(gen); };

    long long n_ok = 0;
    long long n_fallback = 0;
    std::vector<int> histogram(12001, 0);
    long long n_gamma_total = 0;
    double sum_E = 0.0;
    double max_abs_balance_err = 0.0;
    std::map<int, long long> primary_level_counts;

    for (long long i = 0; i < N; ++i) {
        Cascade c;
        if (emitter.Sample(Z, A_target, rand01, c)) {
            ++n_ok;
            double E_total = 0.0;
            const double M = (A_target + 1) * 931494.0;
            for (double E : c.gamma_keV) {
                sum_E += E;
                E_total += E;
                int bin = static_cast<int>(std::floor(E));
                if (bin >= 0 && bin < 12001) {
                    ++histogram[bin];
                }
                ++n_gamma_total;
            }
            E_total += c.conv_keV;
            for (double E : c.gamma_keV) {
                E_total += E * E / (2.0 * M);
            }
            const double Sn = emitter.Sn_keV(Z, A_target);
            max_abs_balance_err = std::max(max_abs_balance_err, std::abs(E_total - Sn));
            ++primary_level_counts[c.primary_level];
        } else {
            ++n_fallback;
        }
    }

    const double Wp_expected = emitter.Wp(Z, A_target);
    const double Wp_observed = static_cast<double>(n_ok) / static_cast<double>(N);

    std::cout << "N=" << N << " n_ok=" << n_ok << " n_fallback=" << n_fallback
              << " Wp_expected=" << std::setprecision(8) << Wp_expected
              << " Wp_observed=" << Wp_observed << "\n";

    const double gammas_per_capture = static_cast<double>(n_gamma_total) / static_cast<double>(n_ok);
    std::cout << "gammas_per_capture=" << std::setprecision(8) << gammas_per_capture << "\n";

    std::cout << "max_abs_balance_err_keV=" << std::setprecision(8) << max_abs_balance_err << "\n";

    std::vector<std::pair<int, int>> bins;
    for (int i = 0; i < 12001; ++i) {
        if (histogram[i] > 0) {
            bins.emplace_back(i, histogram[i]);
        }
    }

    std::sort(bins.begin(), bins.end(), [](const auto& a, const auto& b) {
        return b.second < a.second;  // по убыванию счёта (деление на n_ok лишнее, а n_ok без захвата — ошибка компиляции)
    });

    for (const auto& p : bins) {
        if (static_cast<double>(p.second) / static_cast<double>(n_ok) > 0.001) {
            std::cout << "LINE " << p.first << " per_capture=" << std::setprecision(8)
                      << static_cast<double>(p.second) / static_cast<double>(n_ok) << "\n";
        }
    }

    std::vector<std::pair<int, long long>> levels(primary_level_counts.begin(), primary_level_counts.end());
    std::sort(levels.begin(), levels.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    for (const auto& p : levels) {
        std::cout << "PRIMARY level=" << p.first << " fraction=" << std::setprecision(8)
                  << static_cast<double>(p.second) / static_cast<double>(n_ok) << "\n";
    }

    if (dump) {
        std::ofstream out(dump_prefix + "_hist.csv");
        for (int i = 0; i < 12001; ++i) {
            if (histogram[i] > 0) {
                out << i << "," << histogram[i] << "\n";
            }
        }
    }

    return 0;
}
