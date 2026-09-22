#include "parma_source.hh"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <limits>
#include <iomanip>

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <table_file> <N_events> <radius_cm> <seed>\n";
        return 2;
    }

    std::string table_file = argv[1];
    long long N = std::stoll(argv[2]);
    double radius_cm = std::stod(argv[3]);
    unsigned int seed = static_cast<unsigned int>(std::stoul(argv[4]));

    ParmaSource src;
    std::string err;
    if (!src.Load(table_file, err)) {
        std::cerr << err << "\n";
        return 2;
    }

    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    auto rand01 = [&]() { return u01(gen); };

    const int nebin = src.NEBin();
    const int nabin = src.NABin();
    const double eps_geom = 1e-9;
    const double eps_dist = 1e-6 * radius_cm;
    const double eps_energy = 1e-9;
    const double chi2_ndf_max = 1.6;
    const double z_max = 4.5;
    const double exp_min = 25.0;

    // Подготовка ожидаемых значений
    std::vector<double> expected_energy_groups(nebin / 20, 0.0);
    for (int g = 0; g < nebin / 20; ++g) {
        double sum_prob = 0.0;
        for (int k = 20 * g + 1; k <= 20 * (g + 1); ++k) {
            sum_prob += src.EnergyBinProb(k);
        }
        expected_energy_groups[g] = N * sum_prob;
    }

    std::vector<double> expected_cosine(nabin, 0.0);
    for (int ia = 1; ia <= nabin; ++ia) {
        expected_cosine[ia - 1] = N * src.CosBinProbMarginal(ia);
    }

    // Счетчики
    long long geometry_violations = 0;
    std::vector<long long> observed_energy_groups(nebin / 20, 0);
    std::vector<long long> observed_cosine(nabin, 0);
    long long count_511 = 0;

    double mean_cos_down = 0.0;
    double mean_energy = 0.0;

    for (long long i = 0; i < N; ++i) {
        auto event = src.Sample(radius_cm, rand01);

        // Проверка геометрии
        double d_norm = std::sqrt(event.u * event.u + event.v * event.v + event.w * event.w);
        if (std::abs(d_norm - 1.0) > eps_geom) {
            ++geometry_violations;
        }

        double p_dot_d = event.x_cm * event.u + event.y_cm * event.v + event.z_cm * event.w;
        if (std::abs(p_dot_d + radius_cm) > eps_dist) {
            ++geometry_violations;
        }

        double cross_x = event.y_cm * event.w - event.z_cm * event.v;
        double cross_y = event.z_cm * event.u - event.x_cm * event.w;
        double cross_z = event.x_cm * event.v - event.y_cm * event.u;
        double cross_norm = std::sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
        if (cross_norm > radius_cm * (1.0 + eps_geom)) {
            ++geometry_violations;
        }

        // Энергия
        const auto& edges = src.EdgesMeV();
        int k = static_cast<int>(std::lower_bound(edges.begin() + 1, edges.end(), event.e_MeV) - edges.begin());
        if (k == 0 || k > nebin) {
            k = (event.e_MeV <= edges[0]) ? 1 : nebin;
        }
        int fine_bin = k;

        // Группировка энергии
        int group = (fine_bin - 1) / 20;
        if (group < observed_energy_groups.size()) {
            ++observed_energy_groups[group];
        }

        // Косинус
        double cos_down = event.cos_down;
        int ia = static_cast<int>(std::floor((cos_down + 1.0) * nabin / 2.0)) + 1;
        ia = std::max(1, std::min(nabin, ia));
        ++observed_cosine[ia - 1];

        // Средние значения
        mean_cos_down += cos_down;
        mean_energy += event.e_MeV;

        // 511 кэВ линия
        if (std::abs(event.e_MeV - 0.51099895) < eps_energy) {
            ++count_511;
        }
    }

    mean_cos_down /= N;
    mean_energy /= N;

    // Проверка A
    bool check_a_pass = (geometry_violations == 0);
    std::cout << "CHECK A geometry violations=" << geometry_violations << " " << (check_a_pass ? "PASS" : "FAIL") << "\n";

    // Проверка B
    double chi2_energy = 0.0;
    double max_z_energy = 0.0;
    int ndf_energy = 0;
    bool check_b_pass = true;

    for (size_t g = 0; g < observed_energy_groups.size(); ++g) {
        if (expected_energy_groups[g] >= exp_min) {
            double z = (static_cast<double>(observed_energy_groups[g]) - expected_energy_groups[g]) / std::sqrt(expected_energy_groups[g]);
            chi2_energy += z * z;
            max_z_energy = std::max(max_z_energy, std::abs(z));
            ++ndf_energy;
            if (std::abs(z) > z_max) check_b_pass = false;
        }
    }

    if (ndf_energy > 0) {
        double chi2_ndf = chi2_energy / ndf_energy;
        if (chi2_ndf > chi2_ndf_max) check_b_pass = false;
    } else {
        check_b_pass = false;
    }

    std::cout << "CHECK B energy groups=" << ndf_energy << " max_abs_z=" << max_z_energy <<" chi2_ndf=" << (ndf_energy ? chi2_energy / ndf_energy : 0.0) << " " << (check_b_pass ? "PASS" : "FAIL") << "\n";

    // Проверка C
    double chi2_cosine = 0.0;
    double max_z_cosine = 0.0;
    int ndf_cosine = 0;
    bool check_c_pass = true;

    for (size_t ia = 0; ia < observed_cosine.size(); ++ia) {
        if (expected_cosine[ia] >= exp_min) {
            double z = (static_cast<double>(observed_cosine[ia]) - expected_cosine[ia]) / std::sqrt(expected_cosine[ia]);
            chi2_cosine += z * z;
            max_z_cosine = std::max(max_z_cosine, std::abs(z));
            ++ndf_cosine;
            if (std::abs(z) > z_max) check_c_pass = false;
        }
    }

    if (ndf_cosine > 0) {
        double chi2_ndf = chi2_cosine / ndf_cosine;
        if (chi2_ndf > chi2_ndf_max) check_c_pass = false;
    } else {
        check_c_pass = false;
    }

    std::cout << "CHECK C cosine bins=" << ndf_cosine << " max_abs_z=" << max_z_cosine <<" chi2_ndf=" << (ndf_cosine ? chi2_cosine / ndf_cosine : 0.0) << " " << (check_c_pass ? "PASS" : "FAIL") << "\n";

    // Проверка D
    bool check_d_pass = true;
    if (src.IE511() > 0) {
        double exp_511 = N * src.EnergyBinProb(src.IE511()) * src.LineFrac();
        double z = (count_511 - exp_511) / std::sqrt(exp_511);
        if (exp_511 < exp_min) {
            std::cout << "CHECK D line511 SKIP\n";
        } else {
            if (std::abs(z) > z_max) check_d_pass = false;
            std::cout << "CHECK D line511 obs=" << count_511 << " exp=" << exp_511 << " z=" << z << " " << (check_d_pass ? "PASS" : "FAIL") << "\n";
        }
    } else {
        std::cout << "CHECK D line511 n/a\n";
    }

    // Информация
    std::cout << "INFO mean_cos_down=" << std::setprecision(10) << mean_cos_down << " mean_energy_MeV=" << std::setprecision(10) << mean_energy << "\n";

    // Финальный результат
    bool all_pass = check_a_pass && check_b_pass && check_c_pass && check_d_pass;
    std::cout << "SELFTEST " << (all_pass ? "PASS" : "FAIL") << "\n";

    return all_pass ? 0 : 1;
}
