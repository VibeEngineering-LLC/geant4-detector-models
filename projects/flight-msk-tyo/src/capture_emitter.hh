#pragma once
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <limits>

struct Level {
    double E_level_keV = 0.0;
    std::vector<std::pair<int, double>> transitions; // daughter_idx, prob
    std::vector<double> alphas; // alpha for each transition
    std::vector<double> egs;    // энергия гамма-перехода, кэВ (параллельно transitions; первая версия её выбрасывала)
};

struct Cascade {
    std::vector<double> gamma_keV;
    double conv_keV = 0.0;
    int primary_level = -1;
};

class CaptureEmitter {
private:
    struct Isotope {
        double Sn_keV = 0.0;
        double sigma0_b = 0.0;
        double Wp = 0.0;
        std::vector<std::pair<int, double>> primary; // level_idx, weight
        std::map<int, Level> levels;
    };

    std::unordered_map<int, Isotope> isotopes;

public:
    bool Load(const std::string& path, std::string& err) {
        std::ifstream file(path);
        if (!file.is_open()) {
            err = "Не удалось открыть файл";
            return false;
        }

        isotopes.clear();
        std::string line;
        int lineno = 0;

        while (std::getline(file, line)) {
            ++lineno;
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            std::string tag;
            iss >> tag;

            if (tag == "ISO") {
                int Z, A_target;
                double Sn_keV, sigma0_b, Wp;
                int n_levels, n_primary;
                if (!(iss >> Z >> A_target >> Sn_keV >> sigma0_b >> Wp >> n_levels >> n_primary)) {
                    err = "Ошибка чтения ISO строки";
                    return false;
                }
                if (Wp < 0 || Wp > 1.000001) {
                    err = "Wp вне допустимого диапазона";
                    return false;
                }

                int key = Z * 1000 + A_target;
                auto& iso = isotopes[key];
                iso.Sn_keV = Sn_keV;
                iso.sigma0_b = sigma0_b;
                iso.Wp = Wp;

                iso.primary.reserve(n_primary);
                for (int i = 0; i < n_primary; ++i) {
                    if (!std::getline(file, line)) {
                        err = "Недостаточно строк PRI";
                        return false;
                    }
                    ++lineno;
                    std::istringstream iss2(line);
                    int level_idx;
                    double weight;
                    std::string ptag;  // строка начинается с тега PRI — его надо прочитать (первая версия этого не делала)
                    if (!(iss2 >> ptag >> level_idx >> weight) || ptag != "PRI") {
                        err = "Ошибка чтения PRI строки";
                        return false;
                    }
                    iso.primary.emplace_back(level_idx, weight);
                }

                for (int i = 0; i < n_levels; ++i) {
                    if (!std::getline(file, line)) {
                        err = "Недостаточно строк LEV";
                        return false;
                    }
                    ++lineno;
                    std::istringstream iss2(line);
                    std::string tag2;
                    iss2 >> tag2;
                    if (tag2 != "LEV") {
                        err = "Ожидается LEV";
                        return false;
                    }

                    int idx;
                    double E_level_keV;
                    int n_transitions;
                    if (!(iss2 >> idx >> E_level_keV >> n_transitions)) {
                        err = "Ошибка чтения LEV строки";
                        return false;
                    }

                    auto& level = iso.levels[idx];
                    level.E_level_keV = E_level_keV;

                    for (int j = 0; j < n_transitions; ++j) {
                        if (!std::getline(file, line)) {
                            err = "Недостаточно строк TR";
                            return false;
                        }
                        ++lineno;
                        std::istringstream iss3(line);
                        int daughter_idx;
                        double Eg_keV, prob, alpha;
                        std::string ttag;  // тег TR
                        if (!(iss3 >> ttag >> daughter_idx >> Eg_keV >> prob >> alpha) || ttag != "TR") {
                            err = "Ошибка чтения TR строки";
                            return false;
                        }
                        if (prob < 0 || !std::isfinite(prob)) {
                            err = "Некорректная вероятность TR";
                            return false;
                        }
                        level.transitions.emplace_back(daughter_idx, prob);
                        level.alphas.push_back(alpha);
                        level.egs.push_back(Eg_keV);
                    }
                }

                // Проверка, что все PRI уровни существуют
                for (const auto& p : iso.primary) {
                    if (iso.levels.find(p.first) == iso.levels.end()) {
                        err = "PRI уровень не найден в LEV";
                        return false;
                    }
                }
            } else if (tag == "END") {
                continue;  // конец блока изотопа
            } else {
                err = "Неизвестный тег в строке " + std::to_string(lineno);
                return false;
            }
        }

        if (isotopes.empty()) {
            err = "База данных пуста";
            return false;
        }

        return true;
    }

    bool Has(int Z, int A_target) const {
        return isotopes.find(Z * 1000 + A_target) != isotopes.end();
    }

    double Wp(int Z, int A_target) const {
        auto it = isotopes.find(Z * 1000 + A_target);
        if (it == isotopes.end()) return 0.0;
        return it->second.Wp;
    }

    double Sn_keV(int Z, int A_target) const {
        auto it = isotopes.find(Z * 1000 + A_target);
        if (it == isotopes.end()) return 0.0;
        return it->second.Sn_keV;
    }

    int Count() const {
        return static_cast<int>(isotopes.size());
    }

    template <class Rng>
    bool Sample(int Z, int A_target, Rng&& rand01, Cascade& out) const {
        auto it = isotopes.find(Z * 1000 + A_target);
        if (it == isotopes.end()) return false;

        const auto& iso = it->second;
        double u = rand01();
        if (u >= iso.Wp) return false;

        // Выбор первичного уровня
        double r = u;  // u равномерно на [0, Wp), веса первичных дают в сумме Wp — отдельная выборка не нужна
        double cumsum = 0.0;
        int primary_idx = -1;
        for (const auto& p : iso.primary) {
            cumsum += p.second;
            if (cumsum >= r) {
                primary_idx = p.first;
                break;
            }
        }
        if (primary_idx == -1 && !iso.primary.empty()) {
            primary_idx = iso.primary.back().first;
        }

        out.primary_level = primary_idx;
        out.gamma_keV.clear();
        out.conv_keV = 0.0;

        // Энергия первичного гамма
        double E_level = iso.levels.count(primary_idx) ? iso.levels.at(primary_idx).E_level_keV : 0.0;
        double Eg0 = iso.Sn_keV - E_level;
        if (Eg0 <= 0) return true;

        // Коррекция на релятивистский импульс
        double M = (A_target + 1) * 931494.0;
        Eg0 = Eg0 - Eg0 * Eg0 / (2 * M);

        out.gamma_keV.push_back(Eg0);
        if (E_level <= 0) return true;

        // Каскад
        int j = primary_idx;
        int steps = 0;
        while (j != 0 && steps < 200) {
            auto level_it = iso.levels.find(j);
            if (level_it == iso.levels.end()) {
                break;  // уровня нет в базе: энергию не выдумываем (Sn здесь ставила первая версия)
            }

            const auto& level = level_it->second;
            if (level.transitions.empty()) {
                // Терминальный уровень
                double E_gamma = level.E_level_keV;
                out.gamma_keV.push_back(E_gamma);
                break;
            }

            r = rand01();
            double cumprob = 0.0;
            int trans_idx = -1;
            for (size_t i = 0; i < level.transitions.size(); ++i) {
                cumprob += level.transitions[i].second;
                if (cumprob >= r) {
                    trans_idx = static_cast<int>(i);
                    break;
                }
            }
            if (trans_idx == -1) trans_idx = static_cast<int>(level.transitions.size()) - 1;

            int daughter = level.transitions[trans_idx].first;
            double Eg = level.egs[trans_idx];  // энергия перехода (в первой версии здесь стояла вероятность)
            double alpha = level.alphas[trans_idx];

            // Случай внутреннего преобразования
            if (rand01() < 1.0 / (1.0 + alpha)) {
                out.gamma_keV.push_back(Eg);
            } else {
                out.conv_keV += Eg;
            }

            j = daughter;
            ++steps;
        }

        return true;
    }
};
