#pragma once
// Взвешивающий коэффициент излучения wR по ICRP 103 (2007), таблица 2 и формула для нейтронов (E в МэВ).
#include <cmath>
#include <cstdlib>

inline double WrNeutron(double E) {
    if (E < 1.0) return 2.5 + 18.2 * std::exp(-std::pow(std::log(E), 2) / 6.0);
    if (E <= 50.0) return 5.0 + 17.0 * std::exp(-std::pow(std::log(2.0 * E), 2) / 6.0);
    return 2.5 + 3.25 * std::exp(-std::pow(std::log(0.04 * E), 2) / 6.0);
}

// pdg — код частицы, вошедшей в трубку; E — её кинетическая энергия, МэВ
inline double Wr(int pdg, double E) {
    const int a = std::abs(pdg);
    if (a == 2112) return WrNeutron(E > 1e-12 ? E : 1e-12);
    if (a == 2212 || a == 211) return 2.0;
    if (a >= 1000000000) return 20.0;
    return 1.0;   // фотоны, e±, μ± и прочие
}

// класс источника для таблицы: 0 γ, 1 нейтрон, 2 e±, 3 μ±, 4 протон, 5 прочее
inline int DoseClass(int pdg) {
    const int a = std::abs(pdg);
    return a == 22 ? 0 : a == 2112 ? 1 : a == 11 ? 2 : a == 13 ? 3 : a == 2212 ? 4 : 5;
}
