#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#define _USE_MATH_DEFINES
#ifndef M_PI
constexpr double kPi = 3.14159265358979323846;
#else
constexpr double kPi = M_PI;
#endif

class ParmaSource {
 public:
  struct Event { 
    double e_MeV; 
    double u, v, w; 
    double x_cm, y_cm, z_cm; 
    double cos_down; 
  };

  bool Load(const std::string& path, std::string& err) {
    std::ifstream file(path);
    if (!file.is_open()) {
      err = "Не удалось открыть файл: " + path;
      return false;
    }

    std::string line;
    if (!std::getline(file, line)) {
      err = "Ошибка чтения строки 1";
      return false;
    }

    if (!std::getline(file, line)) {
      err = "Ошибка чтения строки 2";
      return false;
    }

    std::istringstream iss(line);
    int nebin, nabin, ie511;
    double flux511_cont, flux511_line, total_flux;

    if (!(iss >> nebin >> nabin >> ie511 >> flux511_cont >> flux511_line >> total_flux)) {
      err = "Ошибка чтения параметров из строки 2";
      return false;
    }

    if (nebin <= 0 || nabin <= 0) {
      err = "nebin и nabin должны быть положительными";
      return false;
    }

    if (!std::isfinite(total_flux) || total_flux <= 0) {
      err = "total_flux должен быть положительным и конечным";
      return false;
    }

    edges_MeV_.resize(nebin + 1);
    if (!std::getline(file, line)) {
      err = "Ошибка чтения энергетических границ";
      return false;
    }

    std::istringstream iss2(line);
    for (int i = 0; i <= nebin; ++i) {
      if (!(iss2 >> edges_MeV_[i])) {
        err = "Ошибка чтения энергетических границ";
        return false;
      }
    }

    D_.resize(nebin, std::vector<double>(nabin));
    mass_k_.resize(nebin);
    etab_.resize(nebin);

    for (int k = 0; k < nebin; ++k) {
      if (!std::getline(file, line)) {
        err = "Ошибка чтения данных D[k][ia]";
        return false;
      }

      std::istringstream iss3(line);
      double sum = 0.0;
      for (int ia = 0; ia < nabin; ++ia) {
        if (!(iss3 >> D_[k][ia])) {
          err = "Ошибка чтения данных D[k][ia]";
          return false;
        }
        if (!std::isfinite(D_[k][ia]) || D_[k][ia] < 0) {
          err = "Значение D[k][ia] должно быть конечным и неотрицательным";
          return false;
        }
        sum += D_[k][ia];
      }

      double de = edges_MeV_[k+1] - edges_MeV_[k];
      mass_k_[k] = sum * de;

      if (ie511 > 0 && k + 1 == ie511) {
        mass_k_[k] += flux511_line;
      }
    }

    double total_mass = 0.0;
    for (double m : mass_k_) total_mass += m;
    if (total_mass <= 0) {
      err = "Сумма масс равна нулю";
      return false;
    }

    for (int k = 0; k < nebin; ++k) {
      mass_k_[k] /= total_mass;
    }

    etab_[0] = mass_k_[0];
    for (int k = 1; k < nebin; ++k) {
      etab_[k] = etab_[k-1] + mass_k_[k];
    }

    atab_.resize(nebin, std::vector<double>(nabin));
    ahigh_.resize(nabin + 1);  // границы косинуса: 0..nabin (nabin+1 значение)
    for (int ia = 0; ia <= nabin; ++ia) {
      ahigh_[ia] = -1.0 + (2.0 / nabin) * ia;
    }

    for (int k = 0; k < nebin; ++k) {
      atab_[k][0] = D_[k][0];
      for (int ia = 1; ia < nabin; ++ia) {
        atab_[k][ia] = atab_[k][ia-1] + D_[k][ia];
      }

      double sum = atab_[k][nabin-1];
      if (sum > 0) {
        for (int ia = 0; ia < nabin; ++ia) {
          atab_[k][ia] /= sum;
        }
      }
    }

    line_frac_ = 0.0;
    if (ie511 > 0 && mass_k_[ie511-1] > 0) {
      // mass_k_ уже нормирована на total_mass, а flux511_line абсолютный: приводим к одним единицам
      line_frac_ = flux511_line / (mass_k_[ie511-1] * total_mass);
    }

    nebin_ = nebin;
    nabin_ = nabin;
    ie511_ = ie511;
    total_flux_ = total_flux;

    return true;
  }

  double TotalFlux() const { return total_flux_; }
  int NEBin() const { return nebin_; }
  int NABin() const { return nabin_; }
  int IE511() const { return ie511_; }
  double LineFrac() const { return line_frac_; }  // доля линии 511 кэВ внутри ячейки ie511
  const std::vector<double>& EdgesMeV() const { return edges_MeV_; }

  double EnergyBinProb(int k) const {
    if (k < 1 || k > nebin_) return 0.0;
    return mass_k_[k-1];
  }

  double CosBinProbMarginal(int ia) const {
    if (ia < 1 || ia > nabin_) return 0.0;
    double sum = 0.0;
    for (int k = 0; k < nebin_; ++k) {
      // условная вероятность ячейки ia при энергии k = D[k][ia] / сумма по ia
      double rs = 0.0;
      for (int j = 0; j < nabin_; ++j) rs += D_[k][j];
      if (rs > 0.0) sum += mass_k_[k] * D_[k][ia-1] / rs;
    }
    return sum;
  }

  template <class Rng>
  Event Sample(double radius_cm, Rng&& rand01) const {
    double r = rand01();
    int k = static_cast<int>(std::upper_bound(etab_.begin(), etab_.end(), r) - etab_.begin());
    if (k >= nebin_) k = nebin_ - 1;

    double r2 = rand01();
    double e = edges_MeV_[k] * r2 + edges_MeV_[k+1] * (1.0 - r2);

    if (ie511_ > 0 && k + 1 == ie511_) {
      double r3 = rand01();
      if (r3 < line_frac_) {
        e = 0.51099895;
      }
    }

    r = rand01();
    int ia = static_cast<int>(std::upper_bound(atab_[k].begin(), atab_[k].end(), r) - atab_[k].begin());
    if (ia >= nabin_) ia = nabin_ - 1;

    double r2a = rand01();
    double cx = ahigh_[ia] * r2a + ahigh_[ia + 1] * (1.0 - r2a);  // ячейка ia (0-based) = [ahigh[ia], ahigh[ia+1]]

    double phi = 2.0 * kPi * (rand01() - 0.5);
    double xd, yd, zd;
    do {
      xd = (rand01() - 0.5) * 2.0 * radius_cm;
      yd = (rand01() - 0.5) * 2.0 * radius_cm;
    } while (xd*xd + yd*yd > radius_cm * radius_cm);
    zd = radius_cm;

    double sx = std::sqrt(1.0 - cx * cx);

    double x = xd * cx * std::cos(phi) - yd * std::sin(phi) + zd * sx * std::cos(phi);
    double y = xd * cx * std::sin(phi) + yd * std::cos(phi) + zd * sx * std::sin(phi);
    double z = -xd * sx + zd * cx;

    double u = -sx * std::cos(phi);
    double v = -sx * std::sin(phi);
    double w = -cx;

    return {e, u, v, w, x, y, z, cx};
  }

 private:
  int nebin_, nabin_, ie511_;
  double total_flux_, line_frac_;
  std::vector<double> edges_MeV_;
  std::vector<std::vector<double>> D_;
  std::vector<double> mass_k_;
  std::vector<double> etab_;
  std::vector<std::vector<double>> atab_;
  std::vector<double> ahigh_;
};
