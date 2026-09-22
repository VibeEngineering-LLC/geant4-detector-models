#pragma once
#include "psp_record.hh"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

class PspEvents {
 public:
  bool AddFile(const std::string& path) {
    paths_.push_back(path);
    return true;
  }

  bool NextEvent(std::vector<PspRec>& ev) {
    ev.clear();
    if (!Fill()) return false;
    const int32_t id = buf_.evt;
    ev.push_back(buf_);
    have_ = false;

    while (open_ && reader_.Next(buf_)) {
      if (buf_.evt != id) {
        have_ = true;
        break;
      }
      ev.push_back(buf_);
    }

    if (!have_ && !open_) {
      reader_.Close();
      open_ = false;
    }

    ++n_events_;
    return true;
  }

  long long EventsRead() const { return n_events_; }

 private:
  bool Fill() {
    if (have_) return true;

    while (true) {
      if (open_ && reader_.Next(buf_)) {
        have_ = true;
        return true;
      }

      if (open_) {
        reader_.Close();
        open_ = false;
      }

      if (next_file_ >= paths_.size()) return false;

      const std::string& path = paths_[next_file_++];
      if (reader_.Open(path)) {
        open_ = true;
      } else {
        std::fprintf(stderr, "Failed to open file: %s\n", path.c_str());
      }
    }
  }

  std::vector<std::string> paths_;
  size_t next_file_ = 0;
  PspReader reader_;
  bool open_ = false;
  bool have_ = false;
  PspRec buf_;
  long long n_events_ = 0;
};

template <class Rng>
bool ChooseGroup(const std::vector<PspRec>& ev, double ell_cm, double L_cm, Rng&& rand01,
                 std::vector<int>& idx, double& w, double& zd_cm, double half_win_cm = -1.0) {
  // half_win_cm: полуширина окна группы (по умолчанию ell/2); ell — ширина ядра выбора точки прибора; z_d ограничена так, чтобы окно целиком лежало в трубке
  const double h = half_win_cm > 0 ? half_win_cm : ell_cm / 2;
  idx.clear();
  w = 0;
  const int n = static_cast<int>(ev.size());
  if (n == 0) return false;

  const int i = std::min(n - 1, static_cast<int>(rand01() * n));
  zd_cm = ev[i].y_cm + (rand01() - 0.5) * ell_cm;

  if (std::abs(zd_cm) > L_cm / 2 - h) return false;

  const double r = h;
  for (int j = 0; j < n; ++j) {
    if (std::abs(ev[j].y_cm - zd_cm) <= r) {
      idx.push_back(j);
    }
  }

  int mk = 0;   // число записей в ядре выбора (±ell/2) — оно задаёт плотность точки прибора
  for (int j = 0; j < n; ++j) if (std::abs(ev[j].y_cm - zd_cm) <= ell_cm / 2) ++mk;
  w = (mk > 0) ? n * ell_cm / ((L_cm - 2 * h) * mk) : 0.0;
  return true;
}
