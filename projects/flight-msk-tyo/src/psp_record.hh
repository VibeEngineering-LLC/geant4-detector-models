#pragma once
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#pragma pack(push, 1)
struct PspRec {
  int32_t evt;      // event number inside the run
  int32_t pdg;      // PDG code
  float E_MeV;      // kinetic energy at entry
  float phi;        // azimuth of the entry point on the tube side, rad, in (-pi, pi]
  float y_cm;       // axial coordinate of the entry point relative to the tube centre, cm
  float ux, uy, uz; // unit direction of motion (cabin axes)
  float t_ns;       // global time
  float Ep_MeV;     // kinetic energy of the primary of the event
  float wp;         // z component of the primary direction (cosine; +1 = going up, -1 = going down)
  int16_t proc;     // creator process code, see ProcCode
  int16_t vol;      // creator volume code, see VolCode
  int16_t ip;       // source component index 0..6
  int16_t pad;      // always 0
};
#pragma pack(pop)
static_assert(sizeof(PspRec) == 52, "PspRec must be 52 bytes");

inline int ProcCode(const std::string& name) {
  static const char* const names[] = {
    "", "nCapture", "neutronInelastic", "hadElastic", "RadioactiveDecay",
    "annihil", "eBrem", "compt", "phot", "conv", "eIoni", "muIoni",
    "muBrems", "muPairProd", "protonInelastic", "photonNuclear", "ionIoni",
    "msc", "Rayl", "other"
  };
  for (int i = 0; i < 20; ++i)
    if (name == names[i]) return i;
  return 19; // "other"
}

inline int VolCode(const std::string& name) {
  static const char* const names[] = {
    "World", "Skin", "Blanket", "Trim", "Floor", "Pax", "Cargo", "Cabin", "Tube", "other"
  };
  for (int i = 0; i < 10; ++i)
    if (name == names[i]) return i;
  return 9; // "other"
}

inline const char* ProcName(int code) {
  static const char* const names[] = {
    "primary", "nCapture", "neutronInelastic", "hadElastic", "RadioactiveDecay",
    "annihil", "eBrem", "compt", "phot", "conv", "eIoni", "muIoni",
    "muBrems", "muPairProd", "protonInelastic", "photonNuclear", "ionIoni",
    "msc", "Rayl", "other"
  };
  if (code < 0 || code >= 20) return "?";
  return names[code];
}

inline const char* VolName(int code) {
  static const char* const names[] = {
    "World", "Skin", "Blanket", "Trim", "Floor", "Pax", "Cargo", "Cabin", "Tube", "other"
  };
  if (code < 0 || code >= 10) return "?";
  return names[code];
}

class PspWriter {
 public:
  bool Open(const std::string& path) {
    f_ = fopen(path.c_str(), "wb");
    if (!f_) return false;
    setvbuf(f_, nullptr, _IOFBF, 1 << 20);
    n_ = 0;
    return true;
  }
  void Write(const PspRec& r) {
    fwrite(&r, sizeof r, 1, f_);
    ++n_;
  }
  void Close() {
    if (f_) {
      fclose(f_);
      f_ = nullptr;
    }
  }
  long long Count() const { return n_; }
  ~PspWriter() { Close(); }
 private:
  FILE* f_ = nullptr;
  long long n_ = 0;
};

class PspReader {
 public:
  bool Open(const std::string& path) {
    f_ = fopen(path.c_str(), "rb");
    return f_ != nullptr;
  }
  bool Next(PspRec& r) {
    return fread(&r, sizeof r, 1, f_) == 1;
  }
  void Close() {
    if (f_) {
      fclose(f_);
      f_ = nullptr;
    }
  }
  ~PspReader() { Close(); }
 private:
  FILE* f_ = nullptr;
};
