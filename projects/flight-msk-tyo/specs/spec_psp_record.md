You are a C++17 engineer (MSVC, Windows). Write ONE complete header file `psp_record.hh` with NO Geant4 dependency (only the standard library). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short.

## Purpose
Binary record of one particle entering a "tube" in a Geant4 cabin simulation (phase space file, `.psp`), plus a writer, a reader and code tables for creator process and creator volume.

## Interface (exactly these names)
```cpp
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

int ProcCode(const std::string& name);   // process name -> code
int VolCode(const std::string& name);    // logical volume name -> code
const char* ProcName(int code);          // code -> name ("?" if out of range)
const char* VolName(int code);

class PspWriter {
 public:
  bool Open(const std::string& path);    // "wb"; false on failure
  void Write(const PspRec& r);           // buffered; counts records
  void Close();                          // flush + fclose; safe to call twice
  long long Count() const;
  ~PspWriter();
 private:
  FILE* f_ = nullptr;
  long long n_ = 0;
};

class PspReader {
 public:
  bool Open(const std::string& path);    // "rb"; false on failure
  bool Next(PspRec& r);                  // false at end of file
  void Close();
  ~PspReader();
 private:
  FILE* f_ = nullptr;
};
```
Everything must be defined inline in the header (use `inline` for the free functions).

## Code tables (exact)
Process codes: 0 `primary` (empty name = no creator), 1 `nCapture`, 2 `neutronInelastic`, 3 `hadElastic`, 4 `RadioactiveDecay`, 5 `annihil`, 6 `eBrem`, 7 `compt`, 8 `phot`, 9 `conv`, 10 `eIoni`, 11 `muIoni`, 12 `muBrems`, 13 `muPairProd`, 14 `protonInelastic`, 15 `photonNuclear`, 16 `ionIoni`, 17 `msc`, 18 `Rayl`, 19 `other` (any other name; also the fallback of ProcCode).
Volume codes by logical volume name: 0 `World`, 1 `Skin`, 2 `Blanket`, 3 `Trim`, 4 `Floor`, 5 `Pax`, 6 `Cargo`, 7 `Cabin`, 8 `Tube`, 9 `other` (fallback; empty name also 9).
Use two `static const char* const` arrays (inline function returning a reference to a local static array, so that it works in a header) and a linear search for the name -> code direction.

## Requirements
* `Write` buffers with `setvbuf(f_, nullptr, _IOFBF, 1 << 20)` set in Open; use `fwrite(&r, sizeof r, 1, f_)`.
* `Next` returns true only when a full 52-byte record was read (`fread(...) == 1`).
* `Close` must set `f_ = nullptr` after `fclose`; destructors call Close.
* Include every header you use.
