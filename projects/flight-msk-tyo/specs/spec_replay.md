You are a C++17 engineer (MSVC, Windows). Write ONE complete header file `replay.hh` with NO Geant4 dependency. Output ONLY the code (no markdown fences, no prose). Comments in Russian, short.

## Purpose
Stage II of a two-stage simulation replays particles that entered a long tube in stage I (records `PspRec` from `psp_record.hh`, already written: fields `evt`, `y_cm` = axial coordinate of the entry point relative to the tube centre, etc.; classes `PspReader{Open(path), Next(PspRec&), Close()}`). All positions along the tube are equivalent, so for each stage-I event we choose a random detector position along the tube and replay the records that entered near it.

## Interface (exactly)
```cpp
#pragma once
#include "psp_record.hh"
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

// Reads events (groups of consecutive records with equal evt) from several .psp files in the order they were added.
class PspEvents {
 public:
  bool AddFile(const std::string& path);      // false if the file cannot be opened
  bool NextEvent(std::vector<PspRec>& ev);    // clears ev and fills it with the next event; false when all files are exhausted
  long long EventsRead() const;               // number of events returned so far
 private:
  // your members: list of paths, index of the current file, a PspReader for it, a one-record look-ahead buffer (bool + PspRec)
};

// Chooses the replay group of one event. ev must be non-empty.
// n = ev.size(); anchor i = min(n-1, (int)(rand01()*n)); zd = ev[i].y_cm + (rand01() - 0.5) * ell_cm;
// if |zd| > L_cm/2 return false (this draw contributes nothing);
// otherwise idx = indices j with |ev[j].y_cm - zd| <= ell_cm/2 (the anchor is always among them), m = idx.size(),
// w = n * ell_cm / (L_cm * m), return true.
// Exactly TWO rand01() calls per invocation, in this order: anchor, then offset.
template <class Rng>
bool ChooseGroup(const std::vector<PspRec>& ev, double ell_cm, double L_cm, Rng&& rand01,
                 std::vector<int>& idx, double& w, double& zd_cm);
```
Everything inline in the header.

## Behaviour details
* `PspReader::Open` returns `bool` and `PspReader` is default-constructible and not copyable; keep it as a plain member `PspReader reader_` (do NOT use `make_unique`, do NOT assign readers). Members: `std::vector<std::string> paths_; size_t next_file_ = 0; PspReader reader_; bool open_ = false; bool have_ = false; PspRec buf_; long long n_events_ = 0;` (`have_` = the look-ahead buffer `buf_` holds an unread record of the current or a later event).
* Private helper `bool Fill()`: if `have_` return true; loop: if `open_` and `reader_.Next(buf_)` then `have_ = true; return true`; if `open_` then `reader_.Close(); open_ = false;`; if `next_file_ >= paths_.size()` return false; `path = paths_[next_file_++]`; if `reader_.Open(path)` then `open_ = true` else print a message on stderr and continue the loop.
* `NextEvent` algorithm: `ev.clear(); if (!Fill()) return false; const int32_t id = buf_.evt; ev.push_back(buf_); have_ = false;` then loop: `if (!open_) break;` (the file that produced the last record was closed inside Fill — treat as end of the event only when the reader cannot deliver another record of the same file) — implement it exactly like this: `while (open_ && reader_.Next(buf_)) { if (buf_.evt != id) { have_ = true; break; } ev.push_back(buf_); }` — after the loop, if `have_` is false and the loop ended because `reader_.Next` returned false, call `reader_.Close(); open_ = false;` (so the next call opens the next file: events never merge across files); `++n_events_; return true;`.
* `NextEvent`: an event ends when the `evt` field changes OR the current file ends (event numbers restart in every file, so never merge records across files). Use one record look-ahead: read the first record of the next event into the buffer. Files are opened lazily (Open the next file only when the current one is exhausted); a file that fails to open at that moment is skipped with a message on stderr. Empty files are skipped.
* `EventsRead` counts only events actually returned.
* `ChooseGroup` must not modify `ev`; `idx` is cleared first; on `false` return leave `w = 0`.

## Why the weight is this (do not change)
The sampled position density is p(z) = m(z) / (n * ell) (each anchor contributes 1/(n*ell) over its window), so the unbiased weight for the average over positions z uniform on [-L/2, L/2] is w = 1/(L*p(z)) = n*ell/(L*m).

## Requirements
* Include every header you use. No global state. Compile with MSVC `cl /EHsc /std:c++17 /utf-8`.
