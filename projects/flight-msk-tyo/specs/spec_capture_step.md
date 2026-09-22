You are a C++ engineer working with Geant4 11.4.2 (MSVC). Write ONE complete header file `capture_step.hh`. Output ONLY the code (no markdown fences, no prose). Comments in Russian, short.

## Purpose
A Geant4 stepping action that REPLACES the gamma secondaries produced by the standard neutron-capture process (`nCapture`) with a cascade sampled from `CaptureEmitter` (header `capture_emitter.hh`, already written, interface below). When the emitter cannot handle a capture (isotope unknown, or the draw falls into the uncovered fraction), the standard Geant4 secondaries are left untouched.

## Interface (exactly)
```cpp
#pragma once
#include "capture_emitter.hh"
#include "G4UserSteppingAction.hh"
class CaptureStep : public G4UserSteppingAction {
 public:
  explicit CaptureStep(const CaptureEmitter* emitter);
  void UserSteppingAction(const G4Step* step) override;
  long long nReplaced = 0;   // captures whose gammas were replaced
  long long nFallback = 0;   // captures left to Geant4 (emitter returned false)
  long long nNoIon = 0;      // captures where the residual ion was not found among the secondaries
  long long nGammasEmitted = 0;
  double sumEmittedKeV = 0.0;
};
```
Implement everything inline in the header. The emitter API you rely on: `template<class Rng> bool Sample(int Z, int A_target, Rng&& rand01, Cascade& out) const;` with `struct Cascade { std::vector<double> gamma_keV; double conv_keV; int primary_level; };`.

## Behaviour of `UserSteppingAction(const G4Step* step)`
1. Take `const G4StepPoint* post = step->GetPostStepPoint();`. Return immediately unless `post->GetProcessDefinedStep() != nullptr` and its `GetProcessName() == "nCapture"` and `step->GetTrack()->GetDefinition() == G4Neutron::Neutron()`.
2. `G4TrackVector* sec = fpSteppingManager->GetfSecondary();` (`fpSteppingManager` is a protected member of `G4UserSteppingAction`). If `sec == nullptr` return.
3. Find the residual ion among `*sec`: a track whose definition has `GetParticleType() == "nucleus"` (or `GetAtomicNumber() > 0` and `GetAtomicMass() > 0`); `Z = def->GetAtomicNumber()`, `A_product = def->GetAtomicMass()`, target mass number `A_target = A_product - 1`. If none found: `++nNoIon; return;`.
4. `Cascade c; bool ok = emitter->Sample(Z, A_target, []{ return G4UniformRand(); }, c);` — if `!ok`: `++nFallback; return;`.
5. Otherwise delete every gamma secondary from `*sec`: iterate with an index or iterator, for tracks whose definition is `G4Gamma::Gamma()` call `delete` on the track and erase the pointer from the vector (be careful with iterator invalidation).
6. For every `E_keV` in `c.gamma_keV` create a new track: `auto* dp = new G4DynamicParticle(G4Gamma::Gamma(), G4RandomDirection(), E_keV * keV); auto* t = new G4Track(dp, post->GetGlobalTime(), post->GetPosition());` then `t->SetTouchableHandle(post->GetTouchableHandle()); t->SetParentID(step->GetTrack()->GetTrackID()); t->SetCreatorProcess(post->GetProcessDefinedStep());` and `sec->push_back(t);`. Update `nGammasEmitted` and `sumEmittedKeV`.
7. `++nReplaced`.

## Requirements
* Include every Geant4 header you use: `G4UserSteppingAction.hh`, `G4Step.hh`, `G4StepPoint.hh`, `G4Track.hh`, `G4TrackVector.hh`, `G4SteppingManager.hh`, `G4Neutron.hh`, `G4Gamma.hh`, `G4DynamicParticle.hh`, `G4ParticleDefinition.hh`, `G4VProcess.hh`, `G4SystemOfUnits.hh`, `G4ThreeVector.hh`, `Randomize.hh`, `G4RandomDirection.hh`.
* Do not keep pointers across steps; no static state.
* Compile with MSVC `cl /EHsc /std:c++17 /utf-8` against Geant4 11.4.2 headers.
