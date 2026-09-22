You are a C++ engineer working with Geant4 11.4. Write ONE complete header file `parma_primary.hh`. Output ONLY the code (no markdown fences, no prose).

## Purpose
A Geant4 primary generator action that draws primary particles from a `ParmaSource` (header `parma_source.hh`, already written; interface below).

## Interface to write (exactly)
```cpp
#pragma once
#include "parma_source.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include <string>

class ParmaPrimary : public G4VUserPrimaryGeneratorAction {
 public:
  // src: not owned, must outlive this object; particle_name: Geant4 name (e.g. "neutron","proton","mu+","mu-","e-","e+","gamma");
  // radius_cm: radius of the target sphere passed to ParmaSource::Sample
  ParmaPrimary(const ParmaSource* src, const std::string& particle_name, double radius_cm);
  void GeneratePrimaries(G4Event* event) override;
};
```
Implement everything inline in the header.

## Behaviour
* Constructor: create a `G4ParticleGun(1)`; look the particle up with `G4ParticleTable::GetParticleTable()->FindParticle(particle_name)`; if it is null throw `std::runtime_error("unknown particle: " + particle_name)`; set the definition on the gun.
* `GeneratePrimaries`: `ParmaSource::Event ev = src->Sample(radius_cm, []{ return G4UniformRand(); });` then
  `gun.SetParticleEnergy(ev.e_MeV * MeV);`
  `gun.SetParticlePosition(G4ThreeVector(ev.x_cm, ev.y_cm, ev.z_cm) * cm);`
  `gun.SetParticleMomentumDirection(G4ThreeVector(ev.u, ev.v, ev.w));` then `gun.GeneratePrimaryVertex(event);`
* `G4UniformRand()` is declared in `Randomize.hh` — include it. Also include `G4SystemOfUnits.hh`, `G4ThreeVector.hh`, `G4ParticleTable.hh`, `G4Event.hh`, `<stdexcept>`.
* The particle gun must be a member (`G4ParticleGun fGun;` or a `std::unique_ptr`), members initialised in the constructor initializer list.
* Comments in Russian, short. No global state.
