Write ONE C++17 header file `decay_tally.hh`, no prose, no markdown fences, only code. Geant4 11.2 API. Include guards `#pragma once`.

PURPOSE: in a Geant4 radioactive-decay run (one parent ion at rest per event), record per event which decay products are emitted
and kill every non-ion secondary at birth, so there is NO transport. Aggregate identical event signatures in memory, write them at the end.

Provide:
```
namespace DecayTally {
  bool Enabled();                 // true if env GS2020_DECAY_TALLY is set and non-empty (value = output file path)
  const std::string& OutPath();   // that env value
  void BeginEvent();              // clears per-event buffers
  void EndEvent();                // builds signature string, ++counter in std::map<std::string,long long>; ++total events
  void Write(long long seed, int ionZ, int ionA);  // writes file (see format), called once at end of run
  class StackingAction : public G4UserStackingAction { G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*) override; };
}
```
All state: file-scope static variables inside the namespace (single-threaded run manager; no mutex needed).

ClassifyNewTrack(t) rules, in this order:
1. If t->GetParentID() == 0 (primary) -> return fUrgent.
2. Let def = t->GetDefinition(). If def->GetParticleType() == "nucleus" (ion, any excitation):
   record in `std::unordered_map<int,double> gIonExc` : gIonExc[t->GetTrackID()] = excitation energy in keV,
   obtained as `static_cast<const G4Ions*>(def)->GetExcitationEnergy()/keV` (include "G4Ions.hh"); return fUrgent
   (ions must live: excited daughters de-excite through the decay process).
3. Energy E = t->GetKineticEnergy()/keV.
   If def == G4Gamma::Definition(): if E > 20.0 push E to per-event vector `gPhot`. return fKill.
   If def == G4Electron::Definition(): if E > 10.0 AND gIonExc has key t->GetParentID() AND that value > 0.0
     (electron emitted by an EXCITED ion = conversion or Auger electron, not the beta electron) push E to `gElec`. return fKill.
   Any other particle (e+, neutrinos, alpha, etc.): return fKill.
Clear gIonExc in BeginEvent.

EndEvent: sort gPhot and gElec ascending. Signature = "P:" + energies of gPhot joined by ' ' each printed with "%.3f",
then ";E:" + energies of gElec joined by ' ' with "%.2f". Increment map[signature]; ++gTotal.

Write: open OutPath() for writing (std::ofstream); if it fails print "DECAY_TALLY: cannot open <path>" to std::cerr and std::exit(3).
Lines:
`# decay_tally v1`
`ion_Z,<ionZ>`  `ion_A,<ionA>`  `seed,<seed>`  `n_events,<gTotal>`  `n_signatures,<map size>`
`count;signature` header line, then one line per map entry: `<count>;<signature>`.
Also print to std::cout: "DECAY_TALLY: events=<gTotal> signatures=<n> file=<path>".

Includes needed: G4UserStackingAction.hh, G4Track.hh, G4Gamma.hh, G4Electron.hh, G4Ions.hh, G4SystemOfUnits.hh, G4VProcess.hh,
<map>, <unordered_map>, <vector>, <string>, <algorithm>, <fstream>, <iostream>, <cstdio>, <cstdlib>.
Use `static` or `inline` for all non-member functions and variables so the header may be included once in main.cc safely.
