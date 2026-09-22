#pragma once
#include "capture_emitter.hh"
#include <map>
#include <set>
#include <vector>
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4SteppingManager.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

class CaptureStep : public G4UserSteppingAction {
public:
  explicit CaptureStep(const CaptureEmitter* emitter) : emitter(emitter) {}

  void UserSteppingAction(const G4Step* step) override {
    const G4StepPoint* post = step->GetPostStepPoint();
    const G4VProcess* process = post->GetProcessDefinedStep();
    if (!process || process->GetProcessName() != "nCapture" ||
        step->GetTrack()->GetDefinition() != G4Neutron::Neutron()) {
      return;
    }

    G4TrackVector* sec = fpSteppingManager->GetfSecondary();
    if (!sec) return;

    // fpSteppingManager->GetfSecondary() копит вторичные ВСЕГО трека (в т.ч. ядра отдачи упругих рассеяний
    // предыдущих шагов, Fe-56 в железе) — берём только рождённые на ЭТОМ шаге.
    const std::vector<const G4Track*>* cur = step->GetSecondaryInCurrentStep();
    std::set<const G4Track*> thisStep;
    if (cur) thisStep.insert(cur->begin(), cur->end());

    G4Track* residualIon = nullptr;
    for (auto it = sec->begin(); it != sec->end(); ++it) {
      if (!thisStep.count(*it)) continue;
      const G4ParticleDefinition* def = (*it)->GetDefinition();
      if (def->GetParticleType() == "nucleus" ||
          (def->GetAtomicNumber() > 0 && def->GetAtomicMass() > 0)) {
        residualIon = *it;
        break;
      }
    }

    if (!residualIon) {
      ++nNoIon;
      return;
    }

    const G4ParticleDefinition* ionDef = residualIon->GetDefinition();
    int Z = ionDef->GetAtomicNumber();
    int A_product = ionDef->GetAtomicMass();
    int A_target = A_product - 1;

    ++keyCount[Z * 1000 + A_target];
    Cascade c;
    bool ok = emitter->Sample(Z, A_target, []{ return G4UniformRand(); }, c);
    if (!ok) {
      ++nFallback;
      return;
    }

    // Удаляем ВСЕ вторичные штатного захвата, кроме остаточного ядра: не только гамма, но и конверсионные
    // электроны/рентген штатного каскада (там SetICM(true)) — иначе часть энергии считалась бы дважды.
    for (auto it = sec->begin(); it != sec->end();) {
      if (*it != residualIon && thisStep.count(*it)) {
        delete *it;
        it = sec->erase(it);
      } else {
        ++it;
      }
    }

    // Штатная модель оставляет остаточное ядро ВОЗБУЖДЁННЫМ (Fe57[14.4], Fe57[136.5]); каскад эмиттера уже унёс всю
    // энергию до основного состояния, поэтому ядро заменяется основным — иначе его распад посчитал бы её дважды.
    const G4Ions* ion = dynamic_cast<const G4Ions*>(ionDef);
    if (ion && ion->GetExcitationEnergy() > 0.0) {
      auto* gdp = new G4DynamicParticle(G4IonTable::GetIonTable()->GetIon(Z, A_product, 0.0),
                                        residualIon->GetMomentumDirection(), residualIon->GetKineticEnergy());
      auto* gt = new G4Track(gdp, residualIon->GetGlobalTime(), residualIon->GetPosition());
      gt->SetTouchableHandle(residualIon->GetTouchableHandle());
      gt->SetParentID(residualIon->GetParentID());
      gt->SetCreatorProcess(residualIon->GetCreatorProcess());
      for (auto it = sec->begin(); it != sec->end(); ++it) if (*it == residualIon) { *it = gt; break; }
      delete residualIon;
      residualIon = gt;
      ++nIonGround;
    }

    // Создаем новые гамма-треки
    for (double E_keV : c.gamma_keV) {
      auto* dp = new G4DynamicParticle(G4Gamma::Gamma(), G4RandomDirection(), E_keV * keV);
      auto* t = new G4Track(dp, post->GetGlobalTime(), post->GetPosition());
      t->SetTouchableHandle(post->GetTouchableHandle());
      t->SetParentID(step->GetTrack()->GetTrackID());
      t->SetCreatorProcess(post->GetProcessDefinedStep());
      sec->push_back(t);
      ++nGammasEmitted;
      sumEmittedKeV += E_keV;
    }

    ++nReplaced;
  }

  std::map<int, long long> keyCount;   // диагностика: сколько захватов на каждом ядре-мишени (Z*1000+A)
  long long nReplaced = 0;
  long long nFallback = 0;
  long long nNoIon = 0;
  long long nIonGround = 0;   // сколько остаточных ядер заменено с возбуждённого на основное состояние
  long long nGammasEmitted = 0;
  double sumEmittedKeV = 0.0;

private:
  const CaptureEmitter* emitter;
};
