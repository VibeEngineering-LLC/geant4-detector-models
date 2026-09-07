#pragma once

#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"

class NpsmBenchStackingAction : public G4UserStackingAction {
public:
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;
};
