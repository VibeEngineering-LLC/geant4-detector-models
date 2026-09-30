#pragma once

#include "G4UserStackingAction.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Ions.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"

#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <iomanip>

namespace DecayTally {

    static inline bool Enabled() {
        const char* env = std::getenv("GS2020_DECAY_TALLY");
        if (!env || std::string(env).empty()) return false;
        return true;
    }

    static inline const std::string& OutPath() {
        static std::string path;
        const char* env = std::getenv("GS2020_DECAY_TALLY");
        if (env) path = env;
        else path.clear();
        return path;
    }

    // State variables
    static inline std::unordered_map<int, double> gIonExc;
    static inline std::vector<double> gPhot;
    static inline std::vector<double> gElec;
    static inline std::map<std::string, long long> gSignatures;
    static inline long long gTotal = 0;

    static inline void BeginEvent() {
        gIonExc.clear();
        gPhot.clear();
        gElec.clear();
    }

    static inline void EndEvent() {
        // Sort energies ascending
        std::sort(gPhot.begin(), gPhot.end());
        std::sort(gElec.begin(), gElec.end());

        // Build signature string
        std::ostringstream oss;
        oss << "P:";
        for (size_t i = 0; i < gPhot.size(); ++i) {
            if (i > 0) oss << ' ';
            oss << std::fixed << std::setprecision(3) << gPhot[i];
        }
        oss << ";E:";
        for (size_t i = 0; i < gElec.size(); ++i) {
            if (i > 0) oss << ' ';
            oss << std::fixed << std::setprecision(2) << gElec[i];
        }

        std::string sig = oss.str();
        gSignatures[sig]++;
        gTotal++;
    }

    static inline void Write(long long seed, int ionZ, int ionA) {
        const std::string& path = OutPath();
        std::ofstream ofs(path);
        if (!ofs.is_open()) {
            std::cerr << "DECAY_TALLY: cannot open " << path << std::endl;
            std::exit(3);
        }

        ofs << "# decay_tally v1" << std::endl;
        ofs << "ion_Z," << ionZ << std::endl;
        ofs << "ion_A," << ionA << std::endl;
        ofs << "seed," << seed << std::endl;
        ofs << "n_events," << gTotal << std::endl;
        ofs << "n_signatures," << gSignatures.size() << std::endl;
        ofs << "count;signature" << std::endl;

        for (const auto& entry : gSignatures) {
            ofs << entry.second << ";" << entry.first << std::endl;
        }

        ofs.close();

        std::cout << "DECAY_TALLY: events=" << gTotal 
                  << " signatures=" << gSignatures.size() 
                  << " file=" << path << std::endl;
    }

    class StackingAction : public G4UserStackingAction {
    public:
        G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override {
            // Rule 1: Primary
            if (t->GetParentID() == 0) {
                return fUrgent;
            }

            const G4ParticleDefinition* def = t->GetDefinition();

            // Rule 2: Nucleus (Ion)
            if (def->GetParticleType() == "nucleus") {
                double excitation_keV = static_cast<const G4Ions*>(def)->GetExcitationEnergy() / keV;
                gIonExc[t->GetTrackID()] = excitation_keV;
                return fUrgent;
            }

            // Rule 3: Other particles
            double E = t->GetKineticEnergy() / keV;

            if (def == G4Gamma::Definition()) {
                if (E > 20.0) {
                    gPhot.push_back(E);
                }
                return fKill;
            }

            if (def == G4Electron::Definition()) {
                // Check if emitted by an excited ion
                int parentId = t->GetParentID();
                auto it = gIonExc.find(parentId);
                if (it != gIonExc.end() && it->second > 0.0) {
                    if (E > 10.0) {
                        gElec.push_back(E);
                    }
                }
                return fKill;
            }

            // Any other particle: kill
            return fKill;
        }
    };

} // namespace DecayTally
