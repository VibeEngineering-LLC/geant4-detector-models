#ifndef EXT_SOURCE_HH
#define EXT_SOURCE_HH

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>

namespace ExtSource {

struct Config {
    bool on;
    double R;
    double cz;
    int hemi;
};

inline const Config& Get() {
    static Config cfg = []() {
        Config c;
        c.on = false;
        c.R = 0.0;
        c.cz = 0.0;
        c.hemi = 0;

        // GS2020_EXT_R_MM
        const char* r_env = std::getenv("GS2020_EXT_R_MM");
        if (r_env) {
            try {
                double val = std::stod(r_env);
                if (val > 0.0) {
                    c.on = true;
                    c.R = val * mm;
                } else {
                    c.on = false;
                    c.R = 0.0;
                }
            } catch (...) {
                fprintf(stderr, "FATAL: GS2020_EXT_R_MM=%s\n", r_env);   // #GS-44: громкий отказ, не тихое выключение
                std::exit(4);
            }
        }

        // GS2020_EXT_CZ_MM
        const char* cz_env = std::getenv("GS2020_EXT_CZ_MM");
        if (cz_env) {
            try {
                double val = std::stod(cz_env);
                c.cz = val * mm;
            } catch (...) {
                fprintf(stderr, "FATAL: GS2020_EXT_CZ_MM=%s\n", cz_env);
                std::exit(4);
            }
        }

        // GS2020_EXT_HEMI
        const char* hemi_env = std::getenv("GS2020_EXT_HEMI");
        if (hemi_env) {
            std::string s(hemi_env);
            if (s == "all") {
                c.hemi = 0;
            } else if (s == "up") {
                c.hemi = 1;
            } else if (s == "down") {
                c.hemi = -1;
            } else {
                fprintf(stderr, "FATAL: GS2020_EXT_HEMI=%s\n", hemi_env);
                std::exit(4);
            }
        }

        // Print config if on
        if (c.on) {
            const char* hemi_str = "all";
            if (c.hemi == 1) hemi_str = "up";
            else if (c.hemi == -1) hemi_str = "down";
            
            printf("ext_source R_mm=%.6f cz_mm=%.6f hemi=%s\n", c.R/mm, c.cz/mm, hemi_str);
        }

        return c;
    }();
    return cfg;
}

inline long long gN = 0;
inline long long gNear = 0;

inline void Sample(G4ThreeVector& pos, G4ThreeVector& dir) {
    const Config& cfg = Get();
    
    // Point on the sphere
    double u1 = G4UniformRand();
    double w;
    if (cfg.hemi == 0) {
        w = 2.0 * u1 - 1.0;
    } else if (cfg.hemi == 1) {
        w = u1;
    } else { // hemi == -1
        w = -u1;
    }

    double s = std::sqrt(std::max(0.0, 1.0 - w*w));
    double phi = CLHEP::twopi * G4UniformRand();
    
    G4ThreeVector rhat(s * std::cos(phi), s * std::sin(phi), w);
    G4ThreeVector center(0.0, 0.0, cfg.cz);
    pos = center + cfg.R * rhat;

    // Inward normal
    G4ThreeVector n = -rhat;

    // Cosine law w.r.t. n
    double cost = std::sqrt(G4UniformRand());
    double sint = std::sqrt(std::max(0.0, 1.0 - cost*cost));
    double psi = CLHEP::twopi * G4UniformRand();

    // Orthonormal basis around n
    G4ThreeVector e1 = n.orthogonal().unit();
    G4ThreeVector e2 = n.cross(e1).unit();

    dir = (cost * n + sint * (std::cos(psi) * e1 + std::sin(psi) * e2)).unit();

    // Self-check statistic
    double b = cfg.R * sint;
    ++gN;
    if (b < 0.25 * cfg.R) {
        ++gNear;
    }
}

inline void Report() {
    const Config& cfg = Get();
    if (!cfg.on) return;

    double expected = 0.25 * 0.25; // 0.0625
    double frac = (gN > 0) ? static_cast<double>(gNear) / gN : 0.0;
    double sigma = std::sqrt(expected * (1.0 - expected) / std::max(gN, 1LL));
    double pull = (frac - expected) / sigma;

    printf("ext_check n=%lld frac_b_lt_R4=%.6f expected=0.062500 pull=%.2f\n", gN, frac, pull);
}

inline void CheckEnclosure(G4VPhysicalVolume* world) {
    const Config& cfg = Get();
    if (!cfg.on) return;

    // World solid limits
    G4VSolid* worldSolid = world->GetLogicalVolume()->GetSolid();
    G4ThreeVector wmin, wmax;
    worldSolid->BoundingLimits(wmin, wmax);

    // Check if sphere fits in world
    G4ThreeVector center(0.0, 0.0, cfg.cz);
    for (int i = 0; i < 3; ++i) {
        double min_extent = center[i] - cfg.R;
        double max_extent = center[i] + cfg.R;
        if (min_extent <= wmin[i] || max_extent >= wmax[i]) {
            fprintf(stderr, "FATAL: ext sphere R=%.1f mm exceeds world\n", cfg.R/mm);
            std::exit(5);
        }
    }

    // Check daughters
    G4LogicalVolume* worldLV = world->GetLogicalVolume();
    int nDaughters = worldLV->GetNoDaughters();
    
    double dmax = 0.0;
    std::string max_daughter_name = "";

    for (int i = 0; i < nDaughters; ++i) {
        G4VPhysicalVolume* pv = worldLV->GetDaughter(i);
        G4VSolid* solid = pv->GetLogicalVolume()->GetSolid();
        
        G4ThreeVector pmin, pmax;
        solid->BoundingLimits(pmin, pmax);

        // 8 corners of the box
        for (int ix = 0; ix < 2; ++ix) {
            for (int iy = 0; iy < 2; ++iy) {
                for (int iz = 0; iz < 2; ++iz) {
                    G4ThreeVector c(
                        (ix == 0) ? pmin.x() : pmax.x(),
                        (iy == 0) ? pmin.y() : pmax.y(),
                        (iz == 0) ? pmin.z() : pmax.z()
                    );

                    G4ThreeVector p = pv->GetObjectRotationValue() * c + pv->GetObjectTranslation();
                    double d = (p - center).mag();
                    
                    if (d > dmax) {
                        dmax = d;
                        max_daughter_name = pv->GetName();
                    }
                }
            }
        }
    }

    printf("ext_enclosure max_corner_mm=%.2f daughter=%s R_mm=%.6f\n", dmax/mm, max_daughter_name.c_str(), cfg.R/mm);

    if (dmax >= cfg.R) {
        fprintf(stderr, "FATAL: ext sphere R=%.1f mm does not enclose %s (%.1f mm)\n", 
                cfg.R/mm, max_daughter_name.c_str(), dmax/mm);
        std::exit(5);
    }
}

} // namespace ExtSource

#endif
