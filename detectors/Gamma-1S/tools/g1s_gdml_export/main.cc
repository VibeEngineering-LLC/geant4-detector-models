// g1s_gdml_export: Gamma-1S geometry (G1SDetector) -> GDML, no run, no BeamOn. Setup as run_g1s_npsm/main.cc:319-341.
// Usage: g1s_gdml_export.exe out.gdml [shield=1|0] [vessel=marinelli|none|...] [matrix=X] [rho=g_cm3] [sample_cm3=V]
#include "G1SDetector.hh"
#include "G4GDMLParser.hh"
#include <cstdlib>
#include <map>
#include <string>
int main(int argc, char** argv) {
  if (argc < 2) { std::printf("usage: g1s_gdml_export out.gdml [key=value ...]\n"); return 2; }
  std::map<std::string, std::string> a{{"shield", "1"}, {"vessel", "marinelli"}};
  for (int i = 2; i < argc; ++i) { std::string s = argv[i]; auto p = s.find('='); if (p != std::string::npos) a[s.substr(0, p)] = s.substr(p + 1); }
  auto* d = new G1SDetector;
  d->fWithShield = (a["shield"] == "1"); d->fWithVessel = (a["vessel"] != "none");
  d->fShield.variant = "geo1_2026_09_15"; d->fShield.cu = 3.0; d->fShield.cd = 1.5;   // geo1 + liner as in main.cc:326-330
  if (d->fWithVessel) d->fVessel = VesselGeom::Preset(a["vessel"]);
  if (a.count("matrix")) d->fVessel.sampleMatrix = a["matrix"];
  if (a.count("rho")) d->fVessel.sampleDensity = std::atof(a["rho"].c_str());
  if (a.count("sample_cm3")) d->fVessel.sampleCm3 = std::atof(a["sample_cm3"].c_str());
  G4GDMLParser parser; parser.SetRegionExport(false); parser.Write(argv[1], d->Construct(), false);
  std::printf("G1S_GDML_EXPORT: %s shield=%s variant=geo1_2026_09_15 cu=3 cd=1.5 vessel=%s matrix=%s rho=%g sample_cm3=%g\n", argv[1],
              a["shield"].c_str(), a["vessel"].c_str(), d->fVessel.sampleMatrix.c_str(), d->fVessel.sampleDensity, d->fVessel.sampleCm3);
  return 0;
}
