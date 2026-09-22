#pragma once
#include <string>
// Метка нейтронной реакции в приборе: процесс (0 др., 1 упругое, 2 неупругое, 3 захват) x объём (0 др., 1 кристалл, 2 отражатель, 3 корпус, 4 плата)
inline int NtagProc(const std::string& p) { return p == "hadElastic" ? 1 : p == "neutronInelastic" ? 2 : p == "nCapture" ? 3 : 0; }
inline int NtagVol(const std::string& v) {
    if (v == "Crystal") return 1;
    if (v == "PTFE" || v == "AlFoil") return 2;
    if (v == "Body" || v.rfind("Cap", 0) == 0) return 3;
    if (v == "PCB" || v.rfind("Pcb", 0) == 0) return 4;
    return 0;
}
inline int NtagKey(int p, int v) { return p * 5 + v; }
inline const char* NtagProcName(int p) { static const char* n[] = {"other", "elastic recoil", "inelastic", "capture"}; return n[p]; }
inline const char* NtagVolName(int v) { static const char* n[] = {"other", "crystal", "reflector", "housing", "pcb"}; return n[v]; }
