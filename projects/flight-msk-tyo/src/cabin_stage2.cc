#include "psp_record.hh"
#include "replay.hh"
#include "ntag.hh"
#include "phantom_det.hh"
#include "dose_weights.hh"
#include "parma_source.hh"      // только для режима прямого счёта (замыкающая проверка схемы)
#include "capture_step.hh"
#include "capture_emitter.hh"
#include "ASN16Detector.hh"

#include <G4RunManagerFactory.hh>
#include <G4RunManager.hh>
#include <G4PhysListFactory.hh>
#include <G4VModularPhysicsList.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4EmParameters.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VUserPrimaryParticleInformation.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4IonTable.hh>
#include <G4UserTrackingAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserEventAction.hh>
#include <G4Event.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>
#include <G4LogicalVolume.hh>
#include <G4DynamicParticle.hh>
#include <G4UImanager.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

static const double kEll_cm = std::getenv("CABIN_ELL_CM") ? std::atof(std::getenv("CABIN_ELL_CM")) : 30.0;   // окно ℓ; переменная нужна для проверки независимости результата от ℓ
static const double kWinHalf_cm = std::getenv("CABIN_WIN_CM") ? std::atof(std::getenv("CABIN_WIN_CM")) : -1.0;   // полуширина окна группы (по умолчанию ℓ/2); большие значения нужны для косых фотонов вдоль трубки
static const int kRepeat = std::getenv("CABIN_K") ? std::max(1, std::atoi(std::getenv("CABIN_K"))) : 1;   // розыгрышей положения прибора на одно событие этапа I (вес делится на K)
const double kL_cm = 300.0;
const double kTubeR_mm = 139.0;
const double kTrigNs = 1000.0;
const int kNBin = 10000;

struct Ctx {
    PspEvents events;
    vector<PspRec> ev, group;
    vector<int> idx;
    double w;
    double zd;
    bool exhausted;
    long long nDraws, nUsed, nRecs;
    bool direct = false;                           // режим CABIN_DIRECT: первичные из сферы вокруг прибора, без этапа I
    ParmaSource dsrc;
    double directR_cm = 10.0;
    long long directN = 0;
    int directPdg = 22;
    int repeatLeft = 0;                            // сколько розыгрышей осталось у текущего события этапа I
    double sumW;
    vector<int> rootOf;
    vector<double> tot, tot2, lig, lig2;
    struct Deposit { double t, e_keV, light_keV; int root; int tag; };
    vector<Deposit> deposits;                      // энерговыделения в кристалле текущего события
    map<int, vector<double>> catHist;              // категория -> гистограмма весов
    map<int, double> catSum;
    map<int, vector<double>> catHist2;             // категория -> сумма квадратов весов
    vector<int> ntagOf;                            // метка нейтронной реакции у трека (0 — нет); 1 + NtagKey
    vector<char> nPrim;                            // трек — первичный нейтрон
    map<int, vector<double>> ntagHist;             // метка -> гистограмма весов (события с нейтронным корнем)
    long long nTrig = 0, nBelow1keV = 0;
    double sumTrigW = 0;
    double kB, Tflight_s;
    const ASN16Detector* det;
    bool phantom = false;                          // режим CABIN_PHANTOM: фантом тканей вместо прибора, счёт дозы
    PhantomDet* pdet = nullptr;
    double massTot = 0;                            // масса шара, кг
    struct PDep { double t, e_MeV; int lay, root; };
    vector<PDep> pdeps;                            // энерговыделения в теле фантома в текущем событии
    double dabs[4][6] = {}, deq[4][6] = {}, deq2[4][6] = {};   // [слой: кожа, 10 мм, 30 мм, весь шар][класс источника]
} ctx;

class RootInfo : public G4VUserPrimaryParticleInformation {
public:
    explicit RootInfo(int i) : idx(i) {}
    void Print() const override {}
    int idx;
};

class Gen : public G4VUserPrimaryGeneratorAction {
public:
    void GeneratePrimaries(G4Event* e) override {
        ctx.group.clear();
        ctx.w = 0;
        if (ctx.direct) {   // прямой счёт: один первичный на событие, вес 1, кабинная система -> приборная (x,y,z)->(x,z,y)
            if (ctx.nDraws >= ctx.directN) { ctx.exhausted = true; G4RunManager::GetRunManager()->AbortRun(true); return; }
            ++ctx.nDraws; ++ctx.nUsed; ++ctx.nRecs; ctx.w = 1.0; ctx.sumW += 1.0;
            const auto sm = ctx.dsrc.Sample(ctx.directR_cm, []{ return G4UniformRand(); });
            PspRec r{}; r.pdg = ctx.directPdg; r.E_MeV = (float) sm.e_MeV; ctx.group.push_back(r);
            auto* v = new G4PrimaryVertex(sm.x_cm * 10.0, sm.z_cm * 10.0, sm.y_cm * 10.0, 0.0);
            auto* p = new G4PrimaryParticle(G4ParticleTable::GetParticleTable()->FindParticle(ctx.directPdg));
            p->SetMomentumDirection(G4ThreeVector(sm.u, sm.w, sm.v));
            p->SetKineticEnergy(sm.e_MeV * MeV);
            p->SetUserInformation(new RootInfo(0));
            v->SetPrimary(p); e->AddPrimaryVertex(v);
            return;
        }
        while (true) {
            if (ctx.repeatLeft <= 0) {
                if (!ctx.events.NextEvent(ctx.ev)) {
                    ctx.exhausted = true;
                    G4RunManager::GetRunManager()->AbortRun(true);
                    return;
                }
                ctx.repeatLeft = kRepeat;
            }
            --ctx.repeatLeft;
            ++ctx.nDraws;
            bool ok = ChooseGroup(ctx.ev, kEll_cm, kL_cm, []{ return G4UniformRand(); }, ctx.idx, ctx.w, ctx.zd, kWinHalf_cm);
            if (!ok) continue;
            ctx.w /= kRepeat;
            for (int k : ctx.idx) {
                ctx.group.push_back(ctx.ev[k]);
            }
            ++ctx.nUsed;
            ctx.nRecs += ctx.group.size();
            ctx.sumW += ctx.w;
            bool anyGood = false;
            for (size_t i = 0; i < ctx.group.size(); ++i) {
                const PspRec& r = ctx.group[i];
                G4PrimaryVertex* v = new G4PrimaryVertex(
                    kTubeR_mm * cos(r.phi), kTubeR_mm * sin(r.phi), (r.y_cm - ctx.zd) * 10,
                    r.t_ns * ns
                );
                G4ParticleDefinition* def;
                int pdg = r.pdg;
                if (pdg >= 1000000000) {
                    int Z = (pdg / 10000) % 1000;
                    int A = (pdg / 10) % 1000;
                    def = G4IonTable::GetIonTable()->GetIon(Z, A, 0.0);
                } else {
                    def = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
                }
                if (!def) {
                    cerr << "Unknown particle PDG=" << pdg << endl;
                    continue;
                }
                anyGood = true;
                G4PrimaryParticle* p = new G4PrimaryParticle(def);
                p->SetMomentumDirection(G4ThreeVector(r.ux, r.uz, r.uy));   // кабина (x,y,z) -> прибор (X,Y,Z)=(x,z,y)
                p->SetKineticEnergy(r.E_MeV * MeV);
                p->SetUserInformation(new RootInfo(i));
                v->SetPrimary(p);
                e->AddPrimaryVertex(v);
            }
            if (anyGood) return;
        }
    }
};

class TrackingAction : public G4UserTrackingAction {
public:
    void PreUserTrackingAction(const G4Track* t) override {
        int id = t->GetTrackID();
        while (ctx.rootOf.size() <= id) ctx.rootOf.push_back(-1);
        while (ctx.ntagOf.size() <= id) { ctx.ntagOf.push_back(0); ctx.nPrim.push_back(0); }
        const int par = t->GetParentID();
        int tg = 0;
        if (par > 0) {
            tg = ctx.ntagOf[par];
            if (tg == 0 && ctx.nPrim[par]) {   // первое рождение вторичных нейтроном: процесс и объём реакции
                const G4VProcess* cp = t->GetCreatorProcess();
                tg = 1 + NtagKey(NtagProc(cp ? cp->GetProcessName() : ""), NtagVol(t->GetLogicalVolumeAtVertex()->GetName()));
            }
        }
        ctx.ntagOf[id] = tg;
        ctx.nPrim[id] = (par == 0 && t->GetDefinition()->GetPDGEncoding() == 2112);
        if (t->GetParentID() == 0) {
            const RootInfo* info = dynamic_cast<const RootInfo*>(
                t->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation()
            );
            ctx.rootOf[id] = info ? info->idx : -1;
        } else {
            ctx.rootOf[id] = ctx.rootOf[t->GetParentID()];
        }
    }
};

class SteppingAction : public G4UserSteppingAction {
public:
    explicit SteppingAction(CaptureStep* c) : cap(c) {}
    void SetSteppingManagerPointer(G4SteppingManager* p) override {
        G4UserSteppingAction::SetSteppingManagerPointer(p);
        if (cap) cap->SetSteppingManagerPointer(p);
    }
    void UserSteppingAction(const G4Step* s) override {
        if (cap) cap->UserSteppingAction(s);
        const G4StepPoint* pre = s->GetPreStepPoint();
        if (ctx.phantom) {
            const G4LogicalVolume* plv = pre->GetPhysicalVolume() ? pre->GetPhysicalVolume()->GetLogicalVolume() : nullptr;
            for (size_t k = 0; k < ctx.pdet->lv.size(); ++k) {
                if (plv != ctx.pdet->lv[k]) continue;
                const double ed = s->GetTotalEnergyDeposit();
                if (ed > 0) ctx.pdeps.push_back({pre->GetGlobalTime() / ns, ed / MeV, (int) k, ctx.rootOf[s->GetTrack()->GetTrackID()]});
                break;
            }
            return;
        }
        if (!pre->GetPhysicalVolume() || pre->GetPhysicalVolume()->GetLogicalVolume() != ctx.det->fCrystalLV) return;
        double edep = s->GetTotalEnergyDeposit();
        if (edep <= 0) return;
        double e_keV = edep / keV;
        double light_keV = e_keV;
        const G4Track* trk = s->GetTrack();
        const int pdgAbs = abs(trk->GetDefinition()->GetPDGEncoding());
        if (ctx.kB > 0 && trk->GetDefinition()->GetPDGCharge() != 0 && pdgAbs != 11 && pdgAbs != 13) {
            double len = s->GetStepLength();
            if (len > 0) light_keV = e_keV / (1 + ctx.kB * (edep / MeV) / (len / cm));
        }
        ctx.deposits.push_back({pre->GetGlobalTime() / ns, e_keV, light_keV, ctx.rootOf[trk->GetTrackID()], ctx.ntagOf[trk->GetTrackID()]});
    }
private:
    CaptureStep* cap;
};

class EventAction : public G4UserEventAction {
public:
    static void PhantomEvent() {   // доза в слоях фантома от одного воспроизведённого события, взвешенная ctx.w
        double da[4][6] = {}, de[4][6] = {};
        for (const auto& d : ctx.pdeps) {
            int pg = 5; double wr = 1.0;
            if (d.root >= 0 && d.root < (int) ctx.group.size()) { const PspRec& rc = ctx.group[d.root]; pg = DoseClass(rc.pdg); wr = Wr(rc.pdg, rc.E_MeV); }
            const double sat = max(0.0, 1 - d.t / (ctx.Tflight_s * 1e9));
            const double J = d.e_MeV * 1.602176634e-13 * sat;
            const int lay = ctx.pdet->layerOf[d.lay];
            if (lay >= 0) { da[lay][pg] += J / ctx.pdet->massKg[d.lay]; de[lay][pg] += wr * J / ctx.pdet->massKg[d.lay]; }
            da[3][pg] += J / ctx.massTot; de[3][pg] += wr * J / ctx.massTot;
        }
        for (int l = 0; l < 4; ++l) for (int c = 0; c < 6; ++c) {
            ctx.dabs[l][c] += ctx.w * da[l][c]; ctx.deq[l][c] += ctx.w * de[l][c]; ctx.deq2[l][c] += ctx.w * ctx.w * de[l][c] * de[l][c];
        }
    }
    void BeginOfEventAction(const G4Event*) override {
        ctx.rootOf.assign(1000, -1);
        ctx.ntagOf.assign(1000, 0);
        ctx.nPrim.assign(1000, 0);
        ctx.deposits.clear();
        ctx.pdeps.clear();
    }

    void EndOfEventAction(const G4Event*) override {
        if (ctx.phantom) { PhantomEvent(); return; }
        auto& dep = ctx.deposits;
        if (dep.empty()) return;
        sort(dep.begin(), dep.end(), [](const Ctx::Deposit& a, const Ctx::Deposit& b) { return a.t < b.t; });

        for (size_t i = 0; i < dep.size();) {
            size_t j = i;
            while (j + 1 < dep.size() && dep[j + 1].t - dep[j].t <= kTrigNs) ++j;   // разрыв считается до ПРЕДЫДУЩЕГО энерговыделения
            double E = 0, Lg = 0, maxE = 0;
            int root = -1;
            const double t0 = dep[i].t;
            std::map<int, double> byRoot, byTag;
            for (size_t k = i; k <= j; ++k) {
                E += dep[k].e_keV;
                Lg += dep[k].light_keV;
                if (dep[k].root >= 0) byRoot[dep[k].root] += dep[k].e_keV;
                if (dep[k].tag > 0) byTag[dep[k].tag] += dep[k].e_keV;
            }
            for (const auto& kv : byRoot) if (kv.second > maxE) { maxE = kv.second; root = kv.first; }
            i = j + 1;

            const double sat = max(0.0, 1 - t0 / (ctx.Tflight_s * 1e9));
            const double wt = ctx.w * sat;
            if (E < 1) { ++ctx.nBelow1keV; continue; }
            if (wt <= 0) continue;

            const int k = min((int) E, kNBin);
            const int kl = min((int) Lg, kNBin);
            ++ctx.nTrig;
            ctx.sumTrigW += wt;
            ctx.tot[k] += wt;   ctx.tot2[k] += wt * wt;
            if (Lg >= 1) { ctx.lig[kl] += wt; ctx.lig2[kl] += wt * wt; }

            int pg = 5, proc = 19, vol = 9;   // категория «другое»
            if (root >= 0 && root < (int) ctx.group.size()) {
                const PspRec& rec = ctx.group[root];
                const int pdg = abs(rec.pdg);
                pg = pdg == 22 ? 0 : pdg == 2112 ? 1 : pdg == 11 ? 2 : pdg == 13 ? 3 : pdg == 2212 ? 4 : 5;
                proc = rec.proc; vol = rec.vol;
            }
            const int key = (pg * 20 + proc) * 10 + vol;
            auto& h = ctx.catHist[key];
            if (h.empty()) h.assign(kNBin + 1, 0.0);
            h[k] += wt;
            auto& h2 = ctx.catHist2[key];
            if (h2.empty()) h2.assign(kNBin + 1, 0.0);
            h2[k] += wt * wt;
            ctx.catSum[key] += wt;
            if (root >= 0 && root < (int) ctx.group.size() && abs(ctx.group[root].pdg) == 2112 && !byTag.empty()) {   // нейтронные события: по реакции в приборе
                int bt = 0; double bm = -1;
                for (const auto& kv : byTag) if (kv.second > bm) { bm = kv.second; bt = kv.first; }
                auto& hn = ctx.ntagHist[bt];
                if (hn.empty()) hn.assign(kNBin + 1, 0.0);
                hn[k] += wt;
            }
        }
        dep.clear();
    }
};

int main(int argc, char* argv[]) {
    if (argc < 9) {
        cerr << "Usage: cabin_stage2 <out_prefix> <seed> <T_sim_s> <T_flight_s> <table 0|1> "
                "<kB_cm_per_MeV> <capture_db|-> <psp file> [<psp file> ...]" << endl;
        return 2;
    }

    string prefix = argv[1];
    int seed = stoi(argv[2]);
    double T_sim_s = stod(argv[3]);
    double T_flight_s = stod(argv[4]);
    int table = stoi(argv[5]);
    double kB = stod(argv[6]);
    string cap_db = argv[7];
    vector<string> psp_files;
    for (int i = 8; i < argc; ++i) {
        psp_files.push_back(argv[i]);
    }

    G4Random::setTheSeed(seed);

    auto* physList = new G4PhysListFactory();
    auto* pl = physList->GetReferencePhysList(std::getenv("CABIN_PHYS") ? std::getenv("CABIN_PHYS") : "FTFP_BERT_HPT");   // HPT = HP + тепловое рассеяние S(alpha,beta)
    if (!pl) { cerr << "нет физ-листа" << endl; return 2; }
    pl->ReplacePhysics(new G4EmStandardPhysics_option4());
    pl->RegisterPhysics(new G4RadioactiveDecayPhysics());

    G4EmParameters* emp = G4EmParameters::Instance();
    emp->SetFluo(true);
    emp->SetAuger(true);
    emp->SetPixe(std::getenv("CABIN_PIXE") ? std::atoi(std::getenv("CABIN_PIXE")) != 0 : true);   // PIXE включён (полная физика); при падениях — повтор с новым зерном (stage2_run.ps1)
    emp->SetDeexcitationIgnoreCut(true);

    auto* det = new ASN16Detector();
    det->fGeom.table = (table == 1);
    if (const char* cry = std::getenv("CABIN_CRY")) {   // размеры кристалла «X,Y,Z» мм (другой прибор в том же корпусе)
        double a = 0, b = 0, c = 0;
        if (std::sscanf(cry, "%lf,%lf,%lf", &a, &b, &c) == 3) { det->fGeom.cryX = a; det->fGeom.cryY = b; det->fGeom.cryZ = c; }
        cout << "CABIN_CRY crystal mm: " << det->fGeom.cryX << " x " << det->fGeom.cryY << " x " << det->fGeom.cryZ << endl;
    }

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    PhantomDet* pdet = std::getenv("CABIN_PHANTOM") ? new PhantomDet() : nullptr;   // режим дозы: фантом вместо прибора
    if (pdet) runManager->SetUserInitialization(pdet); else runManager->SetUserInitialization(det);
    runManager->SetUserInitialization(pl);

    Gen* gen = new Gen();
    runManager->SetUserAction(gen);

    TrackingAction* tracking = new TrackingAction();
    runManager->SetUserAction(tracking);

    CaptureStep* cap = nullptr;
    static CaptureEmitter emitter;   // должен пережить прогон
    if (cap_db != "-") {
        string err;
        if (!emitter.Load(cap_db, err)) {
            cerr << "Capture load error: " << err << endl;
            return 2;
        }
        cap = new CaptureStep(&emitter);
    }
    runManager->SetUserAction(new SteppingAction(cap));

    EventAction* event = new EventAction();
    runManager->SetUserAction(event);

    runManager->Initialize();

    G4UImanager* ui = G4UImanager::GetUIpointer();
    vector<string> cmds = {
        "/run/setCut 0.05 mm",
        "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+5 s",
        "/run/verbose 0",
        "/event/verbose 0",
        "/tracking/verbose 0"
    };
    for (const auto& cmd : cmds) {
        G4int status = ui->ApplyCommand(cmd);
        cout << "UI " << cmd << " -> " << status << endl;
    }

    for (const string& f : psp_files) {
        if (!ctx.events.AddFile(f)) {
            cerr << "Failed to add file: " << f << endl;
            return 2;
        }
    }

    ctx.det = det;
    ctx.pdet = pdet; ctx.phantom = (pdet != nullptr);
    if (pdet) for (double m : pdet->massKg) ctx.massTot += m;
    ctx.kB = kB;
    ctx.Tflight_s = T_flight_s;
    ctx.tot.assign(kNBin + 1, 0);
    ctx.tot2.assign(kNBin + 1, 0);
    ctx.lig.assign(kNBin + 1, 0);
    ctx.lig2.assign(kNBin + 1, 0);
    ctx.nDraws = 0;
    ctx.nUsed = 0;
    ctx.nRecs = 0;
    ctx.sumW = 0;
    ctx.exhausted = false;
    if (const char* d = std::getenv("CABIN_DIRECT")) {   // формат: <таблица>|<R_см>|<N>|<pdg>
        string spec(d); vector<string> f; size_t a = 0, b;
        while ((b = spec.find('|', a)) != string::npos) { f.push_back(spec.substr(a, b - a)); a = b + 1; }
        f.push_back(spec.substr(a));
        string err;
        if (f.size() != 4 || !ctx.dsrc.Load(f[0], err)) { cerr << "CABIN_DIRECT: <таблица>|<R_см>|<N>|<pdg>; " << err << endl; return 2; }
        ctx.direct = true; ctx.directR_cm = stod(f[1]); ctx.directN = stoll(f[2]); ctx.directPdg = stoi(f[3]);
    }

    runManager->BeamOn(2000000000);

    if (!ctx.exhausted) return 3;

    if (ctx.phantom) {   // доза: [слой 0 кожа 0,07 мм, 1 глубина 10 мм, 2 глубина 30 мм, 3 весь шар] x [класс источника 0 γ, 1 n, 2 e±, 3 μ±, 4 p, 5 прочее]
        ofstream od(prefix + "_dose.csv");
        od.precision(10);
        od << "layer,class,abs_Gy_per_s,eq_Sv_per_s,sigma_eq_Sv_per_s,mass_kg\n";
        const double invTd = 1.0 / T_sim_s;
        const double lm[4] = {ctx.pdet->massKg[0], ctx.pdet->massKg[2], ctx.pdet->massKg[4], ctx.massTot};
        for (int l = 0; l < 4; ++l) for (int c = 0; c < 6; ++c)
            od << l << "," << c << "," << ctx.dabs[l][c] * invTd << "," << ctx.deq[l][c] * invTd << "," << sqrt(ctx.deq2[l][c]) * invTd << "," << lm[l] << "\n";
    }

    // Write output files
    {
        ofstream out(prefix + "_total.csv");
        out << "# rates per second; sumw2 is the sum of squared weights divided by T_sim^2" << endl;
        out.precision(10);
        out << "bin_keV,sumw,sumw2,light_sumw,light_sumw2\n";
        double invT = 1.0 / T_sim_s;
        double invT2 = invT * invT;
        for (int i = 0; i <= kNBin; ++i) {
            out << i << "," << ctx.tot[i] * invT << "," << ctx.tot2[i] * invT2 << ","
                << ctx.lig[i] * invT << "," << ctx.lig2[i] * invT2 << "\n";
        }
    }

    {
        ofstream on(prefix + "_ntag.csv");   // нейтроны в приборе: реакция x объём, скорости с-1 на 1 кэВ
        on << "proc,vol,procname,volname,rate_total";
        for (int i = 0; i <= kNBin; ++i) on << ",bin_" << i;
        on << "\n";
        const double invTn = 1.0 / T_sim_s;
        for (const auto& p : ctx.ntagHist) {
            const int kk = p.first - 1; double sm = 0;
            for (double v : p.second) sm += v;
            on << kk / 5 << "," << kk % 5 << "," << NtagProcName(kk / 5) << "," << NtagVolName(kk % 5) << "," << sm * invTn;
            for (double v : p.second) on << "," << v * invTn;
            on << "\n";
        }
    }

    {
        ofstream out(prefix + "_cat.csv");
        out << "pg,proc,vol,pgname,procname,volname,rate_total";
        for (int i = 0; i <= kNBin; ++i) out << ",bin_" << i;
        out << "\n";

        vector<string> pgNames = {"gamma", "neutron", "e+-", "mu+-", "proton", "other"};
        const double invT = 1.0 / T_sim_s;
        const auto& catHist = ctx.catHist;
        auto& catSum = ctx.catSum;
        for (const auto& p : catHist) {
            int key = p.first;
            int pg = key / 200;
            int proc = (key % 200) / 10;
            int vol = key % 10;
            out << pg << "," << proc << "," << vol << ","
                << pgNames[pg] << "," << ProcName(proc) << "," << VolName(vol) << ","
                << catSum[key] * invT;
            for (int i = 0; i <= kNBin; ++i) {
                out << "," << p.second[i] * invT;
            }
            out << "\n";
        }
        ofstream o2(prefix + "_cat2.csv");   // sumw2 по категориям (для ошибок слоёв)
        o2 << "pg,proc,vol";
        for (int i = 0; i <= kNBin; ++i) o2 << ",bin_" << i;
        o2 << "\n";
        for (const auto& q : ctx.catHist2) {
            o2 << q.first / 200 << "," << (q.first % 200) / 10 << "," << q.first % 10;
            for (double v : q.second) o2 << "," << v * invT * invT;
            o2 << "\n";
        }
    }

    {
        ofstream out(prefix + "_meta.txt");
        out << "seed=" << seed << "\n"
            << "T_sim_s=" << T_sim_s << "\n"
            << "T_flight_s=" << T_flight_s << "\n"
            << "table=" << table << "\n"
            << "kB=" << kB << "\n"
            << "events_read=" << ctx.events.EventsRead() << "\n"
            << "draws=" << ctx.nDraws << "\n"
            << "used=" << ctx.nUsed << "\n"
            << "records_replayed=" << ctx.nRecs << "\n"
            << "sum_weight_replayed=" << ctx.sumW << "\n"
            << "triggers=" << ctx.nTrig << "\n" 
            << "sum_trigger_weight=" << ctx.sumTrigW << "\n"
            << "below_1keV=" << ctx.nBelow1keV << "\n"
            << "exhausted=" << (ctx.exhausted ? 1 : 0) << "\n";
        if (cap) {
            out << "cap_replaced=" << cap->nReplaced << "\n"
                << "cap_fallback=" << cap->nFallback << "\n"
                << "cap_no_ion=" << cap->nNoIon << "\n"
                << "cap_ion_ground=" << cap->nIonGround << "\n";
        }
    }

    const double invT = 1.0 / T_sim_s;
    double rate_total = 0;
    for (int i = 0; i < kNBin; ++i) {
        rate_total += ctx.tot[i];
    }
    cout << "STAGE2 draws=" << ctx.nDraws << " used=" << ctx.nUsed
         << " triggers=" << ctx.nTrig << " rate_total_cps=" << rate_total * invT << endl;

    return 0;
}
