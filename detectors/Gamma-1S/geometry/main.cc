// Отклик спектрометра ГАММА-1С: спектр энерговыделения в кристалле NaI(Tl)
// 63x63 внутри свинцового экрана-защиты.
//
// Запуск:  g1s <макрос> [режим] [плотность]
//   shield — экран собран, крышка закрыта (по умолчанию)
//   open   — экран собран, крышка снята (точечная геометрия 25 см)
//   bare   — устройство детектирования в воздухе, без экрана
//   vessel[:сосуд] — экран закрыт, на головке сосуд комплекта:
//            marinelli (по умолчанию) | denta | petri
//
// Позиционные аргументы: 3 — плотность матрицы, 4 — матрица
// (OISN16 | water | risn379), 5 — объём засыпки мл, 6 — плотность MgO,
// 7 — глубина колодца маринелли мм.
#include "G1SDetector.hh"

// Отпечаток исходников, запечённый в бинарник (provenance.cmake генерирует его
// перед каждой сборкой). Идёт в шапку каждого выходного спектра: без него
// вопрос «этот спектр посчитан ТЕКУЩЕЙ геометрией?» отвечался по mtime, и
// трижды был отвечен неверно.
#if defined(__has_include)
#  if __has_include("g1s_provenance.hh")
#    include "g1s_provenance.hh"
#  endif
#endif
#ifndef G1S_SRC_SHA1
// Сборка мимо CMake — не запрещаем, но помечаем, чтобы такой спектр нельзя
// было принять за прослеженный.
#  define G1S_SRC_SHA1 "БЕЗ-ШТАМПА"
#  define G1S_GIT_DESCRIBE "БЕЗ-ШТАМПА"
#endif

#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4GeneralParticleSource.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4Run.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

class PhysList : public G4VModularPhysicsList {
public:
  PhysList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());   // с флуоресценцией:
    RegisterPhysics(new G4DecayPhysics());                // нужна для ХРИ Pb/Cd/Cu
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    SetDefaultCutValue(0.05 * mm);
  }
};

class Primary : public G4VUserPrimaryGeneratorAction {
  G4GeneralParticleSource fGPS;
public:
  void GeneratePrimaries(G4Event* e) override { fGPS.GeneratePrimaryVertex(e); }

  // Доля телесного угла ФАКТИЧЕСКОГО розыгрыша, (1−cos θmax)/2, спрошенная у
  // самого генератора после исполнения макроса. Прежде это число вычислял
  // драйвер по своей таблице и складывал в отдельный файл рядом с данными:
  // прямой множитель на каждую точечную eps лежал вне цепочки провенанса, и
  // никакая сверка не заметила бы, что прогон шёл с другим углом.
  double SolidAngleFrac() {
    auto* src = fGPS.GetCurrentSource();
    if (!src || !src->GetAngDist()) return 1.0;
    const double th = src->GetAngDist()->GetMaxTheta();
    return 0.5 * (1.0 - std::cos(th));
  }
};

// Время разрешения тракта: энерговыделения, разнесённые больше чем на столько,
// считаются РАЗНЫМИ срабатываниями спектрометра. Подробно — перед EventAct.
constexpr double kResolvingTimeNs = 1000.0;

// --- Накопление спектра -----------------------------------------------------
// 1024 канала — как у самого спектрометра (паспорт, п. 2.12), но шкала здесь
// линейная по энергии и без уширения: разрешение навешивается в постобработке.
class RunAct : public G4UserRunAction {
public:
  // 1 кэВ на канал. Потолок 3700, а не 3200: сетка энергий доходит до
  // 3552,5 кэВ (верхний край паспортной зоны «Денты»), и при потолке 3200
  // пик полного поглощения двух верхних узлов уезжал в канал переполнения.
  // Молчали обе стороны: прогон честно писал файл, а export_curves отбрасывал
  // точку с нулевой площадью без единого слова. Потолок обязан быть выше
  // САМОЙ ВЕРХНЕЙ энергии сетки плюс запас на сумм-пики.
  static constexpr int kBins = 3700;
  static constexpr double kBinKeV = 1.0;

  std::vector<long> fHist;
  // Спектр ИСПУЩЕННЫХ при распаде гамма-квантов: сколько квантов каждой
  // энергии рождается на один распад. Нужен для поправки на каскадное
  // суммирование — чтобы выход линии p_gamma брать из ТОЙ ЖЕ базы
  // PhotonEvaporation, что и транспорт, а не вписывать из справочника руками.
  std::vector<long> fEmit;
  long fWithSignal = 0;
  double fSumEprim = 0;
  G4String fPart = "?";
  G4String fMode = "shield";
  // Фактические параметры прогона — в шапку выхода. Без них штамп отвечает
  // только «из каких исходников собран exe», и два спектра, различающиеся
  // глубиной колодца или матрицей, выглядят одинаково прослеженными.
  G4String fArgs = "?";
  // Доля телесного угла розыгрыша (1−cos θ)/2 при конусе, иначе 1. Прямой
  // множитель на eps конусных сеток; сообщает его тот, кто разыгрывал.
  // Спрашивается в EndOfRunAction, а не при настройке: угол задаётся макросом,
  // то есть ПОСЛЕ создания действий, и опрос до /run/beamOn вернул бы значение
  // по умолчанию — молча и правдоподобно.
  double fSolidAngleFrac = 1.0;
  Primary* fPrimary = nullptr;
  G4String fOut = "spectrum.csv";

  RunAct() : fHist(kBins + 1, 0), fEmit(kBins + 1, 0) {}

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fHist.begin(), fHist.end(), 0L);
    std::fill(fEmit.begin(), fEmit.end(), 0L);
    fWithSignal = 0;
    fSumEprim = 0;
    fPart = "?";
  }

  void FillEmit(double eKeV) {
    if (eKeV <= 0) return;
    int b = static_cast<int>(eKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fEmit[b];
  }

  // Энергия первичной частицы учитывается ОДИН раз на событие Geant4, а Fill()
  // вызывается по разу на каждое РАЗДЕЛЁННОЕ ВО ВРЕМЕНИ срабатывание внутри
  // события (см. EventAct). Раньше это было одной функцией, и при переходе к
  // разделению E_prim_keV множился бы на число срабатываний.
  void FillPrimary(double eprim) { fSumEprim += eprim; }

  void Fill(double edepKeV) {
    if (edepKeV <= 0) return;
    ++fWithSignal;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fHist[b];
  }

  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (N == 0) return;
    if (fPrimary) fSolidAngleFrac = fPrimary->SolidAngleFrac();
    FILE* f = std::fopen(fOut.c_str(), "w");
    if (!f) {
      G4cerr << "!! не открыть " << fOut << G4endl;
      return;
    }
    std::fprintf(f, "# GAMMA-1S, UDS-GC-63x63-USB, NaI(Tl) 63x63 mm\n");
    // src_sha1 — отпечаток main.cc + G1SDetector.cc/.hh, из которых собран ЭТОТ
    // exe; build — время компиляции main.cc (заголовок провенанса меняется при
    // любой правке геометрии, поэтому main.cc заведомо перекомпилируется).
    std::fprintf(f, "# src_sha1 = %s\n", G1S_SRC_SHA1);
    std::fprintf(f, "# git_describe = %s\n", G1S_GIT_DESCRIBE);
    std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
    std::fprintf(f, "# mode = %s\n", fMode.c_str());
    // Фактические параметры прогона, а НЕ только исходники. Отпечаток отвечает
    // «из каких исходников собран exe»; два спектра с одним src_sha1 и разной
    // глубиной колодца или разной матрицей неразличимы, а глубина колодца в
    // этом же файле названа главным подозреваемым. Найдено независимым
    // аудитом: печатался mode и ни один из позиционных аргументов.
    std::fprintf(f, "# run_args = %s\n", fArgs.c_str());
    // Доля телесного угла розыгрыша: конусные сетки приводятся делением на неё,
    // то есть это ПРЯМОЙ множитель на каждую точечную eps. Прежде он лежал в
    // отдельном файле, который писал драйвер, и в цепочку провенанса не входил
    // вовсе. Здесь его сообщает сам exe, разыгравший события.
    std::fprintf(f, "# solid_angle_frac = %.8f\n", fSolidAngleFrac);
    std::fprintf(f, "# particle = %s\n", fPart.c_str());
    std::fprintf(f, "# E_prim_keV = %.4f\n", fSumEprim / N);
    std::fprintf(f, "# N_primaries = %ld\n", N);
    // Срабатываний, а не событий Geant4: в прогоне по цепочке одно первичное
    // ядро даёт несколько разнесённых во времени срабатываний, поэтому
    // N_with_signal может превысить N_primaries. Это не ошибка, а счёт
    // импульсов на распад родителя — именно так задана паспортная активность.
    std::fprintf(f, "# N_with_signal = %ld\n", fWithSignal);
    std::fprintf(f, "# resolving_time_ns = %.0f\n", kResolvingTimeNs);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n", kBinKeV);
    std::fprintf(f, "E_keV,counts\n");
    for (int i = 0; i <= kBins; ++i)
      if (fHist[i]) std::fprintf(f, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fHist[i]);
    std::fclose(f);

    // Спектр испускания — отдельным файлом, если он не пуст (прогон распада)
    long emitted = 0;
    for (long c : fEmit) emitted += c;
    if (emitted > 0) {
      G4String en = fOut;
      const size_t dot = en.rfind('.');
      en = (dot == G4String::npos ? en : en.substr(0, dot)) + "_emit.csv";
      FILE* g = std::fopen(en.c_str(), "w");
      if (g) {
        std::fprintf(g, "# гамма, испущенные при распаде, на %ld распадов\n", N);
        // Штамп и здесь: из _emit.csv берутся выход линии и мера чистоты, то
        // есть он такой же вход пересчёта, как сам спектр.
        std::fprintf(g, "# src_sha1 = %s\n", G1S_SRC_SHA1);
        std::fprintf(g, "# git_describe = %s\n", G1S_GIT_DESCRIBE);
        std::fprintf(g, "# build = %s %s\n", __DATE__, __TIME__);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "E_keV,counts\n");
        for (int i = 0; i <= kBins; ++i)
          if (fEmit[i]) std::fprintf(g, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fEmit[i]);
        std::fclose(g);
      }
      G4cout << "EMIT всего " << emitted << " квантов на " << N
             << " распадов -> " << en << G4endl;
    }

    G4cout << "RESULT E_keV= " << fSumEprim / N << " N= " << N
           << " hits= " << fWithSignal
           << " eff_total= " << double(fWithSignal) / N
           << " file= " << fOut << G4endl;
  }
};

// Событие Geant4 — это НЕ событие спектрометра. Когда розыгрыш идёт по цепочке
// (chain_Ra226, chain_Th232, точка p5_Th228), Geant4 доводит весь ряд до конца
// внутри одного события: порог долгого распада поднят до 1e30 нс, иначе
// долгоживущие звенья вообще не распались бы. В итоге в одном событии
// оказываются кванты Ac-228, Tl-208, Bi-212 — ядер, которые в действительности
// распадаются с разницей в годы. Складывая их энерговыделения, мы получали
// совпадения между РАЗНЫМИ нуклидами: пики обеднялись, а сумма уходила в
// континуум. Спектрометр так себя не ведёт.
//
// Поэтому энерговыделения собираются с отметкой глобального времени и потом
// разбиваются на группы: разрыв больше TAU — это уже другое срабатывание.
// TAU = 1 мкс, обычное время разрешения тракта. Значение некритично: настоящие
// каскады приходят за пикосекунды-наносекунды, а звенья ряда разделены
// секундами и годами, между этими масштабами шесть порядков пустоты.
// Единственные промежуточные — Po-212 (300 нс) и Po-214 (164 мкс), но оба
// альфа-излучатели: до кристалла из пробы они не долетают.
class EventAct : public G4UserEventAction {
  RunAct* fRun;
public:
  // (глобальное время, кэВ) по каждому шагу в кристалле
  std::vector<std::pair<double, double>> fDep;
  explicit EventAct(RunAct* r) : fRun(r) {}
  void BeginOfEventAction(const G4Event*) override { fDep.clear(); }
  void EndOfEventAction(const G4Event* e) override {
    double ep = 0;
    if (e->GetNumberOfPrimaryVertex() > 0) {
      auto* p = e->GetPrimaryVertex(0)->GetPrimary(0);
      ep = p->GetKineticEnergy() / keV;
      if (fRun->fPart == "?" && p->GetParticleDefinition())
        fRun->fPart = p->GetParticleDefinition()->GetParticleName();
    }
    fRun->FillPrimary(ep);
    if (fDep.empty()) return;

    // Шаги приходят в порядке обработки треков, а не по времени, поэтому
    // сортировка обязательна: без неё группировка развалится на первом же
    // треке, начавшемся раньше предыдущего.
    std::sort(fDep.begin(), fDep.end());
    double sum = fDep[0].second, t0 = fDep[0].first;
    for (size_t i = 1; i < fDep.size(); ++i) {
      if (fDep[i].first - t0 > kResolvingTimeNs) {
        fRun->Fill(sum);
        sum = 0;
      }
      t0 = fDep[i].first;
      sum += fDep[i].second;
    }
    fRun->Fill(sum);
  }
};

class Stepping : public G4UserSteppingAction {
  EventAct* fEvt;
  const G4LogicalVolume* fCry;
public:
  Stepping(EventAct* ev, const G4LogicalVolume* c) : fEvt(ev), fCry(c) {}
  void UserSteppingAction(const G4Step* s) override {
    const double e = s->GetTotalEnergyDeposit();
    if (e <= 0) return;
    auto* h = s->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (h && h->GetLogicalVolume() == fCry)
      fEvt->fDep.emplace_back(s->GetPreStepPoint()->GetGlobalTime() / ns,
                              e / keV);
  }
};

// Учёт квантов, РОЖДЁННЫХ радиоактивным распадом (а не рассеянием). Так выход
// линии p_gamma приходит из той же базы PhotonEvaporation, что и транспорт.
class Tracking : public G4UserTrackingAction {
  RunAct* fRun;
public:
  explicit Tracking(RunAct* r) : fRun(r) {}
  void PreUserTrackingAction(const G4Track* t) override {
    if (t->GetDefinition() != G4Gamma::Definition()) return;
    const G4VProcess* p = t->GetCreatorProcess();
    if (!p) return;                       // первичная частица, не распад
    const G4String& n = p->GetProcessName();
    if (n == "RadioactiveDecay" || n == "Radioactivation")
      fRun->FillEmit(t->GetKineticEnergy() / keV);
  }
};

class OutMessenger : public G4UImessenger {
  RunAct* fRun;
  G4UIdirectory* fDir;
  G4UIcmdWithAString* fCmd;
public:
  explicit OutMessenger(RunAct* r) : fRun(r) {
    fDir = new G4UIdirectory("/g1s/");
    fDir->SetGuidance("ГАММА-1С: управление выводом");
    fCmd = new G4UIcmdWithAString("/g1s/outFile", this);
    fCmd->SetGuidance("Файл CSV для спектра следующего прогона");
    fCmd->AvailableForStates(G4State_Idle, G4State_PreInit);
  }
  ~OutMessenger() override { delete fCmd; delete fDir; }
  void SetNewValue(G4UIcommand*, G4String v) override { fRun->fOut = v; }
};

int main(int argc, char** argv) {
  const std::string mode = (argc > 2) ? argv[2] : "shield";

  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* det = new G1SDetector();
  det->fWithShield = (mode != "bare");
  det->fShield.lidClosed = (mode != "open");
  // vessel[:сосуд] — marinelli (по умолчанию) | denta | petri
  det->fWithVessel = (mode.rfind("vessel", 0) == 0);
  if (det->fWithVessel) {
    const size_t c = mode.find(':');
    det->fVessel = VesselGeom::Preset(c == std::string::npos ? "marinelli"
                                                             : mode.substr(c + 1));
  }
  if (argc > 3) det->fVessel.sampleDensity = std::atof(argv[3]);
  if (argc > 4) det->fVessel.sampleMatrix = argv[4];   // OISN16|OISN06|water|risn379
  // Объём засыпки можно переопределить: у источников комплекта разные массы
  // при одном номинальном объёме кюветы, см. kit_inventory.py.
  if (argc > 5) det->fVessel.sampleCm3 = std::atof(argv[5]);
  // Два параметра, вынесенные по итогам сверки с ЛСРМ (см. STATUS.md):
  //   6-й — насыпная плотность отражателя MgO (мягкий край кривой);
  //   7-й — глубина колодца маринелли (главный подозреваемый в превышении
  //         расчёта над измерением на объёмном источнике).
  if (argc > 6) det->fHead.mgoDensity = std::atof(argv[6]);
  if (argc > 7) det->fVessel.wellDepth = std::atof(argv[7]);

  rm->SetUserInitialization(det);
  rm->SetUserInitialization(new PhysList());
  auto* primary = new Primary();
  rm->SetUserAction(primary);

  auto* runAct = new RunAct();
  runAct->fMode = mode;
  // Все аргументы как есть, плюс env-флаг: перечислять поля по одному значит
  // забыть новое поле при следующей правке.
  {
    std::string a;
    for (int i = 1; i < argc; ++i) {
      if (i > 1) a += " ";
      // ТОЛЬКО ИМЯ ФАЙЛА, без каталогов. Драйверы передают макрос абсолютным
      // путём, и полный argv занёс бы путь конкретной машины в шапку каждого
      // спектра — а спектр может быть скопирован в репозиторий. Запрет на
      // локальные пути в публичном дереве действует и для служебных полей.
      std::string v(argv[i]);
      const size_t s = v.find_last_of("/\\");
      a += (s == std::string::npos) ? v : v.substr(s + 1);
    }
    const char* cgv = std::getenv("G1S_CORRELATED_GAMMA");
    a += "; G1S_CORRELATED_GAMMA=";
    a += (cgv ? cgv : "unset");
    runAct->fArgs = a;
  }
  runAct->fPrimary = primary;
  rm->SetUserAction(runAct);
  auto* evtAct = new EventAct(runAct);
  rm->SetUserAction(evtAct);
  auto* mess = new OutMessenger(runAct);

  // Угловые корреляции гамма-квантов каскада. По умолчанию в Geant4 выключены
  // (G4DeexPrecoParameters: fCorrelatedGamma = false), и включить их из
  // обычного макроса НЕЛЬЗЯ: это параметр деэксцитации, он принимается только
  // до инициализации, а макрос выполняется после неё. Команда в макросе
  // отвергается с «Illegal application state», причём КОД ВОЗВРАТА ОСТАЁТСЯ
  // НУЛЕВЫМ — прогон выглядит успешным, корреляции молча не включаются, и
  // сравнение «с флагом против без флага» показывает отсутствие эффекта по
  // причине, не имеющей отношения к физике.
  //
  // Отсюда переменная окружения, а не аргумент: позиционные аргументы заняты
  // и разобраны драйверами, а флаг нужен разово, под одну проверку.
  if (const char* cg = std::getenv("G1S_CORRELATED_GAMMA")) {
    if (std::string(cg) == "1") {
      G4UImanager::GetUIpointer()->ApplyCommand(
          "/process/had/deex/correlatedGamma true");
      G4cout << "SETUP correlatedGamma = true (до инициализации)" << G4endl;
    }
  }

  rm->Initialize();
  det->ReportMasses();
  rm->SetUserAction(new Stepping(evtAct, det->fCrystalLV));
  rm->SetUserAction(new Tracking(runAct));

  auto* ui = G4UImanager::GetUIpointer();
  if (argc > 1) ui->ApplyCommand(G4String("/control/execute ") + argv[1]);

  delete mess;
  delete rm;
  return 0;
}
