// Отклик AtomSpectra Nano 16 PRO: спектр энерговыделения в бруске CsI(Tl)
// 18 x 15 x 57 мм внутри алюминиевой экструзии.
//
// Запуск:  asn16 <макрос>
//
// Положение источника задаётся макросом через GPS. Опорные плоскости, от
// которых отмеряется расстояние, печатает сама программа при старте
// (ReportPlanes): для точечной геометрии «10 см от торца» это НАРУЖНАЯ
// поверхность передней крышки, а не грань кристалла.
#include "ASN16Detector.hh"

// Отпечаток исходников, запечённый в бинарник (provenance.cmake генерирует его
// перед каждой сборкой). Идёт в шапку каждого выходного спектра: без него
// вопрос «этот спектр посчитан ТЕКУЩЕЙ геометрией?» отвечался бы по mtime.
#if defined(__has_include)
#  if __has_include("asn16_provenance.hh")
#    include "asn16_provenance.hh"
#  endif
#endif
#ifndef ASN16_SRC_SHA1
#  define ASN16_SRC_SHA1 "БЕЗ-ШТАМПА"
#  define ASN16_GIT_DESCRIBE "БЕЗ-ШТАМПА"
#endif

#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4Ions.hh"
#include "G4UserStackingAction.hh"
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
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

class PhysList : public G4VModularPhysicsList {
public:
  PhysList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    // Порог рождения вторичных. У Гамма-1С стоит 0,05 мм и проверен на мягком
    // крае (задача 100). Здесь лицевой стек тоньше (0,571 г/см² против
    // торцевой колонны Гамма-1С), поэтому тот же порог тем более достаточен.
    SetDefaultCutValue(0.05 * mm);
  }
};

class Primary : public G4VUserPrimaryGeneratorAction {
  G4GeneralParticleSource fGPS;
public:
  void GeneratePrimaries(G4Event* e) override { fGPS.GeneratePrimaryVertex(e); }

  // Доля телесного угла ФАКТИЧЕСКОГО розыгрыша, (1−cos θmax)/2, спрошенная у
  // самого генератора после исполнения макроса: прямой множитель на точечную
  // eps при конусном розыгрыше. Сообщает её тот, кто разыгрывал, а не
  // сторонняя таблица драйвера.
  double SolidAngleFrac() {
    auto* src = fGPS.GetCurrentSource();
    if (!src || !src->GetAngDist()) return 1.0;
    const double th = src->GetAngDist()->GetMaxTheta();
    return 0.5 * (1.0 - std::cos(th));
  }
};

// Время разрешения тракта: энерговыделения, разнесённые больше чем на столько,
// считаются РАЗНЫМИ срабатываниями спектрометра (обоснование — как у Гамма-1С:
// каскад приходит за наносекунды, звенья ряда разделены секундами и годами).
constexpr double kResolvingTimeNs = 1000.0;

// РАЗЛОЖЕНИЕ ОТКЛИКА ПО КАНАЛАМ (задание автора библиотеки BecqMoni,
// 06.08.2026: «чтобы он в функции отклика раскладывал фотоны: комптон, пары,
// вылеты, тормозные — на все составляющие обязательно»).
//
// Каналы ВЗАИМОИСКЛЮЧАЮЩИЕ и в сумме дают полный спектр: событие попадает
// ровно в один канал. Правило приоритета — по тому, ЧТО унесло энергию из
// кристалла, потому что именно вылет определяет, куда событие уходит из пика:
//
//   1) было рождение пар -> канал по числу вылетевших аннигиляционных квантов
//      (пары часто сопровождаются и тормозным, и рентгеном; отдельные каналы
//      для этих сочетаний дробят статистику, не добавляя смысла);
//   2) иначе ушёл сам первичный квант после комптона -> одно- или
//      многократное рассеяние;
//   3) иначе вылетел характеристический рентген;
//   4) иначе вылетел тормозной квант;
//   5) иначе ничего не вылетело: фотоэффект или комптон с последующим
//      поглощением.
//
// Приоритет объявлен здесь, а не выводится из данных: любое другое правило
// даёт другое разложение при той же физике, и сравнивать разложения можно
// только при одинаковом правиле.
enum Chan {
  kChPhoto = 0,     // фотоэффект, ничего не вылетело
  kChComptFull,     // комптон(ы) и поглощение, ничего не вылетело
  kChComptEsc1,     // однократный комптон, квант ушёл
  kChComptEscN,     // многократный комптон, квант ушёл
  kChXrayEsc,       // вылет характеристического рентгена
  kChBremsEsc,      // вылет тормозного кванта
  kChPairFull,      // пары, оба аннигиляционных кванта поглощены
  kChPairEsc1,      // пары, вылетел один квант 511 кэВ
  kChPairEsc2,      // пары, вылетели оба
  kChExternal,      // первичный квант в кристалле не взаимодействовал:
                    // энергию принесла вторичная частица из корпуса, обёртки
                    // или крышки. На жёстких узлах это не мелочь — 8,7 % на
                    // 3000 кэВ, и в модели без корпуса такого канала нет вовсе
  kChOther,         // остаточный: не должен населяться, служит сторожем
  kNChan
};

static const char* const kChanName[kNChan] = {
  "photo", "compt_full", "compt_esc1", "compt_escN", "xray_esc",
  "brems_esc", "pair_full", "pair_esc1", "pair_esc2", "external", "other"
};

class RunAct : public G4UserRunAction {
public:
  // 1 кэВ на канал. Потолок поднят 3700 -> 4200 (06.08.2026): при розыгрыше
  // РАСПАДА в каскаде Tl-208 складываются 2614,5 + 583,2 + 510,8 + 277,4 =
  // 3985,9 кэВ, и при потолке 3700 весь этот хвост уезжал бы в канал
  // переполнения — то есть сумм-пики, ради которых каскад и считается,
  // пропадали бы из спектра. Шкала линейная и без уширения: приборное
  // разрешение навешивается в постобработке.
  static constexpr int kBins = 4200;
  static constexpr double kBinKeV = 1.0;

  std::vector<long> fHist;
  std::vector<long> fEmit;   // гамма, ИСПУЩЕННЫЕ при распаде: выход линии
  // ВКЛАД БЕТЫ. Отдельная гистограмма событий, в которых в кристалл ВОШЛА
  // заряженная частица извне (электрон или позитрон, рождённый не в
  // кристалле). Ответ «доходит ли бета до кристалла» даёт прогон, а не
  // рассуждение о пробеге: путь от стержня до кристалла складывается из акрила
  // пенала, дна корпуса и платы, и считать его на бумаге значит повторить
  // модель. Гистограмма ведётся отдельно от каналов, чтобы не менять правило
  // приоритета, по которому уже посчитано разложение отклика.
  std::vector<long> fBeta;
  long fBetaEvents = 0;
  // Разложение отклика по каналам. Канал ставится В МОМЕНТ СОБЫТИЯ по истории
  // процессов: из готового спектра его восстановить нельзя, форма к тому
  // времени уже сложена. Каналы взаимоисключающие и в сумме дают fHist —
  // это проверяется при записи.
  std::vector<std::vector<long>> fChan;
  long fWithSignal = 0;
  double fSumEprim = 0;
  G4String fPart = "?";
  G4String fArgs = "?";
  double fSolidAngleFrac = 1.0;
  Primary* fPrimary = nullptr;
  ASN16Detector* fDet = nullptr;
  G4String fOut = "spectrum.csv";

  RunAct() : fHist(kBins + 1, 0), fEmit(kBins + 1, 0), fBeta(kBins + 1, 0),
             fChan(kNChan, std::vector<long>(kBins + 1, 0)) {}

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fHist.begin(), fHist.end(), 0L);
    std::fill(fEmit.begin(), fEmit.end(), 0L);
    std::fill(fBeta.begin(), fBeta.end(), 0L);
    fBetaEvents = 0;
    for (auto& c : fChan) std::fill(c.begin(), c.end(), 0L);
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

  void FillPrimary(double eprim) { fSumEprim += eprim; }

  void Fill(double edepKeV, int chan = -1, bool betaIn = false) {
    if (edepKeV <= 0) return;
    ++fWithSignal;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fHist[b];
    if (chan >= 0 && chan < kNChan) ++fChan[chan][b];
    if (betaIn) { ++fBetaEvents; ++fBeta[b]; }
  }

  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (N == 0) return;
    if (fPrimary) fSolidAngleFrac = fPrimary->SolidAngleFrac();
    // Каталог вывода создаётся здесь, а не считается существующим: макросы
    // пишут в подкаталоги (spectra/, spectra_face/), и без этого прогон на
    // 25 узлов отрабатывал ЦЕЛИКОМ, возвращал ноль и не оставлял ни одного
    // файла — тихий отказ, найденный независимым аудитом 05.08.2026.
    {
      const std::size_t s = fOut.find_last_of("/\\");
      if (s != G4String::npos) {
        std::error_code ec;
        std::filesystem::create_directories(fOut.substr(0, s), ec);
      }
    }
    FILE* f = std::fopen(fOut.c_str(), "w");
    if (!f) {
      // Не предупреждение, а аварийный останов: прежде здесь стоял return, и
      // прогон завершался с кодом 0 без единого файла.
      G4Exception("RunAct::EndOfRunAction", "ASN16_OUT", FatalException,
                  ("не открыть файл вывода " + fOut).c_str());
      return;
    }
    // Габарит кристалла берётся ИЗ ГЕОМЕТРИИ, а не пишется строкой: зашитое
    // «18x15x57» пережило правку cryZ 57 -> 60 и попало неверным в шапку всех
    // 27 спектров прогона 93972a90d30c. Тот же класс дефекта, что уже ловили
    // в рисовальных скриптах, — число, пережившее условия, при которых было
    // верно (06.08.2026).
    if (fDet) {
      const Nano16Geom& gm = fDet->fGeom;
      std::fprintf(f, "# ATOMSPECTRA NANO 16 PRO, CsI(Tl) %gx%gx%g mm\n",
                   gm.cryX, gm.cryY, gm.cryZ);
    } else {
      std::fprintf(f, "# ATOMSPECTRA NANO 16 PRO, CsI(Tl) (геометрия не задана)\n");
    }
    std::fprintf(f, "# src_sha1 = %s\n", ASN16_SRC_SHA1);
    std::fprintf(f, "# git_describe = %s\n", ASN16_GIT_DESCRIBE);
    std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
    std::fprintf(f, "# run_args = %s\n", fArgs.c_str());
    // Опорные плоскости — в шапку спектра. Прибор анизотропен, и «10 см» без
    // указания, от какой поверхности и к какой грани, величины не задаёт.
    if (fDet) {
      std::fprintf(f, "# front_face_z_mm = %.3f  (наружная поверхность крышки)\n",
                   fDet->FrontFaceZ());
      std::fprintf(f, "# crystal_front_z_mm = %.3f\n", fDet->CrystalFrontZ());
      std::fprintf(f, "# work_face_y_mm = %.3f  (наружная поверхность стенки)\n",
                   fDet->WorkFaceY());
      std::fprintf(f, "# crystal_top_y_mm = %.3f\n", fDet->CrystalTopY());
      // Состояние фрезеровки — В ШАПКУ. Две кривые одной ревизии различаются
      // только этим, а опираться на имя каталога нельзя: оно задаётся в
      // макросе и переживает любую путаницу с копированием.
      std::fprintf(f, "# cap_window = %s\n",
                   fDet->fGeom.capWindow ? "on" : "off");
      std::fprintf(f, "# cap_in_beam_mm = %.2f  (алюминий крышки в пучке)\n",
                   fDet->fGeom.capWindow ? fDet->fGeom.wCapWin
                                         : fDet->fGeom.wCap);
      // Источник — в шапку наравне с детектором: шаблон нуклида без геометрии
      // источника величины не задаёт, а имя каталога переживает любую путаницу
      // с копированием.
      std::fprintf(f, "# wt20_pack = %s\n", fDet->fGeom.wt20 ? "on" : "off");
      if (fDet->fGeom.wt20) {
        const Nano16Geom& gp = fDet->fGeom;
        std::fprintf(f, "# wt20_rods = %d x %.2f x %.1f mm, pitch %.2f mm\n",
                     gp.wt20N, gp.wt20D, gp.wt20L, gp.wt20Pitch);
        std::fprintf(f, "# wt20_mass_g = %.2f  (вольфрамовый сплав, %.1f %% ThO2)\n",
                     fDet->PackMassG(), gp.wt20ThO2);
        std::fprintf(f, "# bottom_face_y_mm = %.3f  (дно корпуса — прибор лежит "
                        "на пенале)\n", fDet->BottomFaceY());
      }
    }
    std::fprintf(f, "# solid_angle_frac = %.8f\n", fSolidAngleFrac);
    std::fprintf(f, "# particle = %s\n", fPart.c_str());
    std::fprintf(f, "# E_prim_keV = %.4f\n", fSumEprim / N);
    std::fprintf(f, "# N_primaries = %ld\n", N);
    std::fprintf(f, "# N_with_signal = %ld\n", fWithSignal);
    std::fprintf(f, "# N_charged_entered = %ld  (событий, где в кристалл вошла "
                    "заряженная частица извне)\n", fBetaEvents);
    std::fprintf(f, "# resolving_time_ns = %.0f\n", kResolvingTimeNs);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n",
                 kBinKeV);
    std::fprintf(f, "E_keV,counts\n");
    for (int i = 0; i <= kBins; ++i)
      if (fHist[i]) std::fprintf(f, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fHist[i]);
    std::fclose(f);

    // --- разложение отклика по каналам, отдельным файлом ---
    // Отдельным, а не колонками основного спектра: формат «E_keV,counts»
    // читают все разборные скрипты дерева, и менять его ради нового
    // содержимого значит чинить их все.
    {
      G4String cn = fOut;
      const size_t dot = cn.rfind('.');
      cn = (dot == G4String::npos ? cn : cn.substr(0, dot)) + "_chan.csv";
      FILE* g = std::fopen(cn.c_str(), "w");
      if (g) {
        std::fprintf(g, "# разложение отклика по каналам, канал ставится по "
                        "истории процессов события\n");
        std::fprintf(g, "# src_sha1 = %s\n", ASN16_SRC_SHA1);
        std::fprintf(g, "# git_describe = %s\n", ASN16_GIT_DESCRIBE);
        std::fprintf(g, "# E_prim_keV = %.4f\n", fSumEprim / N);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "# solid_angle_frac = %.8f\n", fSolidAngleFrac);
        std::fprintf(g, "# bin_keV = %.3f\n", kBinKeV);
        std::fprintf(g, "# правило приоритета: пары -> вылет первичного после "
                        "комптона -> рентген -> тормозное -> без вылета\n");
        std::fprintf(g, "E_keV");
        for (int c = 0; c < kNChan; ++c) std::fprintf(g, ",%s", kChanName[c]);
        std::fprintf(g, "\n");
        for (int i = 0; i <= kBins; ++i) {
          if (!fHist[i]) continue;
          std::fprintf(g, "%.1f", (i + 0.5) * kBinKeV);
          for (int c = 0; c < kNChan; ++c)
            std::fprintf(g, ",%ld", fChan[c][i]);
          std::fprintf(g, "\n");
        }
        std::fclose(g);
      }
      // Каналы обязаны в сумме давать полный спектр. Если не дают — правило
      // приоритета пропускает случай, и это дефект разложения, а не мелочь:
      // молча потерянные события выглядели бы как отсутствие канала.
      long sumChan = 0, sumHist = 0;
      for (int i = 0; i <= kBins; ++i) {
        sumHist += fHist[i];
        for (int c = 0; c < kNChan; ++c) sumChan += fChan[c][i];
      }
      if (sumChan != sumHist)
        G4cerr << "ВНИМАНИЕ: сумма каналов " << sumChan
               << " не равна спектру " << sumHist
               << " — правило приоритета неполно" << G4endl;
    }

    // --- вклад беты, отдельным файлом ---------------------------------------
    if (fBetaEvents > 0) {
      G4String bn = fOut;
      const size_t dot = bn.rfind('.');
      bn = (dot == G4String::npos ? bn : bn.substr(0, dot)) + "_beta.csv";
      FILE* g = std::fopen(bn.c_str(), "w");
      if (g) {
        std::fprintf(g, "# события, в которых в кристалл вошла заряженная "
                        "частица извне (бета и вторичные электроны)\n");
        std::fprintf(g, "# src_sha1 = %s\n", ASN16_SRC_SHA1);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "# N_with_signal = %ld\n", fWithSignal);
        std::fprintf(g, "# N_charged_entered = %ld\n", fBetaEvents);
        std::fprintf(g, "# bin_keV = %.3f\n", kBinKeV);
        std::fprintf(g, "E_keV,counts_total,counts_charged_in\n");
        for (int i = 0; i <= kBins; ++i)
          if (fHist[i] || fBeta[i])
            std::fprintf(g, "%.1f,%ld,%ld\n", (i + 0.5) * kBinKeV, fHist[i],
                         fBeta[i]);
        std::fclose(g);
      }
      G4cout << "БЕТА: заряженная частица вошла в кристалл в " << fBetaEvents
             << " событиях из " << fWithSignal << " со сигналом ("
             << 100.0 * fBetaEvents / std::max(1L, fWithSignal) << " %)"
             << G4endl;
    }

    long emitted = 0;
    for (long c : fEmit) emitted += c;
    if (emitted > 0) {
      G4String en = fOut;
      const size_t dot = en.rfind('.');
      en = (dot == G4String::npos ? en : en.substr(0, dot)) + "_emit.csv";
      FILE* g = std::fopen(en.c_str(), "w");
      if (g) {
        std::fprintf(g, "# гамма, испущенные при распаде, на %ld распадов\n", N);
        std::fprintf(g, "# src_sha1 = %s\n", ASN16_SRC_SHA1);
        std::fprintf(g, "# git_describe = %s\n", ASN16_GIT_DESCRIBE);
        std::fprintf(g, "# build = %s %s\n", __DATE__, __TIME__);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "E_keV,counts\n");
        for (int i = 0; i <= kBins; ++i)
          if (fEmit[i]) std::fprintf(g, "%.1f,%ld\n", (i + 0.5) * kBinKeV,
                                     fEmit[i]);
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

// Событие Geant4 — это НЕ событие спектрометра: энерговыделения собираются с
// отметкой глобального времени и разбиваются на группы по времени разрешения.
class EventAct : public G4UserEventAction {
  RunAct* fRun;
public:
  std::vector<std::pair<double, double>> fDep;   // (нс, кэВ)

  // Признаки истории события для разложения по каналам. Ставятся в Stepping.
  int fFirst = 0;          // 1 phot, 2 compt, 3 conv, 0 ничего неупругого
  bool fHadRayl = false;   // было упругое рассеяние (энергии не оставляет)
  int fNCompt = 0;         // сколько раз первичный квант рассеялся в кристалле
  bool fHadConv = false;   // было рождение пар в кристалле
  int fNAnnihEsc = 0;      // сколько аннигиляционных квантов покинуло кристалл
  double fEBremEsc = 0;    // энергия вылетевшего тормозного, кэВ
  double fEXrayEsc = 0;    // энергия вылетевших прочих вторичных гамма, кэВ
  bool fPrimEsc = false;   // сам первичный квант вышел из кристалла
  bool fChargedIn = false; // в кристалл ВОШЛА заряженная частица извне

  explicit EventAct(RunAct* r) : fRun(r) {}
  void BeginOfEventAction(const G4Event*) override {
    fDep.clear();
    fFirst = 0; fHadRayl = false; fNCompt = 0; fHadConv = false;
    fNAnnihEsc = 0; fEBremEsc = 0; fEXrayEsc = 0; fPrimEsc = false;
    fChargedIn = false;
  }

  // Канал по правилу приоритета (см. enum Chan). Возвращает индекс канала.
  int Channel() const {
    if (fHadConv)
      return fNAnnihEsc == 0 ? kChPairFull
           : (fNAnnihEsc == 1 ? kChPairEsc1 : kChPairEsc2);
    if (fPrimEsc && fNCompt > 0)
      return fNCompt == 1 ? kChComptEsc1 : kChComptEscN;
    if (fEXrayEsc > 0) return kChXrayEsc;
    if (fEBremEsc > 0) return kChBremsEsc;
    if (fFirst == 1) return kChPhoto;
    if (fNCompt > 0) return kChComptFull;
    if (fFirst == 0) return kChExternal;
    return kChOther;
  }
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

    // Канал один на событие Geant4. При моноэнергетическом источнике группа
    // энерговыделений одна, и вопроса не возникает; при розыгрыше распада
    // событие может разбиться на несколько срабатываний, и тогда все они
    // получают канал первичного кванта — это огрубление, и оно объявлено.
    const int ch = Channel();

    std::sort(fDep.begin(), fDep.end());
    double sum = fDep[0].second, t0 = fDep[0].first;
    for (size_t i = 1; i < fDep.size(); ++i) {
      if (fDep[i].first - t0 > kResolvingTimeNs) {
        fRun->Fill(sum, ch, fChargedIn);
        sum = 0;
      }
      t0 = fDep[i].first;
      sum += fDep[i].second;
    }
    // Признак «вошла заряженная» — на СОБЫТИЕ, а не на срабатывание: если
    // событие распалось на несколько срабатываний, помечаются все. Огрубление
    // объявлено; при розыгрыше одиночного распада таких событий единицы.
    fRun->Fill(sum, ch, fChargedIn);
  }
};

class Stepping : public G4UserSteppingAction {
  EventAct* fEvt;
  const ASN16Detector* fDet;
public:
  // Держит ДЕТЕКТОР, а не указатель на его логический объём. Кэшированный
  // указатель переживает перестройку геометрии (/asn16/capWindow), но после
  // неё указывает на удалённый объём: сравнение не совпадает никогда, счёт
  // выходит НУЛЕВЫМ, а прогон завершается успешно и пишет пустые спектры.
  // Поймано прямым прогоном 06.08.2026: 0 событий против 504 на том же узле.
  Stepping(EventAct* ev, const ASN16Detector* d) : fEvt(ev), fDet(d) {}
  void UserSteppingAction(const G4Step* s) override {
    auto* pre = s->GetPreStepPoint();
    auto* post = s->GetPostStepPoint();
    auto* h = pre->GetTouchableHandle()->GetVolume();
    const bool inCry = h && fDet &&
                       h->GetLogicalVolume() == fDet->fCrystalLV;
    if (!inCry) return;

    // --- ВХОД ЗАРЯЖЕННОЙ ЧАСТИЦЫ ИЗВНЕ --------------------------------------
    // Признак ставится по вершине трека: частица рождена НЕ в кристалле и
    // пересекла его границу. Так отделяется бета из источника (и вторичный
    // электрон, выбитый ею в корпусе) от электронов, рождённых гамма-квантом
    // уже внутри кристалла, — иначе «вклад беты» получился бы равным почти
    // всему спектру, потому что энергию в CsI в любом случае несёт электрон.
    {
      const G4Track* t0 = s->GetTrack();
      const double q = t0->GetDefinition()->GetPDGCharge();
      if (q != 0.0 && pre->GetStepStatus() == fGeomBoundary &&
          t0->GetLogicalVolumeAtVertex() != fDet->fCrystalLV)
        fEvt->fChargedIn = true;
    }

    // --- пометка канала: что произошло с квантом ВНУТРИ кристалла ---
    const G4Track* trk = s->GetTrack();
    const bool isGamma = trk->GetDefinition() == G4Gamma::Gamma();
    if (isGamma) {
      const G4VProcess* pr = post->GetProcessDefinedStep();
      const G4String pn = pr ? pr->GetProcessName() : G4String();
      if (trk->GetParentID() == 0) {
        if (pn == "compt") {
          ++fEvt->fNCompt;
          if (!fEvt->fFirst) fEvt->fFirst = 2;
        } else if (pn == "phot") {
          if (!fEvt->fFirst) fEvt->fFirst = 1;
        } else if (pn == "conv") {
          fEvt->fHadConv = true;
          if (!fEvt->fFirst) fEvt->fFirst = 3;
        } else if (pn == "Rayl") {
          // Рэлеевское рассеяние УПРУГОЕ: энергии не оставляет и первым
          // взаимодействием не считается. Прежде оно помечалось первым, и
          // событие «рэлей, затем фотоэффект» уходило в остаточный канал —
          // 5,5 % на 180 кэВ, при том что рэлей сам по себе отсчёта не даёт.
          // Поймано разложением: канал, который физически обязан быть пуст,
          // оказался населён.
          fEvt->fHadRayl = true;
        }
      }
      // Выход кванта из кристалла. Квант, вернувшийся обратно из корпуса,
      // здесь уже посчитан вылетевшим — огрубление в пользу каналов вылета;
      // на жёстких узлах доля таких возвратов мала, но она НЕ измерена.
      if (post->GetStepStatus() == fGeomBoundary) {
        auto* hp = post->GetTouchableHandle()->GetVolume();
        const bool outCry = !hp ||
                            hp->GetLogicalVolume() != fDet->fCrystalLV;
        const double ek = post->GetKineticEnergy() / keV;
        if (outCry && ek > 0) {
          const G4VProcess* cp = trk->GetCreatorProcess();
          const G4String cn = cp ? cp->GetProcessName() : G4String();
          if (cn == "annihil") ++fEvt->fNAnnihEsc;
          else if (cn == "eBrem") fEvt->fEBremEsc += ek;
          else if (trk->GetParentID() == 0) fEvt->fPrimEsc = true;
          else fEvt->fEXrayEsc += ek;
        }
      }
    }

    // --- энерговыделение ---
    const double e = s->GetTotalEnergyDeposit();
    if (e > 0)
      fEvt->fDep.emplace_back(pre->GetGlobalTime() / ns, e / keV);
  }
};

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

// ОДНО ЗВЕНО РЯДА НА ПРОГОН. Шаблон отклика строится ОТДЕЛЬНО для каждого
// нуклида цепочки, поэтому розыгрыш обязан остановиться на первом распаде:
// иначе прогон «Pb-212» дал бы заодно спектр Bi-212 и Tl-208, и разложение
// измеренного спектра по нуклидам потеряло бы смысл.
//
// Убивается ДОЧЕРНЕЕ ЯДРО В ОСНОВНОМ СОСТОЯНИИ. Ядро в ВОЗБУЖДЁННОМ состоянии
// не трогается: именно его снятие даёт гамма-каскад распада, ради которого всё
// и считается. Убить всё подряд по признаку «ион» значило бы получить пустой
// спектр при исправном на вид прогоне — отказ того же класса, что уже ловили
// в этом дереве дважды.
//
// Проверка встроена: файл <спектр>_emit.csv считает ИСПУЩЕННЫЕ кванты на
// распад, и их выходы сверяются с библиотечными интенсивностями линий.
class Stacking : public G4UserStackingAction {
public:
  long fKilled = 0;
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override {
    if (t->GetParentID() == 0) return fUrgent;
    const auto* ion = dynamic_cast<const G4Ions*>(t->GetDefinition());
    if (ion && ion->GetExcitationEnergy() <= 0.0) { ++fKilled; return fKill; }
    return fUrgent;
  }
};

class OutMessenger : public G4UImessenger {
  RunAct* fRun;
  ASN16Detector* fDet;
  G4UIdirectory* fDir;
  G4UIcmdWithAString* fCmd;
  G4UIcmdWithAString* fWinCmd;
  G4UIcmdWithAString* fPackCmd;
  G4UIcmdWithAString* fTabCmd;
  G4UIcmdWithAString* fCryYCmd;
  G4UIcmdWithAString* fPcbCuCmd;
  G4UIcmdWithAString* fPcbSolCmd;
  G4UIcmdWithAString* fPcbCmpCmd;
  G4UIcmdWithAString* fPcbPbCmd;
  G4UIcmdWithAString* fCaseOffCmd;
  G4UIcmdWithAString* fRoomCmd;
  G4UIcmdWithAString* fRodDCmd;
public:
  OutMessenger(RunAct* r, ASN16Detector* d) : fRun(r), fDet(d) {
    fDir = new G4UIdirectory("/asn16/");
    fDir->SetGuidance("AtomSpectra Nano 16 PRO: управление выводом");
    fCmd = new G4UIcmdWithAString("/asn16/outFile", this);
    fCmd->SetGuidance("Файл CSV для спектра следующего прогона");
    fCmd->AvailableForStates(G4State_Idle, G4State_PreInit);
    // ФРЕЗЕРОВКА КАК ПАРАМЕТР. Заведена командой, а не второй сборкой, чтобы
    // кривая «с окном» и кривая «без окна» шли от ОДНОЙ ревизии исходников:
    // иначе у них разные src_sha1 и разность двух кривых нельзя приписать
    // одному только окну. Команда обязана стоять ДО /run/initialize.
    fWinCmd = new G4UIcmdWithAString("/asn16/capWindow", this);
    fWinCmd->SetGuidance("on|off — фрезеровка передней крышки напротив "
                         "кристалла (off: сплошная крышка wCap)");
    fWinCmd->SetCandidates("on off");
    // Доступна и в Idle: геометрия к моменту чтения макроса уже построена
    // (rm->Initialize() стоит в main до ApplyCommand), поэтому команда сама
    // просит перестроить её — см. SetNewValue.
    fWinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    // ПАЧКА WT-20 ПОД ПРИБОРОМ. Тем же приёмом, что фрезеровка: геометрия
    // источника — параметр одной ревизии, а не вторая сборка. Иначе спектр «с
    // пачкой» и опорная кривая «без пачки» шли бы от разных штампов и разность
    // нельзя было бы приписать одному только источнику.
    fPackCmd = new G4UIcmdWithAString("/asn16/wt20", this);
    fPackCmd->SetGuidance("on|off — пачка электродов WT-20 под дном прибора");
    fPackCmd->SetCandidates("on off");
    fPackCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    // Столешница под пачкой — отдельным переключателем: разность двух прогонов
    // одной ревизии показывает вклад рассеяния от подложки и ничего больше.
    fTabCmd = new G4UIcmdWithAString("/asn16/table", this);
    fTabCmd->SetGuidance("on|off — столешница под пеналом (рассеиватель)");
    fTabCmd->SetCandidates("on off");
    fTabCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    // ТОЛЩИНА КРИСТАЛЛА ПО ПУЧКУ — только для ДИАГНОСТИКИ. В постановке с
    // пачкой излучение входит через грань 18 x 60 и проходит cryY, тогда как
    // опорный замер Cs-137 просвечивал кристалл вдоль оси Z на 60 мм. Значит
    // размер 15 мм этой сверкой НИКОГДА не проверялся, а именно он задаёт, как
    // быстро падает эффективность с энергией. Команда позволяет измерить эту
    // чувствительность, а не рассуждать о ней. Публикуемые числа считаются
    // ТОЛЬКО при паспортных 15,00 мм.
    fCryYCmd = new G4UIcmdWithAString("/asn16/cryY", this);
    fCryYCmd->SetGuidance("толщина кристалла по Y, мм (диагностика!)");
    fCryYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    // НАПОЛНЕНИЕ ПЛАТЫ — тем же приёмом параметра одной ревизии. Реальная
    // разводка неизвестна, поэтому медь, припой и корпуса компонентов заданы
    // эффективной сплошной толщиной, а команды позволяют по ним СКАНИРОВАТЬ:
    // вклад свинца припоя в горб 75-95 кэВ меряется разностью прогонов, а не
    // назначается. Ноль отключает слой.
    fPcbCuCmd = new G4UIcmdWithAString("/asn16/pcbCu", this);
    fPcbCuCmd->SetGuidance("эффективная сплошная медь дорожек, мм");
    fPcbCuCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fPcbSolCmd = new G4UIcmdWithAString("/asn16/pcbSolder", this);
    fPcbSolCmd->SetGuidance("эффективный сплошной припой, мм");
    fPcbSolCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fPcbCmpCmd = new G4UIcmdWithAString("/asn16/pcbComp", this);
    fPcbCmpCmd->SetGuidance("эффективные корпуса компонентов (Al2O3), мм");
    fPcbCmpCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fPcbPbCmd = new G4UIcmdWithAString("/asn16/pcbSolderPb", this);
    fPcbPbCmd->SetGuidance("on|off — припой Sn63Pb37 (on) или SAC305 (off)");
    fPcbPbCmd->SetCandidates("on off");
    fPcbPbCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    // ДИАМЕТР СТЕРЖНЯ — только для ДИАГНОСТИКИ, как и /asn16/cryY. Этикетка
    // даёт 3,2 мм, фотография другой пачки того же типа — 2,4 мм, и разница
    // бьёт по мягкому концу спектра сильнее всего: тоньше пруток — меньше
    // самопоглощение, выше пик 238,63 и выход K-серии дочерних наружу.
    // Публикуемые числа считаются ТОЛЬКО при подтверждённом замером диаметре.
    // ВНИМАНИЕ: ось стержней зависит от диаметра (yRod = yInT − зазор − r),
    // поэтому область розыгрыша GPS обязана пересчитываться вместе с ним —
    // ловушка, на которой уже был забракован скан по толщине кристалла.
    fRodDCmd = new G4UIcmdWithAString("/asn16/rodD", this);
    fRodDCmd->SetGuidance("диаметр стержня WT-20, мм (ДИАГНОСТИКА!)");
    fRodDCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fRoomCmd = new G4UIcmdWithAString("/asn16/room", this);
    fRoomCmd->SetGuidance("on|off — бетонный пол под столом (рассеиватель)");
    fRoomCmd->SetCandidates("on off");
    fRoomCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fCaseOffCmd = new G4UIcmdWithAString("/asn16/caseOff", this);
    fCaseOffCmd->SetGuidance("on|off — корпус и крышки вакуумом (ДИАГНОСТИКА!)");
    fCaseOffCmd->SetCandidates("on off");
    fCaseOffCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  }
  ~OutMessenger() override {
    delete fRodDCmd; delete fRoomCmd; delete fCaseOffCmd;
    delete fPcbPbCmd; delete fPcbCmpCmd; delete fPcbSolCmd; delete fPcbCuCmd;
    delete fCryYCmd; delete fTabCmd; delete fPackCmd; delete fWinCmd;
    delete fCmd; delete fDir;
  }
  // Общая часть трёх команд по слоям платы: проверка, запись, перестройка.
  // Вынесена, потому что три копии одного кода — это три места, где потом
  // забудут поправить проверку.
  void SetPcbLayer(double& field, double want, const char* what) {
    if (want < 0.0 || want > 1.0) {
      std::fprintf(stderr, "ОТКАЗ: слой платы %s = %g вне (0…1 мм)\n",
                   what, want);
      return;
    }
    if (want == field) return;
    field = want;
    if (auto* rm = G4RunManager::GetRunManager()) {
      rm->ReinitializeGeometry(true);
      rm->GeometryHasBeenModified();
      std::printf("[плата] %s = %.4f мм, геометрия перестроена\n", what, want);
    }
  }
  void SetNewValue(G4UIcommand* c, G4String v) override {
    if (c == fCmd) { fRun->fOut = v; return; }
    if (fDet) {
      if (c == fPcbCuCmd) {
        SetPcbLayer(fDet->fGeom.pcbCuT, std::atof(v.c_str()), "медь"); return;
      }
      if (c == fPcbSolCmd) {
        SetPcbLayer(fDet->fGeom.pcbSnPbT, std::atof(v.c_str()), "припой");
        return;
      }
      if (c == fPcbCmpCmd) {
        SetPcbLayer(fDet->fGeom.pcbCompT, std::atof(v.c_str()), "корпуса");
        return;
      }
      if (c == fRodDCmd) {
        const double want = std::atof(v.c_str());
        if (want <= 0.5 || want > 6.0) {
          std::fprintf(stderr, "ОТКАЗ: /asn16/rodD %g вне разумного "
                               "(0,5…6 мм)\n", want);
          return;
        }
        if (want == fDet->fGeom.wt20D) return;
        fDet->fGeom.wt20D = want;
        if (auto* rm = G4RunManager::GetRunManager()) {
          rm->ReinitializeGeometry(true);
          rm->GeometryHasBeenModified();
          std::printf("[rodD] диаметр стержня %.2f мм, геометрия перестроена "
                      "— ЭТО ДИАГНОСТИКА, не этикетка. Ось стержней сместилась,"
                      " область розыгрыша GPS пересчитать!\n", want);
        }
        return;
      }
      if (c == fRoomCmd) {
        const bool want = (v == "on");
        if (want == fDet->fGeom.room) return;
        fDet->fGeom.room = want;
        if (auto* rm = G4RunManager::GetRunManager()) {
          rm->ReinitializeGeometry(true);
          rm->GeometryHasBeenModified();
          std::printf("[room] бетонный пол %s, мир %s\n",
                      want ? "ПОСТРОЕН" : "УБРАН",
                      want ? "расширен" : "прежний");
        }
        return;
      }
      if (c == fCaseOffCmd) {
        const bool want = (v == "on");
        if (want == fDet->fGeom.caseOff) return;
        fDet->fGeom.caseOff = want;
        if (auto* rm = G4RunManager::GetRunManager()) {
          rm->ReinitializeGeometry(true);
          rm->GeometryHasBeenModified();
          std::printf("[caseOff] корпус и крышки %s — ЭТО ДИАГНОСТИКА\n",
                      want ? "ВАКУУМ" : "алюминий");
        }
        return;
      }
      if (c == fPcbPbCmd) {
        const bool want = (v == "on");
        if (want == fDet->fGeom.pcbSolderPb) return;
        fDet->fGeom.pcbSolderPb = want;
        if (auto* rm = G4RunManager::GetRunManager()) {
          rm->ReinitializeGeometry(true);
          rm->GeometryHasBeenModified();
          std::printf("[плата] припой %s, геометрия перестроена\n",
                      want ? "Sn63Pb37 (свинцовый)" : "SAC305 (бессвинцовый)");
        }
        return;
      }
    }
    if (c == fCryYCmd && fDet) {
      const double want = std::atof(v.c_str());
      if (want <= 0.0 || want > 24.0) {
        std::fprintf(stderr, "ОТКАЗ: /asn16/cryY %g вне разумного (0…24 мм)\n",
                     want);
        return;
      }
      if (want == fDet->fGeom.cryY) return;
      fDet->fGeom.cryY = want;
      if (auto* rm = G4RunManager::GetRunManager()) {
        rm->ReinitializeGeometry(true);
        rm->GeometryHasBeenModified();
        std::printf("[cryY] толщина кристалла по пучку %.2f мм, геометрия "
                    "перестроена — ЭТО ДИАГНОСТИКА, не паспорт\n", want);
      }
      return;
    }
    if (c == fTabCmd && fDet) {
      const bool want = (v == "on");
      if (want == fDet->fGeom.table) return;
      fDet->fGeom.table = want;
      if (auto* rm = G4RunManager::GetRunManager()) {
        rm->ReinitializeGeometry(true);
        rm->GeometryHasBeenModified();
        std::printf("[table] столешница %s, геометрия перестроена\n",
                    want ? "ПОСТРОЕНА" : "УБРАНА");
      }
      return;
    }
    if (c == fPackCmd && fDet) {
      const bool want = (v == "on");
      if (want == fDet->fGeom.wt20) return;
      fDet->fGeom.wt20 = want;
      if (auto* rm = G4RunManager::GetRunManager()) {
        rm->ReinitializeGeometry(true);
        rm->GeometryHasBeenModified();
        std::printf("[wt20] пачка %s, геометрия перестроена\n",
                    want ? "ПОСТРОЕНА" : "УБРАНА");
      }
      return;
    }
    if (c != fWinCmd || !fDet) return;
    const bool want = (v == "on");
    if (want == fDet->fGeom.capWindow) return;
    fDet->fGeom.capWindow = want;
    // ГЕОМЕТРИЮ НАДО ПЕРЕСТРОИТЬ ЯВНО. Флага мало: Construct() к этому моменту
    // уже отработал, и без ReinitializeGeometry прогон пошёл бы со СТАРОЙ
    // геометрией, честно напечатав новое значение флага в шапку. Это ровно тот
    // класс отказа, который здесь ловили дважды: число в шапке верное, а
    // считано другое.
    if (auto* rm = G4RunManager::GetRunManager()) {
      rm->ReinitializeGeometry(true);
      rm->GeometryHasBeenModified();
      std::printf("[capWindow] фрезеровка %s, геометрия перестроена\n",
                  want ? "ВКЛЮЧЕНА" : "ОТКЛЮЧЕНА");
    }
  }
};

int main(int argc, char** argv) {
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* det = new ASN16Detector();
  rm->SetUserInitialization(det);
  rm->SetUserInitialization(new PhysList());

  auto* primary = new Primary();
  rm->SetUserAction(primary);

  auto* runAct = new RunAct();
  runAct->fPrimary = primary;
  runAct->fDet = det;
  {
    std::string a;
    for (int i = 1; i < argc; ++i) {
      if (i > 1) a += " ";
      // ТОЛЬКО ИМЯ ФАЙЛА, без каталогов: полный argv занёс бы путь конкретной
      // машины в шапку каждого спектра, а спектр может попасть в репозиторий.
      std::string v(argv[i]);
      const size_t s = v.find_last_of("/\\");
      a += (s == std::string::npos) ? v : v.substr(s + 1);
    }
    runAct->fArgs = a;
  }
  rm->SetUserAction(runAct);
  auto* evtAct = new EventAct(runAct);
  rm->SetUserAction(evtAct);
  auto* mess = new OutMessenger(runAct, det);

  rm->Initialize();
  det->ReportPlanes();
  det->ReportMasses();
  rm->SetUserAction(new Stepping(evtAct, det));
  rm->SetUserAction(new Tracking(runAct));
  rm->SetUserAction(new Stacking());

  auto* ui = G4UImanager::GetUIpointer();
  int rc = 0;
  if (argc > 1) {
    // Результат ApplyCommand ПРОВЕРЯЕТСЯ. Без проверки прогон с несуществующим
    // макросом или с макросом в UTF-8 с BOM печатал «Command aborted (400)» /
    // «Batch is interrupted!!» и завершался с кодом 0, не создав ни одного
    // файла: обёртка считала такой прогон успешным, а каталог спектров
    // оставался от прошлой ревизии. Отказ был трижды описан в комментариях как
    // известный и не чинился (Ж2 аудита кода 05.08.2026).
    const G4int st = ui->ApplyCommand(G4String("/control/execute ") + argv[1]);
    if (st != 0) {
      std::fprintf(stderr,
                   "ОТКАЗ: макрос «%s» не выполнен, код G4UImanager %d.\n"
                   "Частые причины: файла нет; файл в UTF-8 С BOM (первая\n"
                   "команда не распознаётся); опечатка в команде.\n",
                   argv[1], st);
      rc = 2;
    }
  }

  delete mess;
  delete rm;
  return rc;
}
