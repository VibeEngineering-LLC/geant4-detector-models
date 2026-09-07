#include "Rc103FieldPhysicsList.hh"

#include "G4EmParameters.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"

#include <cstdio>

// Распад ядер и угловые корреляции гамма-каскада. Оба ВЫКЛЮЧЕНЫ по умолчанию:
// прежние прогоны контура стреляют одиночным квантом, и их воспроизводимость
// не должна зависеть от появления этих ключей. Включаются из main.
bool Rc103FieldPhysicsList::gDecay = false;
bool Rc103FieldPhysicsList::gCorrGamma = false;
std::string Rc103FieldPhysicsList::gEmOption = "opt4";
std::string Rc103FieldPhysicsList::gDeexRegion = "";
bool Rc103FieldPhysicsList::gScintLight = false;

double Rc103FieldPhysicsList::gCutMm = 0.05;
std::string Rc103FieldPhysicsList::gDeexMode = "std";
// Пусто = умолчание Geant4. Непустое значение применяется только при deex=max.
std::string Rc103FieldPhysicsList::gPixeModel = "";
std::string Rc103FieldPhysicsList::gPixeElecModel = "";

Rc103FieldPhysicsList::Rc103FieldPhysicsList(double cutMm,
                                             const std::string& deexMode) {
  SetDefaultCutValue(cutMm * mm);
  // Умолчание — option4 (#CFG-2). opt0 регистрируется только по явному
  // ключу и служит суррогатом «без доплеровского уширения» (линия Д).
  if (gEmOption == "opt0") {
    RegisterPhysics(new G4EmStandardPhysics());
  } else {
    RegisterPhysics(new G4EmStandardPhysics_option4());
  }

  // Строго ПОСЛЕ конструктора option4 — см. предупреждение в заголовке.
  auto* p = G4EmParameters::Instance();
  // Режим "fluo" — флуоресценция БЕЗ Оже. Введён 06.09.2026 для разбора W-066:
  // в Geant4 11.2.1 SetAugerCascade сведён к SetAuger (G4EmParameters.hh:139),
  // поэтому Оже нельзя включить без каскада, а каскад ведёт очередь вакансий
  // G4UAtomicDeexcitation::vacancyArray, которая в этой версии не очищается в
  // начале GenerateParticles. Режим разделяет «деэкситацию вообще» и «каскадную
  // очередь»: работы он даёт много, а очередь не наполняет.
  if (deexMode == "fluo" || deexMode == "pixefluo") {
    p->SetFluo(true);
    p->SetDeexcitationIgnoreCut(true);
  }
  // "pixefluo" — это "max" без Оже: PIXE создаёт столько же вакансий, работы
  // столько же, но каскадная очередь не наполняется. Пара к "max" ровно в одну
  // переменную — ею и проверяется, дело в очереди или в объёме работы.
  if (deexMode == "pixefluo") {
    p->SetPixe(true);
  }
  if (deexMode == "deex" || deexMode == "max") {
    p->SetFluo(true);
    p->SetAuger(true);
    p->SetDeexcitationIgnoreCut(true);
  }
  if (deexMode == "max") {
    p->SetPixe(true);
    // Модели сечений PIXE выбираются ТОЛЬКО если заданы явно; пустая строка
    // означает «оставить умолчание Geant4», чтобы прежние прогоны контура
    // воспроизводились байт в байт. Ключи добавлены 05.09 для разбора W-055:
    // порча памяти воспроизводится ровно при включённом PIXE (6 падений из 8
    // против 0 из 5 без него), и требуется отделить дефект КОНКРЕТНОЙ модели
    // от дефекта всей ветви PIXE.
    if (!gPixeModel.empty()) {
      p->SetPIXECrossSectionModel(gPixeModel);
    }
    if (!gPixeElecModel.empty()) {
      p->SetPIXEElectronCrossSectionModel(gPixeElecModel);
    }
  }

  // Адресная деэкситация: флуоресценция, Оже и PIXE включаются ТОЛЬКО в
  // названном регионе, глобально остаётся то, что задано deexMode. Смысл —
  // не платить за PIXE в свинце и корпусе, где сигнал не считается.
  // Работает начиная с 11.4.2: до неё региональная настройка игнорировалась
  // при глобально включённом PIXE (баг 2650, `G4VAtomDeexcitation.cc`).
  // Пустая строка = регион не назначен, поведение прежнее.
  if (!gDeexRegion.empty()) {
    p->SetDeexActiveRegion(gDeexRegion, true, true, true);
  }

  // Штатная сцинтилляция — для НЕЗАВИСИМОЙ сверки модели света (линия Г).
  // Оптический транспорт не нужен: фотоны убиваются при постановке в стек
  // (NpsmBenchStackingAction), считается только их число. Отключаем всё,
  // кроме сцинтилляции, чтобы не платить за границы и отражения.
  if (gScintLight) {
    auto* op = new G4OpticalPhysics();
    RegisterPhysics(op);
    auto* opar = G4OpticalParameters::Instance();
    opar->SetProcessActivation("Cerenkov", false);
    opar->SetProcessActivation("OpAbsorption", false);
    opar->SetProcessActivation("OpRayleigh", false);
    opar->SetProcessActivation("OpMieHG", false);
    opar->SetProcessActivation("OpBoundary", false);
    opar->SetProcessActivation("OpWLS", false);
    opar->SetProcessActivation("OpWLS2", false);
    opar->SetScintTrackSecondariesFirst(false);
    // БЕЗ этого флага Geant4 берёт скалярный SCINTILLATIONYIELD и полностью
    // игнорирует вектор ELECTRONSCINTILLATIONYIELD: первый прогон дал 2181
    // фотон на 2159 МэВ, то есть ровно 1 фотон/МэВ вместо кривой (07.09.2026).
    // Флаг включает выход, зависящий от ТИПА частицы, — тогда читается
    // именно наша таблица L(E).
    opar->SetScintByParticleType(true);
    std::fprintf(stdout, "ScintLight: G4OpticalPhysics зарегистрирован, "
                         "активна только сцинтилляция\n");
  }

  // Распад первичного ядра. Нужен только в режиме primary=ion; без него ядро
  // просто стоит на месте и ни одного кванта не рождается.
  if (gDecay) {
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    // Угловые корреляции гамма-каскада. В Geant4 они выключены по умолчанию —
    // кванты каскада испускаются независимо и изотропно. Для сверки формы
    // спектра с ИЗМЕРЕНИЕМ это допущение: у реального Co-60 направления
    // 1173 и 1332 кэВ статистически связаны, и от этого зависит доля событий,
    // где зарегистрированы оба, то есть величина каскадного суммирования.
    auto* dp = G4NuclearLevelData::GetInstance()->GetParameters();
    if (dp) dp->SetCorrelatedGamma(gCorrGamma);
  }

  gCutMm = cutMm;
  gDeexMode = deexMode;

  // Диагностика ASCII-строкой: читается из лога любой кодировкой консоли и
  // служит проверяемым артефактом того, что флаги ДЕЙСТВИТЕЛЬНО выставлены.
  // Модели печатаются ИЗ САМОГО G4EmParameters, а не из наших переменных:
  // иначе строка подтверждала бы лишь то, что мы просили, а не то, что
  // применилось. На этом уже обжигались (W-052: мутация не применилась,
  // а тест считал её внесённой).
  std::fprintf(stdout, "Rc103FieldPhysicsList: em_option=%s\n",
               gEmOption.c_str());
  std::fprintf(stdout,
               "Rc103FieldPhysicsList: cut_mm=%.4f deex=%s fluo=%d auger=%d "
               "pixe=%d ignore_cut=%d pixe_model=%s pixe_e_model=%s\n",
               cutMm, deexMode.c_str(), static_cast<int>(p->Fluo()),
               static_cast<int>(p->Auger()), static_cast<int>(p->Pixe()),
               static_cast<int>(p->DeexcitationIgnoreCut()),
               p->PIXECrossSectionModel().c_str(),
               p->PIXEElectronCrossSectionModel().c_str());
  // Состояние распада и корреляций — тоже из САМОГО Geant4, а не из наших
  // переменных (W-052). Печатается всегда, в том числе при выключенном распаде:
  // отсутствие строки в старом логе иначе не отличить от выключенной настройки.
  {
    auto* dp = G4NuclearLevelData::GetInstance()->GetParameters();
    std::fprintf(stdout, "Rc103FieldPhysicsList: decay=%d correlated_gamma=%d\n",
                 static_cast<int>(gDecay),
                 dp ? static_cast<int>(dp->CorrelatedGamma()) : -1);
  }
}