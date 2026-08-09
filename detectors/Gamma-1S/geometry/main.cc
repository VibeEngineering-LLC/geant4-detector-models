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

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeneralParticleSource.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Navigator.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
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
#include "G4PhysicsModelCatalog.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4VProcess.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
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

  // Тело розыгрыша и ограничение по объёму — как их видит САМ генератор после
  // исполнения макроса. Спрашивать надо у него, а не разбирать макрос: макрос
  // может быть собран драйвером на лету, склеен из нескольких файлов или
  // переопределён следующей командой.
  struct Region {
    G4String vol;          // имя тома в /gps/pos/confine, пусто — ограничения нет
    bool active = false;   // ограничение ДЕЙСТВУЕТ (том найден в этом режиме)
    G4String type, shape;  // Volume / Cylinder и т.п.
    G4ThreeVector centre;
    double radius = 0, halfz = 0;
  };

  Region Confinement() {
    Region r;
    auto* src = fGPS.GetCurrentSource();
    if (!src || !src->GetPosDist()) return r;
    auto* p = src->GetPosDist();
    r.vol = p->GetConfineVolume();
    if (r.vol == "NULL") r.vol = "";
    r.active = p->GetConfined();
    r.type = p->GetPosDisType();
    r.shape = p->GetPosDisShape();
    r.centre = p->GetCentreCoords();
    r.radius = p->GetRadius();
    r.halfz = p->GetHalfZ();
    return r;
  }
};

// Сторож области розыгрыша: тело /gps/pos против confine-тома.
//
// Ловит ДВА дефекта одного семейства, каждый из которых уже стоил недель.
//
//   R61. Том назван, но в текущем режиме не построен. Geant4 печатает
//   «Volume <...> does not exist / Ignoring confine condition» и МОЛЧА
//   продолжает: ограничение снято, розыгрыш идёт по всему телу /gps/pos,
//   включая детектор. Восемь шаблонов разложения по нуклидам так и были
//   посчитаны — 12,5 % распадов родились внутри кристалла NaI.
//
//   R75. Том построен, но тело розыгрыша МЕНЬШЕ него. Ограничение по объёму
//   лишние точки отбрасывает, недостающих не добавляет: часть пробы не
//   облучается вовсе. После перехода сосуда Маринелли на чертёж изготовителя
//   (проба до r = 75 и до z = +65,2) прежнее тело r = 73, z ∈ [−29; +61]
//   теряло наружное кольцо 2 мм и верхние 4 мм — ровно самые дальние от
//   кристалла области, отчего эффективность выходила завышенной.
//
// Оба дефекта беззвучны: имя выходного файла, его формат и порядок величин
// не меняются. Поэтому проверка стоит В МОДЕЛИ и роняет прогон, а не живёт
// отдельным скриптом, мимо которого можно написать следующий макрос.
//
// Возвращает строку для шапки спектра: чем ограничен розыгрыш и каким телом.
// Тома, названные в /gps/pos/confine ТЕКСТОМ макроса. Спрашивать об этом
// генератор бесполезно: G4SPSPosDistribution::ConfineSourceToVolume, не найдя
// тома, ставит VolName = "NULL" и Confine = false, то есть стирает сам факт,
// что ограничение запрашивалось. После этого «ограничения нет» и «ограничение
// молча снято» с точки зрения генератора неразличимы — а это ровно те два
// состояния, которые надо развести (R61).
static std::vector<G4String> MacroConfineVolumes(const std::string& path,
                                                 int depth = 0) {
  std::vector<G4String> out;
  if (depth > 4) return out;                 // защита от кольцевых include
  std::ifstream f(path);
  if (!f) return out;
  std::string ln;
  while (std::getline(f, ln)) {
    std::istringstream is(ln);
    std::string cmd, arg;
    if (!(is >> cmd)) continue;
    if (cmd == "#" || cmd[0] == '#') continue;
    if (cmd == "/gps/pos/confine" && (is >> arg) && arg != "NULL")
      out.push_back(arg);
    else if (cmd == "/control/execute" && (is >> arg))
      for (const auto& v : MacroConfineVolumes(arg, depth + 1))
        out.push_back(v);
  }
  return out;
}

static G4String CheckConfinement(const Primary::Region& r,
                                 const std::vector<G4String>& asked) {
  // Сначала — по тексту макроса: том, который в этом режиме не построен.
  for (const auto& v : asked) {
    if (G4PhysicalVolumeStore::GetInstance()->GetVolume(v, false)) continue;
    G4ExceptionDescription d;
    d << "/gps/pos/confine " << v << ": тома с таким именем в построенной "
         "геометрии НЕТ. Geant4 в этом случае печатает «Volume <" << v
      << "> does not exist / Ignoring confine condition» и продолжает счёт: "
         "ограничение снято, розыгрыш идёт по всему телу /gps/pos, включая "
         "детектор.\nПроверить режим (аргумент 2): том Sample строится только "
         "в режимах vessel*.";
    G4Exception("CheckConfinement", "G1S_CONFINE_MISSING", FatalException, d);
  }

  if (r.vol.empty()) return asked.empty() ? "нет" : "снято";

  if (!r.active) {
    G4ExceptionDescription d;
    d << "/gps/pos/confine " << r.vol << ": ограничение по объёму не действует.";
    G4Exception("CheckConfinement", "G1S_CONFINE_MISSING", FatalException, d);
    return "снято";
  }

  auto* pv = G4PhysicalVolumeStore::GetInstance()->GetVolume(r.vol, false);
  if (!pv || !pv->GetLogicalVolume()) return r.vol + " (габарит не снят)";
  auto* sol = pv->GetLogicalVolume()->GetSolid();
  G4ThreeVector lo, hi;
  sol->BoundingLimits(lo, hi);
  const G4ThreeVector t = pv->GetObjectTranslation();
  const double zlo = (lo.z() + t.z()) / mm, zhi = (hi.z() + t.z()) / mm;

  // Радиус тела берётся у самого тела, а не из габаритной коробки: у коробки
  // угол лежит на r·√2, и цилиндрический розыгрыш требовалось бы делать в
  // полтора раза шире без всякой нужды.
  double rmax = 0;
  if (auto* tb = dynamic_cast<G4Tubs*>(sol)) {
    rmax = tb->GetOuterRadius();
  } else if (auto* pc = dynamic_cast<G4Polycone*>(sol)) {
    const auto* op = pc->GetOriginalParameters();
    for (G4int i = 0; i < op->Num_z_planes; ++i)
      rmax = std::max(rmax, op->Rmax[i]);
  } else {
    rmax = std::max(std::max(std::abs(lo.x()), std::abs(hi.x())),
                    std::max(std::abs(lo.y()), std::abs(hi.y())));
  }
  rmax = rmax / mm + std::hypot(t.x(), t.y()) / mm;

  char buf[320];
  if (r.type != "Volume" || r.shape != "Cylinder") {
    // Судить не берёмся, но и молчать нельзя: строка уйдёт в шапку спектра.
    std::snprintf(buf, sizeof(buf), "%s (розыгрыш %s/%s, покрытие не проверено)",
                  r.vol.c_str(), r.type.c_str(), r.shape.c_str());
    return buf;
  }

  const double cz = r.centre.z() / mm;
  const double off = std::hypot(r.centre.x(), r.centre.y()) / mm;
  const double rr = r.radius / mm, hz = r.halfz / mm;
  const bool covered = (rr >= rmax + off) && (cz - hz <= zlo) && (cz + hz >= zhi);

  std::snprintf(buf, sizeof(buf),
                "%s: том r<=%.1f z %.1f..%.1f; розыгрыш r<=%.1f z %.1f..%.1f%s",
                r.vol.c_str(), rmax, zlo, zhi, rr, cz - hz, cz + hz,
                covered ? "" : " — НЕ ПОКРЫВАЕТ");
  if (!covered && !std::getenv("G1S_ALLOW_PARTIAL_CONFINE")) {
    G4ExceptionDescription d;
    d << "Тело розыгрыша /gps/pos не покрывает том " << r.vol << ".\n" << buf
      << "\nОграничение по объёму лишние точки отбрасывает, но недостающих не "
         "добавляет: часть тома не будет облучена вовсе, и эффективность "
         "выйдет смещённой.\nРасширить /gps/pos/radius и /gps/pos/halfz до "
         "габарита тома с запасом.\nОсознанный прогон по части тома: "
         "G1S_ALLOW_PARTIAL_CONFINE=1.";
    G4Exception("CheckConfinement", "G1S_CONFINE_PARTIAL", FatalException, d);
  }
  return buf;
}

// Время разрешения тракта: энерговыделения, разнесённые больше чем на столько,
// считаются РАЗНЫМИ срабатываниями спектрометра. Подробно — перед EventAct.
constexpr double kResolvingTimeNs = 1000.0;

// --- Разложение отклика по каналам взаимодействия (директива оператора,
// R45 сессии 08.08.2026, по образцу AtomSpectra Nano 16 PRO geometry/main.cc)
//
// Каналы ВЗАИМОИСКЛЮЧАЮЩИЕ и в сумме дают полный спектр: событие попадает
// ровно в один канал. Правило приоритета — по тому, ЧТО унесло энергию из
// кристалла, потому что именно вылет определяет, куда событие уходит из
// пика (тот же порядок и то же обоснование, что в ASN16):
//
//   1) было рождение пар -> канал по числу вылетевших аннигиляционных
//      квантов;
//   2) иначе ушёл сам первичный квант после комптона -> одно- или
//      многократное рассеяние;
//   3) иначе вылетел характеристический рентген;
//   4) иначе вылетел тормозной квант;
//   5) иначе ничего не вылетело: фотоэффект или комптон с последующим
//      поглощением.
//
// "external" здесь — вторичные кванты из свинцовой защиты, MgO-отражателя,
// алюминиевого корпуса/банки или резины-амортизатора, вошедшие в кристалл
// без взаимодействия самого первичного кванта в нём (в ASN16 — аналог
// "вторичные из корпуса и обёртки").
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
                    // энергию принесла вторичная частица извне (защита,
                    // MgO, корпус, резина)
  kChOther,         // остаточный: не должен населяться, служит сторожем
  kNChan
};

static const char* const kChanName[kNChan] = {
  "photo", "compt_full", "compt_esc1", "compt_escN", "xray_esc",
  "brems_esc", "pair_full", "pair_esc1", "pair_esc2", "external", "other"
};

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
  // Та часть fEmit, что родилась атомной релаксацией (K-, L-серии дочернего
  // атома), а не ядерной деэксцитацией. Подмножество, а не отдельный счёт:
  // вычитанием получается чисто ядерная часть, и сумма сходится по построению.
  std::vector<long> fEmitX;
  // Имена моделей, породивших испущенные кванты, — для проверки самого
  // признака: при G1S_DUMP_MODELS=1 печатаются в конце прогона. Без такой
  // печати классификация опиралась бы на память об именах моделей Geant4,
  // а не на то, что каталог отдаёт в этой сборке.
  std::set<G4String> fEmitModels, fEmitXModels;
  bool fDumpModels = std::getenv("G1S_DUMP_MODELS") != nullptr;
  // Разложение отклика по каналам взаимодействия (см. enum Chan выше).
  // Канал ставится В МОМЕНТ СОБЫТИЯ по истории процессов: из готового
  // спектра его восстановить нельзя, форма к тому времени уже сложена.
  // Каналы взаимоисключающие и в сумме дают fHist — проверяется при записи.
  std::vector<std::vector<long>> fChan;
  long fWithSignal = 0;
  double fSumEprim = 0;
  G4String fPart = "?";
  G4String fMode = "shield";
  // Фактические параметры прогона — в шапку выхода. Без них штамп отвечает
  // только «из каких исходников собран exe», и два спектра, различающиеся
  // глубиной колодца или матрицей, выглядят одинаково прослеженными.
  G4String fArgs = "?";
  // Чем ограничен розыгрыш и покрывает ли его тело этот том — строкой в шапку
  // спектра. Заполняется сторожем CheckConfinement на старте каждого прогона.
  G4String fConfine = "?";
  // Тома, названные в /gps/pos/confine текстом макроса, — снимаются до его
  // исполнения (см. MacroConfineVolumes).
  std::vector<G4String> fMacroConfine;
  // Доля телесного угла розыгрыша (1−cos θ)/2 при конусе, иначе 1. Прямой
  // множитель на eps конусных сеток; сообщает его тот, кто разыгрывал.
  // Спрашивается в EndOfRunAction, а не при настройке: угол задаётся макросом,
  // то есть ПОСЛЕ создания действий, и опрос до /run/beamOn вернул бы значение
  // по умолчанию — молча и правдоподобно.
  double fSolidAngleFrac = 1.0;
  // Первичные вершины, попавшие ВНУТРЬ кристалла. Всегда обязан быть нулём:
  // источник в этой модели — проба или внешний пучок, но не сам детектор.
  // Сторож заведён по разбору R61: макрос запирал источник в объём пробы
  // (/gps/pos/confine Sample), а прогон шёл в режиме без сосуда, где тома
  // «Sample» нет вовсе. Geant4 на несуществующий том в confine печатает
  // предупреждение и МОЛЧА снимает ограничение — розыгрыш пошёл по всему
  // сырому цилиндру, и 12,5 % распадов родились в NaI. Восемь шаблонов
  // разложения по нуклидам оказались негодны, и заметили это только по
  // «лишнему» горбу на графике вебки, спустя недели. Предупреждение в
  // потоке вывода такой прогон не остановило: его никто не читает.
  long fSrcInCrystal = 0;
  // Срабатывания, вызванные вторичным квантом, РОДИВШИМСЯ В СВИНЦЕ ЗАЩИТЫ, —
  // характеристическое рентгеновское излучение (ХРИ) свинца, K-серия 72,8 /
  // 75,0 / 84,9 / 87,3 кэВ. Это НЕ излучение пробы: свинец флуоресцирует под
  // квантами пробы и фона, и его линии лежат ровно там же, где K-серия
  // дочерних ряда. Разделять их по энергии нельзя, по МЕСТУ РОЖДЕНИЯ — можно,
  // и это делает МК, а не модель постфактум (разбор R69: прежняя сущность
  // рентгена отбиралась окном 60-110 кэВ и забирала ядерные гамма-линии
  // Th-228 84,4 и Th-232 63,8, а вычитание резалось по нулю).
  //
  // Разбиение по СРАБАТЫВАНИЯМ, не по энергии: срабатывание целиком
  // относится к ХРИ защиты, если хоть часть его энергии принёс квант
  // свинцового происхождения. Смешанные срабатывания редки (квант 75-87 кэВ
  // поглощается целиком), доля объявлена полем в шапке.
  std::vector<long> fPbX;    // флуоресценция свинца защиты
  std::vector<long> fShX;    // обратное рассеяние в защите
  std::vector<long> fSrcX;   // рентген атомной релаксации самой пробы
  long fPbXHits = 0;
  long fShXHits = 0;
  long fSrcXHits = 0;

  // Разбиение срабатываний по ЧИСЛУ И ТИПУ первичных частиц распада, чьи
  // ветви принесли энергию (R78). Отвечает на прямой вопрос: чего именно
  // недостаёт методу 2, который складывает отклики отдельных гамма-линий и
  // потому моделирует РОВНО одноквантовые срабатывания.
  //   fCo1    — энергию принесла ветвь ОДНОГО гамма-кванта распада: это и
  //             есть то, что метод 2 умеет;
  //   fCoN    — сложились ветви ДВУХ и более гамма-квантов одного распада:
  //             истинное совпадение (true coincidence summing), приход;
  //   fCoBeta — в срабатывании участвовала ветвь беты (сама бета либо её
  //             тормозное): в свёртке гамма-библиотеки такого вклада нет
  //             вовсе;
  //   fCoOth  — прочее (альфа, отдача ядра, конверсионный электрон).
  // Разбиение полное и взаимоисключающее: сумма четырёх равна fHist.
  std::vector<long> fCo1, fCoN, fCoBeta, fCoOth;
  Primary* fPrimary = nullptr;
  G4String fOut = "spectrum.csv";

  RunAct() : fHist(kBins + 1, 0), fEmit(kBins + 1, 0),
             fEmitX(kBins + 1, 0),
             fChan(kNChan, std::vector<long>(kBins + 1, 0)),
             fPbX(kBins + 1, 0), fShX(kBins + 1, 0),
             fSrcX(kBins + 1, 0),
             fCo1(kBins + 1, 0), fCoN(kBins + 1, 0),
             fCoBeta(kBins + 1, 0), fCoOth(kBins + 1, 0) {}

  void BeginOfRunAction(const G4Run*) override {
    // Проверка ДО первого события: смысла считать час, чтобы потом узнать, что
    // источник разыгран не по пробе, нет. Спрашивается каждый раз, а не один
    // раз за прогон: в макросе несколько /run/beamOn, и между ними том или
    // тело розыгрыша могут смениться.
    if (fPrimary) fConfine = CheckConfinement(fPrimary->Confinement(),
                                              fMacroConfine);
    std::fill(fHist.begin(), fHist.end(), 0L);
    std::fill(fEmit.begin(), fEmit.end(), 0L);
    std::fill(fEmitX.begin(), fEmitX.end(), 0L);
    fEmitModels.clear();
    fEmitXModels.clear();
    std::fill(fPbX.begin(), fPbX.end(), 0L);
    std::fill(fShX.begin(), fShX.end(), 0L);
    std::fill(fSrcX.begin(), fSrcX.end(), 0L);
    std::fill(fCo1.begin(), fCo1.end(), 0L);
    std::fill(fCoN.begin(), fCoN.end(), 0L);
    std::fill(fCoBeta.begin(), fCoBeta.end(), 0L);
    std::fill(fCoOth.begin(), fCoOth.end(), 0L);
    fSrcInCrystal = 0;
    fPbXHits = 0;
    fShXHits = 0;
    fSrcXHits = 0;
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

  // Испущенное при распаде, РАЗДЕЛЁННОЕ ПО ПРОИСХОЖДЕНИЮ: рентген атомной
  // релаксации против гаммы ядерной деэксцитации. Разделение принципиально:
  // в полосе 60-110 кэВ лежат и K-серии Z = 80…83, и ядерные линии Th-228
  // 84,373 и Th-232 63,8. Отбор по энергетическому окну их не различает и
  // отдаёт ядерные линии рентгену — то есть лишает Th-228 единственной
  // наблюдаемой линии (разбор R69). Здесь признак берётся у самого Geant4:
  // трек несёт идентификатор породившей его МОДЕЛИ.
  void FillEmitX(double eKeV) {
    if (eKeV <= 0) return;
    int b = static_cast<int>(eKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fEmitX[b];
  }

  // Энергия первичной частицы учитывается ОДИН раз на событие Geant4, а Fill()
  // вызывается по разу на каждое РАЗДЕЛЁННОЕ ВО ВРЕМЕНИ срабатывание внутри
  // события (см. EventAct). Раньше это было одной функцией, и при переходе к
  // разделению E_prim_keV множился бы на число срабатываний.
  void FillPrimary(double eprim) { fSumEprim += eprim; }

  // Классификация срабатывания по защите ВЗАИМОИСКЛЮЧАЮЩАЯ и приоритетная:
  // флуоресценция сильнее рассеяния. Иначе сумма трёх вкладов (проба,
  // рассеяние, флуоресценция) перестала бы равняться полному спектру.
  void Fill(double edepKeV, int chan = -1, bool fromPb = false,
            bool fromShield = false, bool fromSrcX = false,
            int nGammaAnc = 0, bool anyBetaAnc = false,
            bool anyOtherAnc = false) {
    if (edepKeV <= 0) return;
    ++fWithSignal;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fHist[b];
    if (fromPb)            { ++fPbX[b]; ++fPbXHits; }
    else if (fromShield)   { ++fShX[b]; ++fShXHits; }
    // Рентген пробы — ОТДЕЛЬНЫЙ признак, не входящий в приоритетную цепочку
    // «свинец → защита»: он относится к самой пробе и вычитается из шаблона
    // её нуклида, тогда как первые два вычитаются из отклика как чужой вклад.
    if (fromSrcX)          { ++fSrcX[b]; ++fSrcXHits; }
    // Разбиение по происхождению вклада (R78). Приоритет: бета сильнее
    // числа гамма-ветвей, потому что вклад тормозного метод 2 не моделирует
    // ВООБЩЕ, а совпадение — моделирует частично (список сумм-пиков); смешать
    // их в одну корзину значило бы приписать методу 2 то, чего у него нет.
    if (anyBetaAnc)        ++fCoBeta[b];
    else if (nGammaAnc >= 2) ++fCoN[b];
    else if (nGammaAnc == 1) ++fCo1[b];
    else if (anyOtherAnc)  ++fCoOth[b];
    if (chan >= 0 && chan < kNChan) ++fChan[chan][b];
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
    // Чем ограничен розыгрыш и покрывает ли его тело этот том. В шапке, а не
    // только в консоли: спектр переживает лог, а вопрос «по всей ли пробе
    // разыгран источник» задаётся к файлу спустя недели (разбор R61 и R75).
    std::fprintf(f, "# confine = %s\n", fConfine.c_str());
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
    // Сторож «источник в детекторе» в шапке, а не только в консоли: спектр
    // переживает лог, и вопрос «а этот файл прогнан со сторожем?» должен
    // отвечаться по самому файлу. Ноль — прогон проверен; поля нет — файл
    // старше сторожа и доверия не заслуживает (разбор R61).
    std::fprintf(f, "# src_in_crystal = %ld\n", fSrcInCrystal);
    // Срабатываний от ХРИ свинца защиты — подмножество N_with_signal.
    // Их спектр лежит отдельным файлом *_shield.csv (колонка pb_fluor) и
    // ВЫЧИТАЕТСЯ из шаблона нуклида читающей стороной: разбиение точное, по
    // месту рождения кванта, без подгонки формы (см. разбор R69).
    std::fprintf(f, "# pbx_hits = %ld\n", fPbXHits);
    std::fprintf(f, "# shx_hits = %ld\n", fShXHits);
    std::fprintf(f, "# resolving_time_ns = %.0f\n", kResolvingTimeNs);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n", kBinKeV);
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
        std::fprintf(g, "# src_sha1 = %s\n", G1S_SRC_SHA1);
        std::fprintf(g, "# git_describe = %s\n", G1S_GIT_DESCRIBE);
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
      // приоритета пропускает случай, и это дефект разложения, а не мелочь.
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

    // --- вклад защиты, отдельным файлом ---
    // Две колонки: флуоресценция свинца и обратное рассеяние в защите.
    // Это ПОДМНОЖЕСТВА полного спектра, отобранные по истории трека, а не
    // модельные формы. Читающая сторона вычитает их из шаблона нуклида —
    // разбиение точное, обрезки по нулю быть не может (в отличие от
    // прежней сущности рентгена по энергетическому окну, разбор R69).
    if (fPbXHits || fShXHits || fSrcXHits) {
      G4String sn = fOut;
      const size_t dot = sn.rfind('.');
      sn = (dot == G4String::npos ? sn : sn.substr(0, dot)) + "_shield.csv";
      FILE* g = std::fopen(sn.c_str(), "w");
      if (g) {
        std::fprintf(g, "# вклад защиты, отобран по истории трека\n");
        std::fprintf(g, "#   pb_fluor — квант родился в свинце (K-серия Pb "
                        "72,8/75,0/84,9/87,3 кэВ)\n");
        std::fprintf(g, "#   sh_scat  — квант побывал в защите и вернулся "
                        "(обратное рассеяние)\n");
        std::fprintf(g, "#   src_xray — энергию принёс рентген атомной "
                        "релаксации САМОЙ пробы (K-, L-серии дочерних атомов "
                        "ряда)\n");
        std::fprintf(g, "# первые два признака взаимоисключающие, "
                        "флуоресценция приоритетнее рассеяния\n");
        std::fprintf(g, "# src_xray — ОТДЕЛЬНЫЙ признак, к первым двум не "
                        "относится: это вклад самой пробы, и вычитается он "
                        "из шаблона её нуклида\n");
        std::fprintf(g, "# отбор по модели-родителю трека, а НЕ по "
                        "энергетическому окну: окно 60-110 кэВ забирало "
                        "ядерные линии Th-228 84,4 и Th-232 63,8 (R69)\n");
        std::fprintf(g, "# src_sha1 = %s\n", G1S_SRC_SHA1);
        std::fprintf(g, "# mode = %s\n", fMode.c_str());
        std::fprintf(g, "# run_args = %s\n", fArgs.c_str());
        std::fprintf(g, "# src_in_crystal = %ld\n", fSrcInCrystal);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "# pbx_hits = %ld\n", fPbXHits);
        std::fprintf(g, "# shx_hits = %ld\n", fShXHits);
        std::fprintf(g, "# srcx_hits = %ld\n", fSrcXHits);
        std::fprintf(g, "E_keV,pb_fluor,sh_scat,src_xray\n");
        for (int i = 0; i <= kBins; ++i)
          if (fPbX[i] || fShX[i] || fSrcX[i])
            std::fprintf(g, "%.1f,%ld,%ld,%ld\n", (i + 0.5) * kBinKeV,
                         fPbX[i], fShX[i], fSrcX[i]);
        std::fclose(g);
      }
      G4cout << "SHIELD флуоресценция " << fPbXHits << ", рассеяние "
             << fShXHits << ", рентген пробы " << fSrcXHits << " из "
             << fWithSignal << " срабатываний" << G4endl;
    }

    // Разбиение срабатываний по происхождению вклада (R78). Пишется только
    // для прогонов распада: у моноэнергетического источника предок один по
    // построению, и файл был бы копией спектра.
    long n1 = 0, nN = 0, nB = 0, nO = 0;
    for (int i = 0; i <= kBins; ++i) {
      n1 += fCo1[i]; nN += fCoN[i]; nB += fCoBeta[i]; nO += fCoOth[i];
    }
    if (nN + nB + nO > 0) {
      G4String cn = fOut;
      const size_t dot = cn.rfind('.');
      cn = (dot == G4String::npos ? cn : cn.substr(0, dot)) + "_coinc.csv";
      FILE* g = std::fopen(cn.c_str(), "w");
      if (g) {
        std::fprintf(g, "# срабатывания по происхождению вклада: сколько "
                        "РАЗНЫХ первичных частиц распада принесли энергию\n");
        std::fprintf(g, "#   one_gamma  — ровно одна гамма-ветвь: ровно то, "
                        "что моделирует свёртка библиотеки линий\n");
        std::fprintf(g, "#   coinc      — две и более гамма-ветви одного "
                        "распада: истинное совпадение (true coincidence "
                        "summing)\n");
        std::fprintf(g, "#   beta       — участвовала бета-ветвь (сама бета "
                        "или её тормозное); в свёртке гамма-библиотеки "
                        "такого вклада нет\n");
        std::fprintf(g, "#   other      — альфа, отдача ядра, конверсионный "
                        "электрон\n");
        std::fprintf(g, "# разбиение полное и взаимоисключающее: сумма "
                        "четырёх колонок равна counts из основного файла\n");
        std::fprintf(g, "# приоритет: beta > coinc > one_gamma > other\n");
        std::fprintf(g, "# предок трека определяется по цепочке родителей до "
                        "процесса RadioactiveDecay, а не по энергии\n");
        std::fprintf(g, "# src_sha1 = %s\n", G1S_SRC_SHA1);
        std::fprintf(g, "# mode = %s\n", fMode.c_str());
        std::fprintf(g, "# run_args = %s\n", fArgs.c_str());
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "E_keV,one_gamma,coinc,beta,other\n");
        for (int i = 0; i <= kBins; ++i)
          if (fCo1[i] || fCoN[i] || fCoBeta[i] || fCoOth[i])
            std::fprintf(g, "%.1f,%ld,%ld,%ld,%ld\n", (i + 0.5) * kBinKeV,
                         fCo1[i], fCoN[i], fCoBeta[i], fCoOth[i]);
        std::fclose(g);
      }
      G4cout << "COINC одна гамма " << n1 << ", совпадение " << nN
             << ", бета " << nB << ", прочее " << nO << " из " << fWithSignal
             << " срабатываний" << G4endl;
    }

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
        // Режим и сторож — и здесь, теми же именами полей. Спектр испускания
        // от положения источника не зависит (проверено при разборе R61:
        // дефект, менявший депозитные спектры в 3-8 раз, сдвинул выход
        // K-серии на 0,29 %), и проверка ему физически не нужна. Поля всё
        // равно пишутся: читающая сторона не должна знать, для какого файла
        // проверку можно пропустить, — исключение «этот не проверяем» и есть
        // тот механизм, которым R61 прожил недели.
        std::fprintf(g, "# mode = %s\n", fMode.c_str());
        std::fprintf(g, "# run_args = %s\n", fArgs.c_str());
        std::fprintf(g, "# src_in_crystal = %ld\n", fSrcInCrystal);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "E_keV,counts\n");
        for (int i = 0; i <= kBins; ++i)
          if (fEmit[i]) std::fprintf(g, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fEmit[i]);
        std::fclose(g);
      }
      G4cout << "EMIT всего " << emitted << " квантов на " << N
             << " распадов -> " << en << G4endl;

      // --- эмиссия, разделённая по происхождению кванта ---
      // Отдельным файлом, как _chan.csv и _shield.csv: формат «E_keV,counts»
      // основного файла читают все скрипты дерева.
      long emx = 0;
      for (long c : fEmitX) emx += c;
      G4String xn = (dot == G4String::npos ? fOut : fOut.substr(0, dot))
                  + "_emitx.csv";
      FILE* h = std::fopen(xn.c_str(), "w");
      if (h) {
        std::fprintf(h, "# испущенное при распаде, разделённое по "
                        "происхождению кванта, на %ld распадов\n", N);
        std::fprintf(h, "#   x_atomic  — рентген атомной релаксации "
                        "(K-, L-серии дочернего атома)\n");
        std::fprintf(h, "#   g_nuclear — гамма ядерной деэксцитации\n");
        std::fprintf(h, "# признак — модель Geant4, породившая трек "
                        "(G4Track::GetCreatorModelID), не энергетическое "
                        "окно (разбор R69)\n");
        std::fprintf(h, "# сумма колонок равна counts из _emit.csv "
                        "по построению\n");
        std::fprintf(h, "# src_sha1 = %s\n", G1S_SRC_SHA1);
        std::fprintf(h, "# git_describe = %s\n", G1S_GIT_DESCRIBE);
        std::fprintf(h, "# mode = %s\n", fMode.c_str());
        std::fprintf(h, "# run_args = %s\n", fArgs.c_str());
        std::fprintf(h, "# src_in_crystal = %ld\n", fSrcInCrystal);
        std::fprintf(h, "# N_primaries = %ld\n", N);
        std::fprintf(h, "# x_atomic_total = %ld\n", emx);
        std::fprintf(h, "E_keV,x_atomic,g_nuclear\n");
        for (int i = 0; i <= kBins; ++i)
          if (fEmit[i])
            std::fprintf(h, "%.1f,%ld,%ld\n", (i + 0.5) * kBinKeV,
                         fEmitX[i], fEmit[i] - fEmitX[i]);
        std::fclose(h);
      }
      G4cout << "EMITX рентген атомной релаксации " << emx << " из "
             << emitted << G4endl;
      if (fDumpModels) {
        G4cout << "MODELS ядерные:";
        for (const auto& s : fEmitModels) G4cout << " [" << s << "]";
        G4cout << G4endl << "MODELS рентген:";
        for (const auto& s : fEmitXModels) G4cout << " [" << s << "]";
        G4cout << G4endl;
      }
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
  // (глобальное время, кэВ, «принесено квантом свинцового происхождения»)
  // по каждому шагу в кристалле
  struct Dep {
    double t, e;
    bool pb;    // принесено квантом, родившимся в свинце (флуоресценция)
    bool sh;    // принесено квантом, побывавшим в защите (рассеяние)
    bool xr;    // принесено рентгеном атомной релаксации САМОЙ пробы
    int anc;    // предок: трек ПЕРВИЧНОЙ частицы распада, от которой ветвь
    int apdg;   // тип этой первичной частицы (PDG)
    bool operator<(const Dep& o) const { return t < o.t; }
  };
  std::vector<Dep> fDep;

  // Предок вклада (R78). Каждый трек возводится к той ПЕРВИЧНОЙ частице
  // распада, от которой пошла его ветвь: гамма ядерной деэксцитации, рентген
  // атомной релаксации, бета, альфа, отдача ядра. Признак наследуется так же,
  // как fPbBorn, и позволяет ответить на два вопроса, которые иначе неотличимы
  // по спектру:
  //   * сколько РАЗНЫХ квантов одного распада сложилось в срабатывании —
  //     это и есть истинное совпадение (true coincidence summing);
  //   * пришла ли энергия по бета-ветви (сама бета или её тормозное) —
  //     этого вклада нет в свёртке библиотеки гамма-линий вовсе.
  // Метод 2 складывает отклики отдельных линий и потому моделирует РОВНО
  // одноквантовые срабатывания; разделение даёт прямую меру того, чего ему
  // недостаёт, вместо оценки по разности полных спектров.
  std::map<int, int> fAnc;      // трек -> предок
  std::map<int, int> fAncPdg;   // предок -> PDG первичной частицы

  // Защита участвует в отклике ДВУМЯ разными способами, и оба до сих пор
  // сидели внутри шаблонов нуклидов неразличимо:
  //
  //   fPbBorn      — квант РОДИЛСЯ в свинце: характеристическое излучение
  //                  свинца, K-серия 72,8 / 75,0 / 84,9 / 87,3 кэВ.
  //                  Тормозное на этих энергиях пренебрежимо, рассеянный
  //                  квант новым треком не становится, поэтому «рождение
  //                  в свинце» и есть флуоресценция.
  //   fShieldSeen  — квант ПОБЫВАЛ в любом теле защиты (Pb, Cd, Cu, сталь)
  //                  и вернулся в кристалл: обратное рассеяние. Даёт
  //                  характерный горб около 200-250 кэВ и подложку.
  //
  // Оба признака наследуются потомками. Разделять эти вклады по ЭНЕРГИИ
  // нельзя: ХРИ свинца ложится ровно на K-серию дочерних ряда, а обратное
  // рассеяние — на комптоновскую подложку. По ИСТОРИИ ТРЕКА можно, и знает
  // её только МК. Приоритет при классификации срабатывания: флуоресценция
  // сильнее рассеяния, рассеяние сильнее «чистой пробы».
  std::set<int> fPbBorn;
  std::set<int> fShieldSeen;
  // Кванты рентгена атомной релаксации САМОЙ пробы (K-, L-серии дочерних
  // атомов ряда) и их потомки. Тот же приём, что fPbBorn, но признак берётся
  // не по материалу рождения, а по модели-родителю: релаксация идёт в том же
  // объёме, что и ядерный распад, и по месту их не различить. Нужен, чтобы
  // рентген выделялся отдельной сущностью ТОЧНО, вычитанием подмножества, а
  // не приближённой формой по энергетическому окну (разбор R69).
  std::set<int> fXrayBorn;

  // Признаки истории события для разложения по каналам. Ставятся в Stepping.
  // Без fChargedIn/fBeta (были у ASN16 — там источник бета-активен через
  // акриловый пенал; здесь источник в свинцовом экране, вклад беты в
  // кристалл не рассматривается, см. память beta-must-be-checked — вопрос
  // решён прогоном ДЛЯ ASN16, для Гамма-1С отдельно не вставал).
  int fFirst = 0;          // 1 phot, 2 compt, 3 conv, 0 ничего неупругого
  bool fHadRayl = false;   // было упругое рассеяние (энергии не оставляет)
  int fNCompt = 0;         // сколько раз первичный квант рассеялся в кристалле
  bool fHadConv = false;   // было рождение пар в кристалле
  int fNAnnihEsc = 0;      // сколько аннигиляционных квантов покинуло кристалл
  double fEBremEsc = 0;    // энергия вылетевшего тормозного, кэВ
  double fEXrayEsc = 0;    // энергия вылетевших прочих вторичных гамма, кэВ
  bool fPrimEsc = false;   // сам первичный квант вышел из кристалла

  // Сторож «источник в детекторе» (см. RunAct::fSrcInCrystal). Проставляется
  // из main() ПОСЛЕ rm->Initialize(): до неё геометрии ещё нет и fCrystalLV
  // равен nullptr. Навигатор свой, отдельный от трекингового: тот в момент
  // конца события хранит состояние последнего шага, и запрос точки по нему
  // сбил бы транспорт.
  const G4LogicalVolume* fCry = nullptr;
  G4Navigator fNav;
  bool fNavReady = false;

  explicit EventAct(RunAct* r) : fRun(r) {}
  void BeginOfEventAction(const G4Event*) override {
    fDep.clear();
    // ОБА множества чистить обязательно: идентификаторы треков в каждом
    // событии начинаются заново с единицы, и несброшенное множество через
    // десяток событий содержит почти все возможные идентификаторы. Поймано
    // прогоном: доля «рассеяния в защите» вышла 99,8 % вместо единиц
    // процентов. Сначала был обратный промах — признак рассеяния не
    // заполнялся в Stepping вовсе и давал строгий ноль.
    fPbBorn.clear();
    fShieldSeen.clear();
    fXrayBorn.clear();
    fAnc.clear();
    fAncPdg.clear();
    fFirst = 0; fHadRayl = false; fNCompt = 0; fHadConv = false;
    fNAnnihEsc = 0; fEBremEsc = 0; fEXrayEsc = 0; fPrimEsc = false;
  }

  // Где родилась первичная частица. Ноль событий в кристалле — условие
  // осмысленности прогона; одно такое событие роняет прогон немедленно,
  // не через час счёта и не молча в предупреждении (разбор R61).
  void CheckSourcePlacement(const G4Event* e) {
    if (!fCry || e->GetNumberOfPrimaryVertex() == 0) return;
    if (!fNavReady) {
      auto* tm = G4TransportationManager::GetTransportationManager();
      fNav.SetWorldVolume(tm->GetNavigatorForTracking()->GetWorldVolume());
      fNavReady = true;
    }
    const G4ThreeVector pos = e->GetPrimaryVertex(0)->GetPosition();
    const G4VPhysicalVolume* pv =
        fNav.LocateGlobalPointAndSetup(pos, nullptr, false, true);
    if (!pv || pv->GetLogicalVolume() != fCry) return;
    ++fRun->fSrcInCrystal;
    if (std::getenv("G1S_ALLOW_SOURCE_IN_CRYSTAL")) return;
    G4ExceptionDescription d;
    d << "Первичная частица рождена ВНУТРИ кристалла NaI: ("
      << pos.x() / mm << ", " << pos.y() / mm << ", " << pos.z() / mm
      << ") мм. Источник в этой модели — проба или внешний пучок, но не сам "
         "детектор.\nТипичная причина: /gps/pos/confine <том> назван в макросе, "
         "а том в текущем режиме геометрии не построен — Geant4 снимает "
         "ограничение молча, и розыгрыш идёт по всему заданному телу.\n"
         "Проверить: режим прогона (аргумент 2) строит нужный том? Для "
         "confine Sample нужен режим vessel.\n"
         "Осознанный прогон с источником в кристалле: G1S_ALLOW_SOURCE_IN_"
         "CRYSTAL=1.";
    G4Exception("EventAct::CheckSourcePlacement", "G1S_SRC_IN_CRYSTAL",
                FatalException, d);
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
    CheckSourcePlacement(e);
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
    // энерговыделений одна; при розыгрыше распада событие может разбиться
    // на несколько срабатываний — все получают канал первичного кванта,
    // это огрубление, и оно объявлено (как в ASN16).
    const int ch = Channel();

    // Шаги приходят в порядке обработки треков, а не по времени, поэтому
    // сортировка обязательна: без неё группировка развалится на первом же
    // треке, начавшемся раньше предыдущего.
    std::sort(fDep.begin(), fDep.end());
    // Признак «свинцовое» ставится на СРАБАТЫВАНИЕ целиком, если хоть часть
    // его энергии принёс квант, родившийся в свинце. Смешанные срабатывания
    // редки: квант K-серии свинца 73-88 кэВ поглощается в NaI целиком и с
    // квантом пробы в одно срабатывание попадает лишь при наложении.
    // Огрубление объявлено; его цена меряется полем pbx_hits в шапке.
    // Признак рентгена ставится по ДОЛЕ ЭНЕРГИИ, а не по правилу «хоть часть».
    // Правило «хоть часть» (как у pb/sh) для рентгена даёт ложную структуру на
    // жёстком крае: рентген 75 кэВ и ядерная гамма 2614 кэВ, попавшие в одно
    // окно разрешения 1000 нс, помечали ВСЁ срабатывание 2689 кэВ как
    // рентгеновское. Таких срабатываний мало (72 из 4403 по ветви, 1,6 %), но
    // при нормировке доли нуклида на канал они размазываются и дают слой,
    // повторяющий форму полного отклика с множителем ~0,0025 — на лог-шкале
    // это читается как пики Tl-208 в сущности, где выше 110,5 кэВ физически
    // не может быть ничего (K-серия дочерних кончается на 91 кэВ). Поймано
    // оператором по графику вебки 09.08.2026.
    //
    // Порог 0,999, а не строгое равенство: суммирование double по десяткам
    // шагов даёт последний разряд. Смешанные срабатывания уходят НУКЛИДУ —
    // их энергию определяет гамма, а не рентген, и разбиение остаётся точным
    // (каждое срабатывание ровно в одной корзине).
    //
    // Для pb/sh правило «хоть часть» оставлено сознательно: там сущность
    // вычитается как ЧУЖОЙ вклад защиты, и срабатывание, задетое свинцовым
    // квантом, к пробе уже не относится целиком. Цена огрубления объявлена
    // полем pbx_hits в шапке.
    // Предки собираются ПО СРАБАТЫВАНИЮ, а не по событию: окно разрешения
    // может разбить один распад на два срабатывания, и совпадением считается
    // только то, что сложилось ВНУТРИ окна — ровно как у прибора.
    // Ключ — предок, значение — тип его частицы: считать надо РАЗНЫЕ ветви,
    // а не вклады. Одна ветвь даёт десятки шагов, и счёт по шагам объявил бы
    // совпадением любое многошаговое поглощение одного кванта.
    std::map<int, int> anc;
    auto flush = [&](double s, bool p, bool h, bool x) {
      int ng = 0;
      bool beta = false, other = false;
      for (const auto& kv : anc) {
        if (kv.second == 22) ++ng;
        else if (kv.second == 11 || kv.second == -11) beta = true;
        else if (kv.second != 0) other = true;
      }
      fRun->Fill(s, ch, p, h, x, ng, beta, other);
      anc.clear();
    };
    double sum = fDep[0].e, t0 = fDep[0].t;
    double sumx = fDep[0].xr ? fDep[0].e : 0.0;
    bool pb = fDep[0].pb, sh = fDep[0].sh;
    anc[fDep[0].anc] = fDep[0].apdg;
    for (size_t i = 1; i < fDep.size(); ++i) {
      if (fDep[i].t - t0 > kResolvingTimeNs) {
        flush(sum, pb, sh, sumx >= 0.999 * sum);
        sum = 0; sumx = 0; pb = false; sh = false;
      }
      t0 = fDep[i].t;
      sum += fDep[i].e;
      if (fDep[i].xr) sumx += fDep[i].e;
      pb = pb || fDep[i].pb;
      sh = sh || fDep[i].sh;
      anc[fDep[i].anc] = fDep[i].apdg;
    }
    flush(sum, pb, sh, sumx >= 0.999 * sum);
  }
};

class Stepping : public G4UserSteppingAction {
  EventAct* fEvt;
  const G4LogicalVolume* fCry;
public:
  Stepping(EventAct* ev, const G4LogicalVolume* c) : fEvt(ev), fCry(c) {}
  void UserSteppingAction(const G4Step* s) override {
    auto* pre = s->GetPreStepPoint();
    auto* post = s->GetPostStepPoint();
    auto* h = pre->GetTouchableHandle()->GetVolume();
    const bool inCry = h && h->GetLogicalVolume() == fCry;
    if (!inCry) {
      // Шаг вне кристалла: если он в теле защиты, помечаем трек как
      // побывавший в ней. Квант, вернувшийся оттуда в кристалл, принесёт
      // обратно рассеянную энергию — это отдельный вклад защиты, и
      // отличить его от излучения пробы можно только по истории трека.
      // Материалы защиты в этой геометрии больше нигде не встречаются:
      // головка собрана из Al, MgO, резины, стекла и «Electronics».
      const G4Material* m =
          h ? h->GetLogicalVolume()->GetMaterial() : nullptr;
      if (m) {
        const G4String& mn = m->GetName();
        if (mn == "G4_Pb" || mn == "G4_Cd" || mn == "G4_Cu" ||
            mn == "G4_STAINLESS-STEEL")
          fEvt->fShieldSeen.insert(s->GetTrack()->GetTrackID());
      }
      return;
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
          // взаимодействием не считается (см. ASN16: без этой оговорки
          // канал «рэлей+фотоэффект» уходил в остаточный).
          fEvt->fHadRayl = true;
        }
      }
      // Выход кванта из кристалла. Квант, вернувшийся обратно (из защиты,
      // MgO, корпуса), здесь уже посчитан вылетевшим — огрубление в пользу
      // каналов вылета, доля возвратов не измерена.
      if (post->GetStepStatus() == fGeomBoundary) {
        auto* hp = post->GetTouchableHandle()->GetVolume();
        const bool outCry = !hp || hp->GetLogicalVolume() != fCry;
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
    if (e > 0) {
      const G4int id = trk->GetTrackID();
      const bool pb = fEvt->fPbBorn.count(id) > 0;
      const bool sh = pb || fEvt->fShieldSeen.count(id) > 0;
      const bool xr = fEvt->fXrayBorn.count(id) > 0;
      // Предок вклада (R78). Если трека нет в карте — он и есть свой предок:
      // так бывает у самой первичной частицы, стартовавшей до заполнения.
      auto ia = fEvt->fAnc.find(id);
      const int anc = (ia == fEvt->fAnc.end()) ? id : ia->second;
      auto ip = fEvt->fAncPdg.find(anc);
      const int apdg = (ip == fEvt->fAncPdg.end()) ? 0 : ip->second;
      fEvt->fDep.push_back({pre->GetGlobalTime() / ns, e / keV, pb, sh, xr,
                            anc, apdg});
    }
  }
};

// Учёт квантов, РОЖДЁННЫХ радиоактивным распадом (а не рассеянием). Так выход
// линии p_gamma приходит из той же базы PhotonEvaporation, что и транспорт.
class Tracking : public G4UserTrackingAction {
  RunAct* fRun;
  EventAct* fEvt;
public:
  Tracking(RunAct* r, EventAct* ev) : fRun(r), fEvt(ev) {}
  void PreUserTrackingAction(const G4Track* t) override {
    // --- происхождение из свинца защиты (см. EventAct::fPbBorn) ---
    // Geant4 отдаёт вторичные ПОСЛЕ родителя, поэтому к моменту старта
    // потомка идентификатор родителя в множестве уже есть, и признак
    // наследуется по цепочке без отдельного обхода.
    const G4int id = t->GetTrackID(), pid = t->GetParentID();
    bool pb = pid > 0 && fEvt->fPbBorn.count(pid) > 0;
    if (!pb && pid > 0) {
      // Рождение НОВОГО трека в свинце — это флуоресценция (тормозное на
      // этих энергиях пренебрежимо). Рассеянный квант новым треком не
      // становится, поэтому комптон в свинце сюда не попадает.
      const G4VPhysicalVolume* v = t->GetVolume();
      const G4Material* m = v ? v->GetLogicalVolume()->GetMaterial() : nullptr;
      pb = m && m->GetName() == "G4_Pb";
    }
    if (pb) fEvt->fPbBorn.insert(id);
    // Рентген атомной релаксации пробы: признак наследуется так же, как
    // «родился в свинце». Проставляется ниже, при разборе кванта распада;
    // здесь — только передача потомкам.
    if (pid > 0 && fEvt->fXrayBorn.count(pid) > 0) fEvt->fXrayBorn.insert(id);

    // --- предок: первичная частица распада (см. EventAct::fAnc, R78) ---
    // Ставится ДО отсева по типу частицы: бета и альфа тоже начинают ветви,
    // и именно по ним отличается вклад тормозного от вклада гамма-линий.
    const G4VProcess* cp0 = t->GetCreatorProcess();
    const G4String cn0 = cp0 ? cp0->GetProcessName() : G4String("");
    const bool from_decay = (cn0 == "RadioactiveDecay" ||
                             cn0 == "Radioactivation");
    if (pid == 0 || from_decay) {
      fEvt->fAnc[id] = id;
      fEvt->fAncPdg[id] = t->GetDefinition()->GetPDGEncoding();
    } else {
      auto it = fEvt->fAnc.find(pid);
      fEvt->fAnc[id] = (it == fEvt->fAnc.end()) ? id : it->second;
    }

    if (t->GetDefinition() != G4Gamma::Definition()) return;
    const G4VProcess* p = cp0;
    if (!p) return;                       // первичная частица, не распад
    const G4String& n = cn0;
    if (n != "RadioactiveDecay" && n != "Radioactivation") return;
    const double e = t->GetKineticEnergy() / keV;
    fRun->FillEmit(e);

    // Чем именно квант порождён — атомной релаксацией или ядерной
    // деэксцитацией. Признак берётся у Geant4: трек несёт идентификатор
    // породившей его модели, каталог переводит его в имя. Разбор R69:
    // прежний отбор рентгена по окну 60-110 кэВ отдавал рентгену ядерные
    // линии Th-228 84,373 и Th-232 63,8, то есть лишал Th-228 единственной
    // наблюдаемой линии — той самой, на которой держится проверка возраста
    // ряда. Имена моделей не угадываются: они печатаются прогоном с
    // G1S_DUMP_MODELS=1 и сверяются глазами.
    //
    // СНЯТО ПРОГОНОМ 09.08.2026 (Pb-212, 20 тыс. распадов, Geant4 11.2.1),
    // каталог вернул ровно два имени:
    //     MODELS ядерные: [model_RDM_IT]
    //     MODELS рентген: [model_RDM_AtomicRelaxation]
    // Контроль на том же прогоне: K-серия висмута 75,5/77,5/87,5 кэВ ушла
    // целиком в x_atomic (64/70/39 отсчётов из 67/70/39), ядерная линия
    // Pb-212 238,6 кэВ — целиком в g_nuclear (8660 из 8660), её выход
    // 43,3 % на распад против 43,6 % ENSDF. Признак различает то, что
    // должен, и на энергетическое окно не опирается.
    const G4int mid = t->GetCreatorModelID();
    const G4String mn = (mid < 0) ? G4String("?")
                                  : G4PhysicsModelCatalog::GetModelNameFromID(mid);
    const bool xray = mn.find("luo") != std::string::npos     // Fluo/fluo
                   || mn.find("uger") != std::string::npos    // Auger
                   || mn.find("elax") != std::string::npos    // Relaxation
                   || mn.find("ARM") != std::string::npos;    // atomic rearr.
    if (xray) {
      fRun->FillEmitX(e);
      // Признак наследуется потомками (см. инициализацию xr выше): рентген
      // поглощается в кристалле фотоэффектом, и энергию в счёт вносит уже
      // фотоэлектрон, а не сам квант.
      fEvt->fXrayBorn.insert(id);
    }
    if (fRun->fDumpModels)
      (xray ? fRun->fEmitXModels : fRun->fEmitModels).insert(mn);
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

// Выгрузка ПОСТРОЕННОЙ геометрии в CSV: имя, материал, границы по r и z.
// Смысл — рисовать разрез по тому, что Geant4 действительно собрал, а не по
// тому, что задумано в исходнике. Обходится хранилище физических объёмов,
// каждое тело раскладывается на кольцевые сегменты (r_in, r_out, z0, z1).
// Годится потому, что вся геометрия ГАММА-1С осесимметрична и собрана из
// G4Tubs и G4Polycone; тело другого класса выводится строкой «?» и в разрез
// не попадёт молча — это видно в файле.
static void DumpGeometry(const std::string& path, const std::string& args) {
  std::ofstream f(path);
  if (!f) {
    G4cerr << "DUMPGEOM: не открыть " << path << G4endl;
    return;
  }
  f << "# геометрия ГАММА-1С, снята с построенного дерева\n";
  f << "# src_sha1 " << G1S_SRC_SHA1 << "\n";
  // Аргументы запуска — в том же виде, что в шапке спектров. Нужны потому,
  // что геометрия и материалы зависят от них (режим, плотность, матрица), и
  // потребитель выгрузки обязан иметь возможность СВЕРИТЬ, той ли сборкой и
  // теми ли аргументами она снята, что и шаблоны, которые он ей объясняет.
  f << "# run_args = " << args << "\n";
  f << "name,material,rin_mm,rout_mm,z0_mm,z1_mm\n";
  auto* store = G4PhysicalVolumeStore::GetInstance();
  int unknown = 0;
  for (auto* pv : *store) {
    auto* lv = pv->GetLogicalVolume();
    if (!lv) continue;
    const G4String nm = pv->GetName();
    if (nm == "World") continue;
    const G4String mt = lv->GetMaterial() ? lv->GetMaterial()->GetName() : "?";
    const double zc = pv->GetObjectTranslation().z() / mm;
    auto* sol = lv->GetSolid();
    if (auto* t = dynamic_cast<G4Tubs*>(sol)) {
      const double dz = t->GetZHalfLength() / mm;
      f << nm << "," << mt << "," << t->GetInnerRadius() / mm << ","
        << t->GetOuterRadius() / mm << "," << zc - dz << "," << zc + dz << "\n";
    } else if (auto* p = dynamic_cast<G4Polycone*>(sol)) {
      const auto* op = p->GetOriginalParameters();
      for (G4int i = 0; i + 1 < op->Num_z_planes; ++i) {
        const double z0 = op->Z_values[i] / mm, z1 = op->Z_values[i + 1] / mm;
        if (z1 <= z0) continue;                 // вертикальная ступень
        f << nm << "," << mt << "," << op->Rmin[i] / mm << ","
          << op->Rmax[i] / mm << "," << zc + z0 << "," << zc + z1 << "\n";
      }
    } else {
      ++unknown;
      f << nm << "," << mt << ",?,?,?,?\n";
    }
  }
  // Составы материалов — сюда же, отдельным блоком. Нужны потому, что имя
  // «ОИСН-16» состав не определяет: под ним в комплекте ходят две разные
  // рецептуры (см. G1SDetector::MakeMatrix), и страница обязана называть ту,
  // которой реально посчитана, а не переписанную руками в шаблон.
  f << "\n# состав материалов: доли по МАССЕ, плотность г/см3\n";
  f << "MAT,material,density_g_cm3,element,mass_fraction\n";
  for (auto* m : *G4Material::GetMaterialTable()) {
    const G4double rho = m->GetDensity() / (g / cm3);
    const G4double* fr = m->GetFractionVector();
    for (size_t i = 0; i < m->GetNumberOfElements(); ++i) {
      f << "MAT," << m->GetName() << "," << rho << ","
        << m->GetElement(static_cast<G4int>(i))->GetSymbol() << ","
        << fr[i] << "\n";
    }
  }
  G4cout << "DUMPGEOM " << path << " (тел неизвестного класса: " << unknown
         << ")" << G4endl;
}

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
  if (const char* dg = std::getenv("G1S_DUMP_GEOM")) DumpGeometry(dg, runAct->fArgs);
  // Кристалл известен только после Initialize(): до неё Construct() не звался.
  evtAct->fCry = det->fCrystalLV;
  rm->SetUserAction(new Stepping(evtAct, det->fCrystalLV));
  rm->SetUserAction(new Tracking(runAct, evtAct));

  auto* ui = G4UImanager::GetUIpointer();
  // Текст макроса читается ДО исполнения: генератор, не найдя тома, стирает
  // само имя (см. MacroConfineVolumes), и после исполнения спросить уже не у
  // кого.
  if (argc > 1) {
    runAct->fMacroConfine = MacroConfineVolumes(argv[1]);
    ui->ApplyCommand(G4String("/control/execute ") + argv[1]);
  }

  delete mess;
  delete rm;
  return 0;
}
