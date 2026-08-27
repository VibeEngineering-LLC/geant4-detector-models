#include "Rc103RoomFieldPrimaryGeneratorAction.hh"

#include "Rc103RoomFieldGeometry.hh"

#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>

#ifdef ROOMFIELD_BIRTH_STATS
#include <vector>
namespace {
std::vector<long long> gBirthCount;
}
void RoomFieldBirthStatsCount(std::size_t idx) {
  if (gBirthCount.size() != gSlabs.size()) gBirthCount.assign(gSlabs.size(), 0);
  ++gBirthCount[idx];
}
void RoomFieldBirthStatsPrint() {
  long long tot = 0;
  for (long long c : gBirthCount) tot += c;
  if (tot <= 0) {
    std::fprintf(stdout, "BIRTHSTATS: нет ни одного рождения\n");
    return;
  }
  // Ожидание считается по ВЫБРАННЫМ плитам (src=), иначе в режимах
  // brick/concrete диагностика сравнивала бы с заведомо чужим знаменателем.
  double brickBorn = 0, brickMass = 0, brickVol = 0, totMass = 0, totVol = 0;
  for (std::size_t i = 0; i < gSlabs.size(); ++i) {
    std::fprintf(stdout, "BIRTHSTATS slab %-16s born=%.6f m_kg=%.3f\n",
                 gSlabs[i].name.c_str(), double(gBirthCount[i]) / double(tot),
                 gSlabs[i].massG / 1000.0);
    if (!RoomSlabSelected(gSlabs[i])) continue;
    if (gSlabs[i].brick) {
      brickBorn += double(gBirthCount[i]) / double(tot);
      brickMass += gSlabs[i].massG;
      brickVol += gSlabs[i].volumeCm3;
    }
    totMass += gSlabs[i].massG;
    totVol += gSlabs[i].volumeCm3;
  }
  std::fprintf(stdout,
               "BIRTHSTATS ИТОГ: N=%lld  доля_кирпича_рождённых=%.6f  "
               "ожидание_ПО_МАССЕ=%.6f  (по объёму было бы %.6f)\n",
               tot, brickBorn, brickMass / totMass, brickVol / totVol);
}
#endif

Rc103RoomFieldPrimaryGeneratorAction::Rc103RoomFieldPrimaryGeneratorAction(
    int ionZ, int ionA)
    : fIonZ(ionZ), fIonA(ionA), fGun(1) {
  if (gSlabs.empty()) {
    std::fprintf(stderr,
                 "Rc103RoomFieldPrimaryGeneratorAction: FATAL gSlabs пуст.\n");
    std::abort();
  }
  double total = 0.0;
  for (std::size_t i = 0; i < gSlabs.size(); ++i) {
    const RoomSlab& s = gSlabs[i];
    if (!(s.massG > 0.0)) {
      std::fprintf(stderr,
                   "Rc103RoomFieldPrimaryGeneratorAction: FATAL нулевая масса "
                   "плиты '%s' — AssignRoomDensities() не вызван.\n",
                   s.name.c_str());
      std::abort();
    }
    // src=brick|concrete: невыбранный материал в розыгрыш НЕ входит вовсе
    // (геометрически он остаётся и ослабляет излучение — это правильно).
    if (!RoomSlabSelected(s)) continue;
    total += s.massG;
    fSlabIndex.push_back(i);
    fCumMass.push_back(total);
  }
  if (fCumMass.empty()) {
    std::fprintf(stderr,
                 "Rc103RoomFieldPrimaryGeneratorAction: FATAL режим src=%s не "
                 "оставил ни одной плиты для розыгрыша.\n",
                 RoomSrcModeName());
    std::abort();
  }
  for (auto& c : fCumMass) c /= total;
  fCumMass.back() = 1.0;  // защита от накопленной ошибки округления
}

void Rc103RoomFieldPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // 1) плита — с вероятностью, пропорциональной МАССЕ
  const double u = G4UniformRand();
  const auto it = std::lower_bound(fCumMass.begin(), fCumMass.end(), u);
  const std::size_t k = (it == fCumMass.end())
                            ? fCumMass.size() - 1
                            : static_cast<std::size_t>(it - fCumMass.begin());
  const std::size_t idx = fSlabIndex[k];
  const RoomSlab& s = gSlabs[idx];

#ifdef ROOMFIELD_BIRTH_STATS
  // Диагностическая сборка (отдельный каталог build/RadiaCode-103-roomfield-check):
  // считает, в какой плите родилось ядро, и печатает доли в конце. Проверка
  // умеет КРАСНЕТЬ: при розыгрыше по объёму доля кирпича вышла бы 0,7109, при
  // правильном розыгрыше по массе — 0,6580 (разница 8 %). Боевой код-путь при
  // выключенном define побайтно тот же: define добавляет только счётчик.
  RoomFieldBirthStatsCount(idx);
#endif

  // 2) точка равномерно по объёму выбранной плиты
  const double x = s.cxMm + (2.0 * G4UniformRand() - 1.0) * s.hxMm;
  const double y = s.cyMm + (2.0 * G4UniformRand() - 1.0) * s.hyMm;
  const double z = s.czMm + (2.0 * G4UniformRand() - 1.0) * s.hzMm;
  fGun.SetParticlePosition(G4ThreeVector(x * mm, y * mm, z * mm));

  // 3) само ядро в покое; всё остальное сделает RDM
  if (!fIon) {
    fIon = G4IonTable::GetIonTable()->GetIon(fIonZ, fIonA, 0.0);
    if (!fIon) {
      std::fprintf(stderr,
                   "Rc103RoomFieldPrimaryGeneratorAction: FATAL ион Z=%d A=%d "
                   "не найден в таблице ионов.\n",
                   fIonZ, fIonA);
      std::abort();
    }
  }
  fGun.SetParticleDefinition(fIon);
  fGun.SetParticleEnergy(0.0);
  fGun.SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));  // покоится
  fGun.GeneratePrimaryVertex(event);
}
