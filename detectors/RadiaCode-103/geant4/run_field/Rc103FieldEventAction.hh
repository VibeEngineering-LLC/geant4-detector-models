// Копит депонированную в кристалле энергию за одно событие
// (SteppingAction::AddEdep), в конце события сдаёт сумму в RunAction.
#pragma once

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <array>
#include <map>

class Rc103FieldRunAction;
class G4Event;

class Rc103FieldEventAction : public G4UserEventAction {
 public:
  // Категории происхождения счёта (29.08.2026): 0 — фотон дошёл до кристалла,
  // не заходя в свинец домика; 1 — прошёл свинец, рассеявшись в нём; 2 —
  // родился в свинце (флуоресценция K-серии, тормозное, вторичные).
  static constexpr int kCatDirect = 0;
  static constexpr int kCatPbScat = 1;
  static constexpr int kCatPbBorn = 2;
  static constexpr int kNCat = 3;

  explicit Rc103FieldEventAction(Rc103FieldRunAction* runAction);
  ~Rc103FieldEventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void AddEdep(G4double edep) { fEdep += edep; }
  void AddEdepByOrigin(G4double edep, int category) {
    if (category >= 0 && category < kNCat) fEdepCat[category] += edep;
  }
  // Категория ТРЕКА, с наследованием от родителя. Так надо потому, что энергию
  // в кристалле оставляет не фотон, а рождённый им ЗДЕСЬ ЖЕ электрон: его
  // собственная история со свинцом не связана, и классификация по
  // депонирующему треку всегда давала бы «мимо свинца» (проверено прогоном
  // 29.08.2026 — обе свинцовые категории выходили нулевыми).
  int GetCategory(G4int trackID) const {
    auto it = fTrackCat.find(trackID);
    return (it == fTrackCat.end()) ? -1 : it->second;
  }
  void SetCategory(G4int trackID, int category) { fTrackCat[trackID] = category; }

 private:
  Rc103FieldRunAction* fRunAction;
  G4double fEdep = 0.0;
  std::array<G4double, kNCat> fEdepCat{};
  // ⚠ Очищается в начале КАЖДОГО события: иначе таблица растёт на весь прогон
  // и при 3e8 историй съедает память.
  std::map<G4int, int> fTrackCat;
};
