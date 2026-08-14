// Реализация цилиндрической свинцовой защиты. Разбор решений — в PbShield.hh.
#include "PbShield.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace {

// Габариты прибора и сосуда в мировой системе (центр кристалла = 0).
// Выводятся ИЗ ТЕХ ЖЕ полей, которыми строится геометрия, а не набираются
// руками: иначе проверка вместимости проверяла бы не то, что построено.
struct Extent { double r, zLo, zHi; };

Extent AssemblyExtent(const DeviceGeom& d, const VesselGeom& v, bool withVessel) {
  // прибор: сечение caseX x caseY, нос в -crystalZ0, длина caseZ
  double r = 0.5 * std::sqrt(d.caseX * d.caseX + d.caseY * d.caseY);
  double zLo = -d.crystalZ0;
  double zHi = -d.crystalZ0 + d.caseZ;
  if (withVessel) {
    const double zSlot = -d.crystalZ0 + v.wellTip - v.seatGap;
    const double zRim  = zSlot - v.barrelH;
    // радиус: наружный обвод стакана либо юбка крышки, что больше
    double rv = (v.outerR > 0) ? v.outerR : 0.0;
    for (double x : v.profR) rv = std::max(rv, x);
    rv = std::max(rv, v.capSkirtR);
    r = std::max(r, rv);
    zLo = std::min(zLo, zRim - v.capT);          // низ диска крышки
    zHi = std::max(zHi, zSlot);
  }
  return {r, zLo, zHi};
}

G4Material* Nist(const G4String& n) {
  return G4NistManager::Instance()->FindOrBuildMaterial(n);
}

}  // namespace

// ---------------------------------------------------------------------------
RCShieldDetector::RCShieldDetector(bool withVessel, bool withShield,
                                   bool withDevice)
    : RCDetector(withVessel), fWithShield(withShield),
      fWithDevice(withDevice) {}

// ---------------------------------------------------------------------------
void RCShieldDetector::BuildShield(G4LogicalVolume* world) {
  const ShieldGeom& s = fSh;

  // --- проверка вместимости: отказ, а не молчаливо усечённая геометрия ------
  // Для стадии 1 (fWithDevice=false) полость пуста по построению — проверять
  // нечего, а AssemblyExtent на голом DeviceGeom дала бы ложный отказ.
  const double zLoCav = s.zCav - s.hzCav, zHiCav = s.zCav + s.hzCav;
  Extent ex{0, zLoCav, zHiCav};
  if (fWithDevice) {
    ex = AssemblyExtent(fDev, fVes, fWithVessel);
    if (ex.r > s.rCav || ex.zLo < zLoCav || ex.zHi > zHiCav) {
      std::ostringstream m;
      m << "прибор с сосудом не помещается в полость защиты. "
        << "занято: r до " << ex.r << ", z от " << ex.zLo << " до " << ex.zHi
        << " мм; полость: r " << s.rCav << ", z от " << zLoCav << " до "
        << zHiCav << " мм";
      G4Exception("RCShieldDetector::BuildShield", "SHIELD_TOO_SMALL",
                  FatalException, m.str().c_str());
    }
  }

  // --- список слоёв от полости наружу --------------------------------------
  struct Layer { double d; G4Material* mat; char tag; G4Colour col; };
  std::vector<Layer> layers;
  if (s.cu > 0)
    layers.push_back({s.cu, Nist("G4_Cu"), 'u', G4Colour(0.72, 0.45, 0.20)});
  if (s.cd > 0)
    layers.push_back({s.cd, Nist("G4_Cd"), 'd', G4Colour(0.55, 0.55, 0.60)});
  if (s.pb > 0) {
    const int n = (s.nShellPb > 0) ? s.nShellPb : 1;
    const double d = s.pb / n;
    for (int k = 0; k < n; ++k)
      layers.push_back({d, Nist("G4_Pb"), 'p', G4Colour(0.35, 0.35, 0.40)});
  }

  double r = s.rCav, h = s.hzCav;
  double vPb = 0, vCd = 0, vCu = 0;   // см³, считаются по построенным телам
  double vLidSkipped = 0;             // см³ непостроенных крышек (fWithLid=false)

  for (size_t k = 0; k < layers.size(); ++k) {
    const Layer& L = layers[k];
    const double rOut = r + L.d, hOut = h + L.d;
    const int depth = static_cast<int>(k);

    auto place = [&](G4VSolid* sol, const G4String& nm, double z) {
      auto* lv = new G4LogicalVolume(sol, L.mat, nm);
      lv->SetVisAttributes(new G4VisAttributes(L.col));
      auto* pv = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z * mm), lv, nm,
                                   world, false, 0, true);
      fLayerPV.push_back(pv);
      fLayerDepth.push_back(depth);
      return sol->GetCubicVolume() / cm3;
    };

    const char* tag = (L.tag == 'p') ? "Pb" : (L.tag == 'd' ? "Cd" : "Cu");
    std::ostringstream side, bot, top;
    side << "sh" << k << "_" << tag << "_side";
    bot  << "sh" << k << "_" << tag << "_bot";
    top  << "sh" << k << "_" << tag << "_top";

    double v = 0;
    v += place(new G4Tubs(side.str(), r * mm, rOut * mm, h * mm, 0., twopi),
               side.str(), s.zCav);
    v += place(new G4Tubs(bot.str(), 0., rOut * mm, 0.5 * L.d * mm, 0., twopi),
               bot.str(), s.zCav - h - 0.5 * L.d);
    // Крышка — опционально. Реальный домик оператора собран БЕЗ КРЫШКИ
    // (сообщено 13.08.2026), и это не мелочь: поле помещения входит в полость
    // сверху, не проходя свинец, а космические мюоны и так идут
    // преимущественно сверху (cos²θ) — через открытый верх они попадают в
    // кристалл вообще без свинца. Единственная имеющаяся экспериментальная
    // проверка всей цепочки (измерение «Фон домик 23 дня») относится именно
    // к открытой сверху геометрии, поэтому её надо уметь построить.
    if (fWithLid) {
      v += place(new G4Tubs(top.str(), 0., rOut * mm, 0.5 * L.d * mm, 0., twopi),
                 top.str(), s.zCav + h + 0.5 * L.d);
    } else {
      vLidSkipped += pi * rOut * rOut * L.d / 1000.0;   // мм³ -> см³
    }

    if (L.tag == 'p') vPb += v;
    else if (L.tag == 'd') vCd += v;
    else vCu += v;

    r = rOut;
    h = hOut;
  }

  fOuterR = r;
  fOuterHz = h;
  fNDepth = static_cast<int>(layers.size());

  // Массы — по плотностям материалов Geant4, не по справочным числам
  fMassPb = vPb * (Nist("G4_Pb")->GetDensity() / (g / cm3)) / 1000.0;
  fMassCd = vCd * (Nist("G4_Cd")->GetDensity() / (g / cm3)) / 1000.0;
  fMassCu = vCu * (Nist("G4_Cu")->GetDensity() / (g / cm3)) / 1000.0;

  // САМОПРОВЕРКА ОБЪЁМА. Сумма трёх кусков каждого слоя обязана совпасть с
  // разностью объёмов двух цилиндров — иначе где-то щель или перекрытие.
  // Проверка дешёвая и ловит именно ту ошибку, которую легче всего сделать
  // в арифметике стыков.
  const double vTotSolid = vPb + vCd + vCu;
  // Без крышки из аналитического кольца вычитается объём непостроенных
  // крышек — иначе проверка даёт ложный отказ на открытой геометрии.
  const double vTotAnal =
      (pi * r * r * 2 * h - pi * s.rCav * s.rCav * 2 * s.hzCav) / 1000.0
      - vLidSkipped;
  const double rel = (vTotAnal > 0) ? std::abs(vTotSolid / vTotAnal - 1) : 0;
  if (rel > 1e-9) {
    std::ostringstream m;
    m << "объём защиты по телам " << vTotSolid << " см³ против аналитического "
      << vTotAnal << " см³, расхождение " << rel;
    G4Exception("RCShieldDetector::BuildShield", "SHIELD_VOLUME_MISMATCH",
                FatalException, m.str().c_str());
  }

  G4cout << "=== защита ===" << G4endl;
  G4cout << "  полость: r " << s.rCav << ", z " << zLoCav << ".." << zHiCav
         << " мм (занято прибором: r до " << ex.r << ", z " << ex.zLo << ".."
         << ex.zHi << ")" << G4endl;
  G4cout << "  слои от полости наружу: Cu " << s.cu << " + Cd " << s.cd
         << " + Pb " << s.pb << " мм (свинец нарезан на " << s.nShellPb
         << " слоёв), полная стенка " << (s.cu + s.cd + s.pb) << " мм" << G4endl;
  G4cout << "  наружный габарит: r " << fOuterR << ", полувысота " << fOuterHz
         << " мм" << G4endl;
  G4cout << "  МАССЫ: Pb " << fMassPb << " кг, Cd " << fMassCd << " кг, Cu "
         << fMassCu << " кг" << G4endl;
  G4cout << "  объём защиты: " << vTotSolid << " см³ (аналитически "
         << vTotAnal << ")" << G4endl;
}

// ---------------------------------------------------------------------------
G4VPhysicalVolume* RCShieldDetector::Construct() {
  // Мир обязан вместить защиту и поверхность розыгрыша снаружи неё.
  // Запас 30 мм: поверхность источника ставится вплотную к габариту защиты,
  // и ей нужно место, чтобы не совпасть с гранью мира.
  if (fWithShield) {
    const double wall = fSh.cu + fSh.cd + fSh.pb;
    fWorldHalfXY = fSh.rCav + wall + 30.0;
    fWorldHalfZ  = std::abs(fSh.zCav) + fSh.hzCav + wall + 30.0;
  }
  // Мир не может быть УЖЕ поверхности розыгрыша (см. fWorldMinHalfXY в .hh):
  // иначе часть источника молча оказывается снаружи мира.
  if (fWorldMinHalfXY > fWorldHalfXY) fWorldHalfXY = fWorldMinHalfXY;

  DefineMaterials();

  auto* worldS = new G4Box("world", fWorldHalfXY * mm, fWorldHalfXY * mm,
                           fWorldHalfZ * mm);
  auto* worldLV = new G4LogicalVolume(worldS, Nist("G4_AIR"), "world");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPV =
      new G4PVPlacement(nullptr, {}, worldLV, "world", nullptr, false, 0, true);

  // withDevice=false — стадия 1: полость остаётся пустой. Сосуд без прибора
  // внутри не строится ни при каких fWithVessel (см. пояснение в .hh).
  if (fWithDevice) {
    BuildDevice(worldLV);
    if (fWithVessel) BuildVessel(worldLV);
  } else {
    // Именованная "cavity" — см. разбор в .hh. Воздух, точно тех же
    // размеров, что полость (r=rCav, halfz=hzCav, центр zCav); касается
    // внутренней грани sh0 (rmin=rCav) вплотную, не перекрывается —
    // тот же паттерн стыка, что у слоёв самой защиты.
    auto* cavS = new G4Tubs("cavity", 0., fSh.rCav * mm, fSh.hzCav * mm, 0., twopi);
    auto* cavLV = new G4LogicalVolume(cavS, Nist("G4_AIR"), "cavity");
    cavLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    fCavityPV = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fSh.zCav * mm),
                                  cavLV, "cavity", worldLV, false, 0, true);
  }
  if (fWithShield) BuildShield(worldLV);

  G4cout << "\n=== RadiaCode + свинцовая защита ===" << G4endl;
  G4cout << "  мир: " << fWorldHalfXY << " x " << fWorldHalfXY << " x "
         << fWorldHalfZ << " мм (полуразмеры)" << G4endl;
  if (!fWithDevice) G4cout << "  прибор не построен: полость пуста (стадия 1)"
                           << G4endl;
  if (fWithDevice && fWithVessel)
    G4cout << "  сосуд " << fVes.name << ", проба " << fVes.sampleMatrix << " "
           << fVes.sampleDensity << " г/см³, объём " << fSampleVolumeCm3
           << " см³" << G4endl;
  if (!fWithShield) G4cout << "  защита не построена" << G4endl;
  G4cout << G4endl;
  return worldPV;
}
