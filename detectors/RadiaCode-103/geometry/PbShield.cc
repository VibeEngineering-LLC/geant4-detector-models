// Реализация ПРЯМОУГОЛЬНОЙ свинцовой защиты. Разбор решений — в PbShield.hh.
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
// Единственное место, где живёт формула посадки полости на дно. И Construct(),
// и внешний код (мюонный блок shieldrun.cc) обязаны спрашивать её здесь, а не
// повторять у себя: расхождение двух копий такой формулы не упало бы ошибкой,
// а тихо сместило бы источник относительно геометрии.
double RCShieldDetector::PlannedZCav() const {
  if (!fSh.seatOnFloor || !fWithDevice) return fSh.zCav;
  const Extent as = AssemblyExtent(fDev, fVes, fWithVessel);
  return as.zLo + fSh.hzCav;
}

double RCShieldDetector::PlannedOuterR() const {
  const double hx = PlannedOuterHx(), hy = PlannedOuterHy();
  return std::sqrt(hx * hx + hy * hy);
}

// ---------------------------------------------------------------------------
void RCShieldDetector::BuildShield(G4LogicalVolume* world) {
  const ShieldGeom& s = fSh;

  // Посадка полости на дно (seatOnFloor) уже сделана в Construct() — там она
  // нужна раньше, чтобы по правильному zCav посчитать размеры мира.
  //
  // --- проверка вместимости: отказ, а не молчаливо усечённая геометрия ------
  // Для стадии 1 (fWithDevice=false) полость пуста по построению — проверять
  // нечего, а AssemblyExtent на голом DeviceGeom дала бы ложный отказ.
  const double zLoCav = s.zCav - s.hzCav, zHiCav = s.zCav + s.hzCav;
  Extent ex{0, zLoCav, zHiCav};
  if (fWithDevice) {
    ex = AssemblyExtent(fDev, fVes, fWithVessel);
    // Сборка в плане КРУГЛАЯ (радиус ex.r), полость ПРЯМОУГОЛЬНАЯ: круг
    // вписывается, если радиус не больше меньшей полуширины.
    const double rFit = std::min(s.hxCav, s.hyCav);
    // Допуск. При seatOnFloor пол полости ВЫЧИСЛЯЕТСЯ как (zLo + hzCav) - hzCav,
    // а такая пара операций в double не возвращает zLo побитово. Без допуска
    // геометрия, где сборка стоит ровно на дне — то есть штатная, — падала бы
    // с SHIELD_TOO_SMALL на последнем бите мантиссы. 1 нм заведомо меньше любого
    // физически осмысленного зазора и заведомо больше ошибки округления на
    // числах порядка сотен миллиметров.
    const double eps = 1e-6;   // мм
    if (ex.r > rFit + eps || ex.zLo < zLoCav - eps || ex.zHi > zHiCav + eps) {
      std::ostringstream m;
      m << "прибор с сосудом не помещается в полость защиты. "
        << "занято: r до " << ex.r << ", z от " << ex.zLo << " до " << ex.zHi
        << " мм; полость: " << 2 * s.hxCav << " x " << 2 * s.hyCav
        << " мм в плане, z от " << zLoCav << " до " << zHiCav << " мм";
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

  // Границы слоя по оси ведём отдельными координатами, а не полувысотой.
  // ПОЧЕМУ. Наращивать h += d с каждым слоем правильно только для ЗАМКНУТОЙ
  // оболочки. При открытом верхе расти вверх нечему: у реального домика
  // стенки кончаются вровень с полостью, сверху ничего нет. Полувысота дала бы
  // симметричный рост в обе стороны, и стенки поднялись бы над полостью на
  // полную толщину свинца — модель выросла бы до 485 мм наружной высоты и
  // 221 кг вместо 435 мм и 210 кг. Поймано сверкой с массой, названной
  // оператором.
  double hx = s.hxCav, hy = s.hyCav;
  double zLo = s.zCav - s.hzCav;      // низ полости, растёт ВНИЗ с каждым слоем
  double zHi = s.zCav + s.hzCav;      // верх: растёт вверх ТОЛЬКО при крышке
  double vPb = 0, vCd = 0, vCu = 0;   // см³, считаются по построенным телам

  for (size_t k = 0; k < layers.size(); ++k) {
    const Layer& L = layers[k];
    const double d = L.d;
    const double hxOut = hx + d, hyOut = hy + d;
    const double zLoOut = zLo - d;
    const double zHiOut = fWithLid ? (zHi + d) : zHi;
    // Стенки занимают по оси РОВНО [zLo, zHi] — высоту полости этого слоя.
    // Дно ложится снаружи снизу, крышка снаружи сверху, поэтому тянуть стенку
    // до zHiOut нельзя: она перекроется с крышкой (поймано самопроверкой
    // объёма на закрытой геометрии — 21900 см³ по телам против 21650
    // аналитических; на открытой ошибка не проявлялась, крышки-то нет).
    const double zc = 0.5 * (zLo + zHi);
    const double h = 0.5 * (zHi - zLo);
    const int depth = static_cast<int>(k);

    auto place = [&](G4VSolid* sol, const G4String& nm, const G4ThreeVector& p) {
      auto* lv = new G4LogicalVolume(sol, L.mat, nm);
      lv->SetVisAttributes(new G4VisAttributes(L.col));
      auto* pv = new G4PVPlacement(nullptr, p, lv, nm, world, false, 0, true);
      fLayerPV.push_back(pv);
      fLayerDepth.push_back(depth);
      return sol->GetCubicVolume() / cm3;
    };

    const char* tag = (L.tag == 'p') ? "Pb" : (L.tag == 'd' ? "Cd" : "Cu");
    auto nm = [&](const char* part) {
      std::ostringstream o;
      o << "sh" << k << "_" << tag << "_" << part;
      return o.str();
    };

    // Стенки по X берутся во всю ширину слоя (hy + d), стенки по Y — только по
    // внутренней (hx): углы тогда закрыты ровно один раз. Симметрию по знаку
    // разворачиваем циклом, чтобы четыре тела не разошлись опечаткой.
    double v = 0;
    for (int sgn = -1; sgn <= 1; sgn += 2) {
      const char* px = (sgn > 0) ? "xhi" : "xlo";
      const char* py = (sgn > 0) ? "yhi" : "ylo";
      v += place(new G4Box(nm(px), 0.5 * d * mm, hyOut * mm, h * mm), nm(px),
                 G4ThreeVector(sgn * (hx + 0.5 * d) * mm, 0, zc * mm));
      v += place(new G4Box(nm(py), hx * mm, 0.5 * d * mm, h * mm), nm(py),
                 G4ThreeVector(0, sgn * (hy + 0.5 * d) * mm, zc * mm));
    }
    v += place(new G4Box(nm("bot"), hxOut * mm, hyOut * mm, 0.5 * d * mm),
               nm("bot"), G4ThreeVector(0, 0, (zLo - 0.5 * d) * mm));
    // Крышка — опционально. Реальный домик оператора собран БЕЗ КРЫШКИ
    // (сообщено 13.08.2026), и это не мелочь: поле помещения входит в полость
    // сверху, не проходя свинец, а космические мюоны и так идут
    // преимущественно сверху (cos²θ) — через открытый верх они попадают в
    // кристалл вообще без свинца. Единственная имеющаяся экспериментальная
    // проверка всей цепочки (измерение «Фон домик 23 дня») относится именно
    // к открытой сверху геометрии, поэтому её надо уметь построить.
    if (fWithLid) {
      v += place(new G4Box(nm("top"), hxOut * mm, hyOut * mm, 0.5 * d * mm),
                 nm("top"), G4ThreeVector(0, 0, (zHi + 0.5 * d) * mm));
    }

    if (L.tag == 'p') vPb += v;
    else if (L.tag == 'd') vCd += v;
    else vCu += v;

    hx = hxOut;
    hy = hyOut;
    zLo = zLoOut;
    zHi = zHiOut;
  }

  fOuterHx = hx;
  fOuterHy = hy;
  fOuterHz = 0.5 * (zHi - zLo);
  fOuterZc = 0.5 * (zHi + zLo);   // центр наружного габарита, НЕ центр полости
  // Описанная окружность — по ней ставится мюонный диск и полуразмер мира.
  // Полуширины здесь НЕДОСТАТОЧНО: углы короба выходят за неё в sqrt(2) раз.
  fOuterR = std::sqrt(hx * hx + hy * hy);
  fNDepth = static_cast<int>(layers.size());

  // Массы — по плотностям материалов Geant4, не по справочным числам
  fMassPb = vPb * (Nist("G4_Pb")->GetDensity() / (g / cm3)) / 1000.0;
  fMassCd = vCd * (Nist("G4_Cd")->GetDensity() / (g / cm3)) / 1000.0;
  fMassCu = vCu * (Nist("G4_Cu")->GetDensity() / (g / cm3)) / 1000.0;

  // САМОПРОВЕРКА ОБЪЁМА. Сумма кусков всех слоёв обязана совпасть с разностью
  // объёмов двух коробов — иначе где-то щель или перекрытие. Проверка дешёвая
  // и ловит именно ту ошибку, которую легче всего сделать в арифметике стыков,
  // а на коробе таких стыков вдвое больше, чем было на цилиндре.
  const double vTotSolid = vPb + vCd + vCu;
  // Наружный габарит уже учитывает открытый верх (вверх слои не росли), а
  // полость вычитается как есть. Отдельной поправки «на непостроенные крышки»
  // больше не нужно — и это правильнее: раньше она латала следствие того, что
  // габарит считался по замкнутой оболочке.
  const double vTotAnal =
      (4 * hx * hy * (zHi - zLo) - 8 * s.hxCav * s.hyCav * s.hzCav) / 1000.0;
  const double rel = (vTotAnal > 0) ? std::abs(vTotSolid / vTotAnal - 1) : 0;
  if (rel > 1e-9) {
    std::ostringstream m;
    m << "объём защиты по телам " << vTotSolid << " см³ против аналитического "
      << vTotAnal << " см³, расхождение " << rel;
    G4Exception("RCShieldDetector::BuildShield", "SHIELD_VOLUME_MISMATCH",
                FatalException, m.str().c_str());
  }

  G4cout << "=== защита ===" << G4endl;
  G4cout << "  полость: " << 2 * s.hxCav << " x " << 2 * s.hyCav << " мм в плане, z "
         << zLoCav << ".." << zHiCav << " мм (высота " << 2 * s.hzCav
         << "; занято прибором: r до " << ex.r << ", z " << ex.zLo << ".."
         << ex.zHi << ")" << G4endl;
  G4cout << "  верх: " << (fWithLid ? "закрыт крышкой" : "ОТКРЫТ") << G4endl;
  G4cout << "  слои от полости наружу: Cu " << s.cu << " + Cd " << s.cd
         << " + Pb " << s.pb << " мм (свинец нарезан на " << s.nShellPb
         << " слоёв), полная стенка " << (s.cu + s.cd + s.pb) << " мм" << G4endl;
  G4cout << "  наружный габарит: " << 2 * fOuterHx << " x " << 2 * fOuterHy
         << " мм в плане, высота " << 2 * fOuterHz << " мм, z "
         << (fOuterZc - fOuterHz) << ".." << (fOuterZc + fOuterHz)
         << " (описанная окружность r " << fOuterR << ")" << G4endl;
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
  // Посадка полости на дно — ДО расчёта мира: при seatOnFloor центр полости
  // уезжает вверх на сотни миллиметров, и мир, посчитанный по старому zCav,
  // оказался бы уже самой защиты. Формула — в PlannedZCav(), не здесь.
  fSh.zCav = PlannedZCav();

  if (fWithShield) {
    const double wall = fSh.cu + fSh.cd + fSh.pb;
    // Мир квадратный в плане, защита тоже — хватает полуширины, углы короба
    // внутрь мира попадают автоматически.
    fWorldHalfXY = std::max(fSh.hxCav, fSh.hyCav) + wall + 30.0;
    // По оси мир симметричен нулю (центру кристалла), а защита — нет. Значит
    // полуразмер надо брать по ДАЛЬШЕЙ из двух её границ, а не по «центр плюс
    // половина»: при посадке на дно защита уходит вверх на ~400 мм и вниз лишь
    // на ~87, и симметричная оценка от центра габарита обрезала бы верх.
    fWorldHalfZ = std::max(std::abs(PlannedOuterZLo()),
                           std::abs(PlannedOuterZHi())) + 30.0;
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
    // размеров, что полость (hxCav x hyCav x hzCav, центр zCav); касается
    // внутренних граней стенок sh0 вплотную, не перекрывается — тот же
    // паттерн стыка, что у слоёв самой защиты.
    auto* cavS = new G4Box("cavity", fSh.hxCav * mm, fSh.hyCav * mm,
                           fSh.hzCav * mm);
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
