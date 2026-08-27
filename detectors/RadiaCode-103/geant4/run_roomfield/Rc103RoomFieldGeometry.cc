#include "Rc103RoomFieldGeometry.hh"

#include <algorithm>

RoomConfig gRoom;
std::vector<RoomSlab> gSlabs;

namespace {
// Плита задаётся диапазонами по трём осям в мировых координатах (мм) —
// так разбиение читается глазом и проверяется на непересечение арифметикой.
void PushSlab(const char* name, bool brick, double xmin, double xmax,
              double ymin, double ymax, double zmin, double zmax) {
  RoomSlab s;
  s.name = name;
  s.brick = brick;
  s.cxMm = 0.5 * (xmin + xmax);
  s.cyMm = 0.5 * (ymin + ymax);
  s.czMm = 0.5 * (zmin + zmax);
  s.hxMm = 0.5 * (xmax - xmin);
  s.hyMm = 0.5 * (ymax - ymin);
  s.hzMm = 0.5 * (zmax - zmin);
  s.volumeCm3 = 8.0 * s.hxMm * s.hyMm * s.hzMm / 1000.0;  // мм3 -> см3
  gSlabs.push_back(s);
}
}  // namespace

double RoomInnerHalfXMm() { return 0.5 * gRoom.innerXMm; }
double RoomInnerHalfYMm() { return 0.5 * gRoom.innerYMm; }
double RoomInnerHalfZMm() { return 0.5 * gRoom.innerZMm; }

double RoomEffWallXmMm() { return gRoom.wallXmMm + gRoom.extendMm; }
double RoomEffWallXpMm() { return gRoom.wallXpMm + gRoom.extendMm; }
double RoomEffWallYmMm() { return gRoom.wallYmMm + gRoom.extendMm; }
double RoomEffWallYpMm() { return gRoom.wallYpMm + gRoom.extendMm; }
double RoomEffFloorMm() { return gRoom.floorMm + gRoom.extendMm; }
double RoomEffCeilMm() { return gRoom.ceilMm + gRoom.extendMm; }

double RoomWorldHalfXMm() {
  return RoomInnerHalfXMm() +
         std::max(RoomEffWallXmMm(), RoomEffWallXpMm()) + gRoom.worldMarginMm;
}
double RoomWorldHalfYMm() {
  return RoomInnerHalfYMm() +
         std::max(RoomEffWallYmMm(), RoomEffWallYpMm()) + gRoom.worldMarginMm;
}
double RoomWorldHalfZMm() {
  return RoomInnerHalfZMm() + std::max(RoomEffFloorMm(), RoomEffCeilMm()) +
         gRoom.worldMarginMm;
}

double RoomObsXMm() { return -RoomInnerHalfXMm() + gRoom.obsDxMm; }
double RoomObsYMm() { return -RoomInnerHalfYMm() + gRoom.obsDyMm; }
double RoomObsZMm() { return -RoomInnerHalfZMm() + gRoom.obsDzMm; }

void BuildRoomSlabs() {
  gSlabs.clear();

  const double hx = RoomInnerHalfXMm();
  const double hy = RoomInnerHalfYMm();
  const double hz = RoomInnerHalfZMm();

  // Наружные границы следа здания в плане. Толщины ЭФФЕКТИВНЫЕ: ключ
  // extend= добавляет снаружи каждой стены и каждого перекрытия слой того же
  // материала, что и есть «продолжение дома за ограждением».
  const double outXmin = -hx - RoomEffWallXmMm();
  const double outXmax = +hx + RoomEffWallXpMm();
  const double outYmin = -hy - RoomEffWallYmMm();
  const double outYmax = +hy + RoomEffWallYpMm();

  // Перекрытия — на полный наружный след (плита проходит и под стенами).
  PushSlab("floor_concrete", false, outXmin, outXmax, outYmin, outYmax,
           -hz - RoomEffFloorMm(), -hz);
  PushSlab("ceiling_concrete", false, outXmin, outXmax, outYmin, outYmax, +hz,
           +hz + RoomEffCeilMm());

  // Стены по X — на полную наружную ширину по Y, то есть углы плана
  // принадлежат им (материал в углах всё равно кирпич, неоднозначности нет).
  PushSlab("wall_xm_brick", true, outXmin, -hx, outYmin, outYmax, -hz, +hz);
  PushSlab("wall_xp_brick", true, +hx, outXmax, outYmin, outYmax, -hz, +hz);

  // Стены по Y — только на внутреннюю ширину комнаты: углы уже отданы выше.
  PushSlab("wall_ym_brick", true, -hx, +hx, outYmin, -hy, -hz, +hz);
  PushSlab("wall_yp_brick", true, -hx, +hx, +hy, outYmax, -hz, +hz);
}

bool RoomSlabSelected(const RoomSlab& s) {
  switch (gRoom.srcMode) {
    case RoomConfig::kSrcBrick:
      return s.brick;
    case RoomConfig::kSrcConcrete:
      return !s.brick;
    default:
      return true;
  }
}

double RoomSelectedMassG() {
  double m = 0;
  for (const auto& s : gSlabs)
    if (RoomSlabSelected(s)) m += s.massG;
  return m;
}

const char* RoomSrcModeName() {
  switch (gRoom.srcMode) {
    case RoomConfig::kSrcBrick:
      return "brick";
    case RoomConfig::kSrcConcrete:
      return "concrete";
    default:
      return "all";
  }
}

void AssignRoomDensities(double rhoBrickGCm3, double rhoConcreteGCm3) {
  for (auto& s : gSlabs) {
    s.densityGCm3 = s.brick ? rhoBrickGCm3 : rhoConcreteGCm3;
    s.massG = s.volumeCm3 * s.densityGCm3;
  }
}

double RoomBrickVolumeCm3() {
  double v = 0;
  for (const auto& s : gSlabs)
    if (s.brick) v += s.volumeCm3;
  return v;
}
double RoomConcreteVolumeCm3() {
  double v = 0;
  for (const auto& s : gSlabs)
    if (!s.brick) v += s.volumeCm3;
  return v;
}
double RoomBrickMassG() {
  double m = 0;
  for (const auto& s : gSlabs)
    if (s.brick) m += s.massG;
  return m;
}
double RoomConcreteMassG() {
  double m = 0;
  for (const auto& s : gSlabs)
    if (!s.brick) m += s.massG;
  return m;
}
