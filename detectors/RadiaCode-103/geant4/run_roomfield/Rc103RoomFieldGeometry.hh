// Единый источник правды о геометрии комнаты: набор непересекающихся плит
// ограждающих конструкций. Один и тот же список используют и построение
// геометрии (Rc103RoomFieldDetectorConstruction), и розыгрыш точки рождения
// ядра (Rc103RoomFieldPrimaryGeneratorAction), и нормировка (main + RunAction).
// Дублировать размеры в трёх местах нельзя: они разъедутся молча.
//
// ⚠ ПРИВАТНОСТЬ: реальные размеры комнаты оператора в исходники НЕ зашиты.
// Все дефолты ниже — НЕЙТРАЛЬНЫЕ (4000x4000x2800 мм, стены 500 мм), реальные
// значения передаются ключами командной строки при запуске.
//
// Разбиение на плиты (непересекающееся, углы отданы стенам по X):
//   floor   : полный наружный след x[-Xm..+Xp], y[-Ym..+Yp], z ниже пола
//   ceiling : то же, выше потолка
//   wall_xm : x слева от комнаты, y — полный наружный след, z в высоту комнаты
//   wall_xp : x справа, y — полный наружный след, z в высоту комнаты
//   wall_ym : y снизу, x — ТОЛЬКО внутренняя ширина комнаты (углы уже заняты)
//   wall_yp : y сверху, x — ТОЛЬКО внутренняя ширина комнаты
#pragma once

#include <string>
#include <vector>

struct RoomSlab {
  std::string name;
  double cxMm = 0, cyMm = 0, czMm = 0;  // центр в мировых координатах, мм
  double hxMm = 0, hyMm = 0, hzMm = 0;  // полуразмеры, мм
  bool brick = true;                    // true = кирпич, false = бетон
  double volumeCm3 = 0;                 // 8*hx*hy*hz / 1000
  double densityGCm3 = 0;               // назначается AssignRoomDensities()
  double massG = 0;                     // volumeCm3 * densityGCm3
};

struct RoomConfig {
  // Внутренние («в свету») размеры комнаты, мм. НЕЙТРАЛЬНЫЙ дефолт.
  double innerXMm = 4000.0;
  double innerYMm = 4000.0;
  double innerZMm = 2800.0;
  // Толщины ограждений, мм. НЕЙТРАЛЬНЫЕ дефолты.
  double wallXmMm = 500.0;
  double wallXpMm = 500.0;
  double wallYmMm = 500.0;
  double wallYpMm = 500.0;
  double floorMm = 200.0;
  double ceilMm = 200.0;
  // Точка наблюдения: от ВНУТРЕННЕЙ грани стены X-, от внутренней грани Y-,
  // от верха пола. НЕЙТРАЛЬНЫЙ дефолт — по 1000 мм от угла и 1000 мм высоты.
  double obsDxMm = 1000.0;
  double obsDyMm = 1000.0;
  double obsDzMm = 1000.0;
  // Радиус воздушного шара-скорера.
  double ballRMm = 300.0;
  // Плотность кирпича, г/см3 (керамический полнотелый, типично 1,6..1,9).
  double rhoBrickGCm3 = 1.8;
  // Эффективная плотность перекрытий, г/см3. 0 = взять плотность NIST
  // G4_CONCRETE как есть. Для железобетонных КРУГЛОПУСТОТНЫХ плит (ГОСТ 9561)
  // задаётся ГОМОГЕНИЗИРОВАННОЕ значение: паспортная масса, делённая на
  // габаритный объём. Явные цилиндрические пустоты не моделируются.
  double rhoSlabGCm3 = 0.0;
  // Продолжение среды наружу: к КАЖДОЙ стене и КАЖДОМУ перекрытию снаружи
  // добавляется слой ТОГО ЖЕ материала такой толщины (мм). 0 = изолированная
  // комната. Смысл — верхняя оценка вклада остального здания одним прогоном.
  double extendMm = 0.0;
  // Ограничение розыгрыша точки рождения по материалу.
  enum SrcMode { kSrcAll = 0, kSrcBrick = 1, kSrcConcrete = 2 };
  int srcMode = kSrcAll;
  // Порог продукции, мм.
  double cutMm = 1.0;
  // Запас воздуха вокруг наружной границы конструкций, мм.
  double worldMarginMm = 100.0;
};

extern RoomConfig gRoom;
extern std::vector<RoomSlab> gSlabs;

// Пересчитывает gSlabs из gRoom. Идемпотентна.
void BuildRoomSlabs();
// Проставляет плотности и массы. rhoConcrete берётся ФАКТОМ у материала NIST,
// а не переписывается числом.
void AssignRoomDensities(double rhoBrickGCm3, double rhoConcreteGCm3);

// Внутренние полуразмеры комнаты (воздух), мм.
double RoomInnerHalfXMm();
double RoomInnerHalfYMm();
double RoomInnerHalfZMm();
// Полуразмеры мира, мм.
double RoomWorldHalfXMm();
double RoomWorldHalfYMm();
double RoomWorldHalfZMm();
// Точка наблюдения в мировых координатах, мм.
double RoomObsXMm();
double RoomObsYMm();
double RoomObsZMm();

double RoomBrickVolumeCm3();
double RoomConcreteVolumeCm3();
double RoomBrickMassG();
double RoomConcreteMassG();

// Эффективные толщины с учётом gRoom.extendMm.
double RoomEffWallXmMm();
double RoomEffWallXpMm();
double RoomEffWallYmMm();
double RoomEffWallYpMm();
double RoomEffFloorMm();
double RoomEffCeilMm();

// Плита участвует в розыгрыше при текущем srcMode?
bool RoomSlabSelected(const RoomSlab& s);
// Масса ТОЛЬКО выбранного материала — по ней и только по ней считается
// нормировка при src=brick|concrete. Смешать её с полной массой значило бы
// выдать шаблоны в разных единицах, и складывать их было бы нельзя.
double RoomSelectedMassG();
const char* RoomSrcModeName();
