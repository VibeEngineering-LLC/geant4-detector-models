#include "Rc103RunPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4UImanager.hh"

// --- Геометрия источника (пересчитано из RC103_detector.gdml, 26.08.2026,
//     САМОСТОЯТЕЛЬНО — не взято из текста задания) ---
//
// Кристалл (Crystal_log, куб 10x10x10 мм, RC103_detector.gdml:130,175-178):
// мировые координаты его центра = (-49.5, 0, -0.55) мм. Цепочка вложения без
// поворотов, каждый physvol в p_origin=(0,0,0) кроме одного:
//   World -> pv_rc103 @ p_origin -> RC103_device_log
//   RC103_device_log -> pv_case_interior @ p_origin -> Case_interior_log
//   Case_interior_log -> pv_detector @ p_detector=(-49.5,0,-0.55) -> DetectorModule_log
//   DetectorModule_log -> pv_capsule_cavity @ p_origin -> Capsule_cavity_log
//   Capsule_cavity_log -> pv_crystal @ p_origin -> Crystal_log
// (GDML: <position name="p_detector" x="-49.5" y="0" z="-0.55"/>, строка 17;
//  <physvol name="pv_crystal"><positionref ref="p_origin"/>, строки 202-205).
// Итог: мировые координаты кристалла = p_detector дословно, без накопленных
// поворотов/смещений — совпадение с обозначением в GDML-комментарии
// "Кристалл CsI(Tl) 10^3 @ (-49.5,0,-0.55)" (строка 9) подтверждает расчёт.
//
// Наружная грань корпуса по -Z: Case_outer x=123,y=34,z=17.5 мм
// (RC103_detector.gdml:115), Case_shell_log стоит в RC103_device_log в
// p_origin, тот — в World тоже в p_origin => центр корпуса = мировой
// (0,0,0) => полутолщина по Z = 17.5/2 = 8.75 мм => наружная грань
// z = -8.75 мм.
//
// Источник — на оси кристалла (x=-49.5 мм, y=0 мм, совпадает с осью X
// кристалла, так что линия источник->центр кристалла параллельна Z), в
// 100 мм от наружной грани корпуса НАРУЖУ:
//   z_src = -8.75 - 100 = -108.75 мм
// Расстояние источник -> центр кристалла: |-0.55 - (-108.75)| = 108.2 мм.
//
// Узкий конус вместо полного 4π iso — чтобы не тратить статистику на
// кванты, летящие мимо прибора. Требуемый полураствол — угол к САМОЙ
// удалённой от оси точке куба, видимой из источника; для точечного
// источника угол максимален у БЛИЖНЕЙ (по Z) грани куба (меньше расстояние
// => больше угол при том же радиальном сдвиге):
//   ближняя грань кристалла: z = -0.55 - 5 = -5.55 мм,
//   расстояние источник->ближняя грань: |-5.55-(-108.75)| = 103.2 мм,
//   полудиагональ ближней грани (куб 10х10): sqrt(5^2+5^2) = 7.0711 мм,
//   требуемый угол: atan(7.0711/103.2) = 3.92°.
// Берём конус полураствора 5° (запас ~27% на неточности осевого
// совмещения) — гарантированно покрывает весь куб (дальняя грань, z=+4.45,
// расстояние 113.2 мм, тот же радиальный сдвиг => ещё меньший угол 3.58°,
// то есть ближняя грань — действительно худший случай).
//
// GPS-конвенция направления (G4SPSAngDistribution::GenerateOne, без
// поворота Rot1/Rot2 по умолчанию — локальная система = мировая):
//   Px=-sin(theta)cos(phi), Py=-sin(theta)sin(phi), Pz=-cos(theta).
// При theta=0 => направление (0,0,-1); при theta=180° => Pz=-cos(180°)=+1,
// направление (0,0,+1). Источник должен светить в +Z (к кристаллу, который
// находится при БОЛЕЕ БОЛЬШОМ z, чем источник) => берём
// mintheta=175°, maxtheta=180° (конус полураствора 5° вокруг +Z).
//
// Телесный угол конуса: Omega_cone = 2π(1-cos5°) = 0.023911 ср
// (Omega_cone/4π = 0.0019027, т.е. 0.19027% от полной сферы) — используется
// в main.cc/RunAction для пересчёта эффективности «внутри конуса» (что даёт
// прямой BeamOn-счётчик) в АБСОЛЮТНУЮ эффективность (относительно
// изотропного источника 4π), нужную для сравнения с независимой
// аналитической оценкой (см. отчёт).

Rc103RunPrimaryGeneratorAction::Rc103RunPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fGPS(nullptr) {
  fGPS = new G4GeneralParticleSource();

  G4UImanager* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/gps/particle gamma");
  ui->ApplyCommand("/gps/pos/type Point");
  ui->ApplyCommand("/gps/pos/centre -49.5 0 -108.75 mm");
  ui->ApplyCommand("/gps/ang/type iso");
  ui->ApplyCommand("/gps/ang/mintheta 175 deg");
  ui->ApplyCommand("/gps/ang/maxtheta 180 deg");
  ui->ApplyCommand("/gps/ene/type Mono");
  ui->ApplyCommand("/gps/ene/mono 661.657 keV");
}

Rc103RunPrimaryGeneratorAction::~Rc103RunPrimaryGeneratorAction() { delete fGPS; }

void Rc103RunPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  fGPS->GeneratePrimaryVertex(event);
}
