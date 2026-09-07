cmake_minimum_required(VERSION 3.16)
project(g1s_npsm CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

# Этап 5: полная геометрия Гамма-1С + взвешивание светового выхода стенда.
# Все три источника подключены ВКЛЮЧЕНИЕМ, а не копией (§33): правка донора
# расходится на всех потребителей, дублирующей реализации не заводится.
set(GEOM_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/../../geometry")
set(BENCH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../run_npsm_bench")
set(PHYS_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/../../../RadiaCode-103/geant4/run_field")

add_executable(g1s_npsm
  main.cc
  G1sNpsmPrimaryGeneratorAction.cc
  G1sNpsmActionInitialization.cc
  "${GEOM_DIR}/G1SDetector.cc"
  "${BENCH_DIR}/NpsmBenchEventAction.cc"
  "${BENCH_DIR}/NpsmBenchRunAction.cc"
  "${BENCH_DIR}/NpsmBenchRun.cc"
  "${BENCH_DIR}/NpsmBenchStackingAction.cc"
  "${BENCH_DIR}/NpsmBenchActionInitialization.cc"
  "${BENCH_DIR}/NpsmBenchSteppingAction.cc"
  "${BENCH_DIR}/NpsmBenchDetectorConstruction.cc"
  "${BENCH_DIR}/NpsmBenchPrimaryGeneratorAction.cc"
  "${BENCH_DIR}/NpsmLightYield.cc"
  "${PHYS_DIR}/Rc103FieldPhysicsList.cc"
)
# Два файла стенда входят в сборку не ради своей работы, а ради символов:
#  * NpsmBenchDetectorConstruction.cc — NpsmBenchSteppingAction обращается к его
#    статическому методу как к ЗАПАСНОМУ источнику кристалла (здесь кристалл
#    задаётся явно, поэтому запасной путь не используется);
#  * NpsmBenchPrimaryGeneratorAction.cc — NpsmBenchRunAction::WriteCSV печатает
#    в шапку CSV его статики gParticle/gPencilBeam; без этого файла линк падает
#    с LNK2019 (проверено 06.09).
# Ни геометрия стенда, ни его генератор при этом не строятся и в RunManager не
# передаются: объекты не создаются.

target_include_directories(g1s_npsm PRIVATE
  "${G4ROOT}/include/Geant4" "${GEOM_DIR}" "${BENCH_DIR}" "${PHYS_DIR}")
target_link_libraries(g1s_npsm PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(g1s_npsm PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if (MSVC)
  target_compile_options(g1s_npsm PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()
