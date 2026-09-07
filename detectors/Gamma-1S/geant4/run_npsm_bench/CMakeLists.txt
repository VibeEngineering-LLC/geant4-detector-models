cmake_minimum_required(VERSION 3.16)
project(npsm_bench CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")
set(DONOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../RadiaCode-103/geant4/run_field")
add_executable(npsm_bench
  main.cc
  NpsmBenchDetectorConstruction.cc
  NpsmBenchPrimaryGeneratorAction.cc
  NpsmBenchEventAction.cc
  NpsmBenchRunAction.cc
  NpsmBenchRun.cc
  NpsmBenchStackingAction.cc
  NpsmBenchActionInitialization.cc
  NpsmBenchSteppingAction.cc
  NpsmLightYield.cc
  "${DONOR_DIR}/Rc103FieldPhysicsList.cc"
)
target_include_directories(npsm_bench PRIVATE "${G4ROOT}/include/Geant4" "${DONOR_DIR}")
target_link_libraries(npsm_bench PRIVATE ${G4_IMPORT_LIBS})
