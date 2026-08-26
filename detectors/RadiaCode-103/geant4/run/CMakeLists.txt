cmake_minimum_required(VERSION 3.16)
project(rc103_run CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Тот же рецепт direct-link, что common/tools/CMakeLists.txt (эталон) —
# не find_package(Geant4). Собственный CMakeLists.txt, common/tools не тронут.
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR
    "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

# G4GDMLParser.hh -> G4GDMLRead.hh -> xercesc/... — заголовки Xerces-C++ не
# копируются в ${G4ROOT}/include/Geant4, нужен отдельный include path.
if(DEFINED ENV{XERCESC_ROOT} AND NOT DEFINED XERCESC_ROOT)
  set(XERCESC_ROOT "$ENV{XERCESC_ROOT}" CACHE PATH "Xerces-C++ install root (нужен для GDML-заголовков)")
endif()
if(NOT DEFINED XERCESC_ROOT)
  message(FATAL_ERROR
    "Не задан корень Xerces-C++ (нужен для GDML-заголовков). Укажите XERCESC_ROOT или -DXERCESC_ROOT=<путь>.")
endif()
file(GLOB XERCESC_IMPORT_LIBS "${XERCESC_ROOT}/lib/*.lib")

add_executable(rc103_run
  main.cc
  Rc103RunDetectorConstruction.cc
  Rc103RunPhysicsList.cc
  Rc103RunPrimaryGeneratorAction.cc
  Rc103RunRunAction.cc
  Rc103RunEventAction.cc
  Rc103RunSteppingAction.cc
)
target_include_directories(rc103_run PRIVATE "${G4ROOT}/include/Geant4" "${XERCESC_ROOT}/include")
target_link_libraries(rc103_run PRIVATE ${G4_IMPORT_LIBS} ${XERCESC_IMPORT_LIBS})
target_compile_definitions(rc103_run PRIVATE
  G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES
  RC103_GDML_DEFAULT_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../gdml/detector/RC103_detector.gdml"
  RC103_OUTPUT_DEFAULT_PATH="${CMAKE_CURRENT_SOURCE_DIR}/output/rc103_cs137_100mm_run.csv"
)
if(MSVC)
  target_compile_options(rc103_run PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()
