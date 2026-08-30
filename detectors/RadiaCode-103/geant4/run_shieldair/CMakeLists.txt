cmake_minimum_required(VERSION 3.16)
project(rc103_shieldair CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

if(DEFINED ENV{XERCESC_ROOT} AND NOT DEFINED XERCESC_ROOT)
  set(XERCESC_ROOT "$ENV{XERCESC_ROOT}" CACHE PATH "Xerces-C++ install root")
endif()
if(NOT DEFINED XERCESC_ROOT)
  message(FATAL_ERROR "Не задан корень Xerces-C++. Укажите XERCESC_ROOT или -DXERCESC_ROOT=<путь>.")
endif()
file(GLOB XERCESC_IMPORT_LIBS "${XERCESC_ROOT}/lib/*.lib")

# Геометрия домика (посадка P-005, свинец, кристалл) ПЕРЕИСПОЛЬЗУЕТСЯ из
# run_field/ теми же .cc файлами - НЕ копия (§33): один источник истины на
# домик, меняется только источник частиц (объёмный вместо изотропной сферы).
set(FIELD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../run_field")

add_executable(rc103_shieldair
  main.cc
  Rc103ShieldAirPrimaryGeneratorAction.cc
  Rc103ShieldAirPhysicsList.cc
  "${FIELD_DIR}/Rc103FieldDetectorConstruction.cc"
  "${FIELD_DIR}/Rc103FieldRunAction.cc"
  "${FIELD_DIR}/Rc103FieldEventAction.cc"
  "${FIELD_DIR}/Rc103FieldSteppingAction.cc"
)
target_include_directories(rc103_shieldair PRIVATE
  "${G4ROOT}/include/Geant4" "${XERCESC_ROOT}/include" "${FIELD_DIR}")
target_link_libraries(rc103_shieldair PRIVATE ${G4_IMPORT_LIBS} ${XERCESC_IMPORT_LIBS})
target_compile_definitions(rc103_shieldair PRIVATE
  G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES
  RC103_FIELD_GDML_DEFAULT_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../gdml/detector/RC103_detector.gdml"
)
if(MSVC)
  target_compile_options(rc103_shieldair PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()