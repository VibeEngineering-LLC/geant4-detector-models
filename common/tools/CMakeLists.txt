cmake_minimum_required(VERSION 3.16)
project(gdml_vis CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Тот же рецепт прямой линковки, что у RC-103/geometry/CMakeLists.txt —
# direct-link списком import-lib, не find_package(Geant4) (см. комментарий-
# обоснование там же).
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR
    "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

# gdml_vis.cc включает G4GDMLParser.hh -> G4GDMLRead.hh -> xercesc/... —
# заголовки Xerces-C++ не копируются в ${G4ROOT}/include/Geant4, нужен
# отдельный include path на инсталляцию Xerces (тот же env-var-паттерн,
# что и G4ROOT выше).
if(DEFINED ENV{XERCESC_ROOT} AND NOT DEFINED XERCESC_ROOT)
  set(XERCESC_ROOT "$ENV{XERCESC_ROOT}" CACHE PATH "Xerces-C++ install root (нужен для GDML-заголовков)")
endif()
if(NOT DEFINED XERCESC_ROOT)
  message(FATAL_ERROR
    "Не задан корень Xerces-C++ (нужен для GDML-заголовков). Укажите XERCESC_ROOT или -DXERCESC_ROOT=<путь>.")
endif()

file(GLOB XERCESC_IMPORT_LIBS "${XERCESC_ROOT}/lib/*.lib")

add_executable(gdml_vis gdml_vis.cc)
target_include_directories(gdml_vis PRIVATE "${G4ROOT}/include/Geant4" "${XERCESC_ROOT}/include")
target_link_libraries(gdml_vis PRIVATE ${G4_IMPORT_LIBS} ${XERCESC_IMPORT_LIBS})
target_compile_definitions(gdml_vis PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if(MSVC)
  target_compile_options(gdml_vis PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

# gdml_multiview — тот же рецепт, что gdml_vis выше (GDML-заголовки Xerces-C++
# нужны обоим): загрузка GDML из gdml_vis.cc + 4-видовая командная
# последовательность из detectors/RadiaCode-110/geant4/geometry/vis_multiview.cc,
# общий инструмент для любого детектора.
add_executable(gdml_multiview gdml_multiview.cc)
target_include_directories(gdml_multiview PRIVATE "${G4ROOT}/include/Geant4" "${XERCESC_ROOT}/include")
target_link_libraries(gdml_multiview PRIVATE ${G4_IMPORT_LIBS} ${XERCESC_IMPORT_LIBS})
target_compile_definitions(gdml_multiview PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if(MSVC)
  target_compile_options(gdml_multiview PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()
