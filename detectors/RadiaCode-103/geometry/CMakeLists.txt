cmake_minimum_required(VERSION 3.16)
project(rc_curves CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# find_package(Geant4) недоступен: в prebuilt-сборке CERN нет Geant4Config.cmake,
# линкуем напрямую по сгенерированным import-библиотекам (см. скилл geant4-win-native).
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR
    "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

# Штамп провенанса — тем же общим генератором, что у Гамма-1С
# (common/cmake/provenance.cmake). Прежде его здесь не было вовсе, при том что
# RCDetector.cc правится активно: два набора кривых от разных версий геометрии
# выглядели одинаково прослеженными. Отдельная цель, а не PRE_BUILD: у Ninja и
# Makefile тот молча вырождается в PRE_LINK.
set(RC_PROV_SRC "${CMAKE_CURRENT_SOURCE_DIR}/main.cc"
                "${CMAKE_CURRENT_SOURCE_DIR}/RCDetector.cc"
                "${CMAKE_CURRENT_SOURCE_DIR}/RCDetector.hh")
set(RC_PROV_HH "${CMAKE_BINARY_DIR}/rc_provenance.hh")
add_custom_target(rc_provenance ALL
  COMMAND "${CMAKE_COMMAND}"
          "-DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
          "-DOUT=${RC_PROV_HH}"
          "-DPREFIX=RC"
          "-DSRC_LIST=${RC_PROV_SRC}"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/../../../common/cmake/provenance.cmake"
  BYPRODUCTS "${RC_PROV_HH}"
  COMMENT "Отпечаток исходников -> rc_provenance.hh")

add_executable(rc_curves main.cc RCDetector.cc)
add_dependencies(rc_curves rc_provenance)
target_include_directories(rc_curves PRIVATE "${CMAKE_BINARY_DIR}")
target_include_directories(rc_curves PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(rc_curves PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(rc_curves PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if (MSVC)
  target_compile_options(rc_curves PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

# Коэффициенты ослабления матриц пробы — для приведения поправки на
# самопоглощение к физической форме с одним параметром
add_executable(mucalc mucalc.cc RCDetector.cc)
add_dependencies(mucalc rc_provenance)
target_include_directories(mucalc PRIVATE "${CMAKE_BINARY_DIR}")
target_include_directories(mucalc PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(mucalc PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(mucalc PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if (MSVC)
  target_compile_options(mucalc PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

# Спектр поля ЕРН в помещении — отдельная задача со своей геометрией
add_executable(wallfield wallfield.cc)
target_include_directories(wallfield PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(wallfield PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(wallfield PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if (MSVC)
  target_compile_options(wallfield PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

file(GLOB RC_MACROS "${CMAKE_SOURCE_DIR}/*.mac")
foreach(m ${RC_MACROS})
  get_filename_component(n ${m} NAME)
  configure_file(${m} ${CMAKE_BINARY_DIR}/${n} COPYONLY)
endforeach()
