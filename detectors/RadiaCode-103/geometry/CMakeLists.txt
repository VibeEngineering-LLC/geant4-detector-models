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

add_executable(rc_curves main.cc RCDetector.cc)
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
