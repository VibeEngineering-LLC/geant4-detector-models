cmake_minimum_required(VERSION 3.16)
project(gamma1s CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Корень Geant4 — из переменной окружения GEANT4_ROOT или -DG4ROOT=...
# Абсолютных путей, привязанных к машине, в репозитории нет.
#
# ДВА СПОСОБА СБОРКИ:
#  1) Обычная установка Geant4 с Geant4Config.cmake — работает find_package.
#  2) Prebuilt-сборка CERN под Windows: Geant4Config.cmake в ней НЕТ, поэтому
#     линкуем напрямую по import-библиотекам и обязательно определяем
#     G4LIB_BUILD_DLL (переключает заголовки на __declspec(dllimport)).
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR
    "Не задан корень Geant4. Укажите переменную окружения GEANT4_ROOT "
    "или передайте -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

add_executable(g1s main.cc G1SDetector.cc)
target_include_directories(g1s PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(g1s PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(g1s PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if (MSVC)
  target_compile_options(g1s PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

# Коэффициент ослабления матрицы ОИСН-16 — для подгонки эффективной толщины
add_executable(mucalc mucalc.cc G1SDetector.cc)
target_include_directories(mucalc PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(mucalc PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(mucalc PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if (MSVC)
  target_compile_options(mucalc PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

file(GLOB G1S_MACROS "${CMAKE_SOURCE_DIR}/*.mac")
foreach(m ${G1S_MACROS})
  get_filename_component(n ${m} NAME)
  configure_file(${m} ${CMAKE_BINARY_DIR}/${n} COPYONLY)
endforeach()
