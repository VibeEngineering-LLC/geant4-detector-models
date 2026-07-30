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

# Штамп провенанса: отпечаток исходников запекается в бинарник и попадает в
# шапку КАЖДОГО выходного спектра (см. provenance.cmake — там же зачем).
# Отдельная цель, а не PRE_BUILD: PRE_BUILD честно работает только у генератора
# Visual Studio, у Ninja и Makefile он молча превращается в PRE_LINK.
set(G1S_PROV_SRC "${CMAKE_CURRENT_SOURCE_DIR}/main.cc"
                 "${CMAKE_CURRENT_SOURCE_DIR}/G1SDetector.cc"
                 "${CMAKE_CURRENT_SOURCE_DIR}/G1SDetector.hh")
set(G1S_PROV_HH "${CMAKE_BINARY_DIR}/g1s_provenance.hh")
add_custom_target(g1s_provenance ALL
  COMMAND "${CMAKE_COMMAND}"
          "-DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
          "-DOUT=${G1S_PROV_HH}"
          "-DPREFIX=G1S"
          "-DSRC_LIST=${G1S_PROV_SRC}"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/../../../common/cmake/provenance.cmake"
  BYPRODUCTS "${G1S_PROV_HH}"
  COMMENT "Отпечаток исходников -> g1s_provenance.hh")

add_executable(g1s main.cc G1SDetector.cc)
add_dependencies(g1s g1s_provenance)
target_include_directories(g1s PRIVATE "${CMAKE_BINARY_DIR}")
target_include_directories(g1s PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(g1s PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(g1s PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if (MSVC)
  target_compile_options(g1s PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

# Коэффициент ослабления матрицы ОИСН-16 — для подгонки эффективной толщины
# mucalc тоже компилирует G1SDetector.cc, то есть его выход зависит от геометрии
# и материалов, а читается анализом (compare_xcom.py, selfabs_fit.py). Отсутствие
# у него провенанса было дефектом, не намеренным исключением: два набора μ/ρ от
# разных версий материалов выглядели бы одинаково прослеженными.
add_executable(mucalc mucalc.cc G1SDetector.cc)
add_dependencies(mucalc g1s_provenance)
target_include_directories(mucalc PRIVATE "${CMAKE_BINARY_DIR}")
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
