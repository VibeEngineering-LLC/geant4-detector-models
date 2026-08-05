cmake_minimum_required(VERSION 3.16)
project(asn16 CXX)

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
# шапку КАЖДОГО выходного спектра (см. common/cmake/provenance.cmake).
set(ASN16_PROV_SRC "${CMAKE_CURRENT_SOURCE_DIR}/main.cc"
                   "${CMAKE_CURRENT_SOURCE_DIR}/ASN16Detector.cc"
                   "${CMAKE_CURRENT_SOURCE_DIR}/ASN16Detector.hh")
set(ASN16_PROV_HH "${CMAKE_BINARY_DIR}/asn16_provenance.hh")
add_custom_target(asn16_provenance ALL
  COMMAND "${CMAKE_COMMAND}"
          "-DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
          "-DOUT=${ASN16_PROV_HH}"
          "-DPREFIX=ASN16"
          "-DSRC_LIST=${ASN16_PROV_SRC}"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/../../../common/cmake/provenance.cmake"
  BYPRODUCTS "${ASN16_PROV_HH}"
  COMMENT "Отпечаток исходников -> asn16_provenance.hh")

add_executable(asn16 main.cc ASN16Detector.cc)
add_dependencies(asn16 asn16_provenance)
target_include_directories(asn16 PRIVATE "${CMAKE_BINARY_DIR}")
target_include_directories(asn16 PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(asn16 PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(asn16 PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if (MSVC)
  target_compile_options(asn16 PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

file(GLOB ASN16_MACROS "${CMAKE_SOURCE_DIR}/*.mac")
foreach(m ${ASN16_MACROS})
  get_filename_component(n ${m} NAME)
  configure_file(${m} ${CMAKE_BINARY_DIR}/${n} COPYONLY)
endforeach()
