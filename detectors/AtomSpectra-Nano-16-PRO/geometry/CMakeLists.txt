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
file(GLOB G4_IMPORT_LIBS CONFIGURE_DEPENDS "${G4ROOT}/lib/*.lib")

# Штамп провенанса: отпечаток исходников запекается в бинарник и попадает в
# шапку КАЖДОГО выходного спектра (см. common/cmake/provenance.cmake).
# В отпечаток входят И МАКРОСЫ. Результат зависит не только от исходников:
# положение источника, угол конуса, сетка энергий и число первичных заданы в
# .mac. Пока их не было в списке, правка макроса не меняла src_sha1, а в шапке
# спектра лежало только ИМЯ макроса, не его содержимое — два прогона с разной
# постановкой источника выглядели одинаково прослеженными. Найдено независимым
# аудитом 05.08.2026.
# CONFIGURE_DEPENDS обязателен: без него GLOB не перевычисляется при
# инкрементальной сборке, ДОБАВЛЕННЫЙ макрос в отпечаток не попадает, и
# спектры новой постановки получают ПРЕЖНИЙ src_sha1. Удаление и
# переименование ловились и раньше (provenance.cmake), добавление — нет.
file(GLOB ASN16_PROV_MAC CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/../macros/*.mac")
set(ASN16_PROV_SRC "${CMAKE_CURRENT_SOURCE_DIR}/main.cc"
                   "${CMAKE_CURRENT_SOURCE_DIR}/ASN16Detector.cc"
                   "${CMAKE_CURRENT_SOURCE_DIR}/ASN16Detector.hh"
                   "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
                   ${ASN16_PROV_MAC})
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

# Макросы лежат в СОСЕДНЕМ каталоге ../macros, а CMAKE_SOURCE_DIR указывает на
# geometry/ — прежний глоб по CMAKE_SOURCE_DIR был пуст, и шаг копирования
# молча не делал ничего: рецепт из README с чистого клона не работал. Найдено
# независимым аудитом 05.08.2026.
file(GLOB ASN16_MACROS CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/../macros/*.mac")
if(NOT ASN16_MACROS)
  message(WARNING "Макросы не найдены в ../macros — asn16 нечем запускать.")
endif()
foreach(m ${ASN16_MACROS})
  get_filename_component(n ${m} NAME)
  configure_file(${m} ${CMAKE_BINARY_DIR}/${n} COPYONLY)
endforeach()
