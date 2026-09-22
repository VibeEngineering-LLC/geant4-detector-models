cmake_minimum_required(VERSION 3.16)
project(flight_msk_tyo CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Рецепт direct-link по import-библиотекам prebuilt Geant4 (как у gagg-intrinsic/app).
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR "Не задан корень Geant4: GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

# Цели добавляются построчно по мере появления приложений: (имя, исходник).
foreach(app src_check line_bench cabin_stage1 replay_selftest)
  add_executable(${app} ${app}.cc)
  target_include_directories(${app} PRIVATE "${G4ROOT}/include/Geant4" "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(${app} PRIVATE ${G4_IMPORT_LIBS})
  target_compile_definitions(${app} PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
  if(MSVC)
    target_compile_options(${app} PRIVATE /EHsc /MD /bigobj /utf-8 /wd4005 /wd4267 /wd4244)
  endif()
endforeach()

# Этап II: донорский спектрометр подключается по пути, файлы донора НЕ правятся и не копируются.
# DONOR_DIR — путь к geometry-каталогу репозитория geant4-detector-models/detectors/AtomSpectra-Nano-16-PRO
# на локальной машине; задать через -DDONOR_DIR=<путь> или переменную окружения DONOR_DIR.
if(DEFINED ENV{DONOR_DIR} AND NOT DEFINED DONOR_DIR)
  set(DONOR_DIR "$ENV{DONOR_DIR}" CACHE PATH "Каталог geometry донора AtomSpectra-Nano-16-PRO")
endif()
if(EXISTS "${DONOR_DIR}/ASN16Detector.cc" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cabin_stage2.cc")
  add_executable(cabin_stage2 cabin_stage2.cc "${DONOR_DIR}/ASN16Detector.cc")
  target_include_directories(cabin_stage2 PRIVATE "${G4ROOT}/include/Geant4" "${CMAKE_CURRENT_SOURCE_DIR}" "${DONOR_DIR}")
  target_link_libraries(cabin_stage2 PRIVATE ${G4_IMPORT_LIBS})
  target_compile_definitions(cabin_stage2 PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
  if(MSVC)
    target_compile_options(cabin_stage2 PRIVATE /EHsc /MD /bigobj /utf-8 /wd4005 /wd4267 /wd4244)
  endif()
endif()

# Рендер геометрии ASN16 для вебки — паттерн донора RadiaCode-110/geometry/vis_render.cc (см. шапку .cc).
if(EXISTS "${DONOR_DIR}/ASN16Detector.cc" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/asn16_vis_render.cc")
  add_executable(asn16_vis_render asn16_vis_render.cc "${DONOR_DIR}/ASN16Detector.cc")
  target_include_directories(asn16_vis_render PRIVATE "${G4ROOT}/include/Geant4" "${CMAKE_CURRENT_SOURCE_DIR}" "${DONOR_DIR}")
  target_link_libraries(asn16_vis_render PRIVATE ${G4_IMPORT_LIBS})
  target_compile_definitions(asn16_vis_render PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
  if(MSVC)
    target_compile_options(asn16_vis_render PRIVATE /EHsc /MD /bigobj /utf-8 /wd4005 /wd4267 /wd4244)
  endif()
endif()
