cmake_minimum_required(VERSION 3.16)
project(rc110_vis_render CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Тот же рецепт прямой линковки, что у RadiaCode-103/geometry/CMakeLists.txt
# и common/tools/CMakeLists.txt — direct-link списком import-lib, не
# find_package(Geant4) (см. комментарий-обоснование там же: не зависит от
# того, полный архив у собирающего или урезанный).
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR
    "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

add_executable(rc110_vis_render RC110Detector.cc vis_render.cc)
target_include_directories(rc110_vis_render PRIVATE "${CMAKE_BINARY_DIR}")
target_include_directories(rc110_vis_render PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(rc110_vis_render PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(rc110_vis_render PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if(MSVC)
  target_compile_options(rc110_vis_render PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()

# rc110_vis_multiview — 4 вида (спереди/сбоку/разрез/три четверти) одним
# прогоном, см. шапку vis_multiview.cc. Та же геометрия RC110Detector.cc,
# второй executable в этом же CMakeLists (проще интегрировать, чем отдельный
# подкаталог — донор паттерна: сам этот файл, второй add_executable).
add_executable(rc110_vis_multiview RC110Detector.cc vis_multiview.cc)
target_include_directories(rc110_vis_multiview PRIVATE "${CMAKE_BINARY_DIR}")
target_include_directories(rc110_vis_multiview PRIVATE "${G4ROOT}/include/Geant4")
target_link_libraries(rc110_vis_multiview PRIVATE ${G4_IMPORT_LIBS})
target_compile_definitions(rc110_vis_multiview PRIVATE G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if(MSVC)
  target_compile_options(rc110_vis_multiview PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()