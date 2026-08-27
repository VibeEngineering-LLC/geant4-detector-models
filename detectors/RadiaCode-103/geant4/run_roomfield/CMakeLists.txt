cmake_minimum_required(VERSION 3.16)
project(rc103_roomfield CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Тот же рецепт direct-link, что run_field/CMakeLists.txt (эталон) — не
# find_package(Geant4). Соседние комплекты не тронуты.
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR
    "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")

# GDML здесь не используется, но набор импорт-библиотек Geant4 берётся целиком
# тем же глобом, что и в run_field; каталог заголовков Xerces-C++ подключается
# ровно по тому же проверенному рецепту, чтобы не разъезжаться с эталоном.
if(DEFINED ENV{XERCESC_ROOT} AND NOT DEFINED XERCESC_ROOT)
  set(XERCESC_ROOT "$ENV{XERCESC_ROOT}" CACHE PATH "Xerces-C++ install root")
endif()
if(NOT DEFINED XERCESC_ROOT)
  message(FATAL_ERROR
    "Не задан корень Xerces-C++. Укажите XERCESC_ROOT или -DXERCESC_ROOT=<путь>.")
endif()
file(GLOB XERCESC_IMPORT_LIBS "${XERCESC_ROOT}/lib/*.lib")

add_executable(rc103_roomfield
  main.cc
  Rc103RoomFieldGeometry.cc
  Rc103RoomFieldDetectorConstruction.cc
  Rc103RoomFieldPhysicsList.cc
  Rc103RoomFieldPrimaryGeneratorAction.cc
  Rc103RoomFieldRunAction.cc
  Rc103RoomFieldSteppingAction.cc
)
target_include_directories(rc103_roomfield PRIVATE
  "${G4ROOT}/include/Geant4" "${XERCESC_ROOT}/include")
target_link_libraries(rc103_roomfield PRIVATE
  ${G4_IMPORT_LIBS} ${XERCESC_IMPORT_LIBS})
target_compile_definitions(rc103_roomfield PRIVATE
  G4LIB_BUILD_DLL WIN32 NOMINMAX _USE_MATH_DEFINES)
if(MSVC)
  target_compile_options(rc103_roomfield PRIVATE /EHsc /MD /bigobj /utf-8
                         /wd4005 /wd4267 /wd4244)
endif()
