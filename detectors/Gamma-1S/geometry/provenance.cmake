# Отпечаток исходников геометрии, запекаемый В БИНАРНИК.
#
# ЗАЧЕМ. Трижды за линию Гамма-1С вывод строился на спектрах, посчитанных
# ПРЕДЫДУЩЕЙ геометрией: устаревшие моно-спектры моста, устаревший exe в
# каталоге прогонов, точечные сетки от 28.07 против правок торца 29.07.
# Сверка по mtime от этого не спасает — `git checkout`, клон и синхронизация
# облака сбрасывают времена в произвольную сторону. Отвечать надо на вопрос
# «ИЗ ЧЕГО получено», а не «когда записано», а на него может ответить только
# сам бинарник, считавший спектр.
#
# КАК. Скрипт считает SHA1 по содержимому исходников и пишет
# g1s_provenance.hh. Он запускается ПЕРЕД каждой сборкой (add_custom_target +
# add_dependencies), поэтому отпечаток не может отстать от exe: содержимое
# заголовка меняется -> main.cc перекомпилируется -> exe несёт актуальный
# отпечаток -> каждый выходной спектр несёт его в шапке.
#
# copy_if_different обязателен: без него заголовок переписывался бы на каждой
# сборке и тянул за собой перекомпиляцию main.cc даже когда ничего не менялось.
#
# Вызов (из CMakeLists):
#   cmake -DSRC_DIR=... -DOUT=... -DSRC_LIST="a.cc;b.cc" -P provenance.cmake

set(_acc "")
foreach(f ${SRC_LIST})
  if(NOT EXISTS "${f}")
    message(FATAL_ERROR "provenance: нет исходника ${f}")
  endif()
  file(SHA1 "${f}" _h)
  get_filename_component(_n "${f}" NAME)
  string(APPEND _acc "${_n}:${_h}\n")
endforeach()
string(SHA1 _sha "${_acc}")
string(SUBSTRING "${_sha}" 0 12 _short)

# Коммит — справочно: он говорит, что ЗАКОММИЧЕНО, а не что СОБРАНО. Правда о
# собранном — в _short выше; -dirty здесь ровно затем, чтобы несовпадение этих
# двух ответов было видно.
set(_git "нет-git")
find_program(_GIT git)
if(_GIT)
  execute_process(COMMAND "${_GIT}" -C "${SRC_DIR}" describe --always --dirty
                  OUTPUT_VARIABLE _git OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET RESULT_VARIABLE _rc)
  if(NOT _rc EQUAL 0 OR _git STREQUAL "")
    set(_git "нет-git")
  endif()
endif()

file(WRITE "${OUT}.tmp"
"// Сгенерировано provenance.cmake перед сборкой. Руками не править.\n"
"#pragma once\n"
"#define G1S_SRC_SHA1 \"${_short}\"\n"
"#define G1S_GIT_DESCRIBE \"${_git}\"\n"
"// Разбор отпечатка (какой файл что дал):\n"
"/*\n${_acc}*/\n")
execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${OUT}.tmp" "${OUT}")
file(REMOVE "${OUT}.tmp")
message(STATUS "provenance: src_sha1=${_short} git=${_git}")
