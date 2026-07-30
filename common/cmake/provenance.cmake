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
# ОБЩИЙ НА ВСЕ ДЕТЕКТОРЫ. Лежал в geometry/ Гамма-1С, из-за чего у RadiaCode-103
# провенанса не было вовсе — при том что RCDetector.cc правится активно, а
# правила method-rules объявлены обязательными для ЛЮБОГО детектора. Копировать
# скрипт во второй каталог значило бы завести второе правило в двух местах,
# чем этот репозиторий уже наказан четырежды.
#
# Имя макроса-отпечатка задаётся параметром PREFIX, чтобы заголовки разных
# приборов не конфликтовали при общей сборке.
#
# Вызов (из CMakeLists):
#   cmake -DSRC_DIR=... -DOUT=... -DPREFIX=G1S -DSRC_LIST="a.cc;b.cc"
#         -P <репозиторий>/common/cmake/provenance.cmake

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
#
# `describe --dirty` НЕ видит неотслеживаемых файлов (найдено независимым
# аудитом): новый .cc, ещё не добавленный в git, даёт чистое "-dirty"-показание
# на коммите, где производящего файла не существовало. Поэтому отдельно
# проверяем `git status --porcelain` по каждому файлу SRC_LIST и, если
# какой-то из них НЕ отслеживается (статус "??") или изменён, помечаем это
# явно суффиксом — не полагаемся на то, что --dirty это заметит сам.
set(_git "нет-git")
find_program(_GIT git)
if(_GIT)
  execute_process(COMMAND "${_GIT}" -C "${SRC_DIR}" describe --always --dirty
                  OUTPUT_VARIABLE _git OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET RESULT_VARIABLE _rc)
  if(NOT _rc EQUAL 0 OR _git STREQUAL "")
    set(_git "нет-git")
  endif()
  if(NOT _git STREQUAL "нет-git")
    execute_process(COMMAND "${_GIT}" -C "${SRC_DIR}" status --porcelain -- ${SRC_LIST}
                    OUTPUT_VARIABLE _status OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET RESULT_VARIABLE _rc2)
    if(_rc2 EQUAL 0 AND NOT _status STREQUAL "" AND NOT _git MATCHES "-untracked$")
      string(APPEND _git "-untracked")
    endif()
  endif()
endif()

if(NOT DEFINED PREFIX)
  set(PREFIX "G1S")
endif()
file(WRITE "${OUT}.tmp"
"// Сгенерировано common/cmake/provenance.cmake перед сборкой. Не править.\n"
"#pragma once\n"
"#define ${PREFIX}_SRC_SHA1 \"${_short}\"\n"
"#define ${PREFIX}_GIT_DESCRIBE \"${_git}\"\n"
"// Разбор отпечатка (какой файл что дал):\n"
"/*\n${_acc}*/\n")
execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${OUT}.tmp" "${OUT}")
file(REMOVE "${OUT}.tmp")
message(STATUS "provenance: src_sha1=${_short} git=${_git}")
