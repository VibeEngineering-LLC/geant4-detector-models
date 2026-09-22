# SESSION-STATE — проект «АтомНано 16 на эшелоне» (`flight-msk-tyo`)

**Написано:** 2026-09-22 ~16:50. Срок годности 7 дней.

## Что это
DAD–SGN, 20.09.2026, Airbus A320/A321, 9800 м, АтомНано 16 (CsI 18×15×60 мм), 20 кэВ – 10 МэВ. Geant4 11.4.2 +
PARMA/EXPACS 4.10. Двухэтапная схема: этап I — поле в салоне, этап II — воспроизведение на модели прибора
(K-кратный розыгрыш положения, окно ℓ=120 см).

## Текущее состояние
- Итог расчёта: 71,19 ± 0,20 отсч./с. Опубликовано как Artifact, **версия 29**:
  `https://claude.ai/artifact/QVrFt3hnFW17n8nCf7nuEw` — RU-страница + вкладка EN (`en/spectrum-dad-sgn-en.html`
  как отдельный published-файл через параметр `files` при публикации).
- **Английская версия готова:** `web/flight-dad-sgn/text.en.md` (полный перевод), `data.en.js`,
  `research/B-gamma-lines-airframe.en.md` (206 строк раздела 9, словарный перевод — `scripts/translate_dict.py`
  + `scripts/ru_en_reactions.json`/`ru_en_terms.json`), `research/G-neutron-lines-summary.en.md`,
  `results/final/{bands,field,dose,fuel}.en.md`. `scripts/build_page.py` параметризован флагом `--lang en`
  (второй позиционный аргумент CLI не нужен — флаг `--lang en` в конце команды). Видимая кнопка EN/RU в углу
  шапки (`LANGBTN`, `@@LANGHREF@@` placeholder, разный href для index.html/glued-файла).
- **Демо на GitHub Pages (RU):** `https://vibeengineering-llc.github.io/demo-web-pages/atomnano16-flight-dad-sgn/`
  (публичный репо VibeEngineering-LLC/demo-web-pages, метод §29). EN-версия туда ещё не выложена.
- **Исходники расчёта на GitHub:** `github.com/Verter73/GEANT4` (**приватный**, только что создан). Внутри —
  подпроект `flight-msk-tyo/` (см. `.gitignore` в корне GEANT4 — исключены `results/stage1`, `results/stage2`
  (19 ГБ raw), `build/`, `logs/`, `tables/`, `outputs/`, скачанный сторонний xlsx, большой файл измерений).
  Локальная история переписана в ОДИН коммит (`git init` заново, с явного «да» оператора) — личные пути и
  внутренние служебные имена вычищены точечной заменой перед коммитом (одноразовые скрипты удалены после
  использования); `src/CMakeLists.txt` DONOR_DIR и `scripts/build_flight.ps1` переведены на переменную
  окружения/относительный путь (без хардкода личного пути). **Push шёл в фоне на момент записи — статус не
  подтверждён, проверить `gh api repos/Verter73/GEANT4/commits` при возобновлении.**

## Пересборка (RU и EN)
```
PYTHONIOENCODING=utf-8 python scripts/build_page.py web/flight-dad-sgn results/final/bands.md
PYTHONIOENCODING=utf-8 python scripts/build_page.py web/flight-dad-sgn results/final/bands.md --lang en
PYTHONIOENCODING=utf-8 python scripts/build_report.py web/flight-dad-sgn/text.md results/final REPORT-2026-09-21-dad-sgn.md
```
Предпросмотр: `preview_start flight-page`, порт 8765, кэш браузера — добавлять `?v=N`.

## Открытые вопросы / не сделано
- **Push в Verter73/GEANT4 не подтверждён завершённым** — проверить первым делом.
- EN-версия НЕ выложена на GitHub Pages (только RU). Если понадобится — та же процедура §29,
  папка `atomnano16-flight-dad-sgn/en/`.
- `audit_run.py` не перезапускался после серии правок 22.09 (числа могли устареть в `audit/claims-commands.txt`).
- Средние находки внешнего код-аудита не закрыты: 2.1 (радиус ре-инжекции 139/140 мм — нужен пересчёт, #CFG-1),
  2.4–2.6 (методические вопросы: окно кластеризации, слой фантома 30 мм, wR=2 для пиона) — ждут решения
  оператора, не автоправка.
- Перевод 105 reaction-строк + 97 material-строк сделан СЛОВАРНЫМ методом (не вычитан построчно человеком) —
  Ollama-делегация отказала (VRAM guard, W-006). Спот-чек сделан (137.85/512.65/57.61/58.11 кэВ), не полная
  вычитка. Если найдётся неточность — править `scripts/ru_en_reactions.json`/`ru_en_terms.json` точечно.
- TODO (оператор 22.09): разложить непойманные пики (напр. 96,5 кэВ) по изотопу-эмиттеру через per-event
  трекинг этапа II — не начато, отдельный тяжёлый прогон.

## Конвенции
Код ≥25 строк — Edit небольшими кусками (delegation_guard блокирует Write/Edit длиннее 24 строк кода).
Не-код (.md/.json) — Write/Edit, не Bash-heredoc. Долгие Bash — `run_in_background: true`. Публикация Artifact —
republish тем же путём (URL сохраняется), только после проверки в браузере. Push в приватный/публичный git —
только по явному «да» оператора на конкретную цель; секрет-скан (`push-guard` хук) — обязательный автоматический
гейт, скраб делать ДО коммита, не пытаться обойти. Переписывать git-историю — только с явного «да», и только
если история ещё НЕ опубликована (ни один push не прошёл успешно).
