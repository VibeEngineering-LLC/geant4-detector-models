# Внешние источники

Литература и документация, на которую опираются решения по геометрии, физике и
методике. Ссылка здесь означает, что источником пользовались при построении
модели или при выборе конкретного числа, — а не что он просто относится к теме.

Методика ЛСРМ (обработка спектров, мёртвое время, собственная активность)
лежит отдельно, в `detectors/Gamma-1S/reference/lsrm/references/`: она
привязана к конкретной поверочной цепочке, а не к моделированию вообще.

## Geant4

- **Application Developers Guide** —
  https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/index.html
  Основное руководство: конструирование геометрии и материалов, списки физики,
  первичный розыгрыш (GPS), пользовательские действия, макрокоманды UI. Версия
  `dev`; при расхождении с поведением сборки 11.2.1 верить сборке.

- **Исходный код** — https://github.com/Geant4/geant4
  Последняя инстанция, когда документация неоднозначна или молчит. Что реально
  проверяется по коду, а не по руководству: состав материалов NIST
  (`source/materials/src/G4NistMaterialBuilder.cc` — оттуда `G4_MAGNESIUM_OXIDE`
  и прочие), набор процессов в списках физики (`G4EmStandardPhysics_option4`),
  поведение команд GPS при сочетаниях `/gps/ang/*`. Ветка тега сборки — здесь
  `v11.2.1`, не `master`: расхождение версий и есть обычная причина, по которой
  «в документации написано иначе».

## Плотность отражателя MgO

Отражатель — насыпной прессованный порошок, а не сплошной оксид, и это
решает: его массовая толщина управляет мягким краем кривой эффективности.
Сводка обсуждения — в комментарии `detectors/Gamma-1S/geometry/G1SDetector.hh`
у параметра `mgoDensity`.

- **Mendes B.M. et al.** Monitoring internal contamination from Occupationally
  Exposed Workers of an ¹⁸F-FDG production plant. *Braz. J. Rad. Sci.* 07-03A
  (2019) 01-12. Таблица 1 — состав и плотности узлов MCNP-модели NaI(Tl):
  ρ(MgO reflector) = 2,0 г/см³ (по Mouhti 2017 и Salgado 2012).
- **Appl. Radiat. Isot.** — Detection efficiency evaluation for low energy of a
  NaI(Tl) scintillation detector.
  https://www.sciencedirect.com/science/article/abs/pii/S0969806X22003681
- **Appl. Radiat. Isot.** — A computational modelling of low-energy gamma ray
  detection efficiency of a cylindrical NaI(Tl) detector.
  https://www.sciencedirect.com/science/article/abs/pii/S0969806X21002310
- **Appl. Radiat. Isot.** — Optimization of the Monte Carlo simulation model of
  NaI(Tl) detector by Geant4 code.
  https://www.sciencedirect.com/science/article/abs/pii/S0969804317307479

  Три работы выше описывают в том числе приём КАЛИБРОВКИ плотности отражателя
  по отношению эффективностей двух мягких линий. Приём известен, но в этом
  репозитории **намеренно не применяется**: подогнанная так плотность вбирает
  в себя всё прочее неверное в торце (толщину Al, воздушные зазоры, то, чего
  нет в чертеже) и перестаёт быть плотностью, оставаясь эффективной массовой
  толщиной под чужим именем. После подгонки согласие на тех же линиях
  перестаёт быть проверкой. Берётся физическое значение, результат
  проверяется.

- **Bell S. et al.** MgO reflectance data for Monte Carlo simulation of
  LaBr₃:Ce scintillation crystals. *NIM A* 701 (2013) 44–53.
  https://www.sciencedirect.com/science/article/abs/pii/S0168900212012065
  Про ОПТИЧЕСКОЕ отражение как функцию толщины слоя; плотности там нет.
  Полезна одним: отражение с толщиной выходит на насыщение — значит 6 мм на
  входном торце стоят не ради светосбора.
