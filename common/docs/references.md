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

## Коэффициенты ослабления фотонов

Эталон, против которого проверяются сечения самого тулкита: без этого разбор
мягкого края молча опирается на непроверенное допущение, что ослабление в
слоях торца Geant4 считает верно. Сверка — `analysis/compare_xcom.py`,
результат — `results/mu_geant4_vs_xcom.csv`.

- **NIST XCOM** — https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
  Расчёт μ/ρ по компонентам (когерентное, некогерентное, фотоэффект, рождение
  пар) на произвольной энергии. Даёт две суммы — с когерентным рассеянием и
  без; путать их нельзя, на 59,5 кэВ они расходятся на 12 %. Основной расчёт
  ослабления в `mucalc.cc` идёт БЕЗ когерентного, ему соответствует колонка
  «without coherent».
- **Hubbell J.H., Seltzer S.M.** Tables of X-Ray Mass Attenuation Coefficients
  and Mass Energy-Absorption Coefficients. NISTIR 5632.
  https://physics.nist.gov/PhysRefData/XrayMassCoef/cover.html
  Таблицы на фиксированной сетке: Table 3 — элементы, Table 4 — составы.
  Величина в них — μ/ρ **с** когерентным рассеянием (раздел 2 определяет
  полное сечение суммой, куда σ_coh входит явно; проверено на Al: узел
  2,778·10⁻¹ при 60 кэВ сходится именно с этой колонкой). MgO, NaI и резины
  в Table 4 нет — для них первоисточник только XCOM.
- **NIST STAR, состав материалов** —
  https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=243
  Натуральный каучук (ICRU-37): H 0,118371 / C 0,881629, ρ = 0,92 г/см³.
  Готовой таблицы μ/ρ для резины у NIST нет, поэтому она считается в XCOM
  режимом смеси по этому составу.

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
