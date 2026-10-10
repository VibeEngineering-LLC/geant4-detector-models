# Словарь терминов RU → EN для английской страницы GS2020 (#GS-79)

Составлен 10.10.2026 по указанию оператора: «прямые переводы терминов не допускаются, нужны точные термины на
каждом языке». Применяется: (1) в промпте локальной модели при переводе HTML-прозы, строк JS и данных
(`scripts/translate_gs2020_en.py` вставляет таблицу в промпт целиком); (2) стерильным проходом — для сверки
перевода с этим словарём.

**Источник** — употребление термина в англоязычных учебниках и руководствах из библиотеки оператора
(`skills/gamma-library/_incoming/2026-08-27_telegram-export/Export/`). Число в скобках — число вхождений строки
(`grep -oi`, без учёта переносов строк, замер 10.10.2026): **Kn** — Knoll, *Radiation Detection and Measurement*,
4th ed. (`Knoll4thEdition.md`); **Gi** — Gilmore, *Practical Gamma-ray Spectroscopy*, 2nd ed., 2008
(`pgs-gilmore-2008.md`); **FR** — FRMAC Gamma Spectroscopy Knowledge Guide, 2019 (`FRMAC_GammaSpec_KnowledgeGuide_2019-08_UUR.md`);
**AQ48** — IAEA Analytical Quality in Nuclear Applications No. 48 *Determination and Interpretation of Characteristic
Limits for Radioactivity Measurements* (`AQ-48_web.md`); **G4** — имена классов и параметров Geant4 (код физ-листа
контура). Где в библиотеке термина нет — так и сказано, источник назван отдельно.

## Термины

| Русский | English | Источник / примечание |
|---|---|---|
| пик полного поглощения | full-energy peak (FEP) | Kn 89, Gi 36, FR 75; «full-absorption peak» — 0 во всех трёх, не использовать; «photopeak» допустим, но не основной |
| функция полного поглощения (ПП); «функция ПП» в названии метода 2 | full-energy-peak (FEP) efficiency function; в названии метода — «FEP efficiency function», единообразно везде | FR «FEP» 8; Gi «full-energy peak efficiency» 13; принятая терминология #GS-78 (координатор 10.10) |
| ПШПВ | FWHM | Kn 112, Gi 214 |
| разрешение (энергетическое) | energy resolution | Kn 386, Gi 13 |
| форма пика | peak shape | Gi 34 |
| ширина пика | peak width | Gi 93 |
| центроида (пика) | peak centroid | Kn 4, Gi 10 |
| край Комптона | Compton edge | Kn 12, Gi 16, FR 16 |
| пик обратного рассеяния | backscatter peak | Kn 10, FR 17 |
| пик вылета (иода) | (iodine) X-ray escape peak | «escape peak» Kn 67, Gi 40; «X-ray escape» Kn 16; «iodine escape» Gi 2 |
| аннигиляционный пик | annihilation peak | Gi 23 |
| сумм-пик; каскадное суммирование | sum peak; true coincidence summing | «sum peak» Kn 24, Gi 18; «true coincidence summing» Gi 59 |
| поправка на суммирование | summing correction | Gi, AQ48, FR (есть во всех) |
| эффективность регистрации | detection efficiency | Kn 99, Gi 9 |
| эффективность в пике полного поглощения | full-energy peak efficiency | Kn 14, Gi 13 |
| самопоглощение | self-absorption | Kn 32, Gi 31; FR — «self-attenuation» 14 (синоним) |
| ослабление | attenuation | Kn 113, Gi 75, FR 73 |
| площадь пика (нетто) | net peak area | Gi 16 |
| фон; вычитание фона | background; background subtraction | Kn 2, Gi 1, FR 1 |
| континуум | continuum | Kn 109, Gi 120 |
| отклик (детектора); функция отклика | (detector) response; response function | «response function» Kn 121; «detector response» Kn 30 |
| матрица отклика | response matrix | Kn, AQ48 |
| шаблон (МК) | (Monte Carlo) response template | термина в библиотеке нет как устойчивого; «template» в Kn и Gi есть; принят по указанию координатора 10.10 |
| сетка откликов (моноэнергетическая) | grid of monoenergetic responses | «monoenergetic» Kn 73 |
| монолиния | monoenergetic line | «monoenergetic» Kn 73 |
| свёртка; уширение (размытие) | convolution; (Gaussian) broadening | «convolved» Kn; «Gaussian broadening» Gi 1; принятая терминология #GS-78 |
| ядро свёртки (ядро размытия) | convolution kernel | принятая терминология #GS-78 (координатор 10.10; заменяет «smearing kernel» брифа) |
| набор данных детектора; карточка детектора | detector data set | бриф #GS-79 п.5; #GS-78 |
| референсный спектр тория | thorium reference spectrum | #GS-78 |
| измерение (замер) | measurement | #GS-78 |
| собственная шкала энергии (своя шкала) | own energy calibration (of this spectrum) | #GS-78 |
| физический список | physics list | G4 (`G4VModularPhysicsList`); #GS-78 |
| мягкая зона (ниже ~150 кэВ) | soft region — при первом упоминании с пояснением «(below about 150 keV)», как в русском | #GS-78 |
| систематическая погрешность | systematic uncertainty | Kn, Gi «uncertainty»; #GS-78 |
| пики вылета | escape peaks | Kn 67, Gi 40; #GS-78 |
| собственная ширина (замер в этом спектре) | own width (measured in this spectrum) — НЕ «intrinsic» | «intrinsic» в гамма-спектрометрии — собственное разрешение кристалла; здесь другой смысл |
| приведённый фон | background scaled to live time (scaled background) | смысл: умножен на отношение живых времён; не «reduced» |
| к известной / ожидаемой активности (отношение) | ratio to the known / expected activity | — |
| опорный спектр | reference spectrum | #GS-78, второй раунд (координатор 10.10) |
| опорные (реперные) линии | reference lines | #GS-78 |
| перенос (частиц, излучения) | transport | #GS-78; G4 «transportation» |
| сборка кристалла (пакет кристалла) | crystal assembly | #GS-78 |
| девозбуждение атома (деэкситация) | atomic de-excitation | #GS-78; G4 `G4UAtomicDeexcitation` |
| порог образования вторичных частиц (порог продукции) | production cut | #GS-78; G4 `G4ProductionCuts` |
| проба (неизвестная активность) / образец (известная, эталонный) | sample / reference source (reference sample) | #GS-78: по назначению |
| пик суммирования; истинное каскадное суммирование | sum peak (summing peak); true coincidence summing | Kn «sum peak» 24; Gi 59; #GS-78 |
| поправка r(E) | r(E) correction | #GS-78 |
| (библиотека, метод) Гамма-1С (прежнее «донорский») | Gamma-1S (library, method) — имя собственное; «donor» не употреблять | #GS-78 второй раунд |
| мера A2 (критерий A2) | A2 — weighted least-squares amplitude (fit to the net spectrum with weights that include the template variance) | METHOD.md:233–238; координатор 10.10 |
| мера E1 | E1 — total-variation estimator (minimal total variation of normalised shapes) | METHOD.md:233–238; координатор 10.10 |
| фотопик (устар. на странице) | full-energy peak | #GS-78: «фотопик» заменён «пиком полного поглощения» |
| сокращения (абзац «Сокращения») | Abbreviations | #GS-78 |
| хвост (пика, ядра) | (low-energy) tail | «low-energy tail» Kn 5, Gi 9 |
| калибровка по энергии; шкала энергии | energy calibration; energy scale | Kn 30, Gi 62 |
| живое время; реальное время; мёртвое время | live time; real time; dead time | Gi 54 / 12 / 114 |
| время измерения | counting time | Kn 16, Gi 15 |
| отсчёты; отсчётов на канал | counts; counts per channel | Kn 4, Gi 6 |
| скорость счёта | count rate | Gi 393 |
| остатки; нормированные остатки | residuals; normalized residuals | «residual» Kn 45 |
| χ²/ν | chi-square per degree of freedom (χ²/ν), reduced chi-square | «reduced chi» — N42.42-2020, Gulam Razul 2003; в Kn/Gi «chi-square(d)» |
| подгонка; окно подгонки | fit, least-squares fit; fit window | «fitting» Kn 16, Gi 36; «least-squares fit» Kn 2 |
| мешающий параметр | nuisance parameter | в библиотеке нет как термина (Kn «nuisance» 3); стандартный термин статистики (принят по брифу #GS-79 п.5) |
| поправка Бирге | Birge ratio correction | в библиотеке нет; общепринятое имя (R. T. Birge, 1932) — не сверено, оставлять «Birge» |
| поправка на распад | decay correction | Gi 14, FR 5 |
| цепочка (ряд) распада; звено цепочки | decay series (decay chain); chain member, progeny | «decay series» Gi 36; «decay chain» FR 103; «progeny» Kn 1, Gi 1 |
| вековое равновесие | secular equilibrium | Gi 13, FR 9 |
| выход γ-линии (квантов на распад) | γ-ray emission probability | Gi 136 («emission probabilit»); «intensity» — только для I_γ в обозначении |
| библиотека линий (нуклидов) | nuclide library, line library | «nuclide library» Gi 10 |
| удельная активность пробы (Бк/кг) | activity concentration | AQ48 25, FR 6; «massic activity» AQ48 52 — допустимый синоним, не смешивать |
| удельная активность природного калия (Бк/г K) | specific activity of natural potassium | «specific activity» Kn 14; «natural potassium» Gi 3 |
| образец с известной активностью | reference sample of known activity | «reference source» Gi 15 |
| известная (ожидаемая) активность | known (expected) activity | — |
| проба; матрица пробы | sample; sample matrix | «sample matrix» Gi 1 |
| сосуд Маринелли | Marinelli beaker | Kn 4, Gi 21, FR 1 |
| колодец (сосуда Маринелли) | re-entrant well | принят по описанию Marinelli beaker в Gi |
| насыпная плотность; влажность | bulk density; moisture content | «bulk density» Gi 1; «moisture» Kn 5 |
| тормозное излучение (внешнее) | bremsstrahlung (external) | Kn 57, Gi 29 |
| внутреннее тормозное излучение (IB) | internal bremsstrahlung (IB) | в Kn — «inner bremsstrahlung» 2 (синоним); «internal» — указание координатора 10.10 и обзор Italiano et al. 2024, цитируемый на странице |
| электронный захват | electron capture | Kn 13, Gi 31 |
| граничная энергия β-спектра | β endpoint energy | «endpoint energy» Kn 11 |
| разрешённый / запрещённый переход; уникальный запрет | allowed / forbidden transition; unique forbidden | Kn «forbidden» 13, «unique» |
| флуоресценция; Оже; деэкситация | fluorescence; Auger; atomic de-excitation | Kn 39 / есть; G4: `G4UAtomicDeexcitation`, параметры Fluo/Auger/PIXE |
| порог продукции | production cut (range cut) | G4: `G4ProductionCuts`, `SetCuts` |
| перенос (частиц) в веществе | particle transport | G4 |
| энергия, выделенная в кристалле | energy deposited in the crystal | Kn 47 |
| световыход; непропорциональность | light yield; nonproportionality | Kn 116; Kn 17 (писать слитно, «non-proportionality» 0) |
| ФЭУ; кремниевый фотоумножитель | photomultiplier tube (PMT); silicon photomultiplier (SiPM) | Kn 133; Kn 3 / «SiPM» 17 |
| отражатель; кварцевое окно; инкапсулированный | reflector; quartz window; encapsulated | Kn 11 / 3 / 22 |
| порог прибора | lower-level discriminator (LLD) threshold | Kn 2 |
| наложение пиков | pile-up | Kn 239, Gi 75 |
| неопределённость (статистическая / систематическая) | (statistical / systematic) uncertainty | Kn 6, Gi 11 |
| минимально детектируемая активность | minimum detectable activity (MDA) | Gi 19 |
| область интереса | region of interest (ROI) | Kn 3, Gi 5, FR 11 |
| карточка детектора | detector data set | бриф #GS-79 п.5 |
| стерильный пересчёт | independent recalculation | бриф #GS-79 п.5 |
| МК (метод Монте-Карло) | Monte Carlo (MC) | Kn 27, Gi 15 |
| сцинтилляционный гамма-спектрометр | scintillation gamma-ray spectrometer | «scintillation spectrometer» Kn 5 |
| сушёная черника; хлорид калия | dried blueberries; potassium chloride (KCl) | бытовые названия проб |

## Единая терминология страницы (решения координатора по стерильным проходам GS-79, 10.10.2026)

| Русский | English (единственный вариант на странице) | Не употреблять |
|---|---|---|
| сосуд (Маринелли) | Marinelli beaker / beaker; стенка — beaker wall | container, vessel |
| корпус прибора | (instrument) housing | enclosure |
| закон ПШПВ, закон ширины | FWHM(E) function (при первом упоминании — FWHM(E) power law) | FWHM law, width law |
| ветвь (ряд распада) | decay series / decay chain; шаблон ветви — decay-chain template | branch (кроме «branching ratio» = доля ветвления) |
| γ-квант | γ ray / photon | γ-quantum, quanta |
| третий уникальный запрет | third-forbidden unique transition | third unique forbidden |
| каскадное суммирование | true coincidence summing | true cascade summing |
| ионный источник (Geant4) | ion source (source of K-40 ions) | ionic source |
| слой (заливка вклада на графике) | component; заливка — shaded area | layer, fill |
| точки свёртки | FWHM points used in the convolution | convolution points |
| интерполяция в логарифмах | log–log interpolation | interpolation in logarithms |
| при найденной активности | at the fitted activity | determined activity |
| от среза (сосуда) | below the rim | from the cut |
| кювета | dish (sample dish) | cuvette |
| МК по нуклидам (вкладка) | nuclide-by-nuclide MC fit | MC by nuclides |
| протяжка мышью — приближение, двойной клик — сброс | drag to zoom, double-click to reset | mouse selection — zoom |
| СКО | RMS deviation | SD |
| тела (GDML) | solids | bodies |
| стенд | setup | stand |
| аттестация (образца) | certification; аттестованный — certified | accreditation |
| нелинейность (непропорциональность) световыхода | scintillation light-yield nonproportionality (NPSM) | — |
| K-рентген | K X-rays; пара — K and L X-rays | K-X-ray, x-ray |
| собственные линии детектора | intrinsic detector lines | detector's own lines |
| ожидаемая по массе | expected from the sample mass; метка — expected (from mass) | by mass |
| удельная (Бк/кг) | activity concentration | specific (без существительного) |
| «донорская» (страница/модель Гамма-1С) | the Gamma-1S page / model | donor |
| Даты | 9 October 2026 (ставит сборщик en_num) | 09.10.2026 |
| Тысячи | 17 938 с U+202F, только с 5 знаков (ставит сборщик en_num) | 17,938 |
| Библиография | J. Allison et al., Nucl. Instrum. Methods Phys. Res. A 835 (2016) 186–225 | «//», «Vol.», «pp.» |

## Имена собственные и правила

- Бета-1С → **Beta-1S**; Гамма-1С → **Gamma-1S**; СпектраЛайн → **SpectraLine**; BecqMoni, Atom-Spectra, Impulse — как есть.
- Am6er — только по нику.
- Метод 1 / метод 2 / метод 2-2 → **Method 1 / Method 2 / Method 2-2** (как имена разделов страницы).
- Маринелли 1 / Маринелли 2 (два конкретных сосуда) → **Marinelli 1 / Marinelli 2**.
- Кнапп–Уленбек, Блох, Льюис и Форд, Гебхардт → **Knipp–Uhlenbeck, Bloch, Lewis and Ford, Gebhardt**; KUB → KUB.
- Числа: десятичная **точка**, тысячи — пробел (как в русском варианте, ISO 80000-1); единицы keV, MeV, Bq, Bq/kg,
  g/cm³, mm, s, h; диапазон «150–3000 keV»; «×1,05» → «×1.05».
- Регистр: научный английский статей по гамма-спектрометрии (NIM A, Applied Radiation and Isotopes): безличные
  конструкции, без разговорных оборотов и кальки с русского («is calculated», «is not included», а не «we take it»).
- Запрещённое сокращение (#GS-29) на страницах не употребляется и в перевод не вводится.
