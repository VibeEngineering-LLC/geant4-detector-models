# Непропорциональность CsI(Tl) на низких энергиях: правки к сказанному ранее

Дата: 07.08.2026. Повод — возражение собеседника: «в величинах сечений нет такого
излома, ни в XCOM, ни в G4EMLOW». Возражение по существу справедливо, а мои
формулировки содержали четыре ошибки. Ниже каждая правка с дословной цитатой и
адресом первоисточника.

**Разметка уровня проверки:** ✅ прочитано по содержанию (цитата снята с
первоисточника) · 📗 реквизиты верны, содержимое не читалось · ⚠️ не подтверждено.

---

## 1. Не «излом», а разрыв. И знак был не назван

**Было сказано:** «излом у K-краёв 33,2/36,0 кэВ».

**Верно:** на K-крае у кривой непропорциональности **разрыв** (ступенька), а не
перегиб. При движении вверх по энергии световой выход на единицу поглощённой
энергии падает скачком; минимум («K-dip») лежит на 1–1,5 кэВ **выше** края,
дальше идёт восстановление.

> a linear increase from 111.2% to 115.8% with decrease of E_X from 100 keV to
> 50 keV followed by a drop in the range 30 – 45 keV with **a local minimum of
> 114.1% at 34.5 keV**. Next the photopeak-nPR increases up to 117.2% at 20 keV
> followed by a steep decrease of the response with further decrease of E_X. The
> nPR at 9 keV is 111.5%

✅ Khodyuk I. V., Rodnyi P. A., Dorenbos P. Non-proportional scintillation
response of NaI:Tl to low energy X-ray photons and electrons. *J. Appl. Phys.*
**107**, 113513 (2010). Препринт: https://arxiv.org/abs/1102.3799 (PDF прочитан
локально, §III). Измерение синхротронное, шаг сканирования у края 25 эВ.

Там же — что ступенька видна и в самом положении пика, не только в
пересчитанной кривой:

> In the inset of Fig. 4, the data near E_KI has been plotted on an expanded
> scale. Now, **a clear step can be seen in the N_phe exactly at E_KI**.

✅ тот же источник.

Материал в этой работе — NaI:Tl. Для CsI(Tl) отдельного синхротронного
измерения с мелким шагом у обоих краёв найти не удалось ⚠️; форма кривой
CsI(Tl) известна по Aitken 1967 и по обзору (см. п. 3).

## 2. Величина разрыва — около одного процентного пункта, а не 12 %

**Было сказано:** «подъём светового выхода к 10–20 кэВ (+12 % на 10 кэВ) и излом
у K-краёв» — две разные вещи в одной фразе, из-за чего читается, будто край и
даёт 12 %.

**Верно:** это два независимых явления разного масштаба.

| Что | Величина | Причина |
|---|---|---|
| Подъём отклика к 10–20 кэВ | +12 % на 10 кэВ, до +17 % в максимуме | рост плотности ионизации dE/dx вдоль трека при снижении энергии электрона |
| Разрыв на K-крае | ≈1 процентный пункт (115,2 → 114,1 % у NaI:Tl) | смена разбиения поглощённой энергии между фотоэлектроном и каскадом |

**Собеседник прав в узком смысле:** «излома в 12 %» в сечениях нет и быть не
может. Скачок в сечении на K-крае иода — это ×5,5 по массовому коэффициенту
ослабления (6,553 → 35,82 см²/г), а в световом выходе на том же крае — около
1 п.п. Величины разные, масштабы разные, и в одной фразе их смешивать нельзя.

✅ значения μ/ρ: NIST, Tables of X-Ray Mass Attenuation Coefficients, иод (Z=53),
https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z53.html

Поэтому фраза «он и есть причина» (о K-крае) в буквальном чтении неверна: край
объясняет только ступеньку в районе 33–36 кэВ, а не весь низкоэнергетический
подъём.

## 3. Ссылка была не на тот материал

**Было сказано:** в подтверждение разрыва приведена работа arXiv:1101.4485.

**Верно:** в ней измерены **Lu- и Gd-содержащие оксиды**, иода и цезия там нет
вовсе. Из аннотации:

> Using highly monochromatic X-ray synchrotron irradiation ranging from 9 keV to
> 100 keV, accurate **Lu2SiO5:Ce3+,Ca (LSO), Lu3Al5O12:Pr3+ (LuAG),
> Lu2Si2O7:Ce3+ (LPS) and Gd2SiO5:Ce3+ (GSO)** non-proportional response curves
> were determined.

✅ Khodyuk I. V., de Haas J. T. M., Dorenbos P. Non-proportional response between
0.1–100 keV energy by means of highly monochromatic synchrotron X-rays.
https://arxiv.org/abs/1101.4485 (PDF прочитан локально). Таблица II той же
работы: K-оболочка лютеция 63,314 кэВ, гадолиния 50,238 кэВ.

Кажущееся противоречие «провал против возрастания» снимается направлением
обхода. Там сказано:

> **When moving from high energy towards low energy**, we observe a relatively
> slow decrease of proportionality down to the K-edge energy… at the Lu or Gd
> K-edge energy the non-proportionality curve **increases** with a clear
> discontinuity at the edge.

✅ тот же источник. Знак закреплён числом в той же работе: «the number of
photoelectrons **falls** rapidly from 662 at 63.0 keV to 631 at 63.5 keV» —
63,0 кэВ ниже края, 63,5 выше, то есть выше края выход МЕНЬШЕ. Это то же самое,
что «провал», прочитанное с другой стороны.

**Правильные адреса по иодидам:**

- NaI:Tl, синхротрон, шаг 25 эВ — arXiv:1102.3799 (см. п. 1);
- CsI:Tl, форма кривой и величина на 10 кэВ — Khodyuk I. V., Dorenbos P. Trends
  and patterns of scintillator nonproportionality. https://arxiv.org/abs/1204.4350,
  таблица I, строка `CsI:Tl 5.8-6.2 57 1.78 112 6.78 [2, 5]`; колонки таблицы —
  «R,% at 662 keV | Light yield, photons/keV at 662 keV | Refractive index, n |
  **Photon-nPR,% at 10keV** | Degree of photon-nPR, σ,% | Reference». ✅ прочитано.
  Нормировка объявлена там же: все данные приведены к 100 % на 662 кэВ.
- форма для иодидов, дословно оттуда же:

  > CsI:Tl, NaI:Tl and CsI:Na, show the typical shape of the photon-nPR for
  > iodides [2] **with a maximum near 10 - 20 keV**

  ✅; и про края:

  > The only common feature noted was **a dip in the photon-nPR near K-shell and
  > L-shell absorption edges of iodine or calcium**

  ✅ (пересказ обзором работы Aitken D. W. et al., *IEEE Trans. Nucl. Sci.* **14**
  (1967) 468 — сама работа не читалась 📗).

Оговорка о провенансе числа 112: ссылка в строке таблицы — [2, 5], то есть
Aitken 1967 и Mengesha W., Taulbee T. D., Rooney B. D., Valentine J. D. Light
yield nonproportionality of CsI(Tl), CsI(Na), and YAP. *IEEE Trans. Nucl. Sci.*
**45**(3) 456–461 (1998), DOI 10.1109/23.682426 📗 (реквизиты подтверждены OSTI и
Crossref, полный текст за пейволлом IEEE). Собственная кривая CsI:Tl в обзоре
начинается с 12 кэВ, точки ровно на 10 кэВ там нет.

## 4. Непропорциональность в Geant4 есть штатно

**Было сказано:** «Geant4 считает энерговыделение, а не свет. Непропорциональности
в нём нет ни в option4, ни в G4EMLOW; она навешивается отдельной моделью
поверх».

**Верно:** механизм в тулките штатный, нет только ДАННЫХ по CsI(Tl).

Заголовок локальной установки Geant4 11.2.1,
`include/Geant4/G4Scintillation.hh`, строка 121 — комментарий к методу расчёта
числа фотонов:

> deposited (**includes nonlinear dependendency**) and updates the

✅ прочитано на диске. Там же объявлены `SetScintillationByParticleType` (строка
157) и `GetScintillationByParticleType` (161); переключатель вынесен и в
`G4OpticalParameters.hh` (`SetScintByParticleType`, строка 137).

Как считается: свойство материала `ELECTRONSCINTILLATIONYIELD` задаётся не
числом фотонов на МэВ, а **вектором «полный свет как функция кинетической
энергии»**, и на каждом шаге выдаётся разность L(T) − L(T − E_dep) — то есть
любая нелинейность кривой воспроизводится по построению. Исходник
`source/processes/electromagnetic/xrays/src/G4Scintillation.cc`, тег v11.2.1,
строки 874–889 📗 (в бинарной сборке .cc не поставляется, читан по репозиторию
https://github.com/Geant4/geant4/tree/v11.2.1).

Уточнение по `option4`: он к сцинтилляции отношения не имеет — это конструктор
электромагнитной физики; оптические фотоны рождает `G4OpticalPhysics`.

## 5. G4EMLOW не сводится к сечениям — K-края лежат прямо в нём

**Было сказано:** «XCOM и G4EMLOW содержат сечения; светового выхода в них нет
вообще».

Первая половина неполна. Светового выхода там действительно нет ✅, но G4EMLOW —
это ещё и энергии связи оболочек, вероятности флуоресцентных переходов и
оже-данные. K-края иода и цезия прописаны в нём явно:

```
файл  EMLOW8.5/fluor_XDB_EADL/binding.dat
Z=53 (I)  ->  оболочка K:  3.316900e-02 МэВ
Z=55 (Cs) ->  оболочка K:  3.598500e-02 МэВ
```

✅ прочитано на диске (`C:\g4conda\envs\g4data\share\Geant4\data\EMLOW8.5`).
Источник самих значений объявлен в `README_XDB_EADL`: X-Ray Data Booklet, LBNL.

Те же значения независимо: NIST XCOM для соединения CsI печатает строки
`53 K 3.317E-02` и `55 K 3.598E-02` ✅ (прогон XCOM,
https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html).

**О версии.** G4EMLOW 8.8 поставляется не с 11.2.x, а с **Geant4 11.4.0**;
с 11.2.1 идёт 8.5, с 11.3.0 — 8.6.1 (`cmake/Modules/G4DatasetDefinitions.cmake`
на тегах репозитория Geant4; локально то же подтверждается контрольной суммой в
`lib/cmake/Geant4/Geant4Config.cmake`, строка 227). Значит собеседник работает на
11.4. Для нашего спора это ничего не меняет: архив 8.8 скачан и сверен с
локальным 8.5 — **набор записей верхнего уровня совпадает, 28 против 28**, а
`README`, `fluor/binding.dat` и `fluor/fl-tr-pr-82.dat` совпадают побайтово.
Все изменения 8.5 → 8.8 лежат в Geant4-DNA, MicroElec и одном файле кремния. ✅

Проверка «нет ли в наборе чего-то, что можно назвать откликом», по архиву 8.8
целиком: поиск по именам файлов (`*scint*`, `*lightyield*`, `*optical*`,
`*response*`, `*efficien*`) — **ноль совпадений**; поиск по содержимому всех
11 396 файлов на `scintill|birks|light yield|quantum efficien|optical photon` —
**ноль совпадений**. Слово `response` встречается во всём наборе один раз — в
названии цитируемой статьи в `microelec/README`. ✅

Две ловушки, на которые можно опереться по недосмотру, и почему они не работают:
в `estar` есть колонка «Radiation yield», но README определяет её как долю
кинетической энергии, ушедшую в тормозное излучение, а не в свет; константы
Биркса в Geant4 есть, но не в G4EMLOW — они зашиты в исходник `G4EmSaturation.cc`
для четырёх материалов, и NaI/CsI среди них нет. ✅

Итог по формулировке: утверждение «в G4EMLOW нет светового выхода» — **верно**;
формулировка «там только сечения» — **неточна**. Там же формфакторы, энергии
связи, потенциалы ионизации, вероятности радиационных и оже-переходов, спектры
тормозного излучения, тормозные способности и угловые распределения.

## 6. Механизм был изложен неполно

**Было сказано:** «световой выход нелинеен по энергии электрона — значит при
смене разбиения он меняется скачком».

Это верно, но пропущено главное звено: сама нелинейность по энергии электрона
есть следствие роста плотности ионизации вдоль трека.

> According to the Bethe equation, with decrease in the energy of the primary
> electron, **the ionization density (dE/dx) along the track grows**. This leads
> to larger radiationless electron hole recombination rate which forms the basis
> of increasing nonproportionality

✅ arXiv:1204.4350 (Khodyuk & Dorenbos), раздел о механизмах.

И второе: выше края поглощённая энергия делится не на «электрон плюс что-то», а
на **несколько независимых треков** — медленный фотоэлектрон, характеристический
рентген (поглощается в другом месте кристалла или уходит вовсе, давая
escape-пик), либо оже-каскад. Суммарный свет есть сумма откликов на треки разной
энергии, а отклик нелинеен — отсюда и разрыв.

Формулировку «это кривая шкалы, а не потеря событий» тоже стоит смягчить: у
самого края к сдвигу шкалы добавляются escape-пики и скачок разрешения. Для
NaI:Tl измерено:

> A clear step-like change of almost 0.2% can be seen at E_X around E_KI

✅ arXiv:1102.3799 (о разрешении).

## 7. Первоисточник по CsI(Tl): Aitken 1967 — что там сказано на самом деле

Работа Aitken D. W., Beron B. L., Yenicay G., Zulliger H. R. The Fluorescent
Response of NaI(Tl), CsI(Tl), CsI(Na) and CaF2(Eu) to X-Rays and Low Energy Gamma
Rays. *IEEE Trans. Nucl. Sci.* **14**(1) 468–477 (1967),
DOI 10.1109/TNS.1967.4324457. Полный PDF за пейволлом IEEE; получены реферат
(дословно, из двух независимых записей — страница IEEE и выгрузка OpenAIRE) и
фрагменты страницы 472 (OCR отсканированной страницы журнала).

**Разрыв у CsI(Tl) подтверждён самим первоисточником** — из реферата дословно:

> **Both K-shell and L-shell absorption edge discontinuities are evident in the
> fluorescent response curves for Na(Tl) and CsI(Tl).** Only K-shell absorption
> edge discontinuities could be observed in the response curves for CsI(Na) and
> CaF2(Eu).

✅ (в реквизитах издателя опечатка «Na(Tl)» вместо NaI(Tl)).
https://ieeexplore.ieee.org/document/4324457

**Механизм, как его излагают сами авторы** (с. 472, дословно):

> For photon energies lying just above the K-shell absorption edge of iodine, for
> example, a photoelectric interaction takes place with high probability with the
> K-shell electron. The resultant electron emerges with very low energy, since
> most of the absorbed energy is lost to the atomic binding. Immediately below
> the K-shell energy value, however, photoelectric interactions can only take
> place with the L or higher shell electrons, and the resultant photoelectrons
> have higher energies than those produced by photons with energies just above
> the K-shell edge. As a result of the nonlinearity of the electron response of
> the crystal, **the scintillation pulse amplitude is relatively greater just
> below the absorption edge** as the crystal responds to the higher energy
> photoelectrons.

✅ Это ровно тот знак, что и в п. 1: переход через край ВВЕРХ по энергии
понижает отклик.

Важное уточнение к моей формулировке в п. 6. Авторы **не** говорят, что свет
есть сумма откликов на несколько треков. Они пишут обратное — энергия связи
списывается («most of the absorbed energy is lost to the atomic binding»), а
множественность электронов привлекается только для объяснения скачка
РАЗРЕШЕНИЯ:

> the number of electrons that are released through direct photoelectric
> interactions and secondary Auger transitions also increases, thereby improving
> the statistics of the scintillation process

✅ Значит «сумма откликов на несколько треков» — это позднейшая трактовка, а не
формулировка первоисточника; подавать её как то, что «пишут исследователи»,
нельзя.

**Расхождение источников, которое надо держать в уме.** Aitken (с. 472):

> Our results demonstrate that the fluorescent efficiency functions for NaI(Tl),
> CsI(Tl) and CsI(Na) all reach **maxima for photon energies of about 45 keV**,
> and then begin to decrease.

✅ А обзор arXiv:1204.4350, пересказывая ту же работу, говорит про «maximum near
10 – 20 keV» ✅. Числа разные втрое. ⚠️ Причина расхождения не установлена:
возможно, речь о разных величинах (кривая по фотонам против кривой по
электронам — Aitken тут же поясняет, что фотону 45 кэВ отвечает фотоэлектрон
около 12 кэВ), возможно — о разных участках кривой. **До выяснения ссылаться
надо на конкретную величину и конкретный источник, а не на «максимум при
10–20 кэВ» вообще.**

## Что остаётся неподтверждённым

- ⚠️ Разрешены ли в измерении ДВЕ отдельные ступеньки — на крае иода 33,17 и
  цезия 35,98 кэВ — не установлено. Aitken говорит про разрывы у K- и L-краёв
  для CsI(Tl), но раздельность краёв иода и цезия из доступного текста не видна;
  синхротронной работы по CsI(Tl) с мелким шагом не найдено. Корректная
  формулировка: «K-край иода 33,17 кэВ, следом край цезия 35,98 кэВ», без
  утверждения о двух разрешённых ступеньках.
- ⚠️ Расхождение «максимум при 45 кэВ» (Aitken) против «максимум при 10–20 кэВ»
  (обзор) не разобрано.
- 📗 Aitken 1967 по содержанию прочитан ЧАСТИЧНО: реферат целиком и страница 472;
  полного PDF нет (пейволл IEEE). Утверждать, что в остальных страницах чего-то
  нет, нельзя — отсутствие поиском по фрагментам не доказывается.
- 📗 Mengesha 1998 — полный текст за пейволлом, использован через обзор и через
  перепечатку рис. 1 в J. Appl. Phys. 105 (2009) 044507.

## Сводка источников

| Источник | Что берётся | Уровень |
|---|---|---|
| arXiv:1102.3799 (J. Appl. Phys. 107, 113513, 2010) | форма и знак разрыва у K-края, NaI:Tl, шаг 25 эВ | ✅ |
| arXiv:1101.4485 | тот же эффект в Lu/Gd-оксидах, величины скачков | ✅ |
| arXiv:1204.4350 | таблица I (CsI:Tl 112 % на 10 кэВ), форма кривой иодидов, механизм dE/dx | ✅ |
| Mengesha et al., IEEE TNS 45(3) 456 (1998), DOI 10.1109/23.682426 | первичные данные CsI(Tl) | 📗 |
| Aitken et al., IEEE TNS 14 (1967) 468 | первичные данные, провал у краёв | 📗 |
| NIST XCOM | K-края I и Cs, скачок μ/ρ | ✅ |
| NIST X-Ray Mass Attenuation Coefficients, z53 | величина скачка сечения | ✅ |
| Geant4 11.2.1, G4Scintillation.hh / G4OpticalParameters.hh | наличие штатного механизма нелинейного выхода | ✅ |
| Geant4 v11.2.1, G4Scintillation.cc:874-889 | формула L(T) − L(T − E_dep) | 📗 |
| G4EMLOW8.5, fluor_XDB_EADL/binding.dat | K-края I и Cs в данных Geant4 | ✅ |
