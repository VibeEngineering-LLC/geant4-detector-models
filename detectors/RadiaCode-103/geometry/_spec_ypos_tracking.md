# Спека: разложение спектра по Y-координате взаимодействия (Stage-2, #FIT-1)

## Контекст

Проверка гипотезы D-007 (light-collection tailing) требует не качественной прикидки,
а честного сопоставления: для каждого события в кристалле нужна energy-weighted
средняя Y-координата взаимодействия (та же ось Y, что в `opticalcheck.cc` — направление
к SiPM), чтобы потом в постобработке взвесить каждое событие по карте LCE(Y) из D-007
и получить apparent-energy спектр.

Готового трекинга позиции в `main.cc` НЕТ — проверено `Grep` по `geometry/*.cc`
(§33/ступень −1): `GetPosition` есть только в `shieldrun.cc` (другая задача, buildup),
в `main.cc` `Stepping::UserSteppingAction` копит только суммарный `fEdep`, без координат.

Образец механики — уже существующее разложение по каналам истории (`enum Chan`,
`_chan.csv`, строки main.cc:73-105, 180-228): по той же схеме добавляется РАЗЛОЖЕНИЕ
ПО Y, только гистограмма 2D (E_bin × Y_bin) вместо (E_bin × канал), отдельным файлом
`_ypos.csv`.

## Правка файла `main.cc` (единственный файл)

### 1. В классе `EventAct` (после поля `fEdep`, около строки 241)

Добавить накопитель энерговзвешенной позиции:

```cpp
double fEdepY = 0;   // накопитель Edep_i * Y_mid_i по шагам, для взвешенного среднего
```

В `BeginOfEventAction` сбрасывать `fEdepY = 0;` вместе с `fEdep = 0;`.

### 2. В `Stepping::UserSteppingAction` (после строки `fEvt->fEdep += s->GetTotalEnergyDeposit();`, ~337)

Добавить:

```cpp
{
  const double de = s->GetTotalEnergyDeposit();
  if (de > 0) {
    const double yMid = 0.5 * (pre->GetPosition().y() + post->GetPosition().y()) / mm;
    fEvt->fEdepY += de * yMid;
  }
}
```

(`pre`/`post` уже объявлены выше в той же функции — переиспользовать, не заводить новые.)

### 3. В `EventAct::EndOfEventAction` (около строки 283, вызов `fRun->Fill(...)`)

Посчитать средневзвешенную координату и передать пятым параметром:

```cpp
const double yMean = (fEdep > 0) ? (fEdepY / fEdep) : 0.0;
fRun->Fill(fEdep / keV, ep, Channel(), yMean);
```

### 4. В `RunAct` — новые константы и хранилище (после `kBinKeV`, ~строка 111)

```cpp
static constexpr int    kYBins   = 7;
static constexpr double kYBinMM  = 1.5;     // ширина Y-бина, мм
static constexpr double kYMinMM  = -5.25;   // нижняя граница первого бина
// Центры бинов при этих константах: -4.5,-3.0,-1.5,0.0,1.5,3.0,4.5 мм —
// СОВПАДАЮТ с точками карты LCE(Y) из opticalcheck.cc (D-007), чтобы
// постобработка сопоставляла бины напрямую, без интерполяции.
```

Хранилище (аналогично `fChan`, рядом с ним):

```cpp
std::vector<std::vector<long>> fYHist;   // [kYBins][kBins+1]
```

В конструкторе `RunAct()` инициализировать вместе с `fChan`:

```cpp
fYHist(kYBins, std::vector<long>(kBins + 1, 0))
```

В `BeginOfRunAction` сбрасывать вместе с `fChan`:

```cpp
for (auto& y : fYHist) std::fill(y.begin(), y.end(), 0L);
```

### 5. Метод `RunAct::Fill` — расширить сигнатуру (сохранить обратную совместимость)

Было: `void Fill(double edepKeV, double eprim, int chan = -1)`

Стало:

```cpp
void Fill(double edepKeV, double eprim, int chan = -1, double yMeanMM = 0.0) {
  fSumEprim += eprim;
  if (edepKeV <= 0) return;
  ++fWithSignal;
  int b = static_cast<int>(edepKeV / kBinKeV);
  if (b > kBins) b = kBins;
  ++fHist[b];
  if (chan >= 0 && chan < kNChan) ++fChan[chan][b];

  int iy = static_cast<int>((yMeanMM - kYMinMM) / kYBinMM);
  if (iy < 0) iy = 0;
  if (iy >= kYBins) iy = kYBins - 1;
  ++fYHist[iy][b];
}
```

(Единственный вызывающий — `EventAct::EndOfEventAction`, правка сигнатуры безопасна;
других вызовов `Fill` в файле нет — проверено `Grep`.)

### 6. Вывод `_ypos.csv` в `EndOfRunAction` (по образцу блока `_chan.csv`, ~строки 184-228)

Добавить СРАЗУ ПОСЛЕ блока записи `_chan.csv` (перед `G4cout << "RESULT..."`) отдельный
блок с тем же стилем (шапка с провенансом, проверка суммы):

```cpp
{
  G4String yn = fOut;
  const size_t dot = yn.rfind('.');
  yn = (dot == G4String::npos ? yn : yn.substr(0, dot)) + "_ypos.csv";
  FILE* g = std::fopen(yn.c_str(), "w");
  if (g) {
    std::fprintf(g, "# разложение отклика по Y-координате взаимодействия "
                    "(энерговзвешенное среднее за событие), ось как в opticalcheck.cc\n");
    std::fprintf(g, "# src_sha1 = %s\n", RC_SRC_SHA1);
    std::fprintf(g, "# git_describe = %s\n", RC_GIT_DESCRIBE);
    std::fprintf(g, "# particle = %s\n", fPart.c_str());
    std::fprintf(g, "# E_prim_keV = %.4f\n", eMean);
    std::fprintf(g, "# N_primaries = %ld\n", N);
    std::fprintf(g, "# N_with_signal = %ld\n", fWithSignal);
    std::fprintf(g, "# y_bin_mm = %.2f, y_min_mm = %.2f, y_bins = %d\n",
                 kYBinMM, kYMinMM, kYBins);
    std::fprintf(g, "E_keV");
    for (int y = 0; y < kYBins; ++y)
      std::fprintf(g, ",y%.2f", kYMinMM + (y + 0.5) * kYBinMM);
    std::fprintf(g, "\n");
    for (int i = 0; i <= kBins; ++i) {
      long rowSum = 0;
      for (int y = 0; y < kYBins; ++y) rowSum += fYHist[y][i];
      if (!rowSum) continue;
      std::fprintf(g, "%.1f", (i + 0.5) * kBinKeV);
      for (int y = 0; y < kYBins; ++y) std::fprintf(g, ",%ld", fYHist[y][i]);
      std::fprintf(g, "\n");
    }
    std::fclose(g);

    long sumY = 0, sumHist2 = 0;
    for (int i = 0; i <= kBins; ++i) {
      sumHist2 += fHist[i];
      for (int y = 0; y < kYBins; ++y) sumY += fYHist[y][i];
    }
    if (sumY != sumHist2)
      G4cerr << "ВНИМАНИЕ: сумма Y-бинов " << sumY
             << " не равна спектру " << sumHist2 << " — дефект биннинга Y" << G4endl;
    else
      G4cout << "YPOS_OK sum= " << sumY << " file= " << yn << G4endl;
  }
}
```

## Требования к приёмке (обязательно, не пропускать)

1. Компилируется без предупреждений по новому коду (тот же `build\RadiaCode-103` через
   `vcvars64.bat` + ninja).
2. **Мутационная проверка (#SA-3):** прогнать `rc_curves` коротким макросом (малое N,
   например Cs-137 моноэнергетический пучок) ДВАЖДЫ — с правкой и с намеренно
   испорченной версией (например, `yMid` захардкожен в 0 для всех событий).
   Ожидается: у испорченной версии ВСЕ события падают в один Y-бин (тот, что содержит
   0.0, т.е. `y0.00`), у рабочей — распределение по нескольким бинам. Если оба варианта
   дают одинаковое распределение — тест не показывает разницы, значит правка не
   работает, чинить ДО дальнейших шагов.
3. Проверка `YPOS_OK` (сумма бинов равна общей гистограмме) — если ушло в
   `ВНИМАНИЕ: сумма Y-бинов…`, останавливаться и разбирать, не игнорировать.
4. Существующий baseline-тест (`_chan.csv`, `CHAN_OK`) обязан остаться зелёным —
   не задет правкой, но проверить фактом (правка вставлена ПОСЛЕ существующего блока,
   не должна была его сломать).

## Не делать

- Не трогать `enum Chan`, `Channel()`, механику каналов — они не связаны с этой правкой.
- Не менять формат `spectrum.csv`/`_chan.csv` — только добавить новый файл `_ypos.csv`.
- Не переиспользовать `pre`/`post` из другого места функции неверно — они уже объявлены
  в начале `UserSteppingAction`, использовать именно их.
