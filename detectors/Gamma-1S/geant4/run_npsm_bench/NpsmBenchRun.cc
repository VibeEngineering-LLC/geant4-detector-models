#include "NpsmBenchRun.hh"

void NpsmBenchRun::RecordEvent(int nCompt, int /*nRayl*/, bool phot, bool conv,
                               bool escaped, double edepMeV, double edepLightMeV,
                               double energyKeV,
                               const std::vector<std::pair<double, double>>* pulses) {
    ++fNEvents;

    if (phot) {
        ++fNAbsorbed;
        fSumComptAbsorbed += static_cast<double>(nCompt);
        fSumCompt2Absorbed += static_cast<double>(nCompt * nCompt);
        fHistAbsorbed[std::min(nCompt, kMaxCompt)]++;
    } else if (conv) {
        ++fNConv;
    } else if (escaped) {
        ++fNEscaped;
    } else {
        ++fNOther;
    }

    // Увеличиваем общую статистику по всем событиям
    fHistAll[std::min(nCompt, kMaxCompt)]++;

    // Полное поглощение
    if (std::abs(edepMeV - energyKeV / 1000.0) < 0.001) {
        ++fNFullEdep;
        fSumLightFullMeV += edepLightMeV;
        fSumLight2FullMeV2 += edepLightMeV * edepLightMeV;
    }

    // Спектры: одна запись на событие (pulses == nullptr, прежний путь) либо
    // по записи на каждый импульс (ключ pulse_window_s, >1 импульса). Всё
    // остальное в этой функции — один раз на событие в любом режиме.
    const size_t nRec = pulses ? pulses->size() : 1;
    for (size_t k = 0; k < nRec; ++k) {
        const double e = pulses ? (*pulses)[k].first : edepMeV;
        const double l = pulses ? (*pulses)[k].second : edepLightMeV;

        // Спектр энерговыделения
        int indexEdep = static_cast<int>(e * 1000.0 / kBinKeV);
        if (indexEdep >= 0 && indexEdep < kNBins) {
            ++fSpecEdep[indexEdep];
        } else {
            ++fNOverflowEdep;
        }

        // Спектр света
        int indexLight = static_cast<int>(l * 1000.0 / kBinKeV);
        if (indexLight >= 0 && indexLight < kNBins) {
            ++fSpecLight[indexLight];
        } else {
            ++fNOverflowLight;
        }
    }

    // Суммы для статистики
    if (edepMeV > 0.0) {
        fSumEdepMeV += edepMeV;
        fSumLightMeV += edepLightMeV;
        fSumLight2MeV2 += edepLightMeV * edepLightMeV;
        ++fNWithEdep;
    }
}

void NpsmBenchRun::Merge(const G4Run* run) {
    // dynamic_cast, а не static_cast: static_cast от непустого указателя
    // НИКОГДА не даёт нуль, и проверка ниже ловила бы только run == nullptr,
    // то есть не то, ради чего она написана. Поймать надо именно чужой тип:
    // молча слитый не тот объект — заниженный результат при внешне успешном
    // прогоне (приёмка сгенерированного кода 07.09.2026).
    const NpsmBenchRun* benchRun = dynamic_cast<const NpsmBenchRun*>(run);
    if (!benchRun) {
        std::fprintf(stderr, "FATAL: NpsmBenchRun::Merge получил объект "
                             "другого типа или nullptr — результат потока "
                             "ПОТЕРЯН\n");
        std::fflush(stderr);
        std::abort();
    }

    // Складываем массивы
    for (int i = 0; i <= kMaxCompt; ++i) {
        fHistAbsorbed[i] += benchRun->fHistAbsorbed[i];
        fHistAll[i] += benchRun->fHistAll[i];
    }

    for (int i = 0; i < kNBins; ++i) {
        fSpecEdep[i] += benchRun->fSpecEdep[i];
        fSpecLight[i] += benchRun->fSpecLight[i];
    }

    // Складываем скаляры
    fNEvents += benchRun->fNEvents;
    fNAbsorbed += benchRun->fNAbsorbed;
    fNConv += benchRun->fNConv;
    fNEscaped += benchRun->fNEscaped;
    fNOther += benchRun->fNOther;
    fNFullEdep += benchRun->fNFullEdep;
    fNOverflowEdep += benchRun->fNOverflowEdep;
    fNOverflowLight += benchRun->fNOverflowLight;
    fNWithEdep += benchRun->fNWithEdep;

    fSumComptAbsorbed += benchRun->fSumComptAbsorbed;
    fSumCompt2Absorbed += benchRun->fSumCompt2Absorbed;
    fSumEdepMeV += benchRun->fSumEdepMeV;
    fSumLightMeV += benchRun->fSumLightMeV;
    fSumLight2MeV2 += benchRun->fSumLight2MeV2;
    fSumLightFullMeV += benchRun->fSumLightFullMeV;
    fSumLight2FullMeV2 += benchRun->fSumLight2FullMeV2;
    fNScintPhotons += benchRun->fNScintPhotons;

    G4Run::Merge(run);
}
