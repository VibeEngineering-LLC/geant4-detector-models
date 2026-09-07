#include "NpsmBenchStackingAction.hh"
#include "NpsmBenchRun.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include <cstdio>
#include <cstdlib>

G4ClassificationOfNewTrack
NpsmBenchStackingAction::ClassifyNewTrack(const G4Track* track) {
    // Проверка на нулевой указатель
    if (!track) {
        std::fprintf(stderr, "Error: ClassifyNewTrack called with null track pointer\n");
        std::abort();
    }

    // Оптический фотон ищется ПО ИМЕНИ в таблице частиц, а не через
    // G4OpticalPhoton::OpticalPhotonDefinition(). Причина: этот вызов
    // СОЗДАЁТ определение, если его ещё нет, а в режиме light=npsm
    // оптической физики нет вовсе — и Geant4 роняет прогон предупреждением
    // «G4ParticleDefinition should be created in PreInit state», по одному
    // на каждый рабочий поток (поймано 07.09.2026 на первой же серии линии Г).
    // Поиск по имени определения не создаёт: нет оптики — вернётся nullptr,
    // и ни один трек с ним не совпадёт.
    static const G4ParticleDefinition* kOptical =
        G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
    if (kOptical != nullptr && track->GetDefinition() == kOptical) {
        // Получаем текущий объект прогона
        G4RunManager* runManager = G4RunManager::GetRunManager();
        if (!runManager) {
            std::fprintf(stderr, "Error: No run manager found\n");
            std::abort();
        }

        NpsmBenchRun* run = dynamic_cast<NpsmBenchRun*>(runManager->GetNonConstCurrentRun());
        if (!run) {
            std::fprintf(stderr, "Error: Failed to cast current run to NpsmBenchRun\n");
            std::abort();
        }

        // Увеличиваем счётчик фотонов
        run->AddScintPhoton();

        // Убиваем фотон сразу после подсчёта
        return fKill;
    }

    // Для всех остальных частиц — обычное поведение
    return fUrgent;
}
