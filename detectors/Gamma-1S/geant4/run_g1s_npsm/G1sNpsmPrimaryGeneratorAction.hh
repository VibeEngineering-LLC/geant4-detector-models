#ifndef G1S_NPSM_PRIMARY_GENERATOR_ACTION_HH
#define G1S_NPSM_PRIMARY_GENERATOR_ACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"

#include <string>

class G1sNpsmPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    static double gSourceZmm;
    static double gEnergyKeV;

    // Что вылетает из источника. "gamma" — одиночный квант gEnergyKeV, как было
    // всегда. "ion" — покоящееся ядро (gIonZ, gIonA), которое распадается сам
    // Geant4: тогда каскад испускается целиком и в спектре появляется каскадное
    // суммирование, которого при одиночном кванте нет и быть не может.
    // Режим введён 07.09.2026: расчёты этапа 4 по Co-60 шли одиночным квантом
    // 1173,23 кэВ, а измерялся реальный источник, где каскад 1173+1332 есть
    // всегда — эта разница нигде не была учтена и даже не записана.
    static std::string gPrimary;
    static int gIonZ;
    static int gIonA;

    // Откуда вылетает первичка. "point" — из точки на оси (gSourceZmm), как
    // было; "sample" — равномерно по ОБЪЁМУ пробы (логический объём "Sample",
    // сосуд Маринелли/Дента/Петри). Второй режим нужен для источников
    // поверочного комплекта: они не точечные, а засыпаны в кювету, и
    // самопоглощение в матрице — часть отклика, а не поправка к нему.
    // Введён 11.09.2026 для разбора смеси AmTiCsEu (Am-241 59,5 кэВ, где
    // самопоглощение в матрице наибольшее).
    static std::string gSourceMode;

    // Ограничение ОБЛАСТИ розыгрыша внутри пробы — проверка гипотезы о
    // неоднородном распределении активности (осадок за 14 лет хранения).
    // gSrcZFrac: доля высоты пробы СНИЗУ, в которой разыгрывается источник
    //            (1,0 — вся высота, прежнее поведение; 0,2 — нижняя пятая часть).
    // gSrcRFrac: доля радиального размаха СНАРУЖИ (1,0 — весь, 0,3 — внешняя
    //            треть кольца, то есть осадок у внешней стенки).
    // Оба параметра нужны, чтобы ИЗМЕРИТЬ, объясняет ли неоднородность
    // расхождение по мягкой линии, не ломая жёсткую, — а не подогнать одну.
    static double gSrcZFrac;
    static double gSrcRFrac;

public:
    G1sNpsmPrimaryGeneratorAction();
    ~G1sNpsmPrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event* anEvent) override;

private:
    // Точка равномерно внутри тела пробы. Отбор с отклонением по
    // ограничивающему параллелепипеду самого солида: тело — G4Polycone
    // (кольцо вокруг колодца + слой над ним), аналитического розыгрыша для
    // него нет, а Inside() у солида точен. Объём ищется по имени в
    // G4LogicalVolumeStore при первом событии: в многопоточном режиме
    // хранилище общее, а указатель из детектора мастер-потока брать нельзя.
    G4ThreeVector SamplePointInSample();

    G4ParticleGun fGun;
    G4ParticleDefinition* fIon = nullptr;  // создаётся при первом событии
    class G4VSolid* fSampleSolid = nullptr;
    G4ThreeVector fSampleMin, fSampleMax;  // ограничивающий параллелепипед
};

#endif // G1S_NPSM_PRIMARY_GENERATOR_ACTION_HH
