# Понуклидные шаблоны отклика в свинцовом домике, ПОСАДКА ЗАДАЁТСЯ ЯВНО (P-005).
# Прогоны идут по одному: 3e8 историй ~45 мин каждый, это внутри окна W-040
# (длиннее ~2 ч на этой машине не доживают). Зёрна разные — порции независимы.
# Запуск: run_shield_templates.ps1 <суффикс> <stand|asbuilt> <up|down> [n_events]
param([string]$Tag = "real", [string]$Stand = "25", [string]$Flip = "up",
      [long]$NEvents = 300000000)
. C:\g4work\g4setup.ps1 | Out-Null
$exe = "D:\Claude_files\repos\geant4-detector-models\build\RadiaCode-103-field\rc103_field.exe"
$root = "D:\Claude_files\repos\geant4-detector-models\detectors\RadiaCode-103\geant4"
$nuc = @("K40", "Ra226", "Pb214", "Bi214", "Ac228", "Pb212", "Bi212", "Tl208")
$seed = 1000
foreach ($n in $nuc) {
    $src = "$root\run_roomfield\output\wf_room_$n.csv"
    $out = "$root\run_field\output\rc103_field_room_shield${Tag}_$n.csv"
    $seed += 137
    $t = Measure-Command {
        & $exe $src $NEvents $out shield=on stand=$Stand flip=$Flip seed=$seed `
            *> "$root\run_field\logs\shield${Tag}_$n.log"
    }
    $hits = (Select-String -Path $out -Pattern "^n_hits_in_crystal,").Line
    Write-Output ("{0,-6} {1,-24} {2:N1} мин  seed={3}" -f $n, $hits, $t.TotalMinutes, $seed)
}
