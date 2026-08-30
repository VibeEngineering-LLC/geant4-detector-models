# Понуклидные шаблоны в домике, ДВЕ посадки параллельно (P-005, 30.08.2026).
# Машина 32-ядерная, прогон однопоточный и ест ~60 МБ: 16 процессов идут разом
# примерно за время одного (45 мин на 3e8) вместо шести часов подряд.
# Узел asb (asbuilt/down) повторяет ПРЕЖНЮЮ постановку и обязан воспроизвести
# опубликованные 29.08 шаблоны в пределах статистики — регрессия параметризации.
param([long]$NEvents = 300000000)
. C:\g4work\g4setup.ps1 | Out-Null
$exe = "D:\Claude_files\repos\geant4-detector-models\build\RadiaCode-103-field\rc103_field.exe"
$root = "D:\Claude_files\repos\geant4-detector-models\detectors\RadiaCode-103\geant4"
$nuc = @("K40", "Ra226", "Pb214", "Bi214", "Ac228", "Pb212", "Bi212", "Tl208")
$cfg = @(@{tag = "real"; stand = "25"; flip = "up" },
         @{tag = "asb"; stand = "asbuilt"; flip = "down" })
$seed = 2000
foreach ($c in $cfg) { foreach ($n in $nuc) {
    $seed += 137
    $log = "$root\run_field\logs\shield$($c.tag)_$n.log"
    $a = @("$root\run_roomfield\output\wf_room_$n.csv", $NEvents,
           "$root\run_field\output\rc103_field_room_shield$($c.tag)_$n.csv",
           "shield=on", "stand=$($c.stand)", "flip=$($c.flip)", "seed=$seed")
    Start-Process $exe -ArgumentList $a -RedirectStandardOutput $log `
        -RedirectStandardError "$log.err" -WindowStyle Hidden
} }
Write-Output "запущено $($nuc.Count * $cfg.Count) прогонов по $NEvents историй"
