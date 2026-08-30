# Мюоны с домиком в РЕАЛЬНОЙ посадке (P-005), 5 порций по 1e8 параллельно.
# 5e8 в одном процессе — 335 мин (замер 30.08 на 1e7), за пределами W-040.
# Порции ~67 мин каждая, идут разом (машина 32-ядерная) — итог ~67 мин.
. C:\g4work\g4setup.ps1 | Out-Null
$exe = "D:\Claude_files\repos\geant4-detector-models\build\RadiaCode-103-muon\rc103_muon.exe"
$out = "D:\Claude_files\repos\geant4-detector-models\detectors\RadiaCode-103\geant4\run_muon\output"
$log = "D:\Claude_files\repos\geant4-detector-models\detectors\RadiaCode-103\geant4\run_muon\logs"
$seed = 3000
for ($i = 1; $i -le 5; $i++) {
    $seed += 191
    $a = @(100000000, "$out\rc103_muon_shieldreal_z300_r900_p$i.csv",
           "rdisk=900", "zdisk=300", "shield=on", "stand=25", "flip=up", "seed=$seed")
    Start-Process $exe -ArgumentList $a -RedirectStandardOutput "$log\muonreal_p$i.log" `
        -RedirectStandardError "$log\muonreal_p$i.log.err" -WindowStyle Hidden
}
Write-Output "запущено 5 порций по 1e8, seed от 3191"
