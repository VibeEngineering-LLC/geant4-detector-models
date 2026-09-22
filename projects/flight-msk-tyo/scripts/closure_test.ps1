# ASCII only. Closure test: empty world + tube (stage I) -> replay (stage II) vs direct photon count on the detector. Result: results\stage2\closure\*
. "C:\g4work\g4setup-11.4.2.ps1"
$p = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$tab = Join-Path $p "tables\dadsgn\w76_h98_ip33.tab"; $db = Join-Path $p "tables\capture_db.txt"
$s1 = Join-Path $p "results\stage1\closure"; $s2 = Join-Path $p "results\stage2\closure"
$env:CABIN_EMPTY = "1"
& "C:\g4work\build\flight-msk-tyo\cabin_stage1.exe" $tab gamma 33 8000000 21 (Join-Path $s1 "gamma_0.psp") $db *> (Join-Path $s1 "gamma_0.psp.log")
Remove-Item Env:\CABIN_EMPTY
$tsim = [double]((Get-Content (Join-Path $s1 "gamma_0.psp.meta") | Where-Object { $_ -like "T_sim_s=*" }) -replace "T_sim_s=", "")
$env:CABIN_K = "10"
& "C:\g4work\build\flight-msk-tyo\cabin_stage2.exe" (Join-Path $s2 "replay") 5 $tsim 4800 0 0 $db (Join-Path $s1 "gamma_0.psp") *> (Join-Path $s2 "replay.log")
Remove-Item Env:\CABIN_K
$flux = [double]((Get-Content (Join-Path $s1 "gamma_0.psp.meta") | Where-Object { $_ -like "flux_cm2s=*" }) -replace "flux_cm2s=", "")
$tdir = 20000000 / ($flux * [math]::PI * 100.0)
$env:CABIN_DIRECT = "$tab|10|20000000|22"
& "C:\g4work\build\flight-msk-tyo\cabin_stage2.exe" (Join-Path $s2 "direct") 6 $tdir 4800 0 0 $db (Join-Path $s1 "gamma_0.psp") *> (Join-Path $s2 "direct.log")
"CLOSURE DONE tsim=$tsim"
