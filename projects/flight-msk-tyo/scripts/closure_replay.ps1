# ASCII only. Re-run replay of the closure stage-I file with window ELL (default 120 cm), K=10, into results\stage2\closure\replay*.
param([double]$Ell = 120.0, [string]$Out = "replay")
. "C:\g4work\g4setup-11.4.2.ps1"
$p = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$s1 = Join-Path $p "results\stage1\closure"; $s2 = Join-Path $p "results\stage2\closure"
$tsim = [double]((Get-Content (Join-Path $s1 "gamma_0.psp.meta") | Where-Object { $_ -like "T_sim_s=*" }) -replace "T_sim_s=", "")
$env:CABIN_K = "10"; $env:CABIN_ELL_CM = "$Ell"
& "C:\g4work\build\flight-msk-tyo\cabin_stage2.exe" (Join-Path $s2 $Out) 5 $tsim 4800 0 0 (Join-Path $p "tables\capture_db.txt") (Join-Path $s1 "gamma_0.psp") *> (Join-Path $s2 "$Out.log")
"closure replay rc=$LASTEXITCODE ell=$Ell"
