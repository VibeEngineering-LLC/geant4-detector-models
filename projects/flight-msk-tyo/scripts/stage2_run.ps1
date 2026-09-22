# ASCII only (#PS-1). Stage II: ONE PROCESS PER .psp FILE; output results\stage2\<Tag>\<file>_{total.csv,cat.csv,meta.txt}; merge later with stage2_merge.py.
# Usage: pwsh -NoProfile -File stage2_run.ps1 -Tag prod1 [-Table 1] [-Kb 0] [-Tflight 4800] [-Seed 7]
param([string]$Tag = "pilot1", [int]$Table = 1, [double]$Kb = 0.0, [double]$Tflight = 4800.0, [int]$Seed = 7)
# 22.09: скрипт не выставлял CABIN_ELL_CM -> C++ default 30 см занижает результат на 21-22% (audit/calc-sterile-2.md:62).
if (-not $env:CABIN_ELL_CM) { Write-Warning "CABIN_ELL_CM не задан - беру рабочее значение 120 (см. audit/calc-sterile-2.md:62)"; $env:CABIN_ELL_CM = "120" }
. "C:\g4work\g4setup-11.4.2.ps1"
$p = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$in = Join-Path $p "results\stage1\$Tag"; $out = Join-Path $p "results\stage2\$Tag"; New-Item -ItemType Directory -Force -Path $out | Out-Null
$items = @(); $k = 0
foreach ($f in Get-ChildItem $in -Filter "*.psp") {
  $tsim = [double]((Get-Content "$($f.FullName).meta" | Where-Object { $_ -like "T_sim_s=*" }) -replace "T_sim_s=", "")
  $k++; $items += [pscustomobject]@{ psp=$f.FullName; base=$f.BaseName; tsim=$tsim; seed=($Seed + $k) }
}
# retry on PIXE crash (W-066): new seed per attempt, up to 6 attempts; success = exit 0 and _total.csv present
$items | ForEach-Object -ThrottleLimit 3 -Parallel {
  $it = $_
  . "C:\g4work\g4setup-11.4.2.ps1" | Out-Null
  for ($a = 0; $a -lt 6; $a++) {
    if (Test-Path (Join-Path $using:out ($it.base + "_total.csv"))) { Remove-Item (Join-Path $using:out ($it.base + "_total.csv")) }
    & "C:\g4work\build\flight-msk-tyo\cabin_stage2.exe" (Join-Path $using:out $it.base) ($it.seed + 1000 * $a) $it.tsim $using:Tflight $using:Table $using:Kb (Join-Path $using:p "tables\capture_db.txt") $it.psp *> (Join-Path $using:out ($it.base + ".log"))
    if ($LASTEXITCODE -eq 0 -and (Test-Path (Join-Path $using:out ($it.base + "_total.csv")))) { "$($it.base) ok attempt=$($a+1)"; break }
    "$($it.base) FAILED attempt=$($a+1) rc=$LASTEXITCODE"
  }
}
"STAGE2 RUN DONE tag=$Tag files=$k table=$Table kB=$Kb Tflight=$Tflight"
