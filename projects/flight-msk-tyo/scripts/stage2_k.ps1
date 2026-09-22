# ASCII only (#PS-1). Stage II re-run with larger K for chosen files of stage1\prod3 -> stage2\<Out>.
# Usage: pwsh -NoProfile -File stage2_k.ps1 -Out prod3k100 -K 100 -Files neutron_0,neutron_1 [-Throttle 3]
param([string]$Src = "prod3", [string]$Out = "prod3k", [int]$K = 100, [string[]]$Files = @(), [int]$Throttle = 3, [int]$Seed = 500, [int]$IdxOff = 0, [switch]$Phantom)
$p = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$o = Join-Path $p "results\stage2\$Out"; New-Item -ItemType Directory -Force -Path $o | Out-Null
$Files = @($Files | ForEach-Object { $_ -split ',' }); $items = @(); $ix = 0
foreach ($b in $Files) {
  $f = Join-Path $p "results\stage1\$Src\$b.psp"; $ix++; $ob = $b
  if ($b -match '^(.*)_(\d+)$') { $ob = "$($Matches[1])_$([int]$Matches[2] + $IdxOff)" }
  $tsim = [double]((Get-Content "$f.meta" | Where-Object { $_ -like "T_sim_s=*" }) -replace "T_sim_s=", "")
  $items += [pscustomobject]@{ psp=$f; base=$ob; tsim=$tsim; seed=($Seed + $ix) }
}
$items | ForEach-Object -ThrottleLimit $Throttle -Parallel {
  $it = $_; . "C:\g4work\g4setup-11.4.2.ps1" | Out-Null
  if ($using:Phantom) { $env:CABIN_PHANTOM = "1" }; $env:CABIN_K = "$using:K"; $env:CABIN_ELL_CM = "120"; $t = Join-Path $using:o ($it.base + "_total.csv")
  for ($a = 0; $a -lt 6; $a++) {
    if (Test-Path $t) { Remove-Item $t }
    & "C:\g4work\build\flight-msk-tyo\cabin_stage2.exe" (Join-Path $using:o $it.base) ($it.seed + 1000 * $a) $it.tsim 4800 1 0 (Join-Path $using:p "tables\capture_db.txt") $it.psp *> (Join-Path $using:o ($it.base + ".log"))
    if ($LASTEXITCODE -eq 0 -and (Test-Path $t)) { "$($it.base) ok attempt=$($a+1)"; break }
    "$($it.base) FAILED attempt=$($a+1) rc=$LASTEXITCODE"
  }
}
"STAGE2 K DONE out=$Out K=$K files=$($Files.Count)"
