# ASCII only. Response (effective area) runs: CABIN_DIRECT isotropic mono source R=10 cm on the bare instrument (no table), 3 in parallel, retry on crash.
param([long]$N = 200000)
$p = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Get-Content (Join-Path $p "results\eff\jobs.lst") | ForEach-Object -ThrottleLimit 3 -Parallel {
  . "C:\g4work\g4setup-11.4.2.ps1" | Out-Null
  $name, $pdg, $E, $tab = $_ -split " "; $p = $using:p; $n = $using:N
  if ($name -eq "neutron") { $n = 4 * $n }
  $pre = Join-Path $p ("results\eff\{0}_{1}" -f $name, $E); $tsim = $n / ([math]::PI * 100.0)
  $env:CABIN_DIRECT = "$(Join-Path $p $tab)|10|$n|$pdg"
  for ($a = 0; $a -lt 5; $a++) {
    & "C:\g4work\build\flight-msk-tyo\cabin_stage2.exe" $pre (11 + 100 * $a) $tsim 4800 0 0 (Join-Path $p "tables\capture_db.txt") (Join-Path $p "results\stage1\closure\gamma_0.psp") *> "$pre.log"
    if ($LASTEXITCODE -eq 0 -and (Test-Path "$($pre)_total.csv")) { break }
  }
  "$name $E rc=$LASTEXITCODE"
}
"EFF DONE"
