# ASCII only (#PS-1). Final stage II pass: detector (Out prod3n2) and tissue phantom (Out dose3), all stage1 sets, one batch at a time.
# Usage: pwsh -NoProfile -File final_stage2.ps1 [-Throttle 14] [-Comps "proton,mum,mup,ep,k40"]
param([int]$Throttle = 14, [string]$Comps = "neutron,proton,mum,mup,em,ep,gamma,k40")
$s = Join-Path $PSScriptRoot "stage2_k.ps1"; $want = $Comps -split ","
# set: tag, files per comp, index offset, K for n/p/mu, K for e, K for gamma/k40, seed base, comps present
$sets = @(@("prod3b", 6, 3, 100, 30, 0, 100, "neutron,proton,mum,mup,em,ep"), @("prod3c", 8, 9, 25, 15, 0, 200, "proton,mum,mup,em,ep"), @("prod3", 3, 0, 100, 30, 10, 300, "neutron,proton,mum,mup,em,ep,gamma,k40"))
foreach ($ph in @($false, $true)) {
  $out = if ($ph) { "dose3" } else { "prod3n2" }
  foreach ($x in $sets) {
    $tag, $n, $off, $kn, $ke, $kg, $sd, $have = $x
    if ($ph) { $kn = [math]::Max(5, [int]($kn / 5)); $ke = [math]::Max(3, [int]($ke / 5)); $kg = 5 }
    foreach ($c in ($have -split ",")) {
      if ($want -notcontains $c) { continue }
      $k = if ($c -in @("em", "ep")) { $ke } elseif ($c -in @("gamma", "k40")) { $kg } else { $kn }
      $f = (0..($n - 1) | ForEach-Object { "${c}_$_" }) -join ","
      & pwsh -NoProfile -File $s -Src $tag -Out $out -K $k -Throttle $Throttle -Seed ($sd + [math]::Abs($c.GetHashCode() % 50)) -IdxOff $off -Phantom:$ph -Files $f
    }
  }
}
"FINAL STAGE2 DONE"
