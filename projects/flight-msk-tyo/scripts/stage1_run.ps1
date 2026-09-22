# ASCII only (#PS-1). Stage I: 7 components in parallel, Nproc processes each; files results\stage1\<Tag>\<comp>_<k>.psp (+ .meta, .log).
# Usage: pwsh -NoProfile -File stage1_run.ps1 -Tsim 2.0 -Tag pilot1 [-Nproc 2]
param([double]$Tsim = 2.0, [string]$Tag = "pilot1", [int]$Nproc = 2, [string]$Prefix = "tables\dadsgn\w76_h98_ip", [double]$R = 480.0, [string]$Only = "", [int]$SeedOff = 0)
. "C:\g4work\g4setup-11.4.2.ps1"
$p = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$out = Join-Path $p "results\stage1\$Tag"; New-Item -ItemType Directory -Force -Path $out | Out-Null
$comps = @("neutron:0:0", "proton:1:1", "mu+:2:29", "mu-:3:30", "e-:4:31", "e+:5:32", "gamma:6:33"); $jobs = @()
foreach ($c in $comps) {
  $name, $idx, $ip = $c -split ":"; $tab = Join-Path $p "$Prefix$ip.tab"
  $flux = [double](((Get-Content $tab -TotalCount 2)[1] -split "\s+")[5])
  if ($Only -and ($Only -split ",") -notcontains $name) { continue }
  $n = [long][math]::Ceiling($flux * [math]::PI * $R * $R * $Tsim / $Nproc)
  for ($k = 0; $k -lt $Nproc; $k++) {
    $psp = Join-Path $out ("{0}_{1}.psp" -f $name.Replace("+","p").Replace("-","m"), $k)
    $jobs += Start-Job -ArgumentList $tab,$name,$idx,$n,(1000 * ([int]$idx + 1) + 17 * $k + 1 + $SeedOff),$psp,(Join-Path $p "tables\capture_db.txt") -ScriptBlock {
      param($tab,$name,$idx,$n,$seed,$psp,$db)
      & "C:\g4work\build\flight-msk-tyo\cabin_stage1.exe" $tab $name $idx $n $seed $psp $db *> "$psp.log"
    }
  }
}
$jobs | Wait-Job | Out-Null
"STAGE1 RUN DONE tag=$Tag Tsim=$Tsim Nproc=$Nproc"
