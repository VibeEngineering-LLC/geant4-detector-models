# ASCII only. Estimate of wings+fuel+engines: stage I with source disk R=900 cm, with and without CABIN_FUEL; components n, gamma, e-, e+.
param([double]$Tsim = 15.0)
$p = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$env:CABIN_SPHERE_R = "900"
Remove-Item Env:\CABIN_FUEL -ErrorAction SilentlyContinue
& pwsh -NoProfile -File (Join-Path $p "scripts\stage1_run.ps1") -Tsim $Tsim -Tag r900base -Nproc 1 -R 900 -Only "neutron,gamma,e-,e+"
$env:CABIN_FUEL = "1"
& pwsh -NoProfile -File (Join-Path $p "scripts\stage1_run.ps1") -Tsim $Tsim -Tag r900fuel -Nproc 1 -R 900 -Only "neutron,gamma,e-,e+"
Remove-Item Env:\CABIN_FUEL; Remove-Item Env:\CABIN_SPHERE_R
$env:CABIN_ELL_CM = "120"; $env:CABIN_K = "10"
& pwsh -NoProfile -File (Join-Path $p "scripts\stage2_run.ps1") -Tag r900base -Table 1 -Kb 0 -Tflight 4800 -Seed 7
& pwsh -NoProfile -File (Join-Path $p "scripts\stage2_run.ps1") -Tag r900fuel -Table 1 -Kb 0 -Tflight 4800 -Seed 7
"FUEL VARIANT DONE"
