# 2026-11-06: build g1s_gdml_export (Gamma-1S geometry -> GDML) against Geant4 with GDML (C:\geant4-gdml-vfix), recipe of build_gs2020_tagslce_vfix.ps1. ASCII only (#PS-1).
$ErrorActionPreference = "Stop"
$vs = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
$B = "C:\g4work\build\g1s-gdml-export"
$S = $PSScriptRoot   # folder of this script (CMakeLists.txt, main.cc)
Import-Module (Join-Path $vs "Common7\Tools\Microsoft.VisualStudio.DevShell.dll")
Enter-VsDevShell -VsInstallPath $vs -SkipAutomaticLocation -DevCmdArguments "-arch=x64 -host_arch=x64" | Out-Null
New-Item -ItemType Directory -Force $B | Out-Null
Set-Location $B
$ErrorActionPreference = "Continue"
cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DG4ROOT="C:\geant4-gdml-vfix" -DXERCESC_ROOT="C:\g4work\thirdparty\xerces-install" $S *> (Join-Path $B "configure.log")
if ($LASTEXITCODE -ne 0) { Get-Content (Join-Path $B "configure.log") -Tail 20; Write-Host "CONFIGURE_RC=$LASTEXITCODE"; exit 1 }
cmake --build . --target g1s_gdml_export *> (Join-Path $B "build.log")
$rc = $LASTEXITCODE
Get-Content (Join-Path $B "build.log") -Tail 15
Write-Host "BUILD_RC=$rc"
exit $rc
