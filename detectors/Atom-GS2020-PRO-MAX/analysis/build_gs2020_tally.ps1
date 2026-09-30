# GS-49 (30.09): build with decay_tally.hh (DECAY_TALLY mode) into SEPARATE dir; other exes untouched. ASCII only (#PS-1).
# Same configuration as build\gs2020-marinelli-cfg1 (Ninja, Debug, G4ROOT=C:\geant4-gdml). ASCII only (#PS-1).
$ErrorActionPreference = "Stop"
$vs = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
$B = "C:\g4work\build\gs2020-marinelli-tally"
Import-Module (Join-Path $vs "Common7\Tools\Microsoft.VisualStudio.DevShell.dll")
Enter-VsDevShell -VsInstallPath $vs -SkipAutomaticLocation -DevCmdArguments "-arch=x64 -host_arch=x64" | Out-Null
New-Item -ItemType Directory -Force $B | Out-Null
Set-Location $B
$ErrorActionPreference = "Continue"
cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DG4ROOT="C:\geant4-gdml" -DXERCESC_ROOT="C:\g4work\thirdparty\xerces-install" "C:\g4work\gs2020\run_marinelli" *> (Join-Path $B "configure.log")
if ($LASTEXITCODE -ne 0) { Get-Content (Join-Path $B "configure.log") -Tail 20; Write-Host "CONFIGURE_RC=$LASTEXITCODE"; exit 1 }
cmake --build . --target gs2020_marinelli *> (Join-Path $B "build.log")
$rc = $LASTEXITCODE
Get-Content (Join-Path $B "build.log") -Tail 15
Write-Host "BUILD_RC=$rc"
exit $rc
