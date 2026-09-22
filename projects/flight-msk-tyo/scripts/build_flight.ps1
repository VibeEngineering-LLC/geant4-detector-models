# ASCII only (#PS-1). Build flight-msk-tyo apps with g4build from g4setup-11.4.2.ps1 (VS2022 cl + Ninja).
# Usage: pwsh -NoProfile -File build_flight.ps1 [-Configure]   (first time: -Configure)
param([switch]$Configure)
$ErrorActionPreference = "Continue"
. "C:\g4work\g4setup-11.4.2.ps1"
$src = Join-Path (Split-Path -Parent $PSScriptRoot) "src"
$b = "C:\g4work\build\flight-msk-tyo"
New-Item -ItemType Directory -Force -Path $b | Out-Null
if ($Configure) {
    $tools = "$vs\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin;$vs\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja"
    cmd /c "`"$env:G4MODELS_VCVARS`" >nul 2>&1 && set `"PATH=$tools;$env:G4MODELS_SDK_BIN;%PATH%`" && cd /d `"$b`" && cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=`"$env:G4MODELS_CL`" -DG4ROOT=$env:GEANT4_ROOT `"$src`" 2>&1" | Tee-Object "$b\configure.log" | Select-Object -Last 8
}
g4build $b *> "$b\build.log"
"build rc=$LASTEXITCODE"
Get-Content "$b\build.log" -Tail 30
