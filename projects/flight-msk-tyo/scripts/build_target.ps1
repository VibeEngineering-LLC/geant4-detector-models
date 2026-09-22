# ASCII only (#PS-1). Build ONE target of flight-msk-tyo (other exe may be locked by running jobs). Usage: pwsh -NoProfile -File build_target.ps1 -Target cabin_stage2
param([string]$Target = "cabin_stage2")
. "C:\g4work\g4setup-11.4.2.ps1" | Out-Null
$b = "C:\g4work\build\flight-msk-tyo"
$tools = "$vs\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin;$vs\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja"
cmd /c "`"$env:G4MODELS_VCVARS`" >nul 2>&1 && set `"PATH=$tools;$env:G4MODELS_SDK_BIN;%PATH%`" && cd /d `"$b`" && cmake --build . --target $Target 2>&1" *> "$b\build_$Target.log"
"build $Target rc=$LASTEXITCODE"
Get-Content "$b\build_$Target.log" -Tail 12
