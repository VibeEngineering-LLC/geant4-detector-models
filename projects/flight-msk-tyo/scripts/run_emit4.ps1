# ASCII only. Four line_bench runs with capture emitter (Fe, Al, Cu, Si) in parallel; outputs results\linebench\em2_<mat>.{csv,all}
. "C:\g4work\g4setup-11.4.2.ps1"
$p = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$exe = "C:\g4work\build\flight-msk-tyo\line_bench.exe"
$jobs = @()
foreach ($c in @(@("Fe",20,100000),@("Al",50,200000),@("Cu",20,100000),@("Si",100,200000))) {
  $m=$c[0]; $r=$c[1]; $n=$c[2]
  $jobs += Start-Job -ArgumentList $exe,$m,$r,$n,$p -ScriptBlock {
    param($exe,$m,$r,$n,$p)
    & $exe "G4_$m" 2.53e-8 $n 0 21 "$p\results\linebench\em2_$m.csv" $r "$p\tables\capture_db_try.txt" *> "$p\results\linebench\em2_$m.all"
  }
}
$jobs | Wait-Job | Out-Null
"ALL DONE"
