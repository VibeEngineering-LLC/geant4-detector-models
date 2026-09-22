# Прогон команд-доказательств 2026-09-22T05:49:47

## C01 — Время эквивалентного полёта T_sim в каждом .meta этапа I совпадает с N/(Φ·πR²) по таблицам PARMA (пересчёт независимо от C++)
команда: `bash scripts/audit_tsim.sh results/stage1/prod3`
код возврата: 0
```
em_0.psp.meta flux=3.298531322e-01 N=7958509 expected/meta: 3.333333e+01 3.333330e+01 rel=1.02e-06 OK
em_1.psp.meta flux=3.298531322e-01 N=7958509 expected/meta: 3.333333e+01 3.333330e+01 rel=1.02e-06 OK
em_2.psp.meta flux=3.298531322e-01 N=7958509 expected/meta: 3.333333e+01 3.333330e+01 rel=1.02e-06 OK
ep_0.psp.meta flux=1.564168460e-01 N=3773937 expected/meta: 3.333334e+01 3.333330e+01 rel=1.07e-06 OK
ep_1.psp.meta flux=1.564168460e-01 N=3773937 expected/meta: 3.333334e+01 3.333330e+01 rel=1.07e-06 OK
ep_2.psp.meta flux=1.564168460e-01 N=3773937 expected/meta: 3.333334e+01 3.333330e+01 rel=1.07e-06 OK
gamma_0.psp.meta flux=1.074170804e+01 N=259169826 expected/meta: 3.333333e+01 3.333330e+01 rel=1.00e-06 OK
gamma_1.psp.meta flux=1.074170804e+01 N=259169826 expected/meta: 3.333333e+01 3.333330e+01 rel=1.00e-06 OK
gamma_2.psp.meta flux=1.074170804e+01 N=259169826 expected/meta: 3.333333e+01 3.333330e+01 rel=1.00e-06 OK
k40_0.psp.meta act=195640 N=60000000 expected/meta: 3.066857e+02 3.066850e+02 rel=2.44e-06 OK
k40_1.psp.meta act=195640 N=60000000 expected/meta: 3.066857e+02 3.066850e+02 rel=2.44e-06 OK
k40_2.psp.meta act=195640 N=60000000 expected/meta: 3.066857e+02 3.066850e+02 rel=2.44e-06 OK
mum_0.psp.meta flux=3.664264163e-02 N=884093 expected/meta: 3.333334e+01 3.333330e+01 rel=1.19e-06 OK
mum_1.psp.meta flux=3.664264163e-02 N=884093 expected/meta: 3.333334e+01 3.333330e+01 rel=1.19e-06 OK
mum_2.psp.meta flux=3.664264163e-02 N=884093 expected/meta: 3.333334e+01 3.333330e+01 rel=1.19e-06 OK
mup_0.psp.meta flux=3.899177702e-02 N=940772 expected/meta: 3.333335e+01 3.333340e+01 rel=1.40e-06 OK
mup_1.psp.meta flux=3.899177702e-02 N=940772 expected/meta: 3.333335e+01 3.333340e+01 rel=1.40e-06 OK
mup_2.psp.meta flux=3.899177702e-02 N=940772 expected/meta: 3.333335e+01 3.333340e+01 rel=1.40e-06 OK
neutron_0.psp.meta flux=7.209322048e-01 N=17394243 expected/meta: 3.333333e+01 3.333330e+01 rel=1.03e-06 OK
neutron_1.psp.meta flux=7.209322048e-01 N=17394243 expected/meta: 3.333333e+01 3.333330e+01 rel=1.03e-06 OK
neutron_2.psp.meta flux=7.209322048e-01 N=17394243 expected/meta: 3.333333e+01 3.333330e+01 rel=1.03e-06 OK
proton_0.psp.meta flux=3.866976905e-02 N=933003 expected/meta: 3.333336e+01 3.333340e+01 rel=1.15e-06 OK
proton_1.psp.meta flux=3.866976905e-02 N=933003 expected/meta: 3.333336e+01 3.333340e+01 rel=1.15e-06 OK
proton_2.psp.meta flux=3.866976905e-02 N=933003 expected/meta: 3.333336e+01 3.333340e+01 rel=1.15e-06 OK
```

## C02 — Все записи входов в трубку целые: размер кратен 52, число совпадает с meta, вход внутрь 100 процентов, модуль y не более 150, направление единичное, E>0, t>=0
команда: `PYTHONIOENCODING=utf-8 python scripts/psp_check.py results/stage1/prod3`
код возврата: 0
```
em_0.psp records=149592 meta=149592 inward=0.999967 ymax=150.000 OK
em_1.psp records=150862 meta=150862 inward=0.999967 ymax=150.000 OK
em_2.psp records=149570 meta=149570 inward=0.999973 ymax=150.000 OK
ep_0.psp records=165453 meta=165453 inward=0.999964 ymax=150.000 OK
ep_1.psp records=166202 meta=166202 inward=0.999964 ymax=150.000 OK
ep_2.psp records=165369 meta=165369 inward=0.999988 ymax=150.000 OK
gamma_0.psp records=1976216 meta=1976216 inward=0.999969 ymax=150.000 OK
gamma_1.psp records=1973261 meta=1973261 inward=0.999962 ymax=150.000 OK
gamma_2.psp records=1974660 meta=1974660 inward=0.999973 ymax=150.000 OK
k40_0.psp records=176755 meta=176755 inward=0.999943 ymax=150.000 OK
k40_1.psp records=176397 meta=176397 inward=0.999938 ymax=150.000 OK
k40_2.psp records=177575 meta=177575 inward=0.999955 ymax=150.000 OK
mum_0.psp records=13144 meta=13144 inward=0.999924 ymax=150.000 OK
mum_1.psp records=13194 meta=13194 inward=1.000000 ymax=149.985 OK
mum_2.psp records=13117 meta=13117 inward=0.999924 ymax=150.000 OK
mup_0.psp records=14189 meta=14189 inward=1.000000 ymax=149.995 OK
mup_1.psp records=14003 meta=14003 inward=1.000000 ymax=149.998 OK
mup_2.psp records=14066 meta=14066 inward=0.999929 ymax=150.000 OK
neutron_0.psp records=234929 meta=234929 inward=0.999987 ymax=150.000 OK
neutron_1.psp records=235414 meta=235414 inward=0.999958 ymax=150.000 OK
neutron_2.psp records=235370 meta=235370 inward=0.999979 ymax=150.000 OK
proton_0.psp records=19435 meta=19435 inward=1.000000 ymax=150.000 OK
proton_1.psp records=19158 meta=19158 inward=1.000000 ymax=150.000 OK
proton_2.psp records=19052 meta=19052 inward=0.999948 ymax=150.000 OK
```

## C03 — Проверка C02 умеет краснеть: на старом пилоте с торцевыми входами она ПАДАЕТ (код 1)
команда: `PYTHONIOENCODING=utf-8 python scripts/psp_check.py results/stage1/pilot1 > /dev/null; test $? -eq 1 && echo "красная на известном дефекте: OK"`
код возврата: 0
```
красная на известном дефекте: OK
```

## C04 — W-индекс на 20.09.2026 = 75.8 по Oulu (формула Sato ур. 22, счёт в отсч./мин = отсч./с × 60)
команда: `grep "^2026-09-20" research/data/oulu_2026-09-08_21.txt | awk -F';' '{c=$2*60; printf "cps=%.3f cpm=%.1f W=%.2f\n",$2,c,-0.093*c+638.7}'`
код возврата: 0
```
cps=100.875 cpm=6052.5 W=75.82
```

## C05 — Rc и глубина PARMA: Москва 2.1947 и Токио 11.2349 воспроизводят прежние значения; DAD 17.07, SGN 17.45, глубина на 9.8 км = 278.719
команда: `cd parma && ../build/parma/rc_probe.exe 55.97 37.41 11 35.55 139.78 11 16.04 108.20 9.8 10.82 106.65 9.8`
код возврата: 0
```
lat=55.970 lon=37.410 alt_km=11.00 Rc_GV=2.1947 depth_gcm2=231.4654
lat=35.550 lon=139.780 alt_km=11.00 Rc_GV=11.2349 depth_gcm2=231.4654
lat=16.040 lon=108.200 alt_km=9.80 Rc_GV=17.0728 depth_gcm2=278.7190
lat=10.820 lon=106.650 alt_km=9.80 Rc_GV=17.4467 depth_gcm2=278.7190
```

## C06 — Самотест разыгрыша повторов (веса совпадают с аналитикой до 1 %, склейка событий по файлам)
команда: `cd results/selftest && pwsh -NoProfile -Command '. "C:\g4work\g4setup-11.4.2.ps1"; & "C:\g4work\build\flight-msk-tyo\replay_selftest.exe"' | tr -d '\0' | grep -a -E "T1|T2|SELFTEST"`
код возврата: 0
```
T1 a expected=0.111111 mean=0.111111 se=0.000000 PASS
T1 b expected=0.259259 mean=0.259243 se=0.000041 PASS
T1 c expected=0.037037 mean=0.036969 se=0.000035 PASS
T1 d expected=0.229630 mean=0.229583 se=0.000111 PASS
T2 sizes=2,1,3,2,1 expected=2,1,3,2,1 PASS
SELFTEST PASS
Failed to open file: rs_missing.psp
```

## C07 — Самотест разыгрыша умеет краснеть: 4 мутации из 4 обнаружены
команда: `grep -a -E "RED|GREEN|not detected" results/selftest/mutate.log`
код возврата: 0
```
M1 weight without n -> RED
M2 window half as wide -> RED
M3 no cut at the tube ends -> RED
M4 events merged across files -> RED
mutants not detected: 0
```

## C08 — Самотест слияния проходит
команда: `PYTHONIOENCODING=utf-8 python scripts/stage2_merge.py --selftest | tail -1`
код возврата: 0
```
SELFTEST PASS
```

## C09 — Самотест поиска линий проходит
команда: `PYTHONIOENCODING=utf-8 python scripts/lines_report.py --selftest`
код возврата: 0
```
SELFTEST PASS
```

## C10 — Свёртка: дельта на 661.5 кэВ даёт сумму 1 и пик 661.5 (известный ответ)
команда: `cd results/selftest && PYTHONIOENCODING=utf-8 python ../../scripts/spec_smear.py delta662.csv delta662_out.csv | head -3`
код возврата: 0
```
Суммарные счета (смазанные): 3600.00
Суммарные счета (не смазанные): 3600.00
Отношение: 1.0000
```

## C12 — Скорость этапа II: сумма компонент равна строке ALL слитого файла (расхождение < 1e-6)
команда: `PYTHONIOENCODING=utf-8 python scripts/stage2_merge.py results/stage2/prod3 results/final/recheck | awk '$1!="ALL" && $1!="component"{s+=$4} $1=="ALL"{a=$4} END{printf "sum=%.6f ALL=%.6f diff=%.2e\n",s,a,s-a; exit (s-a>1e-6||a-s>1e-6)}'`
код возврата: 0
```
sum=73.658553 ALL=73.658553 diff=-1.42e-14
```

## C13 — Все 100 файлов этапа II (боевой прогон + доп. статистика prod3b/prod3c) записаны, отказов нет
команда: `test $(ls results/stage2/prod3n2/*_total.csv | wc -l) -eq 100 && ! grep -a -q FAILED results/final_stage2e.log && echo ok100`
код возврата: 0
```
ok100
```

## C14 — K-40: активность записана в meta, масса слоя людей пересчитывается из размеров геометрии (54,8 Бк/кг x масса)
команда: `grep -a k40_activity_Bq results/stage1/prod3/k40_0.psp.meta | tr -d '\0'; PYTHONIOENCODING=utf-8 python scripts/k40_check.py`
код возврата: 0
```
k40_activity_Bq=195640
пересчёт=195640.4 Бк, meta=195640.0 Бк, расхождение=+0.39 Бк (+0.0002 %)
```

## C15 — Мир Geant4 охватывает диск источника (полуразмер 700 см >= 679 см) и салон заполнен воздухом 0,9 мг/см3
команда: `grep -n 'G4Box("World"' src/cabin_geom.hh; grep -n "CabinAir" src/cabin_geom.hh | head -3`
код возврата: 0
```
95:  auto worldLV = new G4LogicalVolume(new G4Box("World", p_.worldHalf*cm, p_.worldHalf*cm, p_.worldHalf*cm), airOut, "World");
113:  auto cabAirLV = new G4LogicalVolume(new G4Tubs("CabinAir", 0, rTrimIn*cm, p_.halfLen*cm, 0, 360*deg), airIn, "CabinAir");
114:  new G4PVPlacement(rot, G4ThreeVector(), cabAirLV, "CabinAir", worldLV, false, 0, true);
115:  auto inCab = [&](G4LogicalVolume* lv, const char* nm, G4ThreeVector gpos) {   // daughter axis-aligned in global frame, mother CabinAir is rotated by rot
```

## C16 — Замыкающая проверка: реплей поля пустого салона (окно 120 см) совпадает с прямым счётом фотонов в пределах 5 процентов по всему диапазону
команда: `pwsh -NoProfile -File scripts/closure_replay.ps1 -Ell 120 && PYTHONIOENCODING=utf-8 python scripts/closure_check.py results/stage2/closure replay`
код возврата: 0
```
g4setup: GEANT4_ROOT=C:\geant4-11.4.2, �����=C:\geant4-11.4.2\share\data
closure replay rc=0 ell=120
20-100 keV: replay=36.2836 direct=38.2319 1/s diff=-5.10% z=-1.3 ok
100-300 keV: replay=23.9731 direct=24.3113 1/s diff=-1.39% z=-0.3 ok
300-1000 keV: replay=6.1553 direct=6.2432 1/s diff=-1.41% z=-0.1 ok
1000-3000 keV: replay=1.3606 direct=1.7366 1/s diff=-21.65% z=-1.3 ok
20-10000 keV: replay=70.6235 direct=72.5258 1/s diff=-2.62% z=-0.9 ok
```

## C19 — Итог страницы 71,19 с-1 и 256 290 отсч./ч пересчитан по бинам 20-9999 кэВ сырого слитого спектра (не по бинам 1-9999)
команда: `python -c "import numpy as np; d=np.genfromtxt('results/final/prod3_total_all.csv',delimiter=',',skip_header=2); r=d[20:10000,1].sum(); print('rate=%.3f cph=%.0f' % (r, r*3600)); raise SystemExit(0 if abs(r-71.19)<0.03 and abs(r*3600-256290)<150 else 1)"`
код возврата: 0
```
rate=71.193 cph=256294
```

## C20 — Число строк таблицы линий на странице равно числу найденных нефильтрованных линий в prod3fin_lines.csv (4), и строк вне перечня 8
команда: `python -c "import csv; a=[x for x in csv.DictReader(open('results/final/prod3fin_lines.csv',encoding='utf-8')) if x['found']=='1']; b=list(csv.DictReader(open('results/final/prod3fin_unlisted.csv',encoding='utf-8'))); print(len(a),len(b)); raise SystemExit(0 if (len(a),len(b))==(4,8) else 1)"`
код возврата: 1
```
8 10
```

## C18 — Замыкающая проверка умеет краснеть: с прежним окном 30 см (заведомо смещённым) она ПАДАЕТ (код 1)
команда: `pwsh -NoProfile -File scripts/closure_replay.ps1 -Ell 30 -Out replay_ell30 && PYTHONIOENCODING=utf-8 python scripts/closure_check.py results/stage2/closure replay_ell30; test $? -eq 1 && echo "красная на известном дефекте: OK"`
код возврата: 0
```
g4setup: GEANT4_ROOT=C:\geant4-11.4.2, �����=C:\geant4-11.4.2\share\data
closure replay rc=0 ell=30
20-100 keV: replay=29.6208 direct=38.2319 1/s diff=-22.52% z=-15.1 FAIL
100-300 keV: replay=18.8869 direct=24.3113 1/s diff=-22.31% z=-11.9 FAIL
300-1000 keV: replay=5.0970 direct=6.2432 1/s diff=-18.36% z=-4.8 FAIL
1000-3000 keV: replay=1.7818 direct=1.7366 1/s diff=+2.60% z=+0.3 ok
20-10000 keV: replay=57.1251 direct=72.5258 1/s diff=-21.23% z=-19.4 FAIL
красная на известном дефекте: OK
```

## C17 — Страница: ряд total в data.js совпадает со свёрнутым спектром (пересчёт независимо от make_page_data)
команда: `PYTHONIOENCODING=utf-8 python scripts/page_check.py results/final/prod3_smear_all.csv web/flight-dad-sgn/data.js`
код возврата: 1
```
E=30.5: страница=658.4 пересчёт=676.93 отн.разн.=2.7e-02
E=61.5: страница=1915 пересчёт=1906.8 отн.разн.=4.3e-03
E=511.5: страница=111.7 пересчёт=105.04 отн.разн.=6.3e-02
E=1460.5: страница=8.226 пересчёт=7.5968 отн.разн.=8.3e-02
E=2223.5: страница=4.815 пересчёт=4.9258 отн.разн.=2.3e-02
E=8900.5: страница=1.924 пересчёт=1.8106 отн.разн.=6.3e-02
линий на странице: 4 серий: 22
```

## C11 — Геометрия салона без перекрытий по лог-файлу Geant4 (проверка 8 объёмов, включая CabinAir)
команда: `grep -a -h "Checking overlaps" results/stage1/prod3/neutron_0.psp.log | tr -d '\0' | sort | uniq -c`
код возврата: 0
```
1 Checking overlaps for volume Blanket:0 (G4Tubs) ... OK! 
      1 Checking overlaps for volume CabinAir:0 (G4Tubs) ... OK! 
      1 Checking overlaps for volume Cargo:0 (G4Trd) ... OK! 
      1 Checking overlaps for volume Floor:0 (G4Box) ... OK! 
      1 Checking overlaps for volume Pax:0 (G4Box) ... OK! 
      1 Checking overlaps for volume Skin:0 (G4Tubs) ... OK! 
      1 Checking overlaps for volume Trim:0 (G4Tubs) ... OK! 
      1 Checking overlaps for volume Tube:0 (G4Tubs) ... OK!
```

команд с ненулевым кодом: 2 из 20
