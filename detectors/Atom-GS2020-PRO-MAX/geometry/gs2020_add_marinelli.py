# -*- coding: utf-8 -*-
r"""Добавляет геометрию сосуда Маринелли с источником Th-232 (эпоксидная крошка) в GDML прибора GS2020.
Источник численных данных: README-референсы прибора раздел 1a, согласовано с оператором как таблица
#CFG-1, 2026-09-25. Входной/выходной GDML — переменные окружения GS2020_DEVICE_GDML/GS2020_GDML_OUT."""

import sys, os
sys.stdout.reconfigure(encoding="utf-8")

if not os.environ.get("GS2020_DEVICE_GDML") or not os.environ.get("GS2020_GDML_OUT"):
    raise RuntimeError("Переменные окружения GS2020_DEVICE_GDML и GS2020_GDML_OUT не установлены (см. README.md)")
SRC = os.environ["GS2020_DEVICE_GDML"]
DST = os.environ["GS2020_GDML_OUT"]

with open(SRC, "r", encoding="utf-8") as f:
    text = f.read()

old_worldbox = '<box name="WorldBox" x="200.0" y="200.0" z="200.0" lunit="mm"/>'
new_worldbox = '<box name="WorldBox" x="500.0" y="500.0" z="500.0" lunit="mm"/>'

if text.count(old_worldbox) == 0:
    print("ОШИБКА: WorldBox-строка не найдена, файл не сохранён", file=sys.stderr)
    sys.exit(1)

text = text.replace(old_worldbox, new_worldbox)

materials_xml = """<material name="Epoxy_crumb" state="solid"><D value="0.7987" unit="g/cm3"/><fraction n="0.7057" ref="G4_C"/><fraction n="0.0705" ref="G4_H"/><fraction n="0.2238" ref="G4_O"/></material>
<material name="Vessel_PP" state="solid"><D value="0.905" unit="g/cm3"/><fraction n="1.0" ref="G4_POLYPROPYLENE"/></material>"""

epoxy_marker = '<material name="Epoxy" state="solid"><D value="1.2" unit="g/cm3"/><fraction n="0.7057" ref="G4_C"/><fraction n="0.0705" ref="G4_H"/><fraction n="0.2238" ref="G4_O"/></material>'

try:
    idx = text.index(epoxy_marker)
except ValueError:
    print("ОШИБКА: Epoxy-маркер не найден, файл не сохранён", file=sys.stderr)
    sys.exit(1)

text = text[:idx + len(epoxy_marker)] + materials_xml + text[idx + len(epoxy_marker):]

positions_xml = ""  # ред. 3: вычитаний нет, позиций в <define> не нужно

try:
    idx_define = text.index("</define>")
except ValueError:
    print("ОШИБКА: </define> не найден, файл не сохранён", file=sys.stderr)
    sys.exit(1)

text = text[:idx_define] + positions_xml + text[idx_define:]

# Ред. 3: поликонусы вместо вложенных вычитаний (рендер ред. 2 зависал: BooleanProcessor too many edges).
# Кадр — центр корпуса, мир z = 40,5. Корпус PP с воздухом колодца r41 до z=+1,0 (мир 41,5 = верх прибора, W-135);
# матрица — дочерний объём корпуса: полость r70.5 z-52.5..52.5 минус вкладыш r43 до z=+3,0 (мир 43,5).
solids_xml = """<polycone name="VesselWall_solid" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
  <zplane z="-54.5" rmin="41" rmax="72.5"/><zplane z="1.0" rmin="41" rmax="72.5"/>
  <zplane z="1.0" rmin="0" rmax="72.5"/><zplane z="54.5" rmin="0" rmax="72.5"/></polycone>
<polycone name="SourceMatrix_solid" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
  <zplane z="-52.5" rmin="43" rmax="70.5"/><zplane z="3.0" rmin="43" rmax="70.5"/>
  <zplane z="3.0" rmin="0" rmax="70.5"/><zplane z="52.5" rmin="0" rmax="70.5"/></polycone>"""

try:
    idx_solids = text.index("</solids>")
except ValueError:
    print("ОШИБКА: </solids> не найден, файл не сохранён", file=sys.stderr)
    sys.exit(1)

text = text[:idx_solids] + solids_xml + text[idx_solids:]

volumes_xml = """<volume name="LV_SourceMatrix"><materialref ref="Epoxy_crumb"/><solidref ref="SourceMatrix_solid"/></volume>
<volume name="LV_VesselWall"><materialref ref="Vessel_PP"/><solidref ref="VesselWall_solid"/>
<physvol name="PV_LV_SourceMatrix"><volumeref ref="LV_SourceMatrix"/><position name="pos_SourceMatrix" x="0" y="0" z="0" unit="mm"/></physvol></volume>"""

try:
    idx_volumes = text.index('<volume name="World">')
except ValueError:
    print("ОШИБКА: <volume name=\"World\"> не найден, файл не сохранён", file=sys.stderr)
    sys.exit(1)

text = text[:idx_volumes] + volumes_xml + text[idx_volumes:]

placements_xml = """<physvol name="PV_LV_VesselWall"><volumeref ref="LV_VesselWall"/><position name="pos_VesselWall" x="0" y="0" z="40.5" unit="mm"/></physvol>"""

try:
    idx_physvol = text.rindex("</physvol>")
except ValueError:
    print("ОШИБКА: </physvol> не найден, файл не сохранён", file=sys.stderr)
    sys.exit(1)

text = text[:idx_physvol + len("</physvol>")] + placements_xml + text[idx_physvol + len("</physvol>"):]

import os
os.makedirs(os.path.dirname(DST), exist_ok=True)

with open(DST, "w", encoding="utf-8") as f:
    f.write(text)

print("записано: " + DST + " (" + str(len(text)) + " байт)")