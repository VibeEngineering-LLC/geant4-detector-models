# Geant4 — RadiaCode-110

Полная документация: **[docs/GEANT4-MODEL.md](docs/GEANT4-MODEL.md)**

Быстрый старт:

```powershell
cd scripts
python stl_to_gdml.py
python verify_align.py
```

| Папка | Содержимое |
|---|---|
| `gdml/detector/` | SSOT: все внутренности + полый корпус-бокс |
| `gdml/full/` | STL-корпус + внутренности (генерируется) |
| `gdml/case_mesh/` | только STL (генерируется) |
| `cadmesh/` | C++ загрузчик STL через CADMesh |
| `scripts/` | конвертер и проверка |
| `verify/` | PNG совмещения |
