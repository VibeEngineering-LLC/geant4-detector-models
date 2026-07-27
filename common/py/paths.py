"""Корни путей — из переменных окружения, без привязки к машине.

В репозитории не должно быть ни одного абсолютного пути, привязанного к
машине: иначе код не запустится ни у кого, кроме автора, а имя пользователя
навсегда останется в истории коммитов.

Переменные окружения (все необязательные, есть разумные значения по умолчанию):

  G4MODELS_BUILD_<ДЕТЕКТОР>
                   каталог сборки и расчётных спектров ОДНОГО прибора,
                   например G4MODELS_BUILD_GAMMA_1S или
                   G4MODELS_BUILD_RADIACODE_103. Это основной способ:
                   он не путает приборы между собой.

  G4MODELS_BUILD   то же для случая, когда прибор один. Если внутри есть
                   подкаталог с именем детектора, берётся он.
                   По умолчанию <репозиторий>/build/<детектор>.
                   В репозиторий НЕ коммитится (см. .gitignore): спектров
                   сотни файлов и десятки мегабайт, вместо них — манифест
                   прогонов и готовые кривые в results/.

  G4MODELS_REF     каталог со скачанными эталонными данными ЛСРМ
                   (кривые .efr/.efa, спектры комплекта поверки).
                   Задавать НЕ обязательно: если его нет, берётся
                   закоммиченный набор detectors/<детектор>/reference/lsrm.
                   Раскладка у наборов разная, поэтому файлы ищутся
                   через efficiency_curve(), kit_dir() и find_data(),
                   а не склеиванием пути.

  SPECTRAVIBE_ROOT корень проекта spectravibe (gamma-spectrum-analysis).
                   Нужен только там, где читаются СЫРЫЕ .spe ЛСРМ: отсчёты
                   в них упакованы двоично, и берётся штатный читатель
                   gamma.io.lsrm_spe оттуда. Без этой переменной такие
                   скрипты сообщают, что данных нет, и пропускают шаг.

  GEANT4_ROOT      prebuilt-сборка Geant4 (для CMake). См. common/cmake.

Использование:
    from paths import build, ref, spectravibe
    p = build("Gamma-1S") / "grid" / "rho1.60_E0661.7.csv"
"""
import os
import re
import sys
from pathlib import Path

# Консоль Windows работает в CP-1251/CP-866, где нет ни греческих букв,
# ни знака «≈». Печать отчёта обрывалась на UnicodeEncodeError посреди
# уже посчитанного результата. Пусть лучше подставит вопросительный знак.
# Модуль путей импортируют все скрипты репозитория, поэтому место здесь.
try:
    sys.stdout.reconfigure(errors="replace")
except (AttributeError, ValueError):   # не консоль или старый python
    pass

# корень репозитория: этот файл лежит в common/py/
REPO = Path(__file__).resolve().parents[2]


def _env(name):
    v = os.environ.get(name)
    return Path(v) if v else None


def _slug(detector):
    """Gamma-1S -> GAMMA_1S: имя детектора в виде суффикса переменной."""
    return re.sub(r"[^A-Za-z0-9]+", "_", detector).upper()


def build(detector="Gamma-1S"):
    """Каталог сборки и расчётных спектров данного детектора.

    Порядок:
      1. G4MODELS_BUILD_<ДЕТЕКТОР>  — каталог именно этого прибора,
         например G4MODELS_BUILD_GAMMA_1S;
      2. G4MODELS_BUILD/<детектор>  — если такой подкаталог есть, то есть
         переменная указывает на ОБЩИЙ корень расчётов;
      3. G4MODELS_BUILD             — как есть, для работы с одним прибором;
      4. <репозиторий>/build/<детектор>.

    Пункты 1 и 2 появились не для красоты: раньше G4MODELS_BUILD возвращалась
    независимо от детектора, и при работе с двумя приборами сразу скрипты
    RadiaCode искали свои спектры в каталоге Гамма-1С и падали с
    FileNotFoundError на пути, которого никогда не существовало.
    """
    own = _env("G4MODELS_BUILD_" + _slug(detector))
    if own:
        return own
    common = _env("G4MODELS_BUILD")
    if common:
        sub = common / detector
        return sub if sub.is_dir() else common
    return REPO / "build" / detector


def ref(detector="Gamma-1S"):
    """Каталог эталонных данных.

    Порядок: переменная G4MODELS_REF → скачанные данные reference/data →
    ЗАКОММИЧЕННЫЙ набор reference/lsrm.

    Последнее звено существенно: в чистом клоне каталога reference/data нет,
    и без отката скрипты не находили бы ничего, хотя эталонный набор лежит
    рядом, в репозитории. Раскладка у него другая, поэтому обращаться к
    файлам следует не по имени, а через efficiency_curve() и kit_dir().
    """
    env = _env("G4MODELS_REF")
    if env:
        return env
    base = REPO / "detectors" / detector / "reference"
    data = base / "data"
    if data.is_dir() and any(data.iterdir()):
        return data
    lsrm = base / "lsrm"
    return lsrm if lsrm.is_dir() else data


# Названия геометрий так, как они записаны в именах файлов ЛСРМ.
GEOMETRIES = ("Маринелли", "Дента", "Петри", "Точечная-25см", "Точечная-5см")


def efficiency_curve(geometry, ext="efr", detector="Gamma-1S"):
    """Измеренная кривая эффективности ЛСРМ для геометрии, или None.

    `geometry` — как в имени файла: «Маринелли», «Дента», «Петри»,
    «Точечная-25см», «Точечная-5см». Расширение: efr (точки) или efa
    (точки с аппроксимацией).

    Копий кривой в наборе несколько (рабочее дерево прибора и выделенный
    каталог efficiency/), поэтому предпочтение отдаётся каталогу
    efficiency/ — там лежит итоговая версия, а не промежуточная.
    """
    root = ref(detector)
    if not root.is_dir():
        return None
    hits = [p for p in root.rglob("*." + ext)
            if geometry.lower() in p.name.lower()]
    if not hits:
        return None
    hits.sort(key=lambda p: (0 if "efficiency" in str(p).lower() else 1,
                             len(str(p))))
    return hits[0]


def kit_dir(geometry, fmt="xml", detector="Gamma-1S"):
    """Каталог спектров комплекта поверки для геометрии, или None.

    `geometry` — имя каталога набора: Marinelli_1L, Denta_120mL, Petri_60mL,
    Point_25cm, Point_5cm.

    Комплект лежит в двух форматах: BecqMoni XML (каталог
    reference_kits_becqmoni) и двоичный .spe ЛСРМ (reference_kits). Разбор в
    этом проекте написан под XML, поэтому по умолчанию возвращается он;
    fmt="spe" даёт исходный приборный набор.
    """
    root = ref(detector)
    if not root.is_dir():
        return None
    want = "becqmoni" if fmt == "xml" else None
    hits = [p for p in root.rglob(geometry) if p.is_dir()]
    if not hits:
        return None
    hits.sort(key=lambda p: (0 if want and want in str(p).lower()
                             else (1 if not want and "becqmoni" not in
                                   str(p).lower() else 2), len(str(p))))
    return hits[0]


def find_data(filename, detector="Gamma-1S"):
    """Найти файл эталонных данных по имени, рекурсивно, или None.

    Раскладка скачанного набора плоская, а закоммиченного — по геометриям и
    нуклидам. Скрипты знают имя файла, но не его место, поэтому ищем.
    """
    root = ref(detector)
    if not root.is_dir():
        return None
    hits = sorted(root.rglob(filename), key=lambda p: len(str(p)))
    return hits[0] if hits else None


def read_text(path):
    """Прочитать текстовый файл ЛСРМ, не гадая о кодировке.

    Приборные форматы (.efr, .efa, .cen, .lib) записаны в CP-1251, а
    загрузчик fetch_efr.py сохраняет скачанное в UTF-8. Один и тот же скрипт
    должен читать оба, иначе в чистом клоне он падает с UnicodeDecodeError на
    первом же байте кириллицы.

    Порядок важен: сначала UTF-8. Текст в UTF-8 «успешно» декодируется и как
    CP-1251, но получается мохибейк; обратное почти никогда не верно.
    """
    raw = open(str(path), "rb").read()
    for enc in ("utf-8", "cp1251"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            pass
    return raw.decode("cp1251", errors="replace")


def tools():
    """Каталог инструментов репозитория (там же читатели .efr и загрузчики)."""
    return REPO / "tools"


def results(detector="Gamma-1S"):
    """Итоговые таблицы и кривые — ЭТО коммитится."""
    return REPO / "detectors" / detector / "results"


def measured(detector="RadiaCode-103"):
    """Каталог ИЗМЕРЕННЫХ спектров прибора (не расчётных).

    Переменная `G4MODELS_MEASURED`. Это ЛИЧНЫЕ измерения оператора, в
    репозитории их нет и не будет.

    Если каталога нет, выходим с понятным сообщением прямо здесь. Иначе
    скрипт падал бы ниже с `FileNotFoundError` на длинном пути, и стороннему
    читателю казалось бы, что сломан репозиторий, тогда как данных просто
    нет и взять их неоткуда.
    """
    p = _env("G4MODELS_MEASURED") or (REPO / "measured" / detector)
    if not p.is_dir():
        raise SystemExit(
            "Нет каталога измеренных спектров %s.\n"
            "Это личные измерения оператора, в репозиторий они не входят.\n"
            "Если они у вас есть, укажите G4MODELS_MEASURED; если нет —\n"
            "этот расчёт воспроизвести нельзя, а его результат лежит готовым\n"
            "в detectors/%s/results/." % (p, detector))
    return p


def spectravibe():
    """Корень spectravibe или None, если переменная не задана."""
    return _env("SPECTRAVIBE_ROOT")


def require_spectravibe(what):
    """Понятное сообщение вместо загадочного ImportError."""
    root = spectravibe()
    if root is None or not root.exists():
        raise SystemExit(
            "Для «%s» нужны сырые .spe ЛСРМ и их штатный читатель.\n"
            "Задайте SPECTRAVIBE_ROOT — корень проекта gamma-spectrum-analysis\n"
            "(в нём scripts/gamma/io/lsrm_spe.py). Свой разбор .spe писать не\n"
            "надо: отсчёты там упакованы двоично." % what)
    return root
