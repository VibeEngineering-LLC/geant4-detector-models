# #DET-1 (оператор 09.10): карточка детектора — единственный источник ширины/формы пика для всех подгонок GS2020.
# Читают: fit_gs2020_th232_m1 (ширина в свёртку), make_detector_card (генератор). Путь — GS_DET_CARD, иначе репозиторий.
import os, json, hashlib
_REL = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dataset", "detector_card.json"))   # из репозитория: analysis/ → dataset/
CARD_PATH = os.environ.get("GS_DET_CARD") or _REL   # в репозитории карточка лежит рядом: analysis/ → dataset/


def effective_points(card):
    """[(E, ПШПВ·множитель)] по возрастанию E — ровно то, что идёт в fwhm_points_gs2020.csv."""
    sc = card["fwhm"]["scale"]
    return [(float(E), float(w) * float(sc.get(str(E), 1.0))) for E, w in sorted(card["fwhm"]["points_sl"], key=lambda p: float(p[0]))]


def points_sha(points):
    return hashlib.sha256("\n".join("%r,%r" % (E, w) for E, w in points).encode()).hexdigest()


def load_card(path=CARD_PATH):
    if not os.path.exists(path):
        raise SystemExit("ОТКАЗ #DET-1: нет карточки детектора " + path + " — python scripts/make_detector_card.py")
    card = json.load(open(path, encoding="utf-8"))
    if points_sha(effective_points(card)) != card["fwhm"]["sha256"]:
        raise SystemExit("ОТКАЗ #DET-1: карточка " + path + " изменена — sha256 точек ширины не совпадает")
    return card
