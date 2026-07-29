"""Дрейф параметров обработки со времени аттестации: слепок против текущего.

ЗАЧЕМ. Весь проект сравнивает наши расчёты с аттестационными числами 2024
года. Молчаливое допущение — что обработка спектра с тех пор не менялась.
Оператор дал резервный слепок конфигурации прибора (архив рабочего каталога
`Гамма-1С_№...`), и допущение проверяется прямо: слепок против нынешнего
`lsrm.cnf`.

ЧТО НАЙДЕНО (29.07.2026): параметры, задающие ПОДЛОЖКУ И ПЛОЩАДЬ ПИКА,
изменились. Степень полинома фона 2 -> 0, учёт комптоновских данных
выкл -> вкл, пики вылета (одиночный, двойной, X-escape) учитывались ->
не учитываются, обратное рассеяние не учитывалось -> учитывается.
Обнулены и паспортные характеристики прибора, входящие в бюджет
неопределённости: интегральная нелинейность 0,3 -> 0, временная
нестабильность 0,03 -> 0.

СЛЕДСТВИЕ. Расхождение площади 2614,5 между аттестацией (.efr, 59668) и
нашими повторными сеансами (54000, разброс 1 %) больше нельзя относить к
«вариативности фита»: это в первую очередь РАЗНЫЕ ПАРАМЕТРЫ ОБРАБОТКИ.
Чтобы воспроизвести аттестацию, надо вернуть параметры слепка.

Пути в выводе не печатаются: они содержат личные имена каталогов.
Расхождения путей помечаются как «путь изменён».

    G1S_LSRM_BACKUP_ZIP=<архив слепка> G1S_LSRM_CNF=<текущий lsrm.cnf> \
        python detectors/Gamma-1S/analysis/config_drift.py
"""
import io
import os
import sys
import zipfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402

# Ключи, влияющие на ПЛОЩАДЬ пика и на бюджет неопределённости, с пояснением
MEANING = {
    "Fitting/BackgrPower": "степень полинома фона под пиком",
    "Fitting/UseComptonData": "учёт комптоновских данных в подложке",
    "Fitting/Step": "ступенька в модели пика",
    "Fitting/StepFromPattern": "ступенька из пика-образа",
    "Detector/UseSingleEscape": "пик одиночного вылета",
    "Detector/UseDoubleEscape": "пик двойного вылета",
    "Detector/XEscCorrection": "поправка на вылет характеристического кванта",
    "Processing/ConsiderBackScattering": "учёт обратного рассеяния",
    "Processing/ConsiderApparatusPeaks": "учёт аппаратурных пиков",
    "Processing/LeftPeakBound": "левая граница пика в долях ПШПВ",
    "Processing/RightPeakBound": "правая граница пика в долях ПШПВ",
    "Processing/Sensitivity": "чувствительность поиска пиков",
    "Processing/MaxZoneLength": "максимальная длина зоны в долях ПШПВ",
    "Processing/ChBegin": "начало диапазона обработки (канал)",
    "Processing/ChEnd": "конец диапазона обработки (канал)",
    "Representation/SyseffErr": "систематическая неопределённость eps в %",
    "Representation/ConsiderChainDecay": "коррекция по цепочке",
    "Representation/ConsiderMeasTime": "распад за время измерения",
    "Representation/UseChiSquare": "хи-квадрат в погрешности активности",
    "Representation/CovellPeakArea": "метод Ковелла для площади",
    "Representation/ErrorUpperLevel(%)": "порог вывода верхним пределом в %",
    "Spectrometer/IntegralNonlinearity": "интегральная нелинейность (бюджет)",
    "Spectrometer/DifferentialNonlinearity": "дифф. нелинейность (бюджет)",
    "Spectrometer/TimeInstability": "временная нестабильность (бюджет)",
    "Spectrometer/SummCorrection": "поправка на суммирование",
    "Spectrometer/HighLoadCorrection": "поправка на высокую загрузку",
    "Default/ActivityCalcMethod": "метод расчёта активности",
}


def parse_cnf(text):
    out, sec = {}, ""
    for ln in text.replace("\r\n", "\n").split("\n"):
        ln = ln.strip()
        if ln.startswith("[") and ln.endswith("]"):
            sec = ln[1:-1]
            continue
        if "=" in ln:
            k, v = ln.split("=", 1)
            out["%s/%s" % (sec, k)] = v
    return out


def is_path(v):
    return "\\" in v or ":" in v


if __name__ == "__main__":
    zp = os.environ.get("G1S_LSRM_BACKUP_ZIP")
    cp = os.environ.get("G1S_LSRM_CNF")
    if not (zp and os.path.exists(zp) and cp and os.path.exists(cp)):
        raise SystemExit(
            "Задайте G1S_LSRM_BACKUP_ZIP (архив резервного слепка"
            " конфигурации прибора)\nи G1S_LSRM_CNF (текущий lsrm.cnf)."
            " В репозитории их нет: личные пути и номер прибора.")

    with zipfile.ZipFile(zp) as z:
        name = [n for n in z.namelist() if n.endswith("lsrm.cnf")][0]
        orig = parse_cnf(z.read(name).decode("cp1251"))
    curr = parse_cnf(open(cp, "rb").read().decode("cp1251"))

    rows = []
    for k in sorted(set(orig) | set(curr)):
        a, b = orig.get(k, ""), curr.get(k, "")
        if a == b:
            continue
        if is_path(a) or is_path(b):
            a = b = "(путь изменён)"
        rows.append((k, a, b, MEANING.get(k, "")))

    key_rows = [r for r in rows if r[3]]
    print("Дрейф конфигурации со времени аттестации: расхождений %d,"
          " из них значимых %d\n" % (len(rows), len(key_rows)))
    print("%-42s %-12s %-12s %s"
          % ("ключ", "слепок", "сейчас", "на что влияет"))
    for k, a, b, m in key_rows:
        print("%-42s %-12s %-12s %s" % (k, a[:12], b[:12], m))

    other = [r for r in rows if not r[3]]
    if other:
        print("\nостальные расхождения (имя профиля, даты, плагины, пути):")
        for k, a, b, _m in other:
            print("   %-40s %s -> %s" % (k, a[:20], b[:20]))

    out = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results",
        "config_drift.csv"))
    # csvio запрещает запятую внутри значения: потребитель читает split(",")
    def sane(v):
        return str(v).replace(",", ";")

    csvio.write(
        out, ["key", "snapshot", "current", "affects"],
        [(sane(k), sane(a), sane(b), sane(m)) for k, a, b, m in rows],
        comments=[
            "Сравнение резервного слепка конфигурации прибора (эпоха"
            " аттестации) с текущим lsrm.cnf.",
            "Значения путей заменены на «(путь изменён)»: содержат личные"
            " имена каталогов.",
            "Строки с заполненным affects меняют ПЛОЩАДЬ пика или бюджет"
            " неопределённости — то есть сопоставимость наших чисел с",
            "  аттестационными. Прочие — имя профиля, даты, плагины.",
        ])
    print("\nтаблица: %s" % out)
