"""Обезличивание спектров и паспортов перед публикацией.

Убирает из ДАННЫХ фамилии людей и серийные номера приборов и источников —
и в содержимом файлов, и в их именах.

    python tools/anonymize.py <каталог данных> [--verify] [--dry-run]

Каталог обязателен и указывается явно. Значения по умолчанию нет намеренно:
см. «Три правила» ниже.

--------------------------------------------------------------------------
КАРТА СООТВЕТСТВИЙ ЛЕЖИТ ВНЕ РЕПОЗИТОРИЯ
--------------------------------------------------------------------------
Карта «подлинная фамилия → псевдоним» и «подлинный номер → псевдоним» — это
ключ обратного опознания. Опубликованная вместе с данными, она полностью
отменяет обезличивание: по ней любой читатель восстановит и фамилии, и номера.
Поэтому в этом файле её НЕТ и быть не должно.

Путь к карте задаётся переменной окружения ANON_MAP (или ключом --map=). Файл
JSON в UTF-8:

    {
      "names":   {"Фамилия": "Опер-01", "Другая": "Проба-01"},
      "serials": {"0000-00": "SN-01", "1111-11": null},
      "replace": {"12345": "SRC-16"},
      "remove":  ["_12345"]
    }

  names   — фамилии; инициалы после фамилии убираются вместе с ней;
  serials — числовые ядра номеров без знака «№»; значение null означает
            «убрать вместе со знаком номера»;
  replace — буквальные замены: номер, встречающийся и без знака «№», и не
            подходящий под узкий шаблон заводского номера;
  remove  — строки, вырезаемые как есть.

Держите карту рядом с рабочим каталогом расчётов, а не в репозитории.

--------------------------------------------------------------------------
ТРИ ПРАВИЛА, КАЖДОЕ КУПЛЕНО СЛОМАННЫМИ ФАЙЛАМИ
--------------------------------------------------------------------------

1. ОБРАБАТЫВАТЬ ТОЛЬКО ДАННЫЕ, НЕ ИСХОДНИКИ И НЕ ДОКУМЕНТАЦИЮ.
   Прогон по корню репозитория испортил: заголовки «## 10.» в markdown
   превратились в «##SRC-A01.», CSS-цвета «#667» — в «#SRC-A20», комментарии
   «# 511 кэВ» — в «#SRC-A03». Причина: знак «#» в шаблоне номера, а в
   markdown это заголовок, в CSS — цвет, в питоне — комментарий.
   Поэтому знак «#» из шаблона убран совсем (в приборных форматах номер
   пишется знаком «№»), а расширения исходников и разметки пропускаются
   (SKIP_EXT).

2. НЕ ЗАПУСКАТЬ ПО КАТАЛОГУ, ГДЕ ЛЕЖИТ САМ ИНСТРУМЕНТ.
   Прошлая версия держала карту фамилий внутри себя, обход дошёл до
   собственного файла и заменил ключи карты на псевдонимы: инструмент
   обезличил сам себя и перестал работать. Теперь карта снаружи, а обход
   каталога с инструментом запрещён проверкой.

3. В ДВОИЧНОМ .spe ПРАВИТЬ ТОЛЬКО ЗАГОЛОВОК.
   Формат (по штатному читателю gamma.io.lsrm_spe): текстовый заголовок
   CP-1251 из строк KEY=VALUE, завершающийся меткой «SPECTR=», сразу за
   которой идёт блок отсчётов — беззнаковые 32-битные целые. Таблицы
   смещений нет, начало данных ищется по метке, поэтому изменение длины
   заголовка разбор не ломает — а вот замена по всему файлу незаметно
   испортила бы отсчёты: в двоичном блоке встречается и байт «№», и любая
   цифровая последовательность. Файл режется по метке, правится префикс,
   хвост дописывается как есть.

   Первая версия таких файлов не касалась вообще: она декодировала файл
   целиком, для .spe это не удаётся никогда, и функция возвращала «двоичный
   файл — не трогаем». Молча пропускались ровно те файлы, ради которых всё
   затевалось: 545 из 658 с полем OPERATOR=.

Целость проверяется штатным читателем: сумма отсчётов до и после правки
обязана совпасть (--verify, нужен SPECTRAVIBE_ROOT). Читатель, писатель и
конвертор форматов ЛСРМ живут в spectravibe-toolkit — свой разбор двоичного
.spe писать не надо.

ЧТО СОЗНАТЕЛЬНО ОСТАВЛЕНО: обозначения вида «420-7-15» — это не заводской
номер экземпляра, а каталожное обозначение набора ОИСН, по которому источник
опознаётся в методиках и в кривых эффективности.
"""
import json
import os
import re
import sys

# Знаки номера. «#» в этих данных используется наравне с «№» («Th-232 #420»,
# «Ba-133 #6649» в именах спектров поверки), поэтому поддерживать его надо.
# Опасен он только вне данных: в markdown это заголовок, в CSS — цвет, в
# питоне — комментарий. Защита от этого — не сам шаблон, а SKIP_EXT и запрет
# запуска по каталогу с исходниками (правила 1 и 2).
TOKEN = re.compile(r"([№#])\s*(0?[\d]{2,4}(?:[.\-]\d{2,4})?)")

# Заводской номер прибора без знака «№». Шаблон узкий — четыре цифры с
# ведущим нулём, дефис, две цифры — и замена делается ТОЛЬКО если такой номер
# перечислен в карте. Широкий шаблон испортил бы числовые данные: в спектрах и
# кривых полно последовательностей вида «0609-32».
BARE = re.compile(r"(?<![\d.])(0\d{3}-\d{2})(?![\d.])")

# Инициалы после фамилии («Фамилия Ю.Д.») — тоже персональные данные.
INITIALS = r"(?:\s*[А-ЯЁ]\.\s*[А-ЯЁ]?\.?)?"

# Исходники, разметка, документация, сборка: в данных их не бывает, а
# испортить их замена может — см. правило 1.
SKIP_EXT = (
    ".py", ".pyc", ".md", ".html", ".htm", ".css", ".js",
    ".cc", ".hh", ".cpp", ".h", ".cmake", ".mac", ".yml", ".yaml", ".toml",
    ".png", ".jpg", ".jpeg", ".gif", ".svg", ".pdf",
    ".zip", ".rar", ".7z", ".xlsx", ".docx", ".mdb",
)

MARK = b"SPECTR="


def load_map(path):
    with open(path, encoding="utf-8") as fh:
        m = json.load(fh)
    return (m.get("names") or {}, m.get("serials") or {},
            m.get("replace") or {}, m.get("remove") or [],
            m.get("replace_re") or {}, m.get("auto") or {})


class Scrubber:
    """Замены по карте. Автопсевдонимы — только для номеров со знаком «№».

    В паспортах источников номеров сертификатов десятки, перечислять их руками
    бессмысленно и легко пропустить, поэтому неизвестный номер получает
    следующий свободный псевдоним SRC-Ann. Соответствие держится только в
    памяти и НЕ публикуется.
    """

    def __init__(self, names, serials, replace, remove, replace_re=None,
                 auto=None):
        self.names = names
        self.serials = serials
        self.replace = replace
        self.remove = remove
        # Регулярные замены с якорем по нуклиду: номер без знака «№» голой
        # подстрокой менять нельзя, он совпадёт с отсчётами в спектрах.
        self.replace_re = [(re.compile(k), v) for k, v in
                           (replace_re or {}).items()]
        # Автопсевдонимы, уже присвоенные прошлыми прогонами. Без этого номер
        # получал НОВЫЙ SRC-Ann при каждом запуске: порядок обхода каталога
        # другой — нумерация другая, и один источник в двух прогонах выглядел
        # как два разных. Новые дописываются в карту ключом "auto".
        self.auto = dict(auto or {})
        self.auto_added = {}

    def _pseudonym(self, key):
        if key in self.serials:
            return self.serials[key]
        if key not in self.auto:
            n = 1 + max([0] + [int(v.rsplit("-A", 1)[1])
                               for v in self.auto.values()
                               if re.fullmatch(r"SRC-A\d+", v)])
            self.auto[key] = "SRC-A%02d" % n
            self.auto_added[key] = self.auto[key]
        return self.auto[key]

    def _sub_token(self, m):
        sign, num = m.group(1), m.group(2)
        key = num.lstrip("0") if num.lstrip("0") in self.serials else num
        rep = self._pseudonym(key)
        return "" if rep is None else sign + rep

    def _sub_bare(self, m):
        """Номер без знака «№»: заменяем только заведомо известные."""
        num = m.group(1)
        if num in self.serials:
            rep = self.serials[num]
            return rep if rep is not None else ""
        return num

    def line(self, ln):
        # Якорные замены — ПЕРВЫМИ, пока номер ещё в исходном виде. Повторно
        # их результат не тронется: TOKEN требует после знака «№» цифру, а там
        # уже «SRC-21».
        for rx, rep in self.replace_re:
            ln = rx.sub(rep, ln)
        ln = TOKEN.sub(self._sub_token, ln)
        ln = BARE.sub(self._sub_bare, ln)
        for s in self.remove:
            ln = ln.replace(s, "")
        for src, dst in self.replace.items():
            ln = ln.replace(src, dst)
        for who, alias in self.names.items():
            ln = re.sub(re.escape(who) + INITIALS, alias, ln)
        return ln

    def text(self, txt):
        out = []
        for ln in txt.split("\n"):
            # поле оператора вычищаем целиком, вместе со значением: фамилия в
            # нём может быть и не перечисленной в карте
            m = re.match(r"(OPERATOR\s*=)", ln)
            if m:
                out.append(m.group(1) + ("\r" if ln.endswith("\r") else ""))
                continue
            out.append(self.line(ln))
        return "\n".join(out)

    def name(self, fname):
        return re.sub(r"\s{2,}", " ", self.line(fname)).strip()


def split_header(raw):
    """(заголовок, хвост) для двоичного .spe; хвост правке не подлежит.

    Метка ищется в первых 64 КБ: в блоке отсчётов последовательность
    «SPECTR=» встретиться теоретически может, и брать последнее вхождение
    нельзя.
    """
    i = raw.find(MARK, 0, 65536)
    if i < 0:
        return None, None
    end = i + len(MARK)
    return raw[:end], raw[end:]


def process(path, sc, dry=False):
    """-> ('текст'|'двоичный'|None, изменён?)"""
    if path.lower().endswith(SKIP_EXT):
        return None, False
    raw = open(path, "rb").read()

    # ВАЖЕН ПОРЯДОК: сначала UTF-8, потом CP-1251. Файл в UTF-8 с кириллицей
    # «успешно» декодируется и как CP-1251 — но получается мохибейк, и знак
    # номера в нём уже не распознаётся. Так были пропущены серийники в
    # secondary_peaks_v2.json. Обратный порядок безопасен: кириллица в CP-1251
    # почти никогда не является корректным UTF-8.
    for enc in ("utf-8", "cp1251"):
        try:
            txt = raw.decode(enc)
        except UnicodeDecodeError:
            continue
        new = sc.text(txt)
        if new == txt:
            return "текст", False
        if not dry:
            open(path, "wb").write(new.encode(enc, "replace"))
        return "текст", True

    # Не декодируется целиком — это .spe: заголовок CP-1251, дальше отсчёты.
    head, tail = split_header(raw)
    if head is None:
        return None, False                # картинка, архив, база — не наше
    txt = head.decode("cp1251", errors="replace")
    new = sc.text(txt)
    if new == txt:
        return "двоичный", False
    if not dry:
        open(path, "wb").write(new.encode("cp1251", "replace") + tail)
    return "двоичный", True


def counts_sum(path, scripts_dir):
    """Сумма отсчётов по ШТАТНОМУ читателю ЛСРМ (не своему разбору)."""
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    from gamma.io.lsrm_spe import read_lsrm_spe
    import numpy as np
    return int(np.asarray(read_lsrm_spe(path).counts).sum())


def walk(root):
    for r, _, fs in os.walk(root):
        if ".git" in r.split(os.sep):
            continue
        for f in fs:
            yield os.path.join(r, f)


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    flags = [a for a in argv if a.startswith("--")]
    dry = "--dry-run" in flags
    mp = next((a.split("=", 1)[1] for a in flags if a.startswith("--map=")),
              os.environ.get("ANON_MAP"))

    if not args:
        print("Укажите каталог ДАННЫХ:\n"
              "    python tools/anonymize.py detectors/Gamma-1S/reference "
              "[--verify] [--dry-run]\n"
              "Карта соответствий — в ANON_MAP или --map=, вне репозитория.")
        return 2
    root = os.path.abspath(args[0])
    if not os.path.isdir(root):
        print("нет такого каталога:", root)
        return 2

    # Правило 2: не обрабатывать каталог, где лежит сам инструмент.
    here = os.path.dirname(os.path.abspath(__file__))
    if os.path.commonpath([root, here]) == root:
        print("Отказ: каталог %s содержит сам инструмент (%s).\n"
              "Прошлый такой прогон заменил ключи карты внутри anonymize.py на\n"
              "псевдонимы — инструмент обезличил сам себя. Укажите каталог\n"
              "ДАННЫХ, например detectors/Gamma-1S/reference." % (root, here))
        return 2

    if not mp or not os.path.isfile(mp):
        print("Не задана карта соответствий. Укажите ANON_MAP=<файл.json> или\n"
              "--map=<файл.json>. Карта НЕ должна лежать в репозитории: это\n"
              "ключ обратного опознания. Формат — в начале этого файла.")
        return 2
    sc = Scrubber(*load_map(mp))

    sv = os.environ.get("SPECTRAVIBE_ROOT")
    scripts = os.path.join(sv, "scripts") if sv else None
    verify = "--verify" in flags
    if verify and not (scripts and os.path.isdir(scripts)):
        print("SPECTRAVIBE_ROOT не задан или в нём нет scripts/ — проверить\n"
              "целость .spe штатным читателем нечем, останавливаюсь.")
        return 2

    # 1. суммы отсчётов до правки
    before = {}
    if verify:
        for p in walk(root):
            if p.lower().endswith(".spe"):
                try:
                    before[p] = counts_sum(p, scripts)
                except Exception as exc:                      # noqa: BLE001
                    # Обычно это ТЕКСТОВЫЙ вариант .spe (в комплекте их три):
                    # двоичный читатель не находит в нём метку SPECTR=. Такой
                    # файл правится как текст и по отсчётам НЕ проверяется —
                    # дыра в контроле, о которой надо знать. Если правка их
                    # когда-нибудь заденет, разбирать их надо lsrm_spe_text.py.
                    print("не читается штатным читателем: %s (%s)" % (p, exc))

    # 2. содержимое
    stat = {"текст": 0, "двоичный": 0}
    for p in walk(root):
        kind, changed = process(p, sc, dry)
        if changed:
            stat[kind] += 1
    print("изменено: текстовых %d, двоичных .spe %d%s"
          % (stat["текст"], stat["двоичный"], " (пробный прогон)" if dry else ""))

    # 3. проверка целости: сумма отсчётов обязана совпасть
    if verify and not dry:
        bad = 0
        for p, s0 in before.items():
            try:
                s1 = counts_sum(p, scripts)
            except Exception as exc:                          # noqa: BLE001
                print("ПОСЛЕ ПРАВКИ НЕ ЧИТАЕТСЯ: %s (%s)" % (p, exc))
                bad += 1
                continue
            if s1 != s0:
                print("СУММА ОТСЧЁТОВ ИЗМЕНИЛАСЬ: %s  %d -> %d" % (p, s0, s1))
                bad += 1
        print("проверено .spe: %d, расхождений: %d" % (len(before), bad))
        if bad:
            return 1

    # 4. имена файлов и каталогов (снизу вверх, чтобы не рвать пути)
    ren = 0
    for r, dirs, fs in os.walk(root, topdown=False):
        if ".git" in r.split(os.sep):
            continue
        for f in fs + dirs:
            n = sc.name(f)
            if n != f:
                ren += 1
                if not dry:
                    os.replace(os.path.join(r, f), os.path.join(r, n))
    print("переименовано:", ren)

    # 5. новые автопсевдонимы — в карту, иначе следующий прогон присвоит тому
    # же номеру другой SRC-Ann и один источник разъедется на два.
    if sc.auto_added and not dry:
        with open(mp, encoding="utf-8") as fh:
            m = json.load(fh)
        m.setdefault("auto", {}).update(sc.auto_added)
        with open(mp, "w", encoding="utf-8") as fh:
            json.dump(m, fh, ensure_ascii=False, indent=2)
        print("дописано в карту автопсевдонимов: %d" % len(sc.auto_added))
    if sc.auto:
        print("автопсевдонимов всего: %d (соответствие не публикуется)"
              % len(sc.auto))

    # 6. Контрольный прогон проверки по обработанному каталогу. Раньше
    # --verify означал только «отсчёты в .spe целы», и обезличивание могло
    # отчитаться успехом, оставив номера в .src и путях: два разных вопроса
    # проверялись одной галочкой.
    if not dry:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import check_paths
        hits, opaque = check_paths.scan(root)
        seen = sorted(set(hits))
        for rel, what, s in seen:
            print("ОСТАЛОСЬ  %-52s %-24s %s" % (rel[:52], what, s))
        if opaque:
            print("непрозрачных файлов: %d (содержимое проверить нечем)"
                  % len(set(opaque)))
        print("проверка после обезличивания: находок %d" % len(seen))
        if seen:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
