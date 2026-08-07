/* AtomSpectra Nano 16 PRO — функция отклика.
   Отрисовка холстов, всплывающие пояснения, управление каналами.
   Числа в прозе подставляются НЕ здесь, а при сборке страницы: сборка
   падает, если подстановка не разрешилась. Здесь — только то, что зависит
   от действий: подписи легенды, показания курсора, сводная таблица. */
(function () {
'use strict';

var DATA = window.ASN16;
var $ = function (s) { return document.querySelector(s); };

/* ── язык ────────────────────────────────────────────────────────────────
   Русский — по умолчанию; переключается кнопкой в тулбаре. Хранится в
   localStorage, поэтому язык переживает перезагрузку. */
var lang = 'ru';
try { var st = localStorage.getItem('asn16_lang'); if (st==='ru'||st==='en') lang = st; }
catch (e) {}
// URL-параметр приоритетнее localStorage: «…/atomspectra-nano16/?lang=en»
// или «…?en» открывает страницу сразу на английском, ссылка шарится как
// English. Поддерживаем оба короткой и полной формы.
try {
  var q = new URLSearchParams(location.search);
  var qLang = q.get('lang');
  if (!qLang && q.has('en')) qLang = 'en';
  if (!qLang && q.has('ru')) qLang = 'ru';
  if (qLang === 'en' || qLang === 'ru') lang = qLang;
} catch (e) {}

// pick — универсальный доступ к паре {ru,en} из data.js; строку возвращает
// как есть (для меток, которые ещё не переведены).
function pick(v) {
  if (v == null) return '';
  return typeof v === 'string' ? v : (v[lang] || v.ru || v.en || '');
}
function tr(ru, en) { return lang === 'en' ? en : ru; }

/* ── формат чисел ────────────────────────────────────────────────────────
   ГОСТ 8.417: десятичная запятая в РУ, точка в EN. Группировка тысяч
   только у счёта событий: на шкале энергии она ломает чтение
   («2 615» вместо 2614,5). */
function num(x, d) {
  var s = x.toFixed(d === undefined ? 1 : d);
  return lang === 'en' ? s : s.replace('.', ',');
}
function cnt(x) { return x.toLocaleString(lang === 'en' ? 'en-US' : 'ru-RU'); }
// Доля канала: одного знака хватает почти везде, но самые редкие каналы при
// одном знаке округляются в ноль — а это не «канала нет», а «канал редок».
function pctf(x) { return num(x, x < 1 ? 2 : 1) + ' %'; }
function keV(x) { return num(x, 1).replace(/[.,]0$/, ''); }
function keVunit() { return tr('кэВ', 'keV'); }
function plural(n, one, few, many) {
  var a = Math.abs(n) % 100, b = a % 10;
  return n + ' ' + ((a > 10 && a < 20) || b > 4 || b === 0 ? many
    : (b === 1 ? one : few));
}
// Показатель степени целиком надстрочный, ВКЛЮЧАЯ минус: обычный «−» рядом с
// надстрочной цифрой даёт «10−⁶» вместо «10⁻⁶».
function sup(n) {
  return String(n).replace(/-/g, '⁻')
    .replace(/[0-9]/g, function (d) { return '⁰¹²³⁴⁵⁶⁷⁸⁹'[+d]; });
}
function sci(x) {
  if (x === 0) return '0';
  var e = Math.floor(Math.log10(Math.abs(x)));
  var m = x / Math.pow(10, e);
  if (m >= 9.995) { m /= 10; e += 1; }      // иначе печаталось бы «10,00·10ⁿ»
  return num(m, 2) + '·10' + sup(e);
}
function css(n) {
  return getComputedStyle(document.documentElement).getPropertyValue(n).trim();
}
// Цвет канала берётся из таблицы стилей: светлая и тёмная тема держат разные
// значения, иначе заливка не даёт нужного контраста к фону сразу в обеих.
function chColor(key) { return css('--ch-' + key) || '#888888'; }
function rgba(col, a) {
  var h = col.replace('#', '');
  if (h.length === 3) h = h[0] + h[0] + h[1] + h[1] + h[2] + h[2];
  return 'rgba(' + parseInt(h.slice(0, 2), 16) + ',' + parseInt(h.slice(2, 4), 16)
    + ',' + parseInt(h.slice(4, 6), 16) + ',' + a + ')';
}
function labelOf(key) { return pick(DATA.labels[key]) || key; }
function tabAt(e0) {
  return DATA.tabs.filter(function (t) { return t.e0 === e0; })[0];
}
function chanAt(e0, key) {
  var t = tabAt(e0);
  return t && t.channels.filter(function (c) { return c.key === key; })[0];
}

/* ── монотонная кубическая интерполяция (Фрич — Карлсон) ─────────────────
   Кривая проходит через КАЖДЫЙ отсчёт и по построению не выходит за их
   пределы между узлами. Обычный кубический сплайн этого не гарантирует и
   дорисовал бы на спектре пики и провалы, которых в данных нет.
   Интерполяция ведётся в той шкале, в какой рисуется: для логарифмической
   оси — по логарифмам, иначе низкие декады выродились бы в прямые.        */
function spline(xs, ys) {
  var n = xs.length, i;
  if (n < 2) return function () { return ys[0] || 0; };
  var dx = [], m = [], c = [];
  for (i = 0; i < n - 1; i++) {
    dx[i] = xs[i + 1] - xs[i];
    m[i] = (ys[i + 1] - ys[i]) / dx[i];
  }
  c[0] = m[0]; c[n - 1] = m[n - 2];
  for (i = 1; i < n - 1; i++) {
    if (m[i - 1] * m[i] <= 0) c[i] = 0;               // экстремум — полка
    else {
      var w1 = 2 * dx[i] + dx[i - 1], w2 = dx[i] + 2 * dx[i - 1];
      c[i] = (w1 + w2) / (w1 / m[i - 1] + w2 / m[i]);
    }
  }
  return function (x) {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    var lo = 0, hi = n - 1, mid;
    while (hi - lo > 1) { mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
    var h = dx[lo], t = (x - xs[lo]) / h, t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * ys[lo] + (t3 - 2 * t2 + t) * h * c[lo]
         + (-2 * t3 + 3 * t2) * ys[lo + 1] + (t3 - t2) * h * c[lo + 1];
  };
}

/* ── общий подгон холста под контейнер ───────────────────────────────── */
function fit(c, h) {
  var dpr = window.devicePixelRatio || 1;
  // Размеры НЕ прибиваются пикселями в стиль: прибитая ширина пережила бы
  // сужение окна и вылезла за страницу. Раскладку делает CSS, сюда берётся
  // фактический размер после неё — тогда внутренняя система координат
  // холста совпадает с экранной ВСЕГДА, включая совсем узкие окна.
  // Высота задаётся только там, где она фиксированная (h задан); у главного
  // графика её определяет высота панели.
  if (h) c.style.height = h + 'px';
  var r = c.getBoundingClientRect();
  var w = Math.max(2, Math.round(r.width));
  var hh = Math.max(2, Math.round(r.height));
  c.width = Math.round(w * dpr);
  c.height = Math.round(hh * dpr);
  var g = c.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  return {g: g, W: w, H: hh};
}

/* ── состояние ───────────────────────────────────────────────────────── */
var tabIdx = 0;
var mode = 'log';
var shown = DATA.tabs.map(function (t) {
  var s = {total: true};
  t.channels.forEach(function (c) { s[c.key] = true; });
  return s;
});

/* ── всплывающие пояснения ───────────────────────────────────────────── */
var pop = $('#pop'), scrim = $('#scrim'), popOpener = null;
function closePop() {
  pop.hidden = true; scrim.hidden = true;
  document.querySelectorAll('[aria-expanded="true"]').forEach(function (b) {
    b.setAttribute('aria-expanded', 'false');
  });
  if (popOpener && popOpener.focus) popOpener.focus();
  popOpener = null;
}
function openPop(node, opener, wide) {
  pop.textContent = '';
  var close = document.createElement('button');
  close.type = 'button'; close.className = 'close';
  close.setAttribute('aria-label', tr('закрыть', 'close'));
  close.textContent = '✕';
  close.onclick = closePop;
  pop.appendChild(close);
  pop.appendChild(node);
  pop.classList.toggle('wide', !!wide);
  pop.style.left = '50%'; pop.style.top = '50%';
  pop.style.transform = 'translate(-50%, -50%)';
  pop.hidden = false; scrim.hidden = false;
  popOpener = opener || null;
  if (opener) opener.setAttribute('aria-expanded', 'true');
  // Применить язык к содержимому popup'а — шаблоны клонированы «сырыми»,
  // атрибуты data-en у них есть, но в isolatiion они ещё не отработали.
  if (typeof applyLangTo === 'function') applyLangTo(pop);
  close.focus();
}
scrim.addEventListener('click', closePop);
document.addEventListener('keydown', function (e) {
  if (e.key === 'Escape' && !pop.hidden) closePop();
});

/* ── переключение языка ──────────────────────────────────────────────────
   Каждый статичный текст в шаблоне размечен парой атрибутов
   `data-ru`/`data-en` (для innerHTML), либо `data-i18n-ru`/`data-i18n-en`
   (для textContent), либо `data-en-title`/`data-en-aria` (для атрибутов).
   applyLangTo(root) прогоняет ВСЕ такие узлы в корне; вызывается на
   документе и на каждом клонированном <template> перед показом. */
function ruTextFor(el, key) {
  var stored = el.getAttribute(key);
  if (stored != null) return stored;
  // Первый переход RU→EN — RU-версии ещё нет; берём текущее содержимое.
  return null;
}
function applyLangTo(root) {
  var docLang = lang;
  root.querySelectorAll('[data-en]').forEach(function (el) {
    if (!el.hasAttribute('data-ru')) el.setAttribute('data-ru', el.innerHTML);
    el.innerHTML = docLang === 'en' ? el.getAttribute('data-en')
                                    : el.getAttribute('data-ru');
  });
  root.querySelectorAll('[data-i18n-en]').forEach(function (el) {
    if (!el.hasAttribute('data-i18n-ru')) el.setAttribute('data-i18n-ru', el.textContent);
    el.textContent = docLang === 'en' ? el.getAttribute('data-i18n-en')
                                      : el.getAttribute('data-i18n-ru');
  });
  root.querySelectorAll('[data-en-aria]').forEach(function (el) {
    if (!el.hasAttribute('data-ru-aria'))
      el.setAttribute('data-ru-aria', el.getAttribute('aria-label') || '');
    el.setAttribute('aria-label',
      docLang === 'en' ? el.getAttribute('data-en-aria')
                       : el.getAttribute('data-ru-aria'));
  });
}
function applyLang(newLang) {
  lang = newLang;
  try { localStorage.setItem('asn16_lang', lang); } catch (e) {}
  // Обновляем URL, чтобы ссылка сразу открывалась на нужном языке.
  // RU — параметр удаляется (это дефолт), EN — ставится «?lang=en».
  try {
    var u = new URL(location.href);
    u.searchParams.delete('en'); u.searchParams.delete('ru');
    if (lang === 'en') u.searchParams.set('lang', 'en');
    else u.searchParams.delete('lang');
    history.replaceState(null, '', u.pathname
      + (u.search || '') + (u.hash || ''));
  } catch (e) {}
  document.documentElement.lang = lang;
  applyLangTo(document);
  var docTitle = document.querySelector('title');
  if (docTitle && docTitle.hasAttribute('data-en')) {
    if (!docTitle.hasAttribute('data-ru'))
      docTitle.setAttribute('data-ru', docTitle.textContent);
    docTitle.textContent = lang === 'en' ? docTitle.getAttribute('data-en')
                                         : docTitle.getAttribute('data-ru');
  }
  var b = $('#btnLang');
  if (b) {
    b.textContent = lang === 'en' ? 'RU' : 'EN';
    b.setAttribute('aria-label', lang === 'en' ? 'switch to Russian' : 'на русский');
  }
  // Легенда, hint и подписи элементов спектра лежат в DOM с прежней
  // языковой версией — перестроим их вместе с холстами. Проверка на
  // определённость нужна, потому что applyLang вызывается ещё до
  // объявления buildLegend/redraw.
  if (typeof buildLegend === 'function' && tabIdx < DATA.tabs.length) buildLegend();
  if (typeof redraw === 'function') redraw();
}
document.addEventListener('DOMContentLoaded', function () {
  var b = $('#btnLang');
  if (b) b.addEventListener('click', function () {
    applyLang(lang === 'en' ? 'ru' : 'en');
  });
  // Первичное применение — если из localStorage взят en, поставим его.
  applyLang(lang);
});
function tplNode(sel) {
  var t = document.querySelector(sel);
  return t ? t.content.cloneNode(true) : null;
}
document.querySelectorAll('[data-pop]').forEach(function (b) {
  b.addEventListener('click', function () {
    var id = b.getAttribute('data-pop');
    var node = tplNode('#' + id);
    if (!node) return;
    openPop(node, b, id === 't-numbers');
    if (id === 't-numbers') buildTable();
  });
});
function openChannel(key, opener) {
  var node = tplNode('[data-ch="' + key + '"]');
  if (!node) return;
  var head = document.createElement('div');
  var h2 = document.createElement('h2');
  h2.id = 'popTitle'; h2.textContent = labelOf(key);
  var em = document.createElement('span');
  em.className = 'key'; em.textContent = key;
  var tb = document.createElement('table');
  tb.className = 'mini';
  tb.style.setProperty('--c', chColor(key));
  var rows = '';
  DATA.tabs.forEach(function (t) {
    var f = chanAt(t.e0, key);
    rows += '<tr><th scope="row">' + num(t.e0, 0) + ' ' + keVunit() + '</th>'
      + (f ? '<td class="bar-cell" style="--w:' + f.pct.toFixed(2) + '">'
           + pctf(f.pct) + '</td>' : '<td>—</td>')
      + '<td>' + (f ? cnt(f.n) : '—') + '</td></tr>';
  });
  tb.innerHTML = '<thead><tr><th scope="col">' + tr('узел','node')
    + '</th><th scope="col">' + tr('доля','fraction') + '</th>'
    + '<th scope="col">' + tr('событий','events') + '</th></tr></thead>'
    + '<tbody>' + rows + '</tbody>';
  head.appendChild(h2); head.appendChild(em); head.appendChild(tb);
  var frag = document.createDocumentFragment();
  frag.appendChild(head); frag.appendChild(node);
  // Внутри клонированного шаблона — тоже применить перевод и подмену чисел.
  if (typeof applyLangTo === 'function') applyLangTo(frag);
  openPop(frag, opener);
}
// Полное имя элемента спектра. На графике подпись короткая (иначе налезает
// на кривые), а здесь — то, что читается в списке и в заголовке пояснения.
var ELEM_NAME = {
  peak:    {ru:'пик полного поглощения',       en:'full-energy peak'},
  edge:    {ru:'комптоновский край',            en:'Compton edge'},
  esc511:  {ru:'вылет одного кванта 511 кэВ',   en:'single 511 keV escape'},
  esc1022: {ru:'вылет обоих квантов 511 кэВ',   en:'double 511 keV escape'},
  xray:    {ru:'вылет K-рентгена',              en:'K X-ray escape'}
};
function elemName(id, fallback) {
  return pick(ELEM_NAME[id]) || pick(fallback) || id;
}
function openElement(id) {
  var node = tplNode('[data-el="' + id + '"]');
  if (node) openPop(node, null);
}

/* ── сводная таблица (внутри всплывающего окна) ──────────────────────── */
function buildTable() {
  var el = $('#tbl');
  if (!el) return;
  var maxPct = 0;
  DATA.tabs.forEach(function (t) {
    t.channels.forEach(function (c) { if (c.pct > maxPct) maxPct = c.pct; });
  });
  var h = '<thead><tr><th scope="col">канал</th>'
    + DATA.tabs.map(function (t) {
        return '<th scope="col">' + num(t.e0, 0) + ' ' + keVunit() + '</th>';
      }).join('') + '</tr></thead><tbody>';
  DATA.order.forEach(function (key) {
    h += '<tr style="--c:' + chColor(key) + '"><th scope="row">'
      + '<span class="sw"></span>' + labelOf(key) + ', %</th>';
    DATA.tabs.forEach(function (t) {
      var f = chanAt(t.e0, key);
      h += f ? '<td class="num" style="--w:' + (100 * f.pct / maxPct).toFixed(1)
        + '">' + num(f.pct, 2) + '</td>' : '<td>—</td>';
    });
    h += '</tr>';
  });
  h += '<tr class="sum"><th scope="row">'
    + tr('событий с сигналом, шт', 'events with signal') + '</th>'
    + DATA.tabs.map(function (t) { return '<td>' + cnt(t.n_signal) + '</td>'; })
      .join('') + '</tr>';
  h += '<tr><th scope="row">'
    + tr('ничего не вылетело из кристалла, %',
         'nothing escaped the crystal, %') + '</th>'
    + DATA.tabs.map(function (t) { return '<td>' + num(t.nofly_pct, 2) + '</td>'; })
      .join('') + '</tr>';
  h += '<tr><th scope="row">'
    + tr('в пике полного поглощения, %', 'in the full-energy peak, %') + '</th>'
    + DATA.tabs.map(function (t) { return '<td>' + num(t.peak_pct, 2) + '</td>'; })
      .join('') + '</tr>';
  el.innerHTML = h + '</tbody>';
}

/* ── полоса состава в шапке ──────────────────────────────────────────── */
var cvComp = $('#cvComp'), roComp = $('#roComp'), compGeom = null;

function drawComp() {
  var f = fit(cvComp, 150), g = f.g, W = f.W, H = f.H;
  var C = DATA.comp, L = 44, R = 10, T = 10, B = 24, i, k;
  g.clearRect(0, 0, W, H);
  var e0 = C.es[0], e1 = C.es[C.es.length - 1];
  var px = function (e) { return L + (e - e0) / (e1 - e0) * (W - L - R); };
  var py = function (v) { return T + (1 - v) * (H - T - B); };
  compGeom = {L: L, R: R, T: T, B: B, W: W, H: H};

  var N = Math.max(2, Math.round(W - L - R));
  var xs = [];
  for (i = 0; i <= N; i++) xs.push(e0 + i / N * (e1 - e0));

  // Интерполируется КАЖДАЯ накопленная граница, а не каждый канал по
  // отдельности: иначе между полосами появились бы щели и наложения.
  var prev = xs.map(function () { return 0; });
  for (k = 0; k < C.keys.length; k++) {
    var acc = C.f.map(function (r) {
      return r.slice(0, k + 1).reduce(function (a, b) { return a + b; }, 0);
    });
    var s = spline(C.es, acc);
    var cur = xs.map(function (x) { return Math.min(1, Math.max(0, s(x))); });
    var any = false;
    for (i = 0; i <= N; i++) if (cur[i] - prev[i] > 1e-4) { any = true; break; }
    if (any) {
      g.beginPath();
      g.moveTo(px(xs[0]), py(prev[0]));
      for (i = 0; i <= N; i++) g.lineTo(px(xs[i]), py(cur[i]));
      for (i = N; i >= 0; i--) g.lineTo(px(xs[i]), py(prev[i]));
      g.closePath();
      g.fillStyle = chColor(C.keys[k]);
      g.fill();
    }
    prev = cur;
  }

  g.strokeStyle = 'rgba(255,255,255,.26)'; g.lineWidth = 1;
  for (var v = 0.25; v < 1; v += 0.25) {
    g.beginPath(); g.moveTo(L, py(v) + .5); g.lineTo(W - R, py(v) + .5); g.stroke();
  }
  g.strokeStyle = css('--ink'); g.setLineDash([3, 3]); g.globalAlpha = .5;
  DATA.tabs.forEach(function (t) {
    g.beginPath(); g.moveTo(px(t.e0) + .5, T); g.lineTo(px(t.e0) + .5, H - B);
    g.stroke();
  });
  g.setLineDash([]); g.globalAlpha = 1;

  g.fillStyle = css('--faint'); g.font = '10px ' + css('--mono');
  g.textAlign = 'right';
  for (var q = 0; q <= 1.0001; q += 0.5)
    g.fillText(Math.round(q * 100) + ' %', L - 6, py(q) + 3.5);
  g.textAlign = 'center';
  for (var e = 500; e <= e1; e += 500) g.fillText(String(e), px(e), H - B + 15);
  g.textAlign = 'left'; g.fillText(keV(e0), px(e0) + 6, H - B + 15);
  g.strokeStyle = css('--rule'); g.lineWidth = 1;
  g.strokeRect(L + .5, T + .5, W - L - R - 1, H - T - B - 1);
  compGeom.px = px;
  drawZones($('#zonesComp'), DATA.comp.zones, W, L, R,
            function (e) { return px(e); }, e0, e1);
}

/* ── зоны спектра под графиком ───────────────────────────────────────────
   Полосы под холстом с реальным текстом: они читаются с экрана и с
   клавиатуры, в отличие от штриховых линий на самом графике. Ширина
   каждой полосы в пикселях считается той же координатной сеткой, что и
   ось X главного холста, — иначе цветной блок и подпись под осью уехали
   бы друг относительно друга при любом изменении размера окна. */
function drawZones(host, zones, W, L, R, px, xlo, xhi) {
  if (!host) return;
  host.textContent = '';
  host.style.paddingLeft = L + 'px';
  host.style.paddingRight = R + 'px';
  var vis = zones.filter(function (z) { return z.hi > xlo && z.lo < xhi; });
  // Плашки равной ширины — единственный ряд; связь с осью читается по
  // цвету и по диапазону в подписи. Пропорциональная цветная линейка над
  // плашками читалась как отдельная сущность, у которой другой масштаб,
  // чем у графика, и путала — снята.
  var row = document.createElement('div');
  row.className = 'zrow';
  vis.forEach(function (z) {
    var a = Math.max(z.lo, xlo), b = Math.min(z.hi, xhi);
    var el = document.createElement('div');
    el.className = 'zone'; el.setAttribute('data-id', z.id);
    var zl = pick(z.label);
    el.innerHTML = '<b>' + zl + '</b><i>' + keV(a) + '–' + keV(b)
      + ' ' + keVunit() + '</i>';
    el.title = zl + ': ' + keV(a) + '–' + keV(b) + ' ' + keVunit();
    row.appendChild(el);
  });
  host.appendChild(row);
}

cvComp.addEventListener('pointermove', function (ev) {
  if (!compGeom) return;
  var C = DATA.comp, r = cvComp.getBoundingClientRect();
  var x = ev.clientX - r.left, y = ev.clientY - r.top;
  if (x < compGeom.L || x > compGeom.W - compGeom.R
      || y < compGeom.T || y > compGeom.H - compGeom.B) {
    roComp.style.opacity = 0; return;
  }
  var e0 = C.es[0], e1 = C.es[C.es.length - 1];
  var e = e0 + (x - compGeom.L) / (compGeom.W - compGeom.L - compGeom.R) * (e1 - e0);
  var i = 0;
  C.es.forEach(function (v, k) {
    if (Math.abs(v - e) < Math.abs(C.es[i] - e)) i = k;
  });
  var rows = C.keys.map(function (k, j) { return [k, C.f[i][j]]; })
    .filter(function (p) { return p[1] >= 0.001; })
    .sort(function (a, b) { return b[1] - a[1]; });
  roComp.innerHTML = '<b>' + keV(C.es[i]) + ' ' + keVunit() + '</b>' + rows.map(function (p) {
    return '<span><i style="--c:' + chColor(p[0]) + '"></i>' + labelOf(p[0])
      + '<span class="v">' + pctf(p[1] * 100) + '</span></span>';
  }).join('');
  roComp.style.opacity = 1;
  var bw = roComp.offsetWidth;
  // #roComp вынесен из .hero в .app; top отсчитывается от .app и ложится
  // сразу ПОД канвас полосы состава — над плашками и подписью, но НЕ
  // налезая на приборный модуль. Горизонталь — по курсору в системе
  // координат канваса, дальше пересчитанным в .app.
  var appEl = cvComp.closest('.app');
  var cvBox = cvComp.getBoundingClientRect();
  var appBox = appEl.getBoundingClientRect();
  var cx = cvBox.left - appBox.left + x;
  var side = cx + bw + 28 > appBox.width ? -1 : 1;
  var lx = side > 0 ? cx + 14 : cx - bw - 14;
  roComp.style.left = Math.max(6,
    Math.min(lx, appBox.width - bw - 6)) + 'px';
  // Оверлей поднимается НАД канвасом полосы состава: между низом canvas и
  // приборным модулем места мало (высота списка каналов), поэтому ложим
  // сверху шапки. Если не влезает — прижимаем к верху .app.
  var ry = cvBox.top - appBox.top - roComp.offsetHeight - 6;
  if (ry < 4) ry = 4;
  roComp.style.top = ry + 'px';
});
cvComp.addEventListener('pointerleave', function () { roComp.style.opacity = 0; });

/* ── вкладки ─────────────────────────────────────────────────────────── */
var tabsEl = $('#tabs');
function addTab(label, i) {
  var b = document.createElement('button');
  b.type = 'button'; b.className = 'btn'; b.setAttribute('role', 'tab');
  b.textContent = label;
  b.setAttribute('aria-controls', 'cv');
  b.onclick = function () { selectTab(i); };
  b.onkeydown = function (ev) {
    var d = ev.key === 'ArrowRight' ? 1 : ev.key === 'ArrowLeft' ? -1 : 0;
    if (!d) return;
    ev.preventDefault();
    var n = tabsEl.children.length, k = (i + d + n) % n;
    tabsEl.children[k].focus(); selectTab(k);
  };
  tabsEl.appendChild(b);
}
DATA.tabs.forEach(function (t, i) {
  addTab(num(t.e0, 0) + ' ' + keVunit(), i);
});
addTab(tr('карта ', 'map ') + keV(DATA.run.e_lo) + '–' + keV(DATA.run.e_hi),
       DATA.tabs.length);

$('#rowSel').max = DATA.matrix.es.length - 1;
$('#rowSel').oninput = function () { drawMap(); drawSlice(); };
$('#cRaw').onchange = function () { drawMap(); drawSlice(); };

function selectTab(i) {
  tabIdx = i;
  Array.prototype.forEach.call(tabsEl.children, function (b, k) {
    b.setAttribute('aria-selected', k === i ? 'true' : 'false');
    b.tabIndex = k === i ? 0 : -1;
  });
  var isMap = i >= DATA.tabs.length;
  $('#viewChart').hidden = isMap;
  $('#viewMap').hidden = !isMap;
  $('#mLog').parentNode.hidden = isMap;
  if (isMap) {
    drawMap(); drawSlice();
    $('#elems').textContent = '';
    $('#zones').textContent = '';                       // на карте зон нет
    var b = $('#btnMap'); if (b) b.hidden = false;
    $('#hint').textContent = 'строка — энергия падающего кванта, столбец — '
      + 'энерговыделение, цвет — вероятность на квант; штриховая — пик '
      + 'полного поглощения, точечные — вылет одного и обоих квантов '
      + '511 кэВ';
  } else {
    var b = $('#btnMap'); if (b) b.hidden = true;
    buildLegend(); draw();
  }
}

/* ── легенда ─────────────────────────────────────────────────────────── */
function buildLegend() {
  var t = DATA.tabs[tabIdx], st = shown[tabIdx], leg = $('#leg');
  leg.textContent = '';
  function row(key, label, color, pct, cls) {
    var l = document.createElement('label');
    l.className = 'ch' + (cls ? ' ' + cls : '') + (st[key] ? '' : ' off');
    l.style.setProperty('--c', color);
    var cb = document.createElement('input');
    cb.type = 'checkbox'; cb.checked = st[key];
    cb.onchange = function () {
      st[key] = cb.checked; l.classList.toggle('off', !cb.checked); draw();
    };
    var sw = document.createElement('span'); sw.className = 'sw';
    var nm = document.createElement('span'); nm.className = 'nm';
    nm.textContent = label;
    var pc = document.createElement('span'); pc.className = 'pc';
    pc.textContent = pct === null ? '' : pctf(pct);
    l.appendChild(cb); l.appendChild(sw); l.appendChild(nm); l.appendChild(pc);
    if (cls !== 'total') {
      var go = document.createElement('button');
      go.type = 'button'; go.className = 'go'; go.textContent = '?';
      go.setAttribute('aria-label', 'что такое «' + label + '»');
      go.onclick = function (e) { e.preventDefault(); e.stopPropagation();
        openChannel(key, go); };
      l.appendChild(go);
    } else { l.appendChild(document.createElement('span')); }
    leg.appendChild(l);
  }
  if (mode === 'log')
    row('total', tr('полный отклик','full response'), css('--ink'), null, 'total');
  DATA.order.forEach(function (key) {
    var c = chanAt(t.e0, key);
    if (c) row(c.key, pick(c.label), chColor(c.key), c.pct);
  });
  $('#sideTitle').textContent = mode === 'log'
    ? tr('каналы','channels') : tr('каналы, доли','channels, fractions');
  $('#hint').textContent = tr(
    'событий с сигналом ' + cnt(t.n_signal) + ' из '
      + cnt(t.n_primaries) + ' выпущенных · канал отображения '
      + keV(t.step) + ' кэВ',
    'events with signal ' + cnt(t.n_signal) + ' out of '
      + cnt(t.n_primaries) + ' emitted · display bin '
      + keV(t.step) + ' keV');
  buildElems();
}

function buildElems() {
  var el = $('#elems');
  el.textContent = '';
  if (tabIdx >= DATA.tabs.length) return;
  var t = DATA.tabs[tabIdx];
  var cap = document.createElement('span');
  cap.textContent = tr('элементы спектра:', 'spectrum features:');
  el.appendChild(cap);
  t.markers.forEach(function (m) {
    var b = document.createElement('button');
    b.type = 'button';
    b.textContent = elemName(m.id, m.short) + ' · ' + keV(m.e)
      + ' ' + keVunit();
    b.onclick = function () { openElement(m.id); };
    el.appendChild(b);
  });
}
$('#bAll').onclick = function () {
  var st = shown[tabIdx];
  Object.keys(st).forEach(function (k) { st[k] = true; });
  buildLegend(); draw();
};
$('#bNone').onclick = function () {
  var st = shown[tabIdx];
  Object.keys(st).forEach(function (k) { st[k] = (k === 'total'); });
  buildLegend(); draw();
};
function setMode(m) {
  mode = m;
  $('#mLog').setAttribute('aria-pressed', m === 'log' ? 'true' : 'false');
  $('#mFrac').setAttribute('aria-pressed', m === 'frac' ? 'true' : 'false');
  buildLegend(); draw();
}
$('#mLog').onclick = function () { setMode('log'); };
$('#mFrac').onclick = function () { setMode('frac'); };

/* ── график отклика ──────────────────────────────────────────────────── */
var cv = $('#cv'), ro = $('#ro'), geom = null, hoverIdx = -1;

function drawMarkers(g, t, px, xmax, top, bot) {
  var acc = css('--accent'), faint = css('--faint'), paper = css('--paper');
  var boxes = [];
  t.markers.forEach(function (m) {
    if (m.e > xmax) return;
    var col = m.kind === 'peak' ? acc : faint;
    var x = px(m.e);
    // Штриховая линия сверху донизу поля, флажок-маркер над ней. Полная
    // подпись — в горизонтальном тултипе при наведении: повёрнутые подписи
    // на жёстких узлах сходились и не читались, а горизонтальные не
    // помещались бы под графиком.
    g.setLineDash([4, 3]); g.lineWidth = 1; g.strokeStyle = col;
    g.beginPath(); g.moveTo(x + .5, top); g.lineTo(x + .5, bot); g.stroke();
    g.setLineDash([]);
    g.beginPath();
    g.moveTo(x - 5.5, top - 2);
    g.lineTo(x + 5.5, top - 2);
    g.lineTo(x, top + 6);
    g.closePath();
    g.fillStyle = col; g.fill();
    g.strokeStyle = paper; g.lineWidth = 1.5; g.stroke();
    boxes.push({id: m.id, x: x, y0: top - 6, y1: bot,
                e: m.e, short: m.short, kind: m.kind});
  });
  return boxes;
}

function draw() {
  var t = DATA.tabs[tabIdx], st = shown[tabIdx];
  var f = fit(cv, 0), g = f.g, W = f.W, H = f.H;
  g.clearRect(0, 0, W, H);
  var L = 62, R = 14, T = 14, B = 40, i;
  var xmax = Math.min(3200, t.e0 * 1.15);
  var px = function (e) { return L + (e / xmax) * (W - L - R); };
  drawZones($('#zones'), t.zones, W, L, R, px, 0, xmax);
  var dim = css('--dim'), faint = css('--faint'), ink = css('--ink'),
      grid = css('--grid'), gridS = css('--grid-soft');

  g.lineWidth = 1;
  for (var e = 0; e <= xmax; e += 50) {
    g.strokeStyle = (e % 250 === 0) ? grid : gridS;
    g.beginPath(); g.moveTo(px(e) + .5, T); g.lineTo(px(e) + .5, H - B); g.stroke();
  }
  g.font = '11px ' + css('--mono'); g.fillStyle = faint; g.textAlign = 'center';
  var xTick = xmax > 2000 ? 500 : (xmax > 700 ? 250 : 50);
  for (var e2 = 0; e2 <= xmax; e2 += xTick) g.fillText(String(e2), px(e2), H - B + 16);
  g.font = '12px ' + css('--sans'); g.fillStyle = dim;
  g.fillText(tr('энерговыделение, кэВ','energy deposition, keV'),
             L + (W - L - R) / 2, H - 7);

  var N = Math.max(4, Math.round(W - L - R));
  var xs = [];
  for (i = 0; i <= N; i++) xs.push(i / N * xmax);
  var boxes = [];

  if (mode === 'log') {
    var top = Math.max.apply(null, t.total);
    // Три с небольшим десятичных порядка вниз: ниже лежит уже не структура
    // отклика, а разброс розыгрыша — единицы событий на канал отображения.
    // Верх с запасом под подписи элементов спектра.
    var ymax = top * 5, ymin = top / 2e3;
    var lgMin = Math.log10(ymin), lgMax = Math.log10(ymax);
    var py = function (v) {
      var q = Math.log10(Math.max(v, ymin * 0.35));
      return T + (1 - (Math.max(q, lgMin - 0.4) - lgMin) / (lgMax - lgMin))
        * (H - T - B);
    };
    geom = {L: L, R: R, T: T, B: B, W: W, H: H, xmax: xmax, px: px, py: py};

    g.textAlign = 'right'; g.font = '11px ' + css('--mono');
    for (var d = Math.ceil(lgMin); d <= lgMax; d++) {
      var y = py(Math.pow(10, d));
      if (y < T || y > H - B) continue;
      g.strokeStyle = grid;
      g.beginPath(); g.moveTo(L, y + .5); g.lineTo(W - R, y + .5); g.stroke();
      g.fillStyle = faint; g.fillText('10' + sup(d), L - 8, y + 3.5);
    }
    g.save(); g.translate(15, T + (H - T - B) / 2); g.rotate(-Math.PI / 2);
    g.textAlign = 'center'; g.font = '12px ' + css('--sans'); g.fillStyle = dim;
    g.fillText(tr('вероятность на квант в 4π','probability per γ in 4π'),
               0, 0); g.restore();

    boxes = drawMarkers(g, t, px, xmax, T, H - B);

    g.save(); g.beginPath(); g.rect(L, T, W - L - R, H - T - B); g.clip();
    // Сперва самые населённые каналы, чтобы редкие не тонули под заливкой.
    // Сперва ВСЕ заливки (самые населённые каналы ниже, чтобы не накрывали
    // редкие), и только потом ВСЕ линии: иначе заливка соседнего канала
    // ложится поверх уже проведённой линии и та пропадает.
    var vis = t.channels.filter(function (c) { return st[c.key]; })
      .slice().sort(function (a, b) { return b.pct - a.pct; });
    var paths = vis.map(function (c) {
      var s = spline(t.xs, c.ys.map(function (v) {
        return Math.log10(Math.max(v, ymin * 0.3));
      }));
      return {c: c, pts: xs.map(function (x) {
        return [px(x), py(Math.pow(10, s(x)))];
      })};
    });
    paths.forEach(function (o) {
      var gr = g.createLinearGradient(0, T, 0, H - B);
      gr.addColorStop(0, rgba(chColor(o.c.key), .22));
      gr.addColorStop(1, rgba(chColor(o.c.key), .015));
      g.beginPath(); g.moveTo(o.pts[0][0], H - B);
      o.pts.forEach(function (p) { g.lineTo(p[0], p[1]); });
      g.lineTo(o.pts[o.pts.length - 1][0], H - B); g.closePath();
      g.fillStyle = gr; g.fill();
    });
    paths.forEach(function (o) {
      g.beginPath();
      o.pts.forEach(function (p, k) { k ? g.lineTo(p[0], p[1]) : g.moveTo(p[0], p[1]); });
      g.strokeStyle = chColor(o.c.key); g.lineWidth = 1.6;
      g.lineJoin = 'round'; g.stroke();
    });
    if (st.total) {
      var st2 = spline(t.xs, t.total.map(function (v) {
        return Math.log10(Math.max(v, ymin * 0.3));
      }));
      g.beginPath();
      xs.forEach(function (x, k) {
        var yy = py(Math.pow(10, st2(x)));
        k ? g.lineTo(px(x), yy) : g.moveTo(px(x), yy);
      });
      g.strokeStyle = ink; g.lineWidth = 2.1; g.stroke();
    }
    g.restore();
  } else {
    var pyf = function (v) { return T + (1 - v) * (H - T - B); };
    geom = {L: L, R: R, T: T, B: B, W: W, H: H, xmax: xmax, px: px, py: pyf,
            frac: true};
    g.strokeStyle = grid; g.textAlign = 'right'; g.font = '11px ' + css('--mono');
    for (var q = 0; q <= 1.0001; q += 0.2) {
      g.beginPath(); g.moveTo(L, pyf(q) + .5); g.lineTo(W - R, pyf(q) + .5);
      g.stroke();
      g.fillStyle = faint; g.fillText(Math.round(q * 100) + ' %', L - 8, pyf(q) + 3.5);
    }
    g.save(); g.translate(15, T + (H - T - B) / 2); g.rotate(-Math.PI / 2);
    g.textAlign = 'center'; g.font = '12px ' + css('--sans'); g.fillStyle = dim;
    g.fillText(tr('доля канала в этой точке шкалы',
                  'channel fraction at this energy'), 0, 0); g.restore();

    g.save(); g.beginPath(); g.rect(L, T, W - L - R, H - T - B); g.clip();
    var prev = xs.map(function () { return 0; });
    var acc2 = t.xs.map(function () { return 0; });
    DATA.order.forEach(function (key) {
      var c = chanAt(t.e0, key);
      if (!c || !st[c.key]) return;
      acc2 = acc2.map(function (v, k) {
        return v + (t.total[k] > 0 ? c.ys[k] / t.total[k] : 0);
      });
      var s = spline(t.xs, acc2);
      var cur = xs.map(function (x) { return Math.min(1, Math.max(0, s(x))); });
      g.beginPath();
      g.moveTo(px(xs[0]), pyf(prev[0]));
      cur.forEach(function (v, k) { g.lineTo(px(xs[k]), pyf(v)); });
      for (var k2 = xs.length - 1; k2 >= 0; k2--) g.lineTo(px(xs[k2]), pyf(prev[k2]));
      g.closePath(); g.fillStyle = rgba(chColor(c.key), .9); g.fill();
      prev = cur;
    });
    g.restore();
    boxes = drawMarkers(g, t, px, xmax, T, H - B);
  }
  geom.boxes = boxes;

  if (hoverIdx >= 0 && hoverIdx < t.xs.length && t.xs[hoverIdx] <= xmax) {
    g.strokeStyle = css('--accent'); g.globalAlpha = .55; g.lineWidth = 1;
    g.beginPath(); g.moveTo(px(t.xs[hoverIdx]) + .5, T);
    g.lineTo(px(t.xs[hoverIdx]) + .5, H - B); g.stroke(); g.globalAlpha = 1;
  }
  g.strokeStyle = grid; g.lineWidth = 1;
  g.strokeRect(L + .5, T + .5, W - L - R - 1, H - T - B - 1);
}

var mkTip = $('#mkTip');
function hideMark() { if (mkTip && !mkTip.hidden) mkTip.hidden = true; }
cv.addEventListener('pointermove', function (ev) {
  if (!geom || tabIdx >= DATA.tabs.length) return;
  var r = cv.getBoundingClientRect();
  var x = ev.clientX - r.left, y = ev.clientY - r.top;
  if (x < geom.L || x > geom.W - geom.R || y < geom.T || y > geom.H - geom.B) {
    ro.style.opacity = 0;
    hideMark();
    if (hoverIdx !== -1) { hoverIdx = -1; draw(); }
    return;
  }
  var t = DATA.tabs[tabIdx], st = shown[tabIdx];
  var e = (x - geom.L) / (geom.W - geom.L - geom.R) * geom.xmax;
  // xs[i] — ЦЕНТР канала шириной step, первый центр = step/2; обратный
  // переход поэтому e/step − 1/2, а не e/step.
  var i = Math.max(0, Math.min(t.xs.length - 1, Math.round(e / t.step - 0.5)));
  // Ловим маркер шире, чем при клике: наведение мыши на штриховую линию
  // должно показывать полное имя маркера. Тултип позиционируется по
  // маркеру, а не по курсору: он не двигается пока курсор идёт вдоль
  // линии одного маркера.
  var near = geom.boxes.filter(function (b) {
    return Math.abs(b.x - x) <= 5 && y >= b.y0 && y <= b.y1;
  })[0];
  cv.style.cursor = near ? 'pointer' : 'crosshair';
  if (near && mkTip) {
    ro.style.opacity = 0;
    var name = elemName(near.id, near.short);
    mkTip.innerHTML = '<b>' + name + '</b><span>' + keV(near.e)
      + ' ' + keVunit() + tr(' · клик — пояснение',
                             ' · click for details') + '</span>';
    mkTip.hidden = false;
    var tw = mkTip.offsetWidth, th = mkTip.offsetHeight;
    // Тултип центрируется по маркеру, зажимается в поле графика.
    // Верх кладём чуть ниже маркера; если не влезает — переносим выше.
    var lx = Math.max(geom.L + 2,
      Math.min(near.x - tw / 2, geom.W - geom.R - tw - 2));
    var ty = geom.T + 14;
    if (ty + th > geom.H - geom.B - 4) ty = geom.T - th - 6;
    mkTip.style.left = lx + 'px';
    mkTip.style.top = ty + 'px';
    return;
  }
  hideMark();
  // Readout — доли каналов в %, а не абсолютные вероятности: сравнение
  // «сколько чего в этой точке шкалы» читается прямее, чем sci. Абсолютный
  // масштаб уходит в заголовок отдельным числом (вероятность полного
  // отклика на квант в 4π).
  var rows = [];
  var tot = t.total[i];
  t.channels.forEach(function (c) {
    if (!st[c.key] || c.ys[i] === 0) return;
    var pct = tot > 0 ? 100 * c.ys[i] / tot : 0;
    rows.push([chColor(c.key), pick(c.label), pctf(pct)]);
  });
  ro.innerHTML = '<b>' + keV(t.xs[i]) + ' ' + keVunit()
    + (mode === 'log' && tot > 0 ? ' · ' + sci(tot) : '')
    + '</b>' + rows.map(function (p) {
    return '<span><i style="--c:' + p[0] + '"></i>' + p[1]
      + '<span class="v">' + p[2] + '</span></span>';
  }).join('');
  ro.style.opacity = 1;
  var bw = ro.offsetWidth, bh = ro.offsetHeight;
  // Оверлей уходит ВЛЕВО, как только курсор пересекает середину поля,
  // и всегда с зазором в 16 px по горизонтали: не накрывает точку под
  // указателем, не выпадает за правый край.
  var side = x > geom.W / 2 ? -1 : 1;
  var lx = side > 0 ? x + 16 : x - bw - 16;
  ro.style.left = Math.max(4, Math.min(lx, geom.W - bw - 6)) + 'px';
  var ly = y - bh - 14;
  if (ly < geom.T) ly = Math.min(geom.H - bh - 6, y + 22);
  ro.style.top = ly + 'px';
  if (i !== hoverIdx) { hoverIdx = i; draw(); }
});
cv.addEventListener('pointerleave', function () {
  ro.style.opacity = 0;
  hideMark();
  if (hoverIdx !== -1) { hoverIdx = -1; draw(); }
});
cv.addEventListener('click', function (ev) {
  if (!geom || !geom.boxes) return;
  var r = cv.getBoundingClientRect();
  var x = ev.clientX - r.left, y = ev.clientY - r.top;
  var hit = geom.boxes.filter(function (b) {
    return Math.abs(b.x - x) <= 4 && y >= b.y0 && y <= b.y1;
  })[0];
  if (hit) openElement(hit.id);
});

/* ── карта отклика ───────────────────────────────────────────────────── */
var MAGMA = [[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],
             [229,80,100],[251,135,97],[254,194,135],[252,253,191]];
function magmaRGB(u) {
  u = Math.max(0, Math.min(1, u));
  var p = u * (MAGMA.length - 1);
  var i = Math.min(MAGMA.length - 2, Math.floor(p)), f = p - i;
  var a = MAGMA[i], b = MAGMA[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f,
          a[2] + (b[2] - a[2]) * f];
}
function magma(u) {
  var c = magmaRGB(u);
  return 'rgb(' + Math.round(c[0]) + ',' + Math.round(c[1]) + ','
    + Math.round(c[2]) + ')';
}
function matKey() { return $('#cRaw') && $('#cRaw').checked ? 'raw' : 'broadened'; }
var matLim = {}, matImg = {};
// Пределы цветовой шкалы СВОИ у каждой матрицы: у неразмытой пик собран в
// один канал и задирает верх на полтора порядка, при общей шкале размытая
// карта теряет весь контраст. Подпись шкалы перерисовывается вместе с ней.
function limits() {
  var key = matKey();
  if (matLim[key]) return matLim[key];
  var lo = Infinity, hi = 0;
  DATA.matrix[key].forEach(function (r) {
    r.forEach(function (v) { if (v > 0) { if (v < lo) lo = v; if (v > hi) hi = v; } });
  });
  if (!(hi > lo)) { lo = 1e-9; hi = 1e-3; }      // защита от вырожденной шкалы
  matLim[key] = [Math.log10(lo), Math.log10(hi)];
  return matLim[key];
}
// Матрица рисуется в буфер размером ровно в данные и растягивается с
// билинейной фильтрацией: попиксельные прямоугольники давали лестницу
// шириной в узел сетки, которой в самих данных нет.
function matrixImage() {
  var key = matKey();
  if (matImg[key]) return matImg[key];
  var rows = DATA.matrix[key], nc = DATA.matrix.cols.length, nr = rows.length;
  var lim = limits(), lo = lim[0], hi = lim[1];
  var off = document.createElement('canvas');
  off.width = nc; off.height = nr;
  var og = off.getContext('2d'), im = og.createImageData(nc, nr);
  for (var i = 0; i < nr; i++) {
    var y = nr - 1 - i;                       // ось энергии снизу вверх
    for (var j = 0; j < nc; j++) {
      var v = rows[i][j];
      // Пустая ячейка и самая редкая — разные вещи: ноль уводится ниже
      // начала шкалы, иначе они неразличимы.
      var c = magmaRGB(v > 0 ? (Math.log10(v) - lo) / (hi - lo) : -0.05);
      var o = (y * nc + j) * 4;
      im.data[o] = c[0]; im.data[o + 1] = c[1]; im.data[o + 2] = c[2];
      im.data[o + 3] = 255;
    }
  }
  og.putImageData(im, 0, 0);
  matImg[key] = off;
  return off;
}

function drawMap() {
  var f = fit($('#cvMap'), 300), g = f.g, W = f.W, H = f.H;
  g.clearRect(0, 0, W, H);
  var L = 62, R = 74, T = 12, B = 38;
  var cols = DATA.matrix.cols, es = DATA.matrix.es;
  var lim = limits(), lo = lim[0], hi = lim[1];
  var xhi = cols[cols.length - 1] + DATA.run.cell_keV / 2;
  var px = function (e) { return L + (e / xhi) * (W - L - R); };
  var py = function (e) {
    return T + (1 - (e - es[0]) / (es[es.length - 1] - es[0])) * (H - T - B);
  };

  g.save();
  g.imageSmoothingEnabled = true; g.imageSmoothingQuality = 'high';
  // Буфер кладётся ровно на ту же ось, по которой подписаны деления: строки
  // равномерны по индексу, а сетка узлов равномерна по энергии всюду, кроме
  // последнего интервала, — поэтому верхняя строка растягивается отдельно.
  var nr = es.length;
  var hFull = (H - T - B) * (nr - 1) / (nr - 1);
  g.drawImage(matrixImage(), px(cols[0] - DATA.run.cell_keV / 2), T,
              px(xhi) - px(cols[0] - DATA.run.cell_keV / 2), hFull);
  g.restore();

  g.lineWidth = 1;
  function diag(dE, dash, col) {
    g.setLineDash(dash); g.strokeStyle = col; g.beginPath();
    var first = true;
    es.forEach(function (e) {
      if (e - dE < 0) return;
      var x = px(e - dE), y = py(e);
      if (first) { g.moveTo(x, y); first = false; } else g.lineTo(x, y);
    });
    g.stroke(); g.setLineDash([]);
  }
  diag(0, [6, 4], 'rgba(180,225,255,.85)');
  diag(511, [2, 4], 'rgba(150,245,190,.8)');
  diag(1022, [2, 4], 'rgba(150,245,190,.5)');

  var sel = +$('#rowSel').value;
  g.strokeStyle = css('--accent'); g.lineWidth = 2;
  g.beginPath(); g.moveTo(L, py(es[sel])); g.lineTo(W - R, py(es[sel])); g.stroke();

  var faint = css('--faint'), dim = css('--dim'), grid = css('--grid');
  g.strokeStyle = grid; g.lineWidth = 1;
  g.strokeRect(L + .5, T + .5, W - L - R - 1, H - T - B - 1);
  g.font = '11px ' + css('--mono'); g.fillStyle = faint; g.textAlign = 'center';
  for (var e = 0; e <= xhi; e += 500) g.fillText(String(e), px(e), H - B + 16);
  g.textAlign = 'right';
  for (var e3 = 500; e3 <= es[es.length - 1]; e3 += 500)
    g.fillText(String(e3), L - 8, py(e3) + 3.5);
  g.font = '12px ' + css('--sans'); g.fillStyle = dim; g.textAlign = 'center';
  g.fillText(tr('энерговыделение, кэВ','energy deposition, keV'),
             L + (W - L - R) / 2, H - 7);
  g.save(); g.translate(15, T + (H - T - B) / 2); g.rotate(-Math.PI / 2);
  g.fillText(tr('энергия падающего кванта, кэВ',
                'incident γ energy, keV'), 0, 0); g.restore();

  var bx = W - R + 16, bw = 14, bh = H - T - B;
  for (var k = 0; k < bh; k++) {
    g.fillStyle = magma(1 - k / bh); g.fillRect(bx, T + k, bw, 1);
  }
  g.strokeStyle = grid; g.strokeRect(bx + .5, T + .5, bw, bh);
  g.fillStyle = faint; g.textAlign = 'left'; g.font = '10px ' + css('--mono');
  for (var d2 = Math.ceil(lo); d2 <= hi; d2++) {
    var yy = T + (1 - (d2 - lo) / (hi - lo)) * bh;
    g.fillText('10' + sup(d2), bx + bw + 4, yy + 3);
  }
}

function drawSlice() {
  var cvs = $('#cvSlice');
  if (!cvs) return;
  var f = fit(cvs, 190), g = f.g, W = f.W, H = f.H;
  g.clearRect(0, 0, W, H);
  var L = 62, R = 74, T = 12, B = 38;
  var sel = +$('#rowSel').value;
  var es = DATA.matrix.es, cols = DATA.matrix.cols;
  var yv = DATA.matrix[matKey()][sel];
  var lim = limits(), lo = lim[0], hi = lim[1];
  var xhi = cols[cols.length - 1] + DATA.run.cell_keV / 2;
  var px = function (e) { return L + (e / xhi) * (W - L - R); };
  var py = function (v) {
    return T + (1 - (Math.log10(Math.max(v, Math.pow(10, lo))) - lo) / (hi - lo))
      * (H - T - B);
  };
  var faint = css('--faint'), dim = css('--dim'),
      grid = css('--grid'), gridS = css('--grid-soft'), acc = css('--accent');
  g.lineWidth = 1;
  for (var e = 0; e <= xhi; e += 250) {
    g.strokeStyle = (e % 500 === 0) ? grid : gridS;
    g.beginPath(); g.moveTo(px(e) + .5, T); g.lineTo(px(e) + .5, H - B); g.stroke();
  }
  g.font = '11px ' + css('--mono'); g.textAlign = 'right';
  for (var d = Math.ceil(lo); d <= hi; d++) {
    var y = py(Math.pow(10, d));
    g.strokeStyle = grid;
    g.beginPath(); g.moveTo(L, y + .5); g.lineTo(W - R, y + .5); g.stroke();
    g.fillStyle = faint; g.fillText('10' + sup(d), L - 8, y + 3.5);
  }
  var N = Math.max(4, Math.round(W - L - R));
  var s = spline(cols, yv.map(function (v) {
    return Math.log10(Math.max(v, Math.pow(10, lo)));
  }));
  var pts = [];
  for (var i = 0; i <= N; i++) {
    var ee = i / N * xhi;
    pts.push([px(ee), py(Math.pow(10, s(Math.min(ee, cols[cols.length - 1]))))]);
  }
  g.save(); g.beginPath(); g.rect(L, T, W - L - R, H - T - B); g.clip();
  var gr = g.createLinearGradient(0, T, 0, H - B);
  gr.addColorStop(0, rgba(acc, .3)); gr.addColorStop(1, rgba(acc, .02));
  g.beginPath(); g.moveTo(pts[0][0], H - B);
  pts.forEach(function (p) { g.lineTo(p[0], p[1]); });
  g.lineTo(pts[pts.length - 1][0], H - B); g.closePath();
  g.fillStyle = gr; g.fill();
  g.beginPath();
  pts.forEach(function (p, k) { k ? g.lineTo(p[0], p[1]) : g.moveTo(p[0], p[1]); });
  g.strokeStyle = acc; g.lineWidth = 1.9; g.lineJoin = 'round'; g.stroke();
  g.restore();
  g.strokeStyle = grid; g.lineWidth = 1;
  g.strokeRect(L + .5, T + .5, W - L - R - 1, H - T - B - 1);
  g.fillStyle = faint; g.textAlign = 'center'; g.font = '11px ' + css('--mono');
  for (var e4 = 0; e4 <= xhi; e4 += 500) g.fillText(String(e4), px(e4), H - B + 16);
  g.font = '12px ' + css('--sans'); g.fillStyle = dim;
  g.fillText(tr('энерговыделение, кэВ','energy deposition, keV'),
             L + (W - L - R) / 2, H - 7);
  g.save(); g.translate(15, T + (H - T - B) / 2); g.rotate(-Math.PI / 2);
  g.fillText('вероятность на квант в 4π', 0, 0); g.restore();
  $('#rowOut').textContent = keV(es[sel]) + ' ' + keVunit();
}

/* ── запуск и перерисовка ────────────────────────────────────────────── */
function redraw() {
  drawComp();
  if (tabIdx >= DATA.tabs.length) { drawMap(); drawSlice(); } else draw();
}
var rt = null;
function redrawSoon() { clearTimeout(rt); rt = setTimeout(redraw, 110); }
addEventListener('resize', redrawSoon);
// Событие resize окна ловит не всё: полоса прокрутки, раскрытие вкладки и
// смена ширины контейнера при той же ширине окна оставляют холст с прежней
// внутренней шириной, и рисунок оказывается растянут относительно осей.
if (window.ResizeObserver) {
  var obs = new ResizeObserver(redrawSoon);
  document.querySelectorAll('.plotwrap, .mapwrap, .hero, .stage').forEach(function (el) {
    obs.observe(el);
  });
}
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', redraw);
new MutationObserver(redraw).observe(document.documentElement,
  {attributes: true, attributeFilter: ['data-theme']});

drawComp();
selectTab(0);
})();
