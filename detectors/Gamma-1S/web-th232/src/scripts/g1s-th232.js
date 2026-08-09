/* Разложение спектра Th-232 в маринелли по нуклидам ветви.
   Данные — window.G1S (собраны build_page.py из export_data.py).
   Ничего не тянет по сети; canvas → HiDPI-корректная отрисовка.

   Числа сводок и таблиц выводятся ЗДЕСЬ, а не подставляются при сборке:
   страница переключает закон ширины линии (по спектру / по цезию) и
   состав библиотеки метода 2, и у каждого числа есть четыре значения
   вместо одного. Источник тот же — выгрузка export_data.py, меняется
   только момент подстановки. */
(function () {
  "use strict";
  var D = window.G1S;
  if (!D) { console.error("Нет window.G1S"); return; }

  /* ── число/формат ───────────────────────────────────────────── */
  function num(x, d) {
    if (d === undefined) d = 1;
    return (Number(x)).toFixed(d).replace(".", ",");
  }
  function cnt(x) {
    var s = String(Math.round(Number(x))), out = "", neg = s[0] === "-";
    if (neg) s = s.slice(1);
    while (s.length > 3) { out = " " + s.slice(-3) + out; s = s.slice(0, -3); }
    return (neg ? "-" : "") + s + out;
  }
  function signedPct(ratio) {
    var s = 100 * (ratio - 1);
    return (s < 0 ? "−" : "+") + num(Math.abs(s), 1) + " %";
  }
  function esc(s) {
    return String(s === undefined || s === null ? "" : s)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

  /* ── состояние ──────────────────────────────────────────────── */
  var ST = {
    on:  {},          // нуклид -> показывать
    log: true,
    sum: true,        // суммарный слой (сумма отмеченных шаблонов)
    fwhmLaw: "lines", // закон ширины линии: lines (по спектру) | cs (цезий)
    lib: "fixed",     // состав библиотеки метода 2: fixed | full
    m2sort: "contrib",// сортировка таблицы метода 2: contrib | nuclide | energy
    cursorE: null,
  };
  D.nuclides.forEach(function (n) { ST.on[n.key] = true; });

  /* ── аксессоры: какой набор чисел действует сейчас ──────────── */
  function SRC()  { return ST.fwhmLaw === "cs" ? D.cs : D; }
  function SPEC() { return ST.fwhmLaw === "cs" ? D.cs.spectrum : D.spectrum; }
  function M1()   { return SRC().method1; }
  function M2()   {
    return ST.lib === "full" ? SRC().method2_full : SRC().method2;
  }
  function STACK1() { return SPEC().stack; }
  function STACK2() {
    return ST.lib === "full" ? SPEC().stack2_full : SPEC().stack2;
  }
  // R45: то же разложение метода 2, но по каналам взаимодействия, а не по
  // нуклидам — прямой ответ на директиву "метод 2 это и есть функция
  // отклика по каналам взаимодействия".
  function STACK2CHAN() {
    return ST.lib === "full" ? SPEC().stack2_chan_full : SPEC().stack2_chan;
  }
  // Сторож статистики (R66): по каналу на нуклид — набран ли МК-шаблон
  // настолько, чтобы доле нуклида в этом канале верить. Маска своя у каждого
  // закона ширины линии: n_eff считается ПОСЛЕ свёртки, и при другой ширине
  // те же отсчёты МК размазываются по другому числу каналов.
  function TRUSTED() { return SPEC().trusted || D.spectrum.trusted || null; }
  function TRUST_MASK(key) {
    var t = TRUSTED();
    var v = t && t[key];
    return (v && v.length === D.spectrum.e_of_ch.length) ? v : null;
  }

  /* ── HiDPI-canvas ───────────────────────────────────────────── */
  function fit(cv) {
    // Стиль canvas НЕ трогаем: правило CSS растягивает его на всю
    // плот-область; заданная в пикселях cv.style.width переопределила бы
    // его и заморозила размер до применения раскладки.
    var dpr = window.devicePixelRatio || 1;
    var r = cv.getBoundingClientRect();
    var w = Math.max(200, Math.floor(r.width));
    var h = Math.max(120, Math.floor(r.height || cv.height || 260));
    cv.width = Math.floor(w * dpr);
    cv.height = Math.floor(h * dpr);
    var g = cv.getContext("2d");
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, w, h);
    return { g: g, w: w, h: h };
  }

  /* ── палитра ────────────────────────────────────────────────── */
  function css(name, fallback) {
    var v = getComputedStyle(document.documentElement).getPropertyValue(name);
    return (v || "").trim() || fallback;
  }
  function pal() {
    return {
      ink:   css("--ink",   "#16140f"),
      rule:  css("--rule",  "#16140f"),
      faint: css("--faint", "#6a6558"),
      grid:  css("--grid",  "#dcd7c8"),
      paper: css("--paper", "#f5f2ea"),
      // Суммарная кривая разложения (R79) — красным, отдельно от черноты
      // измеренного спектра, с которым её и сравнивают.
      sum:   css("--sum-line", "#d21f1f"),
    };
  }

  /* ── шкалы ─────────────────────────────────────────────────── */
  function mapX(v, lo, hi, x0, x1) {
    return x0 + (v - lo) / (hi - lo) * (x1 - x0);
  }
  function makeY(logY, lo, hi) {
    if (logY) {
      lo = Math.max(0.5, lo);
      hi = Math.max(lo * 10, hi);
      var l0 = Math.log10(lo), l1 = Math.log10(hi);
      return {
        lo: lo, hi: hi, log: true,
        map: function (v, y0, y1) {
          var t = (Math.log10(Math.max(v, lo)) - l0) / (l1 - l0);
          return y0 + (1 - t) * (y1 - y0);
        }
      };
    }
    return {
      lo: 0, hi: hi, log: false,
      map: function (v, y0, y1) { return y0 + (1 - v / hi) * (y1 - y0); }
    };
  }
  function supNum(n) {
    var d = {"-":"⁻","0":"⁰","1":"¹","2":"²","3":"³","4":"⁴","5":"⁵",
             "6":"⁶","7":"⁷","8":"⁸","9":"⁹"};
    return String(n).split("").map(function (c) { return d[c] || c; }).join("");
  }

  /* ── стек ───────────────────────────────────────────────────── */
  function stackTotal(stk, i) {
    var acc = 0;
    for (var j = 0; j < D.nuclides.length; j++) {
      var k = D.nuclides[j].key;
      if (!ST.on[k] || !stk[k]) continue;
      acc += stk[k][i];
    }
    return acc;
  }

  /* ── общая отрисовка спектра с разложением ──────────────────── */
  // Один код на три вкладки: меняются только источник стека (метод 1 или
  // метод 2) и набор дополнительных кривых поверх заливок.
  function drawSpectrum(cvId, stk, extra) {
    var cv = document.getElementById(cvId);
    if (!cv) return null;
    var p = pal();
    var f = fit(cv);
    var g = f.g, W = f.w, H = f.h;
    var m = { l: 62, r: 14, t: 12, b: 34 };
    var e = D.spectrum.e_of_ch;
    var yy = D.spectrum.counts;
    var xLo = 0, xHi = e[e.length - 1];

    var vMax = 1;
    for (var i0 = 0; i0 < e.length; i0++) {
      if (e[i0] < 30 || e[i0] > 3500) continue;
      var v0 = Math.max(yy[i0], stackTotal(stk, i0));
      if (v0 > vMax) vMax = v0;
    }
    var Y = makeY(ST.log, ST.log ? 0.5 : 0, vMax * (ST.log ? 2.0 : 1.1));

    // сетка
    g.strokeStyle = p.grid; g.lineWidth = 1; g.beginPath();
    var xTicks = [250, 500, 750, 1000, 1250, 1500, 1750, 2000,
                  2250, 2500, 2750, 3000, 3250, 3500];
    for (var xi = 0; xi < xTicks.length; xi++) {
      if (xTicks[xi] > xHi) break;
      var xx = mapX(xTicks[xi], xLo, xHi, m.l, W - m.r);
      g.moveTo(xx, m.t); g.lineTo(xx, H - m.b);
    }
    var yTicks = [];
    if (ST.log) {
      for (var d = 0; d <= Math.ceil(Math.log10(Y.hi)); d++)
        yTicks.push(Math.pow(10, d));
    } else {
      var stp = Math.pow(10, Math.floor(Math.log10(Y.hi / 4)));
      var s0 = Math.ceil(Y.hi / 4 / stp) * stp;
      for (var v = s0; v < Y.hi; v += s0) yTicks.push(v);
    }
    for (var yi = 0; yi < yTicks.length; yi++) {
      var yg = Y.map(yTicks[yi], m.t, H - m.b);
      g.moveTo(m.l, yg); g.lineTo(W - m.r, yg);
    }
    g.stroke();

    // подписи осей
    g.fillStyle = p.faint;
    g.font = "11px system-ui, sans-serif";
    g.textAlign = "center"; g.textBaseline = "top";
    for (var xj = 0; xj < xTicks.length; xj++) {
      if (xTicks[xj] > xHi) break;
      g.fillText(String(xTicks[xj]),
                 mapX(xTicks[xj], xLo, xHi, m.l, W - m.r), H - m.b + 4);
    }
    g.textAlign = "right"; g.textBaseline = "middle";
    for (var yj = 0; yj < yTicks.length; yj++) {
      var lbl = ST.log ? "10" + supNum(Math.round(Math.log10(yTicks[yj])))
                       : cnt(yTicks[yj]);
      g.fillText(lbl, m.l - 4, Y.map(yTicks[yj], m.t, H - m.b));
    }
    g.textAlign = "center"; g.textBaseline = "bottom";
    g.fillText("энергия, кэВ", (m.l + W - m.r) / 2, H - 2);

    g.strokeStyle = p.rule; g.lineWidth = 2;
    g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);

    // Заливки — НЕЗАВИСИМЫМ НАЛОЖЕНИЕМ, каждая от нижней границы окна до
    // СВОЕГО значения. Прежде здесь был накопительный стек: слой рисовался
    // поверх суммы предыдущих, поэтому по картинке нельзя было прочитать
    // вклад нуклида — он ехал вверх вместе с соседями, а включение второго
    // флажка визуально «поднимало» первый. Двигаться при переключении
    // флажков должна СУММАРНАЯ кривая, а сами шаблоны стоять на месте.
    //
    // Порядок: сначала самые крупные (уходят назад), мелкие рисуются
    // последними и оказываются впереди. При обратном порядке крупный слой
    // накрыл бы мелкие целиком — они заливаются от одного и того же низа.
    // Ранжирование — по интегралу в видимой полосе, а не по порядку в
    // легенде: порядок ряда и порядок по величине вклада не совпадают.
    var order = [];
    for (var ni = 0; ni < D.nuclides.length; ni++) {
      var nk = D.nuclides[ni].key;
      if (!ST.on[nk] || !stk[nk]) continue;
      var s = 0;
      for (var si = 0; si < e.length; si++) {
        if (e[si] < xLo || e[si] > xHi) continue;
        s += stk[nk][si];
      }
      order.push({ nuc: D.nuclides[ni], area: s });
    }
    order.sort(function (a, b) { return b.area - a.area; });

    // Сторож статистики (R66). Долю нуклида в канале определяет его МК-шаблон;
    // там, где шаблон набран единицами отсчётов, доля — пуассоновский шум, а не
    // отклик детектора. Умноженная на большой полный отклик, она повторяет его
    // форму, и на логарифмической шкале это читается как физические пики
    // (тот же класс дефекта, что пойман в слое рентгена, R76). Поэтому ниже
    // порога заливки нет — только тонкий пунктир: значение из данных не
    // выброшено (сумма слоёв обязана сходиться с полным откликом, и суммарная
    // кривая считается по всем каналам), но читателю видно, где числу верить.
    //
    // Отрезок по каналам берётся ПО ГРАНИЦАМ бинов, а не по их центрам: иначе
    // одиночный достоверный канал вырождается в отрезок нулевой ширины и
    // пропадает с рисунка совсем.
    function edgeLo(i) {
      return mapX(i > 0 ? (e[i - 1] + e[i]) / 2 : e[0], xLo, xHi, m.l, W - m.r);
    }
    function edgeHi(i) {
      return mapX(i + 1 < e.length ? (e[i] + e[i + 1]) / 2 : e[e.length - 1],
                  xLo, xHi, m.l, W - m.r);
    }
    // Разбиение канальной оси на отрезки постоянного статуса достоверности.
    // Маски нет (старые данные без сторожа) — весь диапазон считается
    // достоверным, страница ведёт себя как до R66.
    function runs(msk, want) {
      var out = [], i2 = 0, n2 = e.length;
      while (i2 < n2) {
        var ok = msk ? !!msk[i2] : true;
        var j2 = i2;
        while (j2 + 1 < n2 && (msk ? !!msk[j2 + 1] : true) === ok) j2++;
        if (ok === want) out.push([i2, j2]);
        i2 = j2 + 1;
      }
      return out;
    }
    function segPath(vec, a, b) {
      g.beginPath();
      for (var q = a; q <= b; q++) {
        var x = mapX(e[q], xLo, xHi, m.l, W - m.r);
        var y = Y.map(vec[q] > Y.lo ? vec[q] : Y.lo, m.t, H - m.b);
        if (q === a) g.moveTo(edgeLo(a), y);
        g.lineTo(x, y);
        if (q === b) g.lineTo(edgeHi(b), y);
      }
    }

    var yBase = Y.map(Y.lo, m.t, H - m.b);
    function fillRuns(vec, rr, style) {
      g.fillStyle = style;
      for (var ri = 0; ri < rr.length; ri++) {
        segPath(vec, rr[ri][0], rr[ri][1]);
        g.lineTo(edgeHi(rr[ri][1]), yBase);
        g.lineTo(edgeLo(rr[ri][0]), yBase);
        g.closePath();
        g.fill();
      }
    }
    for (var oi = 0; oi < order.length; oi++) {
      var kf = order[oi].nuc.key;
      var vec = stk[kf];
      var mf = TRUST_MASK(kf);
      // Заливка ОДИНАКОВАЯ на всём слое (директива оператора 09.08.2026).
      // Отметка сторожа R66 держится только на контуре: достоверный участок
      // сплошной линией, недостоверный — пунктиром. Пробовали и разрыв
      // заливки, и штриховку: разрыв читался как дефект отрисовки, штриховка
      // на логарифмической шкале занимала треть картинки при вкладе в доли
      // процента. Числовая мера недостоверности осталась в подсказке легенды.
      fillRuns(vec, [[0, e.length - 1]], order[oi].nuc.color);
    }
    // Контур каждого шаблона своим цветом ПОВЕРХ всех заливок: там, где
    // мелкий слой локально выше крупного, заливка крупного перекрыта, и
    // без контура его ход в этом месте не прочитать.
    for (var oj = 0; oj < order.length; oj++) {
      var ks = order[oj].nuc.key;
      var vs = stk[ks];
      var msk = TRUST_MASK(ks);
      g.strokeStyle = order[oj].nuc.color;
      var solid = runs(msk, true);
      g.lineWidth = 1.2; g.setLineDash([]);
      for (var sj = 0; sj < solid.length; sj++) {
        segPath(vs, solid[sj][0], solid[sj][1]); g.stroke();
      }
      var noisy = runs(msk, false);
      g.lineWidth = 0.8; g.setLineDash([2, 3]);
      for (var nj = 0; nj < noisy.length; nj++) {
        segPath(vs, noisy[nj][0], noisy[nj][1]); g.stroke();
      }
      g.setLineDash([]);
    }

    function trace(getV, color, width, dash) {
      g.strokeStyle = color; g.lineWidth = width;
      if (dash) g.setLineDash(dash);
      g.beginPath();
      var st = false;
      for (var k = 0; k < e.length; k++) {
        if (e[k] < xLo || e[k] > xHi) continue;
        var v = getV(k);
        var x = mapX(e[k], xLo, xHi, m.l, W - m.r);
        var y = Y.map(v > Y.lo ? v : Y.lo, m.t, H - m.b);
        if (!st) { g.moveTo(x, y); st = true; } else g.lineTo(x, y);
      }
      g.stroke();
      if (dash) g.setLineDash([]);
    }

    // суммарный слой — верхняя граница стека отмеченных нуклидов
    if (ST.sum) {
      // Считается по ВСЕМ каналам, включая помеченные сторожем (R66):
      // сумма слоёв обязана сходиться с полным откликом, иначе кривая
      // перестанет быть моделью, к которой велась подгонка.
      trace(function (i) { return stackTotal(stk, i); }, p.sum, 1.8, [6, 3]);
    }
    (extra || []).forEach(function (ex) {
      if (!ex.arr) return;
      trace(function (i) { return ex.arr[i]; }, ex.color, ex.width || 1.4,
            ex.dash);
    });

    // измерение поверх всего
    trace(function (i) { return yy[i]; }, p.ink, 0.8);

    if (ST.cursorE !== null) {
      var xC = mapX(ST.cursorE, xLo, xHi, m.l, W - m.r);
      g.strokeStyle = p.rule; g.lineWidth = 1; g.setLineDash([4, 4]);
      g.beginPath(); g.moveTo(xC, m.t); g.lineTo(xC, H - m.b); g.stroke();
      g.setLineDash([]);
    }
    g.strokeStyle = p.rule; g.lineWidth = 2;
    g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);
    return { m: m, W: W, H: H, xLo: xLo, xHi: xHi, Y: Y };
  }

  /* ── легенда: флажки нуклидов + служебные слои ──────────────── */
  // Сколько каналов шаблона набрано выше порога сторожа (R66). null — маски
  // нет (данные посчитаны до введения сторожа), тогда легенда без пометок.
  function nucOf(key) {
    for (var i = 0; i < D.nuclides.length; i++)
      if (D.nuclides[i].key === key) return D.nuclides[i];
    return null;
  }
  // Доля ИНТЕГРАЛА слоя, лежащая в недостоверной области. Считать каналы
  // нельзя: у нуклида с короткой шкалой (Pb-212 обрывается на 479 кэВ) почти
  // все каналы пусты по физике, и счёт каналов объявил бы ненадёжными 89 %
  // шаблона, надёжного на 99,7 % по вкладу. null — маски нет (данные
  // посчитаны до введения сторожа), тогда легенда без пометок.
  function noiseFrac(key) {
    var nf = SPEC().noise_frac || D.spectrum.noise_frac;
    if (!nf || !(key in nf)) return null;
    return nf[key];
  }
  function guardHint(nf) {
    return "шаблон набран статистикой МК: " + num(100 * (1 - nf), 1) + " %"
      + " вклада слоя; остальное — области ниже "
      + (D.spectrum.n_eff_min || 0) + " отсчётов МК на канал, там слой идёт "
      + "пунктиром без заливки (доля нуклида определяется шумом шаблона)";
  }
  function buildLegend(elId) {
    var el = document.getElementById(elId);
    if (!el || el.dataset.built) return;
    el.dataset.built = "1";
    var html = "";
    D.nuclides.forEach(function (nuc) {
      // SECOND (вторичные пики) — временная сущность метода 2 (R33/R37),
      // отдельный чекбокс убран директивой R47: с переходом метода 2 на
      // полную матрицу отклика (R45) вторичные войдут в шаблон каждого
      // нуклида естественно, отдельной сущности не будет вовсе. До тех
      // пор вклад остаётся включённым в общую сумму (ST.on["SECOND"]
      // не трогается, по умолчанию true, см. инициализацию ST.on) — без
      // переключателя, но и без потери из total/«сумма».
      if (nuc.key === "SECOND") return;
      // Сторож статистики (R66) в легенде. Шаблон, нигде не набранный до
      // порога, помечается штриховкой вместо сплошного цвета: на графике
      // такой слой идёт одним пунктиром, без единой заливки, и без метки в
      // легенде читатель принял бы это за «нуклида просто мало».
      var nf = noiseFrac(nuc.key);
      var sw = nuc.color;
      var hint = (nf === null) ? "" : " title='" + esc(guardHint(nf)) + "'";
      html += "<label class='chip' data-nuc='" + nuc.key + "'" + hint + ">"
           + "<input type='checkbox' " + (ST.on[nuc.key] ? "checked" : "") + ">"
           + "<span class='sw' style='background:" + sw + "'></span>"
           + "<span class='nm'>" + esc(nuc.label_ru) + "</span>"
           + "</label>";
    });
    html += "<label class='chip'>"
         + "<input type='checkbox' class='c-sum' " + (ST.sum ? "checked" : "") + ">"
         + "<span class='sw sum'></span>"
         + "<span class='nm'>сумма</span></label>";
    html += "<label class='chip toggle'>"
         + "<input type='checkbox' class='c-log' " + (ST.log ? "checked" : "") + ">"
         + "<span class='sw log'></span><span class='nm'>лог</span></label>";
    el.innerHTML = html;
    el.querySelectorAll("label.chip[data-nuc] input").forEach(function (inp) {
      var key = inp.parentNode.getAttribute("data-nuc");
      inp.addEventListener("change", function (ev) {
        ST.on[key] = ev.target.checked; syncLegends(); redraw();
      });
    });
    el.querySelectorAll(".c-sum").forEach(function (inp) {
      inp.addEventListener("change", function (ev) {
        ST.sum = ev.target.checked; syncLegends(); redraw();
      });
    });
    el.querySelectorAll(".c-log").forEach(function (inp) {
      inp.addEventListener("change", function (ev) {
        ST.log = ev.target.checked; syncLegends(); redraw();
      });
    });
  }

  // Легенд две (по одной на вкладку метода), состояние — одно: после
  // клика на любой из них остальная обязана показывать то же самое.
  function syncLegends() {
    ["legendM1", "legendM2"].forEach(function (id) {
      var el = document.getElementById(id);
      if (!el || !el.dataset.built) return;
      el.querySelectorAll("label.chip[data-nuc]").forEach(function (lab) {
        var k = lab.getAttribute("data-nuc");
        lab.querySelector("input").checked = !!ST.on[k];
        // Пометка сторожа (R66) пересчитывается здесь, а не только при сборке
        // легенды: маска своя у каждого закона ширины линии, и переключатель
        // закона обязан её обновить — иначе штриховка останется от прежнего.
        var nf = noiseFrac(k), nuc = nucOf(k);
        if (nf === null || !nuc) return;
        lab.querySelector(".sw").style.background = nuc.color;
        lab.title = guardHint(nf);
      });
      el.querySelectorAll(".c-sum").forEach(function (i) { i.checked = ST.sum; });
      el.querySelectorAll(".c-log").forEach(function (i) { i.checked = ST.log; });
    });
  }

  /* ── сводки методов ─────────────────────────────────────────── */
  function cell(lab, val, big, hint) {
    return "<div" + (hint ? " title='" + hint + "'" : "") + "><span class='lab'>"
         + lab + "</span><span class='val"
         + (big ? " big-num" : "") + "'>" + val + "</span></div>";
  }
  // Подпись и справка величины bg_amplitude. Модель обоих методов — не
  // только нуклидные шаблоны/линии, но и измеренный фон той же геометрии
  // как ВТОРОЙ, свободный столбец матрицы плана (после приведения к живому
  // времени образца). Если бы шаблоны/линии полностью объясняли континуум
  // и подложку сами, множитель этого столбца выходил бы к единице — «это
  // и есть фон, без добавок». Он выходит заметно больше единицы и РАЗНЫЙ у
  // метода 1 и метода 2 (см. живые числа в карточках) — то есть столбец
  // работает заплаткой под ту часть комптоновского континуума, которую
  // нуклидная часть модели не воспроизводит, а не показывает кратность
  // реального фона в кювете (находка R44). Публичная подпись это и
  // называет — не «фон», а именно поправку/заплатку; название описательное,
  // не термин: устоявшегося слова в словаре контура (523 записи, домен
  // gamma-spec) не нашлось, чеканить нельзя, вердикт «термина нет»
  // выносит Терминолог. Полная справка — попап «как посчитано».
  var CONT_LAB = "поправка континуума, множитель";
  var CONT_HINT = "коэффициент второго, не нуклидного столбца подгонки "
                + "(приведённый фон); заметно больше единицы — заплатка "
                + "под континуум, а не кратность реального фона. Подробнее "
                + "— «как посчитано».";

  function fillSummaries() {
    var m1 = M1(), m2 = M2(), pass = D.passport;
    var s1 = document.getElementById("sumM1");
    if (s1) {
      s1.innerHTML =
        cell("активность ветви", cnt(m1.A_Bq) + " Бк <em>± "
             + cnt(m1.dA_Bq) + " Бк</em>", true)
        + cell("против паспорта", num(m1.A_Bq / pass.A_Bq, 3) + " ("
               + signedPct(m1.A_Bq / pass.A_Bq) + ")")
        + cell("χ²/ν", num(m1.chi2_ndof, 2))
        + cell("каналов в подгонке", cnt(m1.ndof))
        + cell("диапазон", num(m1.E_fit_lo, 0) + "–" + num(m1.E_fit_hi, 0)
               + " кэВ")
        + cell(CONT_LAB, num(m1.bg_amplitude, 2), false, CONT_HINT);
    }
    var s2 = document.getElementById("sumM2");
    if (s2) {
      s2.innerHTML =
        cell("активность ветви", cnt(m2.A_Bq) + " Бк <em>± "
             + cnt(m2.dA_Bq) + " Бк</em>", true)
        + cell("против паспорта", num(m2.A_Bq / pass.A_Bq, 3) + " ("
               + signedPct(m2.A_Bq / pass.A_Bq) + ")")
        + cell("χ²/ν", num(m2.chi2_ndof, 2))
        + cell("линий в модели", cnt(m2.n_lines) + " + " + cnt(m2.n_sum_peaks)
               + " сумм-пиков + K-рентген")
        + cell("каналов в подгонке", cnt(m2.n_channels_fit))
        + cell(CONT_LAB, num(m2.bg_amplitude, 2), false, CONT_HINT);
    }
  }

  /* ── таблицы методов ────────────────────────────────────────── */
  function buildM1() {
    var tbl = document.getElementById("tblM1");
    if (!tbl) return;
    var m1 = M1();
    var stk = STACK1();
    var head = "<thead><tr><th>нуклид</th><th class='num'>амплитуда, Бк</th>"
             + "<th class='num'>к паспорту</th><th class='num'>доля в спектре</th>"
             + "<th>пояснение</th></tr></thead>";
    var body = "<tbody>";
    var grand = 0;
    D.nuclides.forEach(function (nuc) {
      var arr = stk[nuc.key];
      if (!arr) return;
      for (var i = 0; i < arr.length; i++) grand += arr[i];
    });
    D.nuclides.forEach(function (nuc) {
      var arr = stk[nuc.key];
      if (!arr) return;
      var sum = 0;
      for (var i = 0; i < arr.length; i++) sum += arr[i];
      var amp, damp, tag;
      var v = m1.per_nuclide[nuc.key];
      if (!v) return;
      amp = v.A_Bq; damp = v.dA_Bq; tag = "";
      body += "<tr>"
        + "<td><span class='sw' style='background:" + nuc.color + "'></span>"
        + esc(nuc.label_ru) + "</td>"
        + "<td class='num'>" + cnt(amp) + " \u00b1 " + cnt(damp) + "</td>"
        + "<td class='num'>" + num(amp / D.passport.A_Bq, 3) + "</td>"
        + "<td class='num'>" + num(100 * sum / Math.max(grand, 1e-9), 1)
          + " %</td>"
        + "<td>" + esc(nuc.note) + tag + "</td></tr>";
    });
    tbl.innerHTML = head + body + "</tbody>";
  }

  function buildM2Chan() {
    var tbl = document.getElementById("tblM2Chan");
    if (!tbl || !D.channels) return;
    var sc = STACK2CHAN();
    if (!sc) { tbl.innerHTML = ""; return; }
    var rows = D.channels.map(function (ch) {
      var arr = sc[ch.key] || [];
      var sum = 0;
      for (var i = 0; i < arr.length; i++) sum += arr[i];
      return { ch: ch, sum: sum };
    });
    var total = rows.reduce(function (a, r) { return a + r.sum; }, 0);
    rows.sort(function (a, b) { return b.sum - a.sum; });
    var head = "<thead><tr><th>\u043a\u0430\u043d\u0430\u043b</th><th class='num'>\u0434\u043e\u043b\u044f</th>"
             + "<th class='num'>\u0441\u0447\u0451\u0442</th></tr></thead>";
    var body = "<tbody>";
    rows.forEach(function (r) {
      var pct = total > 0 ? 100 * r.sum / total : 0;
      body += "<tr>"
        + "<td><span class='sw' style='background:" + r.ch.color + "'></span>"
        + esc(r.ch.label_ru) + "</td>"
        + "<td class='num'>" + num(pct, 1) + " %</td>"
        + "<td class='num'>" + cnt(r.sum) + "</td></tr>";
    });
    tbl.innerHTML = head + body + "</tbody>";
  }

  function buildM2() {
    var tbl = document.getElementById("tblM2");
    if (!tbl) return;
    var nucCol = {}, nucLab = {}, nucIdx = {};
    D.nuclides.forEach(function (n, i) {
      nucCol[n.key] = n.color; nucLab[n.key] = n.label_ru; nucIdx[n.key] = i;
    });
    function predicted(r) { return r.predicted_net; }
    var sorters = {
      contrib: function (a, b) { return predicted(b) - predicted(a); },
      nuclide: function (a, b) {
        var d = (nucIdx[a.nuclide] || 0) - (nucIdx[b.nuclide] || 0);
        return d !== 0 ? d : a.E_keV - b.E_keV;
      },
      energy: function (a, b) { return a.E_keV - b.E_keV; },
    };
    var rows = M2().lines.slice().sort(sorters[ST.m2sort] || sorters.contrib);
    var head = "<thead><tr>"
      + "<th>\u043b\u0438\u043d\u0438\u044f</th><th>\u043d\u0443\u043a\u043b\u0438\u0434</th>"
      + "<th class='num'>I<sub>\u03b3</sub> \u043d\u0430 \u0440\u0430\u0441\u043f\u0430\u0434 \u043d\u0443\u043a\u043b\u0438\u0434\u0430</th>"
      + "<th class='num'>\u03b5<sub>\u041f\u041f</sub></th>"
      + "<th class='num'>\u043f\u0440\u0435\u0434\u0441\u043a\u0430\u0437\u0430\u043d\u043e</th><th>\u043f\u0440\u0438\u043c\u0435\u0447\u0430\u043d\u0438\u0435</th></tr></thead>";
    var body = "<tbody>";
    for (var i = 0; i < rows.length; i++) {
      var r = rows[i];
      var lineTxt, iTxt;
      if (r.kind === "sum") {
        lineTxt = num(r.E1_keV, 1) + "+" + num(r.E2_keV, 1) + " = "
                + num(r.E_keV, 1) + " \u043a\u044d\u0412";
        iTxt = num(r.I1_pct, 2) + " % \u00d7 " + num(r.I2_pct, 2) + " %";
      } else if (r.kind === "xray") {
        lineTxt = "K-\u0441\u0435\u0440\u0438\u044f, \u0446\u0435\u043d\u0442\u0440 " + num(r.E_keV, 1) + " \u043a\u044d\u0412";
        iTxt = num(r.I_gamma_pct, 1) + " % \u043d\u0430 \u0440\u0430\u0441\u043f\u0430\u0434 \u0432\u0435\u0442\u0432\u0438";
      } else {
        lineTxt = num(r.E_keV, 1) + " \u043a\u044d\u0412";
        // \u0412\u044b\u0445\u043e\u0434 \u043b\u0438\u043d\u0438\u0438 \u2014 \u041d\u0410 \u0420\u0410\u0421\u041f\u0410\u0414 \u0421\u0412\u041e\u0415\u0413\u041e \u041d\u0423\u041a\u041b\u0418\u0414\u0410, \u043a\u0430\u043a \u043e\u043d \u0441\u0442\u043e\u0438\u0442 \u0432 ENSDF:
        // \u0447\u0438\u0442\u0430\u0442\u0435\u043b\u044c \u0441\u0432\u0435\u0440\u044f\u0435\u0442 \u043a\u043e\u043b\u043e\u043d\u043a\u0443 \u0441 \u0431\u0438\u0431\u043b\u0438\u043e\u0442\u0435\u043a\u043e\u0439. \u0412\u0435\u0442\u0432\u043b\u0435\u043d\u0438\u0435 \u043e\u0442 \u0440\u043e\u0434\u0438\u0442\u0435\u043b\u044f \u0440\u044f\u0434\u0430
        // \u2014 \u043e\u0442\u0434\u0435\u043b\u044c\u043d\u044b\u0439 \u043c\u043d\u043e\u0436\u0438\u0442\u0435\u043b\u044c, \u0432 \u043c\u043e\u0434\u0435\u043b\u044c \u043e\u043d\u043e \u0432\u0445\u043e\u0434\u0438\u0442 (export_data.py,
        // w = BR\u00b7I/100) \u0438 \u043f\u043e\u043a\u0430\u0437\u044b\u0432\u0430\u0435\u0442\u0441\u044f \u0437\u0434\u0435\u0441\u044c \u0442\u0430\u043c, \u0433\u0434\u0435 \u043d\u0435 \u0440\u0430\u0432\u043d\u043e \u0435\u0434\u0438\u043d\u0438\u0446\u0435.
        // \u0415\u0434\u0438\u043d\u0441\u0442\u0432\u0435\u043d\u043d\u044b\u0439 \u0442\u0430\u043a\u043e\u0439 \u043d\u0443\u043a\u043b\u0438\u0434 \u0432 \u0432\u0435\u0442\u0432\u0438 \u2014 Tl-208: \u043e\u043d \u043e\u0431\u0440\u0430\u0437\u0443\u0435\u0442\u0441\u044f \u043b\u0438\u0448\u044c \u0432
        // \u03b1-\u0432\u0435\u0442\u043a\u0435 \u0440\u0430\u0441\u043f\u0430\u0434\u0430 Bi-212.
        iTxt = num(r.I_gamma_pct, r.I_gamma_pct < 0.1 ? 4 : 2) + " %";
        if (typeof r.branch === "number" && r.branch < 0.999)
          iTxt += " <em>\u00d7 " + num(100 * r.branch, 2) + " % \u0432\u0435\u0442\u0432\u044c</em>";
      }
      var tag = r.kind === "sum" ? " \u00b7 \u0441\u0443\u043c\u043c-\u043f\u0438\u043a"
              : (r.kind === "xray" ? " \u00b7 \u0440\u0435\u043d\u0442\u0433\u0435\u043d" : "");
      var epsTxt = (typeof r.eps_peak === "number")
        ? r.eps_peak.toExponential(3).replace(".", ",") : "\u2014";
      body += "<tr>"
        + "<td><span class='sw' style='background:"
        + (nucCol[r.nuclide] || "#999") + "'></span>" + lineTxt + "</td>"
        + "<td>" + esc(nucLab[r.nuclide] || r.nuclide) + tag + "</td>"
        + "<td class='num'>" + iTxt + "</td>"
        + "<td class='num'>" + epsTxt + "</td>"
        + "<td class='num'>" + cnt(predicted(r)) + "</td>"
        + "<td>" + esc(r.note || "") + "</td></tr>";
    }
    tbl.innerHTML = head + body + "</tbody>";
  }

  /* ── курсорные строки ───────────────────────────────────────── */
  function cursorText(elId, stk, tipId) {
    var el = document.getElementById(elId);
    if (!el) return;
    if (ST.cursorE === null) {
      el.textContent = "наведи курсор на канал спектра";
      return;
    }
    var e = D.spectrum.e_of_ch;
    var i = 0, best = Infinity;
    for (var k = 0; k < e.length; k++) {
      var dd = Math.abs(e[k] - ST.cursorE);
      if (dd < best) { best = dd; i = k; }
    }
    var meas = D.spectrum.counts[i];
    var contribs = [];
    D.nuclides.forEach(function (nuc) {
      if (!ST.on[nuc.key] || !stk[nuc.key]) return;
      var v = stk[nuc.key][i];
      if (v > 0.5) contribs.push({ nuc: nuc, v: v });
    });
    contribs.sort(function (a, b) { return b.v - a.v; });
    var top = contribs.slice(0, 4).map(function (c) {
      return c.nuc.label_ru + " " + cnt(c.v);
    }).join(" · ");
    var txt = num(e[i], 0) + " кэВ — измерено " + cnt(meas)
            + ", модель " + cnt(stackTotal(stk, i));
    if (top) txt += " — " + top;
    el.textContent = txt;
    var tip = document.getElementById(tipId);
    if (tip) tip.textContent = num(e[i], 0) + " кэВ · " + cnt(meas);
  }

  function attachCursor(cvId, tipId, onMove) {
    var cv = document.getElementById(cvId);
    var tip = document.getElementById(tipId);
    if (!cv) return;
    cv.addEventListener("pointermove", function (ev) {
      var r = cv.getBoundingClientRect();
      var x = ev.clientX - r.left, y = ev.clientY - r.top;
      var m = { l: 62, r: 14 };
      var e = D.spectrum.e_of_ch;
      var xHi = e[e.length - 1];
      if (x < m.l || x > r.width - m.r) ST.cursorE = null;
      else ST.cursorE = ((x - m.l) / (r.width - m.r - m.l)) * xHi;
      onMove();
      if (tip) {
        if (ST.cursorE === null) tip.hidden = true;
        else {
          tip.hidden = false;
          tip.style.left = x + "px";
          tip.style.top = Math.max(0, y) + "px";
        }
      }
    });
    cv.addEventListener("pointerleave", function () {
      ST.cursorE = null; onMove();
      if (tip) tip.hidden = true;
    });
  }

  /* ── перерисовка активной вкладки ───────────────────────────── */
  function redraw() {
    var vis = document.querySelector(".stage:not([hidden])");
    if (!vis) return;
    if (vis.id === "viewM1") {
      drawSpectrum("cvM1", STACK1(), []);
      cursorText("cursorM1", STACK1(), "m1-tip");
    } else if (vis.id === "viewM2") {
      drawSpectrum("cvM2", STACK2(), []);
      cursorText("cursorM2", STACK2(), "m2-tip");
    } else if (vis.id === "viewCmp") {
      drawCmp(); fillCmpTable();
    } else if (vis.id === "viewCal") {
      drawCal(); drawFwhm();
    }
  }

  // Смена закона ширины или состава библиотеки меняет ВСЕ числа страницы,
  // а не только текущую вкладку: сводки и таблицы перестраиваются целиком.
  function refreshAll() {
    // syncLegends здесь не ради флажков (их состояние не менялось), а ради
    // пометок сторожа R66: маска достоверности своя у каждого закона ширины
    // линии, и переключатель закона обязан её обновить в легенде.
    syncLegends();
    fillSummaries();
    buildM1();
    buildM2();
    buildM2Chan();
    fillCmpTable();
    redraw();
  }

  /* ── попапы ─────────────────────────────────────────────────── */
  function openPop(id) {
    var tpl = document.getElementById(id);
    if (!tpl) return;
    var pop = document.getElementById("pop");
    var scrim = document.getElementById("scrim");
    pop.innerHTML = "";
    pop.appendChild(tpl.content.cloneNode(true));
    var close = document.createElement("button");
    close.type = "button"; close.className = "pop-close";
    close.setAttribute("aria-label", "закрыть");
    close.textContent = "×";
    close.addEventListener("click", closePop);
    pop.appendChild(close);
    pop.hidden = false; scrim.hidden = false;
  }
  function closePop() {
    document.getElementById("pop").hidden = true;
    document.getElementById("scrim").hidden = true;
  }

  /* ── переключение вкладок ───────────────────────────────────── */
  function switchTab(name) {
    document.querySelectorAll("#tabs .tab").forEach(function (t) {
      t.setAttribute("aria-selected",
                     t.getAttribute("data-tab") === name ? "true" : "false");
    });
    var map = { m1: "viewM1", m2: "viewM2", cmp: "viewCmp", cal: "viewCal" };
    Object.keys(map).forEach(function (k) {
      var el = document.getElementById(map[k]);
      if (el) el.hidden = (k !== name);
    });
    if (name === "m1") { buildLegend("legendM1");
                         fillSummaries(); buildM1(); }
    if (name === "m2") { buildLegend("legendM2");
                         fillSummaries(); buildM2(); buildM2Chan(); }
    if (name === "cal") buildCal();
    redraw();
  }

  /* ── калибровка образца и фона ──────────────────────────────── */
  function buildCal() {
    var tbl = document.getElementById("tblCal");
    if (!tbl || tbl.dataset.built) return;
    tbl.dataset.built = "1";
    var m = D.meta;
    function coefsHtml(coefs) {
      return coefs.map(function (c, i) {
        var abs = Math.abs(c), s;
        if (abs === 0) s = "0";
        else if (abs >= 0.01 && abs < 10000) s = num(c, 6);
        else s = c.toExponential(4).replace(".", ",");
        return "<span class='mono'>c" + i + " = " + s + "</span>";
      }).join("<br>");
    }
    var head = "<thead><tr><th>параметр</th>"
             + "<th>образец (Th-232)</th><th>фон той же геометрии</th></tr></thead>";
    var body = "<tbody>"
      + "<tr><td>каналов</td><td class='num'>" + m.cal_sample.n_channels
      + "</td><td class='num'>" + m.cal_bg.n_channels + "</td></tr>"
      + "<tr><td>живое время, с</td><td class='num'>" + num(m.live_s, 2)
      + "</td><td class='num'>" + num(m.bg_live_s, 2) + "</td></tr>"
      + "<tr><td>реальное время, с</td><td class='num'>" + num(m.real_s, 2)
      + "</td><td class='num'>" + num(m.bg_real_s, 2) + "</td></tr>"
      + "<tr><td>мёртвое время, %</td><td class='num'>"
      + num(100 * (m.real_s - m.live_s) / m.real_s, 3) + "</td><td class='num'>"
      + num(100 * (m.bg_real_s - m.bg_live_s) / m.bg_real_s, 3) + "</td></tr>"
      + "<tr><td>степень полинома E(канал)</td><td class='num'>"
      + m.cal_sample.order + "</td><td class='num'>" + m.cal_bg.order + "</td></tr>"
      + "<tr><td>коэффициенты</td><td>" + coefsHtml(m.cal_sample.coefs)
      + "</td><td>" + coefsHtml(m.cal_bg.coefs) + "</td></tr>"
      + "<tr><td>масштаб фона (t_обр / t_фон)</td>"
      + "<td class='num' colspan='2'>" + num(m.bg_scale_time, 4) + "</td></tr>"
      + "</tbody>";
    tbl.innerHTML = head + body;
    buildFwhmTable();
  }

  var CAL = { smp: true, bg: true, diff: false, log: true, anch: true,
              zoom: null, drag: null, dragging: false };
  var CAL_wired = false;
  var CAL_MARGIN = { l: 62, r: 14, t: 12, b: 32 };
  var CAL_MARK_H = 7;        // высота маркера-флажка репера, px
  var CAL_MARK_HIT_PX = 5;   // допуск наведения по x на маркер, px

  function drawCal() {
    var cv = document.getElementById("cvCal");
    if (!cv) return;
    if (!CAL_wired) wireCal();
    var p = pal();
    var f = fit(cv);
    var g = f.g, W = f.w, H = f.h;
    var m = CAL_MARGIN;
    var e = D.spectrum.e_of_ch;
    var y = D.spectrum.counts;
    var b = D.spectrum.bg_counts;
    var xLo = CAL.zoom ? CAL.zoom.xLo : 0;
    var xHi = CAL.zoom ? CAL.zoom.xHi : e[e.length - 1];
    var vMax = 1, series = [];
    if (CAL.smp) series.push(y);
    if (CAL.bg)  series.push(b);
    if (CAL.diff) {
      var diff = new Array(e.length);
      for (var i = 0; i < e.length; i++) diff[i] = y[i] - b[i];
      series.push(diff);
    }
    series.forEach(function (arr) {
      for (var i = 0; i < arr.length; i++) {
        if (e[i] < xLo || e[i] > xHi) continue;
        var v = CAL.log ? Math.abs(arr[i]) : arr[i];
        if (v > vMax) vMax = v;
      }
    });
    var Y = makeY(CAL.log, CAL.log ? 0.5 : 0, vMax * (CAL.log ? 2 : 1.1));
    g.strokeStyle = p.grid; g.lineWidth = 1; g.beginPath();
    var xTicks = [250, 500, 750, 1000, 1250, 1500, 1750, 2000,
                  2250, 2500, 2750, 3000];
    for (var xi = 0; xi < xTicks.length; xi++) {
      if (xTicks[xi] > xHi) break;
      var xx = mapX(xTicks[xi], xLo, xHi, m.l, W - m.r);
      g.moveTo(xx, m.t); g.lineTo(xx, H - m.b);
    }
    var yTicks = [];
    if (CAL.log) {
      for (var d = 0; d <= Math.ceil(Math.log10(Y.hi)); d++)
        yTicks.push(Math.pow(10, d));
    } else {
      var stp = Math.pow(10, Math.floor(Math.log10(Y.hi / 4)));
      var s0 = Math.ceil(Y.hi / 4 / stp) * stp;
      for (var v = s0; v < Y.hi; v += s0) yTicks.push(v);
    }
    for (var yi = 0; yi < yTicks.length; yi++) {
      var yg = Y.map(yTicks[yi], m.t, H - m.b);
      g.moveTo(m.l, yg); g.lineTo(W - m.r, yg);
    }
    g.stroke();
    g.fillStyle = p.faint;
    g.font = "11px system-ui, sans-serif";
    g.textAlign = "center"; g.textBaseline = "top";
    for (var xj = 0; xj < xTicks.length; xj++) {
      if (xTicks[xj] > xHi) break;
      g.fillText(String(xTicks[xj]),
                 mapX(xTicks[xj], xLo, xHi, m.l, W - m.r), H - m.b + 4);
    }
    g.textAlign = "right"; g.textBaseline = "middle";
    for (var yj = 0; yj < yTicks.length; yj++) {
      var label = CAL.log ? "10" + supNum(Math.round(Math.log10(yTicks[yj])))
                          : cnt(yTicks[yj]);
      g.fillText(label, m.l - 4, Y.map(yTicks[yj], m.t, H - m.b));
    }
    g.textAlign = "center"; g.textBaseline = "bottom";
    g.fillText("энергия, кэВ", (m.l + W - m.r) / 2, H - 2);
    g.strokeStyle = p.rule; g.lineWidth = 2;
    g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);

    // Реперы: пунктир + маленький маркер-флажок сверху. Подпись НЕ на
    // канве — при густой сетке линий (16 реперов) повёрнутый текст
    // накладывается друг на друга и нечитаем. Подпись — в #cal-tip по
    // наведению на сам маркер (см. wireCal/calRefLineAt).
    if (CAL.anch && D.reference_lines) {
      var nucCol = {};
      D.nuclides.forEach(function (n) { nucCol[n.label_ru] = n.color; });
      D.reference_lines.forEach(function (r) {
        var E = r[0], nuc = r[1];
        if (E < xLo || E > xHi) return;
        var xa = mapX(E, xLo, xHi, m.l, W - m.r);
        var col = nucCol[nuc] || p.faint;
        g.strokeStyle = col; g.lineWidth = 1;
        g.setLineDash([3, 3]);
        g.beginPath(); g.moveTo(xa, m.t + CAL_MARK_H + 2); g.lineTo(xa, H - m.b);
        g.stroke();
        g.setLineDash([]);
        g.fillStyle = col;
        g.beginPath();
        g.moveTo(xa - 4, m.t); g.lineTo(xa + 4, m.t);
        g.lineTo(xa, m.t + CAL_MARK_H); g.closePath();
        g.fill();
        g.strokeStyle = p.paper; g.lineWidth = 1; g.stroke();
      });
    }

    function drawTrace(arr, color, allowNeg) {
      g.strokeStyle = color; g.lineWidth = 1.2;
      g.beginPath();
      var started = false;
      for (var k = 0; k < e.length; k++) {
        if (e[k] < xLo || e[k] > xHi) continue;
        var v = arr[k];
        if (!allowNeg && v < 0) v = 0;
        var vv = CAL.log ? Math.max(v, Y.lo) : v;
        var xt = mapX(e[k], xLo, xHi, m.l, W - m.r);
        var yt = Y.map(vv, m.t, H - m.b);
        if (!started) { g.moveTo(xt, yt); started = true; } else g.lineTo(xt, yt);
      }
      g.stroke();
    }
    if (CAL.bg)   drawTrace(b, "#0f5aa8", false);
    if (CAL.diff) drawTrace((function () {
      var dd = new Array(e.length);
      for (var i = 0; i < e.length; i++) dd[i] = y[i] - b[i];
      return dd;
    })(), "#c8541c", CAL.log ? false : true);
    if (CAL.smp) drawTrace(y, p.ink, false);

    if (CAL.drag) {
      var xa2 = Math.min(CAL.drag.x0, CAL.drag.x1);
      var xb2 = Math.max(CAL.drag.x0, CAL.drag.x1);
      g.fillStyle = "rgba(246,211,28,.22)";
      g.fillRect(xa2, m.t, xb2 - xa2, H - m.b - m.t);
      g.strokeStyle = "#16140f"; g.lineWidth = 1.5;
      g.setLineDash([4, 4]);
      g.strokeRect(xa2, m.t, xb2 - xa2, H - m.b - m.t);
      g.setLineDash([]);
    }
  }

  function calEfromX(x, rectWidth) {
    var e = D.spectrum.e_of_ch;
    var xLo = CAL.zoom ? CAL.zoom.xLo : 0;
    var xHi = CAL.zoom ? CAL.zoom.xHi : e[e.length - 1];
    var m = CAL_MARGIN;
    return xLo + ((x - m.l) / (rectWidth - m.r - m.l)) * (xHi - xLo);
  }

  // Репер под курсором (x,y в px канвы) — ближайший маркер-флажок в
  // пределах CAL_MARK_HIT_PX по x и полосы маркера по y (та же геометрия,
  // что рисует drawCal). null, если курсор не над маркером.
  function calRefLineAt(x, y, rectWidth) {
    if (!CAL.anch || !D.reference_lines) return null;
    if (y < CAL_MARGIN.t - 2 || y > CAL_MARGIN.t + CAL_MARK_H + 3) return null;
    var e = D.spectrum.e_of_ch;
    var xLo = CAL.zoom ? CAL.zoom.xLo : 0;
    var xHi = CAL.zoom ? CAL.zoom.xHi : e[e.length - 1];
    var best = null, bestD = CAL_MARK_HIT_PX;
    D.reference_lines.forEach(function (r) {
      var E = r[0];
      if (E < xLo || E > xHi) return;
      var xa = mapX(E, xLo, xHi, CAL_MARGIN.l, rectWidth - CAL_MARGIN.r);
      var d = Math.abs(xa - x);
      if (d <= bestD) { bestD = d; best = r; }
    });
    return best;
  }

  function wireCal() {
    ["c-smp", "c-bg", "c-diff", "c-log", "c-anch"].forEach(function (id) {
      var el = document.getElementById(id);
      if (!el) return;
      var key = { "c-smp": "smp", "c-bg": "bg", "c-diff": "diff",
                  "c-log": "log", "c-anch": "anch" }[id];
      el.addEventListener("change", function (ev) {
        CAL[key] = ev.target.checked; drawCal();
      });
    });
    var calReset = document.getElementById("cal-reset");
    if (calReset) calReset.addEventListener("click", function () {
      CAL.zoom = null; drawCal();
    });
    var cv = document.getElementById("cvCal");
    var ro = document.getElementById("cal-ro");
    var tip = document.getElementById("cal-tip");
    if (cv && ro) {
      cv.addEventListener("pointermove", function (ev) {
        var r = cv.getBoundingClientRect();
        var x = ev.clientX - r.left, y = ev.clientY - r.top;
        if (CAL.dragging) {
          CAL.drag.x1 = Math.max(CAL_MARGIN.l,
                                 Math.min(r.width - CAL_MARGIN.r, x));
          drawCal();
          return;
        }
        if (x < CAL_MARGIN.l || x > r.width - CAL_MARGIN.r) {
          ro.textContent = "";
          if (tip) tip.hidden = true;
          return;
        }
        var ref = calRefLineAt(x, y, r.width);
        if (ref) {
          var refTxt = ref[1] + " · " + num(ref[0], 1) + " кэВ";
          ro.textContent = refTxt;
          if (tip) {
            tip.hidden = false;
            tip.textContent = refTxt;
            tip.style.left = x + "px";
            tip.style.top = Math.max(0, y) + "px";
          }
          return;
        }
        var e = D.spectrum.e_of_ch;
        var E = calEfromX(x, r.width);
        var i = 0, best = Infinity;
        for (var k = 0; k < e.length; k++) {
          var dd = Math.abs(e[k] - E);
          if (dd < best) { best = dd; i = k; }
        }
        var smp = D.spectrum.counts[i], bgv = D.spectrum.bg_counts[i];
        ro.textContent = num(e[i], 0) + " кэВ · образец " + cnt(smp)
          + " · фон " + cnt(bgv) + " · разность " + cnt(smp - bgv);
        if (tip) {
          tip.hidden = false;
          tip.textContent = num(e[i], 0) + " кэВ · " + cnt(smp);
          tip.style.left = x + "px";
          tip.style.top = Math.max(0, y) + "px";
        }
      });
      cv.addEventListener("pointerleave", function () {
        if (!CAL.dragging) { ro.textContent = ""; if (tip) tip.hidden = true; }
      });
      // Протяжка мышью — выбор диапазона зумом. Отпускание слушается на
      // документе: за пределами канвы оно иначе теряется и drag «залипает».
      cv.addEventListener("mousedown", function (ev) {
        var r = cv.getBoundingClientRect();
        var x = ev.clientX - r.left;
        if (x < CAL_MARGIN.l || x > r.width - CAL_MARGIN.r) return;
        ev.preventDefault();
        CAL.dragging = true;
        CAL.drag = { x0: x, x1: x };
        drawCal();
      });
      document.addEventListener("mouseup", function () {
        if (!CAL.dragging) return;
        CAL.dragging = false;
        var r = cv.getBoundingClientRect();
        var x0 = CAL.drag.x0, x1 = CAL.drag.x1;
        CAL.drag = null;
        if (Math.abs(x1 - x0) < 6) { drawCal(); return; }   // это клик, не зум
        CAL.zoom = { xLo: Math.max(0, calEfromX(Math.min(x0, x1), r.width)),
                     xHi: calEfromX(Math.max(x0, x1), r.width) };
        drawCal();
      });
      cv.addEventListener("dblclick", function () {
        CAL.zoom = null; drawCal();
      });
    }
    CAL_wired = true;
  }

  /* ── калибровка по разрешению: точки и степенной закон ──────── */
  function buildFwhmTable() {
    var tbl = document.getElementById("tblFwhm");
    if (!tbl || !D.fwhm_cal) return;
    var fw = D.fwhm_cal;
    var head = "<thead><tr><th>линия, кэВ</th><th class='num'>центроида</th>"
      + "<th class='num'>ПШПВ, кэВ</th><th class='num'>разрешение</th>"
      + "<th class='num'>линий в окне</th><th class='num'>закон k·E<sup>p</sup></th>"
      + "<th class='num'>отклонение</th><th>статус</th></tr></thead>";
    var body = "<tbody>";
    fw.points.forEach(function (q) {
      if (!q.used) {
        body += "<tr class='row-dirty'><td>" + num(q.E_nominal, 1) + "</td>"
          + "<td class='num'>—</td><td class='num'>—</td><td class='num'>—</td>"
          + "<td class='num'>—</td><td class='num'>—</td><td class='num'>—</td>"
          + "<td>отброшена: " + esc(q.reject) + "</td></tr>";
        return;
      }
      body += "<tr><td>" + num(q.E_nominal, 1) + "</td>"
        + "<td class='num'>" + num(q.E_centroid, 1) + "</td>"
        + "<td class='num'>" + num(q.fwhm_keV, 2) + " ± "
        + num(q.d_fwhm_keV, 2) + "</td>"
        + "<td class='num'>" + num(q.res_pct, 2) + " %</td>"
        + "<td class='num'>" + q.n_lines_window + "</td>"
        + "<td class='num'>" + num(q.fwhm_model_keV, 2) + "</td>"
        + "<td class='num'>" + (q.dev_pct >= 0 ? "+" : "−")
        + num(Math.abs(q.dev_pct), 1) + " %</td>"
        + "<td>в подгонке</td></tr>";
    });
    body += "<tr class='sum'><td>степенной закон</td>"
      + "<td class='num' colspan='2'>ПШПВ = " + num(fw.k, 3) + "·E<sup>"
      + num(fw.p, 4) + "</sup></td>"
      + "<td class='num'>" + num(fw.res662_pct, 2) + " % на 662</td>"
      + "<td class='num'>" + fw.n_used + " из " + fw.n_anchors + "</td>"
      + "<td class='num'>" + num(fw.fwhm662_law, 1) + " кэВ</td>"
      + "<td class='num'>СКО " + num(fw.rms_dev_pct, 1) + " %</td>"
      + "<td>по цезию комплекта " + num(fw.fwhm662_cs, 1) + " кэВ</td></tr>";
    tbl.innerHTML = head + body + "</tbody>";
  }

  function drawFwhm() {
    var cv = document.getElementById("cvFwhm");
    if (!cv || !D.fwhm_cal) return;
    var fw = D.fwhm_cal;
    var p = pal();
    var f = fit(cv);
    var g = f.g, W = f.w, H = f.h;
    var m = { l: 62, r: 16, t: 14, b: 34 };
    var used = fw.points.filter(function (q) { return q.used; });
    if (!used.length) return;
    var xLo = 0, xHi = 2900;
    var vMax = 0;
    used.forEach(function (q) { vMax = Math.max(vMax, q.fwhm_keV); });
    vMax = Math.max(vMax, fw.k * Math.pow(xHi, fw.p)) * 1.15;

    g.strokeStyle = p.grid; g.lineWidth = 1; g.beginPath();
    var xTicks = [500, 1000, 1500, 2000, 2500];
    xTicks.forEach(function (t) {
      var x = mapX(t, xLo, xHi, m.l, W - m.r);
      g.moveTo(x, m.t); g.lineTo(x, H - m.b);
    });
    var yTicks = [25, 50, 75, 100, 125];
    yTicks.forEach(function (t) {
      if (t > vMax) return;
      var y = m.t + (1 - t / vMax) * (H - m.b - m.t);
      g.moveTo(m.l, y); g.lineTo(W - m.r, y);
    });
    g.stroke();
    g.fillStyle = p.faint; g.font = "11px system-ui, sans-serif";
    g.textAlign = "center"; g.textBaseline = "top";
    xTicks.forEach(function (t) {
      g.fillText(String(t), mapX(t, xLo, xHi, m.l, W - m.r), H - m.b + 4);
    });
    g.textAlign = "right"; g.textBaseline = "middle";
    yTicks.forEach(function (t) {
      if (t > vMax) return;
      g.fillText(String(t), m.l - 4, m.t + (1 - t / vMax) * (H - m.b - m.t));
    });
    g.textAlign = "center"; g.textBaseline = "bottom";
    g.fillText("энергия, кэВ", (m.l + W - m.r) / 2, H - 2);
    g.save();
    g.translate(12, (m.t + H - m.b) / 2);
    g.rotate(-Math.PI / 2);
    g.textAlign = "center"; g.textBaseline = "top";
    g.fillText("ПШПВ, кэВ", 0, 0);
    g.restore();
    g.strokeStyle = p.rule; g.lineWidth = 2;
    g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);

    function yOf(v) { return m.t + (1 - v / vMax) * (H - m.b - m.t); }

    g.strokeStyle = "#0f5aa8"; g.lineWidth = 2;
    g.beginPath();
    for (var E = 40; E <= xHi; E += 10) {
      var x = mapX(E, xLo, xHi, m.l, W - m.r);
      var y = yOf(fw.k * Math.pow(E, fw.p));
      if (E === 40) g.moveTo(x, y); else g.lineTo(x, y);
    }
    g.stroke();

    // Закон по записи цезия — для сравнения: одна опорная точка, корень.
    g.strokeStyle = p.faint; g.lineWidth = 1.5; g.setLineDash([5, 4]);
    g.beginPath();
    for (var E2 = 40; E2 <= xHi; E2 += 10) {
      var x2 = mapX(E2, xLo, xHi, m.l, W - m.r);
      var y2 = yOf(fw.fwhm662_cs * Math.sqrt(E2 / 661.657));
      if (E2 === 40) g.moveTo(x2, y2); else g.lineTo(x2, y2);
    }
    g.stroke();
    g.setLineDash([]);

    g.fillStyle = "#c8541c"; g.strokeStyle = "#c8541c"; g.lineWidth = 1.5;
    used.forEach(function (q) {
      var x = mapX(q.E_centroid, xLo, xHi, m.l, W - m.r);
      var y = yOf(q.fwhm_keV);
      g.beginPath();
      g.moveTo(x, yOf(q.fwhm_keV - q.d_fwhm_keV));
      g.lineTo(x, yOf(q.fwhm_keV + q.d_fwhm_keV));
      g.stroke();
      g.beginPath(); g.arc(x, y, 4, 0, 2 * Math.PI); g.fill();
    });

    g.font = "600 11px system-ui, sans-serif";
    g.textAlign = "left"; g.textBaseline = "top";
    g.fillStyle = "#c8541c";
    g.fillText("снято с этого спектра", m.l + 10, m.t + 8);
    g.fillStyle = "#0f5aa8";
    g.fillText("степенной закон k·E^p", m.l + 10, m.t + 24);
    g.fillStyle = p.faint;
    g.fillText("корневой закон по записи цезия", m.l + 10, m.t + 40);
  }

  /* ── сравнение ──────────────────────────────────────────────── */
  function cmpItems() {
    var p = pal();
    return [
      { lab: "паспорт", A: D.passport.A_Bq, dA: D.passport.dA_Bq, col: p.ink },
      { lab: "метод 1", A: M1().A_Bq, dA: M1().dA_Bq, col: "#0f5aa8" },
      { lab: "метод 2", A: M2().A_Bq, dA: M2().dA_Bq, col: "#c8541c" }
    ];
  }

  function fillCmpTable() {
    var el = document.getElementById("cmpTable");
    if (!el) return;
    var pass = D.passport, m1 = M1(), m2 = M2();
    var libTxt = ST.lib === "full" ? "все известные линии" : "отобранная библиотека";
    var modeTxt = ST.fwhmLaw === "cs" ? "ПШПВ по модели (цезий)"
                                      : "ПШПВ по линиям спектра";
    // Удельная активность — та же масса заливки для всех трёх строк
    // (паспорт, метод 1, метод 2 меряют один и тот же образец), поэтому
    // Бк/кг сравнимы напрямую и показаны здесь же, не только в паспорте.
    var massKg = pass.mass_g / 1000;
    function row(cls, lab, A, dA, note) {
      return "<div class='cmp-row " + cls + "'>"
        + "<span class='cmp-lab'>" + lab + "</span>"
        + "<span class='cmp-val big-num'>" + cnt(A) + " Бк <em>± "
        + cnt(dA) + " Бк</em> <em>· " + cnt(A / massKg) + " Бк/кг ± "
        + cnt(dA / massKg) + " Бк/кг</em></span>"
        + "<span class='cmp-note'>" + note + "</span></div>";
    }
    el.innerHTML =
      row("cmp-pass", "паспорт",
          pass.A_Bq, pass.dA_Bq,
          "удельная " + cnt(pass.Bq_per_kg) + " Бк/кг × "
          + num(pass.mass_g, 0) + " г; распад между аттестацией и "
          + "измерением при периоде 1,405·10¹⁰ лет от единицы неотличим")
    + row("cmp-m1", "метод 1: МК по нуклидам, " + num(m1.E_fit_lo, 0) + "–"
          + num(m1.E_fit_hi, 0) + " кэВ",
          m1.A_Bq, m1.dA_Bq,
          num(m1.A_Bq / pass.A_Bq, 3) + " паспорта, "
          + signedPct(m1.A_Bq / pass.A_Bq) + "; χ²/ν = " + num(m1.chi2_ndof, 2)
          + " на " + cnt(m1.ndof) + " каналах; " + modeTxt)
    + row("cmp-m2", "метод 2: функция ПП + библиотека, " + cnt(m2.n_lines)
          + " линий",
          m2.A_Bq, m2.dA_Bq,
          num(m2.A_Bq / pass.A_Bq, 3) + " паспорта, "
          + signedPct(m2.A_Bq / pass.A_Bq) + "; χ²/ν = " + num(m2.chi2_ndof, 2)
          + " на " + cnt(m2.n_channels_fit) + " каналах окон пиков; "
          + libTxt + ", " + modeTxt)
    + "<div class='cmp-row'><span class='cmp-lab'>расхождение методов</span>"
      + "<span class='cmp-val big-num'>" + signedPct(m1.A_Bq / m2.A_Bq)
      + "</span><span class='cmp-note'>метод 1 относительно метода 2 при "
      + "текущих переключателях</span></div>";
  }

  function drawCmp() {
    var cv = document.getElementById("cvCmp");
    if (!cv) return;
    var p = pal();
    var f = fit(cv);
    var g = f.g, W = f.w, H = f.h;
    var m = { l: 26, r: 20, t: 18, b: 34 };
    var items = cmpItems();
    var lo = Infinity, hi = -Infinity;
    items.forEach(function (it) {
      lo = Math.min(lo, it.A - it.dA);
      hi = Math.max(hi, it.A + it.dA);
    });
    var pad = (hi - lo) * 0.15 + 1;
    lo -= pad; hi += pad; lo = Math.max(0, lo);

    g.strokeStyle = p.rule; g.lineWidth = 2;
    g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);

    var range = hi - lo;
    var stp = Math.pow(10, Math.floor(Math.log10(range / 4)));
    var s = Math.max(stp, Math.ceil(range / 5 / stp) * stp);
    g.strokeStyle = p.grid; g.beginPath();
    g.fillStyle = p.faint; g.font = "11px system-ui, sans-serif";
    g.textAlign = "center"; g.textBaseline = "top";
    for (var v = Math.ceil(lo / s) * s; v <= hi; v += s) {
      var x = mapX(v, lo, hi, m.l, W - m.r);
      g.moveTo(x, m.t); g.lineTo(x, H - m.b);
      g.fillText(cnt(v), x, H - m.b + 4);
    }
    g.stroke();
    g.textAlign = "center"; g.textBaseline = "bottom";
    g.fillText("активность, Бк", (m.l + W - m.r) / 2, H - 2);

    var innerH = H - m.b - m.t;
    var rowH = innerH / items.length;
    for (var j = 0; j < items.length; j++) {
      var it = items[j];
      var yc = m.t + rowH * (j + 0.5);
      var xl = mapX(it.A - it.dA, lo, hi, m.l, W - m.r);
      var xr = mapX(it.A + it.dA, lo, hi, m.l, W - m.r);
      var xm = mapX(it.A, lo, hi, m.l, W - m.r);
      g.fillStyle = it.col; g.globalAlpha = 0.28;
      g.fillRect(xl, yc - rowH * 0.28, Math.max(xr - xl, 1), rowH * 0.56);
      g.globalAlpha = 1;
      g.strokeStyle = it.col; g.lineWidth = 3;
      g.beginPath();
      g.moveTo(xm, yc - rowH * 0.36); g.lineTo(xm, yc + rowH * 0.36);
      g.stroke();
      g.fillStyle = it.col;
      g.textAlign = "left"; g.textBaseline = "middle";
      g.font = "bold 13px system-ui, sans-serif";
      g.fillText(it.lab, m.l + 6, yc - rowH * 0.28);
      g.fillStyle = p.ink;
      g.font = "12px ui-monospace, Menlo, monospace";
      g.fillText(cnt(it.A) + " ± " + cnt(it.dA) + " Бк",
                 Math.min(xr + 8, W - m.r - 130), yc);
    }
  }

  /* ── события ────────────────────────────────────────────────── */
  function wire() {
    document.querySelectorAll("#tabs .tab").forEach(function (t) {
      t.addEventListener("click", function () {
        switchTab(t.getAttribute("data-tab"));
      });
    });
    document.querySelectorAll("#fwhmseg .btn").forEach(function (b) {
      b.addEventListener("click", function () {
        ST.fwhmLaw = b.getAttribute("data-fwhmlaw");
        document.querySelectorAll("#fwhmseg .btn").forEach(function (o) {
          o.setAttribute("aria-pressed",
            o.getAttribute("data-fwhmlaw") === ST.fwhmLaw ? "true" : "false");
        });
        refreshAll();
      });
    });
    // Обработчиков переключателя режима разложения больше нет: сам
    // переключатель снят из разметки (директива оператора 09.08.2026).
    document.querySelectorAll("#libseg .btn").forEach(function (b) {
      b.addEventListener("click", function () {
        ST.lib = b.getAttribute("data-lib");
        document.querySelectorAll("#libseg .btn").forEach(function (o) {
          o.setAttribute("aria-pressed",
            o.getAttribute("data-lib") === ST.lib ? "true" : "false");
        });
        refreshAll();
      });
    });
    document.querySelectorAll("#m2sortseg .btn").forEach(function (b) {
      b.addEventListener("click", function () {
        ST.m2sort = b.getAttribute("data-sort");
        document.querySelectorAll("#m2sortseg .btn").forEach(function (o) {
          o.setAttribute("aria-pressed",
            o.getAttribute("data-sort") === ST.m2sort ? "true" : "false");
        });
        buildM2();
      });
    });
    document.querySelectorAll("[data-pop]").forEach(function (b) {
      b.addEventListener("click", function () {
        openPop(b.getAttribute("data-pop"));
      });
    });
    var scrim = document.getElementById("scrim");
    if (scrim) scrim.addEventListener("click", closePop);
    document.addEventListener("keydown", function (ev) {
      if (ev.key === "Escape") closePop();
    });
    window.addEventListener("resize", redraw);
    attachCursor("cvM1", "m1-tip", function () {
      cursorText("cursorM1", STACK1(), "m1-tip");
      drawSpectrum("cvM1", STACK1(), []);
    });
    attachCursor("cvM2", "m2-tip", function () {
      cursorText("cursorM2", STACK2(), "m2-tip");
      drawSpectrum("cvM2", STACK2(), []);
    });
  }

  /* ── старт ──────────────────────────────────────────────────── */
  document.addEventListener("DOMContentLoaded", function () {
    wire();
    switchTab("m1");
  });
})();
