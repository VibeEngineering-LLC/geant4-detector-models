/* Разложение спектра Ra-226 в маринелли — ТОЛЬКО метод 2 (пилот
   обобщённого конвейера, задача #182/#183). Данные — window.RA226
   (export_ra226_data.py). Лёгкая версия g1s-th232.js: нет метода 1,
   нет варианта "cs" по отдельной калибровке цезия, нет масок
   достоверности МК-статистики (нет самих МК-шаблонов). */
(function () {
  "use strict";
  var D = window.RA226;
  if (!D) { console.error("Нет window.RA226"); return; }

  function num(x, d) {
    if (d === undefined) d = 1;
    return (Number(x)).toFixed(d).replace(".", ",");
  }
  function cnt(x) {
    var s = String(Math.round(Number(x))), out = "", neg = s[0] === "-";
    if (neg) s = s.slice(1);
    while (s.length > 3) { out = " " + s.slice(-3) + out; s = s.slice(0, -3); }
    return (neg ? "-" : "") + s + out;
  }
  function esc(s) {
    return String(s === undefined || s === null ? "" : s)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

  var ST = { on: {}, log: true, lib: "sel" };
  D.nuclides.forEach(function (n) { ST.on[n.key] = true; });

  function M2() { return ST.lib === "full" ? D.method2_full : D.method2_sel; }

  function fit(cv) {
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
  function css(name, fallback) {
    var v = getComputedStyle(document.documentElement).getPropertyValue(name);
    return (v || "").trim() || fallback;
  }
  function pal() {
    return {
      ink: css("--ink", "#16140f"), rule: css("--rule", "#16140f"),
      faint: css("--faint", "#6a6558"), grid: css("--grid", "#dcd7c8"),
      paper: css("--paper", "#f5f2ea"), sum: css("--sum-line", "#d21f1f"),
    };
  }
  function mapX(v, lo, hi, x0, x1) { return x0 + (v - lo) / (hi - lo) * (x1 - x0); }
  function makeY(logY, lo, hi) {
    if (logY) {
      lo = Math.max(0.5, lo); hi = Math.max(lo * 10, hi);
      var l0 = Math.log10(lo), l1 = Math.log10(hi);
      return { map: function (v, y0, y1) {
        var t = (Math.log10(Math.max(v, lo)) - l0) / (l1 - l0);
        return y0 + (1 - t) * (y1 - y0);
      } };
    }
    return { map: function (v, y0, y1) { return y0 + (1 - v / hi) * (y1 - y0); } };
  }

  function stackTotal(stk, i) {
    var acc = 0;
    for (var j = 0; j < D.nuclides.length; j++) {
      var k = D.nuclides[j].key;
      if (!ST.on[k] || !stk[k]) continue;
      acc += stk[k][i];
    }
    return acc;
  }

  function drawSpectrum() {
    var cv = document.getElementById("cvM2");
    if (!cv) return;
    var p = pal(), f = fit(cv), g = f.g, W = f.w, H = f.h;
    var m = { l: 60, r: 14, t: 12, b: 34 };
    var e = D.spectrum.e_of_ch;
    // Фон ВЫЧТЕН (замечание оператора №2, 09.08.2026): модель (сумма
    // нуклидов) не несёт фонового столба, поэтому сравнивать её надо с
    // ЧИСТЫМ спектром, не с сырыми отсчётами образца — иначе расхождение
    // в хвосте (там, где сигнал слабый, а фон относительно большой)
    // читается как ошибка модели, а не как невычтенный фон.
    var yy = D.spectrum.counts.map(function (c, i) {
      return c - D.spectrum.bg_counts[i];
    });
    var stk = M2().stack;
    var xLo = 0, xHi = e[e.length - 1];

    var vMax = 1;
    for (var i0 = 0; i0 < e.length; i0++) {
      if (e[i0] < 30 || e[i0] > 2400) continue;
      var v0 = Math.max(yy[i0], stackTotal(stk, i0));
      if (v0 > vMax) vMax = v0;
    }
    var Y = makeY(ST.log, ST.log ? 0.5 : 0, vMax * (ST.log ? 2.0 : 1.1));

    g.strokeStyle = p.grid; g.lineWidth = 1; g.beginPath();
    var xTicks = [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250];
    xTicks.forEach(function (tv) {
      var x = mapX(tv, xLo, xHi, m.l, W - m.r);
      g.moveTo(x, m.t); g.lineTo(x, H - m.b);
    });
    g.stroke();

    g.fillStyle = p.faint; g.font = "11px var(--mono, monospace)";
    g.textAlign = "center"; g.textBaseline = "top";
    xTicks.forEach(function (tv) {
      var x = mapX(tv, xLo, xHi, m.l, W - m.r);
      g.fillText(String(tv), x, H - m.b + 6);
    });

    var x0 = m.l, x1 = W - m.r, y0 = m.t, y1 = H - m.b;
    var n = e.length;

    // Независимое наложение (R65: НЕ накопительный стек) — каждый нуклид
    // заливкой от нуля до СВОЕГО вклада, полупрозрачно; общая красная
    // кривая ниже несёт сумму. Порядок отрисовки (замечание оператора
    // №4, 09.08.2026): КРУПНЕЙШИЙ по площади нуклид рисуется ПЕРВЫМ (на
    // заднем плане), мелкие — поверх, иначе доминирующий Bi-214,
    // нарисованный последним по алфавиту, хоронит все остальные заливки
    // под собой (та же сортировка, что g1s-th232.js:228-239).
    var order = D.nuclides.filter(function (nd) {
      return ST.on[nd.key] && stk[nd.key];
    }).map(function (nd) {
      var s = 0;
      for (var k = 0; k < n; k++) s += stk[nd.key][k];
      return { nd: nd, area: s };
    });
    order.sort(function (a, b) { return b.area - a.area; });

    order.forEach(function (o) {
      var nd = o.nd;
      g.fillStyle = nd.color; g.globalAlpha = 0.55;
      g.beginPath();
      g.moveTo(mapX(e[0], xLo, xHi, x0, x1), Y.map(0, y0, y1));
      for (var i = 0; i < n; i++) {
        var x = mapX(e[i], xLo, xHi, x0, x1);
        g.lineTo(x, Y.map(stk[nd.key][i], y0, y1));
      }
      g.lineTo(mapX(e[n - 1], xLo, xHi, x0, x1), Y.map(0, y0, y1));
      g.closePath(); g.fill();
      g.globalAlpha = 1;
    });

    g.strokeStyle = p.sum; g.lineWidth = 1.4; g.beginPath();
    for (var i2 = 0; i2 < n; i2++) {
      var xs = mapX(e[i2], xLo, xHi, x0, x1);
      var ys = Y.map(stackTotal(stk, i2), y0, y1);
      if (i2 === 0) g.moveTo(xs, ys); else g.lineTo(xs, ys);
    }
    g.stroke();

    g.strokeStyle = p.ink; g.lineWidth = 1; g.globalAlpha = 0.85;
    g.beginPath();
    for (var i3 = 0; i3 < n; i3++) {
      var xm = mapX(e[i3], xLo, xHi, x0, x1);
      var ym = Y.map(yy[i3], y0, y1);
      if (i3 === 0) g.moveTo(xm, ym); else g.lineTo(xm, ym);
    }
    g.stroke(); g.globalAlpha = 1;

    g.strokeStyle = p.rule; g.lineWidth = 1.2;
    g.beginPath(); g.moveTo(m.l, y0); g.lineTo(m.l, y1); g.lineTo(x1, y1); g.stroke();
  }

  // Не у каждого звена цепочки есть линия в модели — Rn-222/Po-218/Po-214
  // гамма практически не дают (IAEA Live Chart: Po-218 не отдал ни одной
  // строки, у Rn-222/Po-214 максимум найденной линии <0,1%). Оставлены в
  // списке для полноты цепочки распада, но их вклад в модель ТОЧНЫЙ НОЛЬ —
  // чекбокс декоративен, и легенда обязана честно это показывать, а не
  // тихо рисовать пустой квадратик рядом с рабочими нуклидами.
  function isInert(key) {
    var v = D.method2_sel.stack[key];
    if (!v) return true;
    for (var i = 0; i < v.length; i++) if (v[i] !== 0) return false;
    return true;
  }

  function buildLegend() {
    var el = document.getElementById("legendM2");
    if (!el) return;
    el.innerHTML = "";
    D.nuclides.forEach(function (nd) {
      var inert = isInert(nd.key);
      var chip = document.createElement("label");
      chip.className = "chip";
      if (inert) { chip.style.opacity = "0.45"; chip.title =
        "γ-линий в модели нет (пренебрежимо малый выход) — оставлен для полноты цепочки"; }
      var cb = document.createElement("input");
      cb.type = "checkbox"; cb.checked = ST.on[nd.key];
      cb.disabled = inert;
      cb.addEventListener("change", function () {
        ST.on[nd.key] = cb.checked; drawSpectrum();
      });
      var sw = document.createElement("span");
      sw.className = "sw"; sw.style.background = nd.color;
      var lb = document.createElement("span");
      lb.className = "nm";
      lb.textContent = nd.label_ru + (inert ? " (γ нет)" : "");
      chip.appendChild(cb); chip.appendChild(sw); chip.appendChild(lb);
      el.appendChild(chip);
    });
  }

  function cell(lab, val, big) {
    return "<div><span class='lab'>" + lab + "</span><span class='val" +
      (big ? " big-num" : "") + "'>" + val + "</span></div>";
  }

  function fillSummary() {
    var el = document.getElementById("sumM2");
    if (!el) return;
    var m2 = M2(), pass = D.passport;
    el.innerHTML =
      cell("активность (метод 2)", cnt(m2.A_Bq) + " Бк <em>± " + cnt(m2.dA_Bq) + " Бк</em>", true) +
      cell("против паспорта", num(m2.A_Bq / pass.A_Bq, 3) + " (" + signedPct(m2.A_Bq / pass.A_Bq) + ")") +
      cell("χ²/ν", num(m2.chi2_ndof, 2)) +
      cell("линий в модели", cnt(m2.n_lines) + " + " + cnt(m2.n_sum_peaks) + " сумм-пиков") +
      cell("амплитуда фона", num(m2.bg_amplitude, 2));
  }
  function signedPct(ratio) {
    var s = 100 * (ratio - 1);
    return (s < 0 ? "−" : "+") + num(Math.abs(s), 1) + " %";
  }

  function labelRu(key) {
    for (var i = 0; i < D.nuclides.length; i++)
      if (D.nuclides[i].key === key) return D.nuclides[i].label_ru;
    return key;
  }

  function fillTable() {
    var tbl = document.getElementById("tblM2");
    if (!tbl) return;
    var m2 = M2();
    var rows = m2.lines.slice().sort(function (a, b) { return a.E_keV - b.E_keV; });
    var html = "<thead><tr><th>E, кэВ</th><th>нуклид</th>" +
      "<th>I<sub>γ</sub>, %</th><th>тип</th><th>примечание</th></tr></thead><tbody>";
    rows.forEach(function (r) {
      html += "<tr><td>" + num(r.E_keV, 3) + "</td><td>" + esc(labelRu(r.nuclide)) +
        "</td><td>" + (r.I_pct === null || r.I_pct === undefined ? "—" : num(r.I_pct, 3)) +
        "</td><td>" + (r.kind === "sum" ? "сумма" : "линия") +
        "</td><td>" + esc(r.note || "") + "</td></tr>";
    });
    html += "</tbody>";
    tbl.innerHTML = html;
  }

  function fillCompare() {
    var cv = document.getElementById("cvCmp");
    var tbl = document.getElementById("cmpTable");
    var pass = D.passport;
    var items = [
      { lab: "паспорт", A: pass.A_Bq, dA: pass.dA_Bq, col: "#6a6558" },
      { lab: "метод 2 (отобр.)", A: D.method2_sel.A_Bq, dA: D.method2_sel.dA_Bq, col: "#0f5aa8" },
      { lab: "метод 2 (полн.)", A: D.method2_full.A_Bq, dA: D.method2_full.dA_Bq, col: "#c8541c" },
    ];
    if (cv) {
      var f = fit(cv), g = f.g, W = f.w, H = f.h;
      var maxA = Math.max.apply(null, items.map(function (it) { return it.A + it.dA; })) * 1.15;
      var m = { l: 130, r: 60, t: 14, b: 14 };
      var rowH = (H - m.t - m.b) / items.length;
      items.forEach(function (it, i) {
        var y = m.t + i * rowH + rowH * 0.22;
        var barH = rowH * 0.56;
        var x0 = m.l, x1 = W - m.r;
        var wA = (it.A / maxA) * (x1 - x0);
        var wD = (it.dA / maxA) * (x1 - x0);
        g.fillStyle = it.col; g.globalAlpha = 0.85;
        g.fillRect(x0, y, wA, barH);
        g.globalAlpha = 0.35;
        g.fillRect(x0 + wA - wD, y, 2 * wD, barH);
        g.globalAlpha = 1;
        g.fillStyle = css("--ink", "#16140f");
        g.font = "12px var(--sans, sans-serif)"; g.textAlign = "right"; g.textBaseline = "middle";
        g.fillText(it.lab, m.l - 10, y + barH / 2);
        g.textAlign = "left";
        g.fillText(cnt(it.A) + " ± " + cnt(it.dA), x0 + wA + 8, y + barH / 2);
      });
    }
    if (tbl) {
      var html = "<table class='big'><thead><tr><th>оценка</th><th>A, Бк</th><th>отношение к паспорту</th></tr></thead><tbody>";
      items.forEach(function (it) {
        var r = it.A / pass.A_Bq;
        html += "<tr><td>" + it.lab + "</td><td>" + cnt(it.A) + " ± " + cnt(it.dA) +
          "</td><td>" + num(r, 3) + "</td></tr>";
      });
      html += "</tbody></table>";
      tbl.innerHTML = html;
    }
    var rn = document.getElementById("radonNote");
    if (rn && D.radon_check) {
      // Честный статус (замечание оператора №3, 09.08.2026): попытка
      // независимой проверки утечки радона НЕ дала надёжного числа —
      // см. export_ra226_data.py.radon_check.reason. Метод 2 в текущем
      // виде (одна амплитуда на всю цепочку) утечку в принципе не видит.
      rn.innerHTML = "<b>Утечка радона (Rn-222, T½=3,82 сут) — не проверена.</b> " +
        esc(D.radon_check.reason);
    }
  }

  /* ── калибровка: спектры образца/фона (порт g1s-th232.js один в один,
     задача — калибровочная вкладка обязана быть идентична референсу
     Th-232: те же переключатели, зум протяжкой, наведение на реперы) ──── */
  function supNum(n) {
    var d = {"-":"⁻","0":"⁰","1":"¹","2":"²","3":"³","4":"⁴","5":"⁵",
             "6":"⁶","7":"⁷","8":"⁸","9":"⁹"};
    return String(n).split("").map(function (c) { return d[c] || c; }).join("");
  }

  function buildCal() {
    var tbl = document.getElementById("tblCal");
    if (tbl && !tbl.dataset.built) {
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
               + "<th>образец (Ra-226)</th><th>фон той же геометрии</th></tr></thead>";
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
    }
    var fw = D.fwhm_cal;
    var fwcs = document.getElementById("cal-fwcs");
    var fw662 = document.getElementById("cal-fw662");
    if (fwcs) fwcs.textContent = num(fw.fwhm662_cs, 1) + " кэВ";
    if (fw662) fw662.textContent = num(fw.fwhm662_law, 1) + " кэВ";
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
        var smp = D.spectrum.counts[i], bgv2 = D.spectrum.bg_counts[i];
        ro.textContent = num(e[i], 0) + " кэВ · образец " + cnt(smp)
          + " · фон " + cnt(bgv2) + " · разность " + cnt(smp - bgv2);
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
        if (Math.abs(x1 - x0) < 6) { drawCal(); return; }
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
        + "<td class='num'>" + num(q.res_pct, 2) + " %</td>"
        + "<td class='num'>" + q.n_lines_window + "</td>"
        + "<td class='num'>" + num(q.fwhm_model_keV, 2) + "</td>"
        + "<td class='num'>" + (q.dev_pct >= 0 ? "+" : "−")
        + num(Math.abs(q.dev_pct), 1) + " %</td>"
        + "<td>в подгонке</td></tr>";
    });
    body += "<tr class='sum'><td>степенной закон</td>"
      + "<td class='num' colspan='2'>ПШПВ = " + num(fw.k, 3) + "·E<sup>"
      + num(fw.p, 4) + "</sup></td>"
      + "<td class='num'>" + num(fw.res662_pct, 2) + " % на 662</td>"
      + "<td class='num'>" + fw.n_used + " из " + fw.n_anchors + "</td>"
      + "<td class='num'>" + num(fw.fwhm662_law, 1) + " кэВ</td>"
      + "<td class='num'>СКО " + num(fw.rms_dev_pct, 1) + " %</td>"
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

  function fillHeader() {
    var pass = D.passport;
    document.getElementById("p-aksp").textContent = cnt(pass.Bq_per_kg) + " Бк/кг";
    document.getElementById("p-uncpct").textContent = "±" + num(pass.unc_pct, 0) + " %";
    document.getElementById("p-datepass").textContent = pass.date_certified;
    document.getElementById("p-mass").textContent = num(pass.mass_g, 0) + " г";
    document.getElementById("p-apass").textContent = cnt(pass.A_Bq) + " Бк";
    document.getElementById("p-apassdev").textContent = "± " + cnt(pass.dA_Bq) + " Бк";
  }

  var VIEW_ID = { m2: "viewM2", cmp: "viewCmp", cal: "viewCal" };
  function switchTab(name) {
    Object.keys(VIEW_ID).forEach(function (t) {
      document.getElementById(VIEW_ID[t]).hidden = t !== name;
    });
    document.querySelectorAll("#tabs .tab").forEach(function (b) {
      b.setAttribute("aria-selected", String(b.dataset.tab === name));
    });
    if (name === "m2") { buildLegend(); fillSummary(); fillTable(); drawSpectrum(); }
    if (name === "cmp") fillCompare();
    if (name === "cal") { buildCal(); drawCal(); drawFwhm(); }
  }

  document.querySelectorAll("#tabs .tab").forEach(function (b) {
    b.addEventListener("click", function () { switchTab(b.dataset.tab); });
  });
  document.querySelectorAll("#libseg .btn").forEach(function (b) {
    b.addEventListener("click", function () {
      ST.lib = b.dataset.lib;
      document.querySelectorAll("#libseg .btn").forEach(function (x) {
        x.setAttribute("aria-pressed", String(x === b));
      });
      fillSummary(); fillTable(); drawSpectrum();
    });
  });

  window.addEventListener("resize", function () {
    if (!document.getElementById("viewM2").hidden) drawSpectrum();
    if (!document.getElementById("viewCal").hidden) { drawCal(); drawFwhm(); }
  });

  fillHeader();
  switchTab("m2");
})();
