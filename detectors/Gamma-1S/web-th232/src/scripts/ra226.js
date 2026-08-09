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

  function drawCalibration() {
    var cv = document.getElementById("cvCal");
    if (!cv) return;
    var p = pal(), f = fit(cv), g = f.g, W = f.w, H = f.h;
    var m = { l: 60, r: 14, t: 12, b: 34 };
    var e = D.spectrum.e_of_ch, yy = D.spectrum.counts, bgv = D.spectrum.bg_counts;
    var xLo = 0, xHi = e[e.length - 1];
    var vMax = 1;
    for (var i0 = 0; i0 < e.length; i0++) {
      if (e[i0] < 30 || e[i0] > 2400) continue;
      if (yy[i0] > vMax) vMax = yy[i0];
    }
    var Y = makeY(true, 0.5, vMax * 2.0);
    var x0 = m.l, x1 = W - m.r, y0 = m.t, y1 = H - m.b, n = e.length;

    g.strokeStyle = p.grid; g.lineWidth = 1; g.beginPath();
    [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250].forEach(function (tv) {
      var x = mapX(tv, xLo, xHi, x0, x1); g.moveTo(x, y0); g.lineTo(x, y1);
    });
    g.stroke();

    function line(vec, color, alpha, width) {
      g.strokeStyle = color; g.globalAlpha = alpha; g.lineWidth = width;
      g.beginPath();
      for (var i = 0; i < n; i++) {
        var x = mapX(e[i], xLo, xHi, x0, x1);
        var y = Y.map(vec[i], y0, y1);
        if (i === 0) g.moveTo(x, y); else g.lineTo(x, y);
      }
      g.stroke(); g.globalAlpha = 1;
    }
    line(bgv, p.faint, 0.9, 1);
    line(yy, p.ink, 0.9, 1);

    // реперные линии — вертикальные метки
    g.strokeStyle = p.sum; g.lineWidth = 1;
    (D.reference_lines || []).forEach(function (rl) {
      var E = rl[0];
      if (E < xLo || E > xHi) return;
      var x = mapX(E, xLo, xHi, x0, x1);
      g.globalAlpha = 0.55;
      g.beginPath(); g.moveTo(x, y0); g.lineTo(x, y1); g.stroke();
      g.globalAlpha = 1;
    });

    g.fillStyle = p.faint; g.font = "11px var(--mono, monospace)";
    g.textAlign = "center"; g.textBaseline = "top";
    [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250].forEach(function (tv) {
      var x = mapX(tv, xLo, xHi, x0, x1);
      g.fillText(String(tv), x, y1 + 6);
    });
    g.strokeStyle = p.rule; g.lineWidth = 1.2;
    g.beginPath(); g.moveTo(x0, y0); g.lineTo(x0, y1); g.lineTo(x1, y1); g.stroke();
  }

  function fillFwhmTable() {
    var tbl = document.getElementById("tblFwhm");
    if (!tbl) return;
    var fw = D.fwhm_cal;
    tbl.innerHTML = "<thead><tr><th>параметр</th><th>значение</th></tr></thead><tbody>" +
      "<tr><td>ПШПВ(E) = k·E<sup>p</sup></td><td>k=" + num(fw.k, 4) +
      ", p=" + num(fw.p, 4) + "</td></tr>" +
      "<tr><td>СКО отклонения</td><td>" + num(fw.rms_dev_pct, 1) + " %</td></tr>" +
      "<tr><td>опорных линий</td><td>" + (fw.n_used || "—") + "</td></tr></tbody>";
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
    if (name === "cal") { fillFwhmTable(); drawCalibration(); }
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
    if (!document.getElementById("viewCal").hidden) drawCalibration();
  });

  fillHeader();
  switchTab("m2");
})();
