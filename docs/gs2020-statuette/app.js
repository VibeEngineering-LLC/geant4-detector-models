 
(function () {
  "use strict";
  var D = window.G1S;
  var ST = { on: {}, log: true, sum: true, cursorE: null, zoom: null, drag: null, win: "main" };
  D.nuclides.forEach(function (n) { ST.on[n.key] = true; });
  function SPEC() { return D.spectrum; }
  function STACK1() { return D.fits[ST.win].stack; }
  function STACK2() { return STACK1(); }
  function TRUST_MASK() { return null; }
  function syncLegends() {}
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
  function numg(x, d) { var p = num(x, d).split(","); return cnt(p[0]) + "," + p[1]; }
  function signedPct(ratio) {
    var s = 100 * (ratio - 1);
    return (s < 0 ? "−" : "+") + num(Math.abs(s), 1) + " %";
  }
  function esc(s) {
    return String(s === undefined || s === null ? "" : s)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }
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
      ink:   css("--ink",   "#16140f"),
      rule:  css("--rule",  "#16140f"),
      faint: css("--faint", "#6a6558"),
      grid:  css("--grid",  "#dcd7c8"),
      paper: css("--paper", "#f5f2ea"),
      sum:   css("--sum-line", "#d21f1f"),
    };
  }
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
  function stackTotal(stk, i) {
    var acc = 0;
    for (var j = 0; j < D.nuclides.length; j++) {
      var k = D.nuclides[j].key;
      if (!ST.on[k] || !stk[k]) continue;
      acc += stk[k][i];
    }
    return acc;
  }
  function drawSpectrum(cvId, stk, extra) {
    var cv = document.getElementById(cvId);
    if (!cv) return null;
    var p = pal();
    var f = fit(cv);
    var g = f.g, W = f.w, H = f.h;
    var m = { l: 62, r: 14, t: 12, b: 34 };
    var e = D.spectrum.e_of_ch;
    var yy = D.spectrum.counts;
    var xLo = ST.zoom ? ST.zoom.xLo : 25, xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];
    var vMax = 1;
    for (var i0 = 0; i0 < e.length; i0++) {
      if (e[i0] < 30 || e[i0] > 3500) continue;
      var v0 = Math.max(yy[i0], stackTotal(stk, i0));
      if (v0 > vMax) vMax = v0;
    }
    var Y = makeY(ST.log, ST.log ? 0.5 : 0, vMax * (ST.log ? 2.0 : 1.1));
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
    g.save(); g.beginPath(); g.rect(m.l, m.t, W - m.r - m.l, H - m.b - m.t); g.clip();
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
    function edgeLo(i) {
      return mapX(i > 0 ? (e[i - 1] + e[i]) / 2 : e[0], xLo, xHi, m.l, W - m.r);
    }
    function edgeHi(i) {
      return mapX(i + 1 < e.length ? (e[i] + e[i + 1]) / 2 : e[e.length - 1],
                  xLo, xHi, m.l, W - m.r);
    }
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
      fillRuns(vec, [[0, e.length - 1]], order[oi].nuc.color);
      strokeLayer(oi);
    }
    function strokeLayer(oj) {
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
    if (ST.sum) {
      trace(function (i) { return stackTotal(stk, i); }, p.sum, 1.8, [6, 3]);
    }
    (extra || []).forEach(function (ex) {
      if (!ex.arr) return;
      trace(function (i) { return ex.arr[i]; }, ex.color, ex.width || 1.4,
            ex.dash);
    });
    trace(function (i) { return yy[i]; }, p.ink, 0.8);
    if (ST.cursorE !== null) {
      var xC = mapX(ST.cursorE, xLo, xHi, m.l, W - m.r);
      g.strokeStyle = p.rule; g.lineWidth = 1; g.setLineDash([4, 4]);
      g.beginPath(); g.moveTo(xC, m.t); g.lineTo(xC, H - m.b); g.stroke();
      g.setLineDash([]);
    }
    if (ST.drag) {
      var xa3 = mapX(Math.min(ST.drag.e0, ST.drag.e1), xLo, xHi, m.l, W - m.r);
      var xb3 = mapX(Math.max(ST.drag.e0, ST.drag.e1), xLo, xHi, m.l, W - m.r);
      g.fillStyle = "rgba(246,211,28,.22)";
      g.fillRect(xa3, m.t, xb3 - xa3, H - m.b - m.t);
      g.strokeStyle = "#16140f"; g.lineWidth = 1.5; g.setLineDash([4, 4]);
      g.strokeRect(xa3, m.t, xb3 - xa3, H - m.b - m.t);
      g.setLineDash([]);
    }
    g.restore();
    g.strokeStyle = p.rule; g.lineWidth = 2;
    g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);
    return { m: m, W: W, H: H, xLo: xLo, xHi: xHi, Y: Y };
  }
  function nucOf(key) {
    for (var i = 0; i < D.nuclides.length; i++)
      if (D.nuclides[i].key === key) return D.nuclides[i];
    return null;
  }
  function noiseFrac(key) {
    var nf = SPEC().noise_frac || D.spectrum.noise_frac;
    if (!nf || !(key in nf)) return null;
    return nf[key];
  }
  function guardHint(nf) {
    return "статистики МК-шаблона достаточно для " + num(100 * (1 - nf), 1) + " %"
      + " вклада слоя; остальное — области ниже "
      + (D.spectrum.n_eff_min || 0) + " отсчётов МК на канал, там слой показан "
      + "пунктиром без заливки (доля нуклида определяется шумом шаблона)";
  }
  function buildLegend(elId) {
    var el = document.getElementById(elId);
    if (!el || el.dataset.built) return;
    el.dataset.built = "1";
    var html = "";
    var stkL = (elId === 'legendM2') ? STACK2() : STACK1();
    var arrL = D.nuclides.map(function (n, i) { var v = stkL && stkL[n.key], s = 0;
      if (v) for (var q = 0; q < v.length; q++) s += v[q]; return { n: n, s: s, i: i }; });
    arrL.sort(function (a, b) { return (b.s - a.s) || (a.i - b.i); });
    arrL.map(function (o) { return o.n; }).forEach(function (nuc) {
      if (nuc.key === "SECOND") return;
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
         + "<span class='sw log'></span><span class='nm'>логарифмическая шкала</span></label>";
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
  function cell(lab, val, big, hint) {
    return "<div" + (hint ? " title='" + hint + "'" : "") + "><span class='lab'>"
         + lab + "</span><span class='val"
         + (big ? " big-num" : "") + "'>" + val + "</span></div>";
  }
  function cursorText(elId, stk, tipId) {
    var el = document.getElementById(elId);
    if (!el) return;
    if (ST.cursorE === null) {
      el.textContent = "наведите курсор на канал спектра";
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
    var txt = num(e[i], 0) + " кэВ — измерено " + cnt(meas)
            + ", модель " + cnt(stackTotal(stk, i));
    if (top) txt += " — " + top;
    el.textContent = txt;
    var tip = document.getElementById(tipId);
    if (tip) tip.textContent = num(e[i], 0) + " кэВ · " + cnt(meas);
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
      var xLo = ST.zoom ? ST.zoom.xLo : 25, xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];
      if (x < m.l || x > r.width - m.r) ST.cursorE = null;
      else ST.cursorE = xLo + ((x - m.l) / (r.width - m.r - m.l)) * (xHi - xLo);
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
function zEfromX(x, rectWidth) {
  var e = D.spectrum.e_of_ch;
  var xLo = ST.zoom ? ST.zoom.xLo : 25;
  var xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];
  var m = { l: 62, r: 14 };
  return xLo + ((x - m.l) / (rectWidth - m.r - m.l)) * (xHi - xLo);
}
function wireZoom(cvId) {
  var cv = document.getElementById(cvId);
  if (!cv) return;
  var dragging = false;
  var dragState = null;
  var m = { l: 62, r: 14 };
  cv.addEventListener("mousedown", function (ev) {
    var r = cv.getBoundingClientRect();
    var x = ev.clientX - r.left;
    if (x < m.l || x > r.width - m.r) return;
    ev.preventDefault();
    dragging = true;
    dragState = { x0: x, x1: x, w: r.width };
    ST.drag = { e0: zEfromX(x, r.width), e1: zEfromX(x, r.width) };
    redraw();
  });
  document.addEventListener("mousemove", function (ev) {
    if (!dragging) return;
    var r = cv.getBoundingClientRect();
    var x = ev.clientX - r.left;
    if (x < m.l) x = m.l;
    if (x > r.width - m.r) x = r.width - m.r;
    dragState.x1 = x;
    ST.drag.e1 = zEfromX(x, r.width);
    redraw();
  });
  document.addEventListener("mouseup", function () {
    if (!dragging) return;
    dragging = false;
    var x0 = dragState.x0, x1 = dragState.x1, w = dragState.w;
    dragState = null;
    ST.drag = null;
    if (Math.abs(x1 - x0) < 6) { redraw(); return; }
    ST.zoom = { xLo: Math.max(0, zEfromX(Math.min(x0, x1), w)),
                xHi: zEfromX(Math.max(x0, x1), w) };
    redraw();
  });
  cv.addEventListener("dblclick", function () {
    ST.zoom = null;
    redraw();
  });
}
  function redraw() { drawSpectrum("cvM1", STACK1(), []); cursorText("cursorM1", STACK1(), "m1-tip"); }
  function refresh() {
    var f = D.fits[ST.win];
    document.getElementById("sumM1").innerHTML = f.sum_html;
    document.getElementById("tblM1").innerHTML = f.tbl_html;
    redraw();
  }
  document.querySelectorAll("#tabs .tab").forEach(function (b) {
    b.addEventListener("click", function () {
      ST.win = b.getAttribute("data-win");
      document.querySelectorAll("#tabs .tab").forEach(function (x) { x.setAttribute("aria-selected", x === b ? "true" : "false"); });
      refresh();
    });
  });
  buildLegend("legendM1");
  document.querySelectorAll("#legendM1 label.chip[data-nuc]").forEach(function (lb) {   
    var k = lb.getAttribute("data-nuc");
    if ((D.undetected || []).indexOf(k) >= 0) { ST.on[k] = false; lb.remove(); }
  });
  document.querySelectorAll("#legendM1 .c-sum").forEach(function (inp) { inp.parentNode.querySelector(".nm").textContent = "сумма отмеченных слоёв"; });
  attachCursor("cvM1", "m1-tip", redraw); wireZoom("cvM1");
  window.addEventListener("resize", redraw); refresh();
  function closeLB() { document.querySelectorAll(".lbox").forEach(function (d) { d.remove(); }); }
  document.querySelectorAll(".dev-fig img").forEach(function (im) { im.addEventListener("click", function () {
    var d = document.createElement("div"); d.className = "lbox"; d.addEventListener("click", closeLB);
    d.innerHTML = "<img src='" + im.src + "' alt=''><p>" + esc(im.alt) + "</p>"; document.body.appendChild(d);
  }); });
  document.addEventListener("keydown", function (e) { if (e.key === "Escape") closeLB(); });
})();
