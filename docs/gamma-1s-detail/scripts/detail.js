(function () { "use strict";

function init() {
  var D = window.G1SD, T = window.G1SD_TEXT;
  if (!D || !T) return;

  /* ---- 3. Number formatting ---- */
  function fmt(x, d) {
    if (x === null || x === undefined || (typeof x === "number" && isNaN(x))) return "\u2014";
    var s = x.toFixed(d);
    s = s.replace(".", ",").replace("-", "\u2212");
    return s;
  }
  function signed(x, d) {
    if (x === null || x === undefined || (typeof x === "number" && isNaN(x))) return "\u2014";
    var raw = x.toFixed(d);
    if (parseFloat(raw.replace("-", "")) === 0) return (0).toFixed(d).replace(".", ",");
    var s = raw.replace(".", ",").replace("-", "\u2212");
    if (x > 0) s = "+" + s;
    return s;
  }
  function rel(r) { return signed((r - 1) * 100, 1); }
  function thin(n) {
    n = Math.round(n);
    var neg = n < 0;
    var s = String(Math.abs(n));
    var parts = [];
    while (s.length > 3) { parts.unshift(s.slice(-3)); s = s.slice(0, -3); }
    parts.unshift(s);
    return (neg ? "\u2212" : "") + parts.join("\u202F");
  }
  function sup(n) {
    var m = {"0":"\u2070","1":"\u00B9","2":"\u00B2","3":"\u00B3","4":"\u2074","5":"\u2075","6":"\u2076","7":"\u2077","8":"\u2078","9":"\u2079","-":"\u207B"};
    var s = String(n), out = "";
    for (var i = 0; i < s.length; i++) out += m[s[i]] || s[i];
    return out;
  }

  /* ---- 4. Path resolver ---- */
  function get(path) {
    var parts = path.split("."), obj = window.G1SD;
    for (var i = 0; i < parts.length; i++) {
      if (obj === null || obj === undefined) return undefined;
      var p = parts[i];
      if (obj instanceof Array) {
        var idx = parseInt(p, 10);
        obj = (idx === idx && idx >= 0) ? obj[idx] : obj[p];
      } else {
        obj = obj[p];
      }
    }
    return obj;
  }

  /* ---- 1. Header and toolbar ---- */
  document.title = T.title;
  var hdr = document.getElementById("hdr");
  hdr.innerHTML = '<p class="eyebrow">' + T.eyebrow + '</p><h1>' + T.h1 + '</h1>' + T.stand;

  var tb = document.getElementById("toolbar"), th = "";
  for (var i = 0; i < T.sections.length; i++) {
    th += '<a class="btn" href="#' + T.sections[i].id + '">' + T.sections[i].h + '</a>';
  }
  th += '<span class="sp"></span>';
  th += '<a class="btn" href="https://github.com/VibeEngineering-LLC/geant4-detector-models/tree/main/detectors/Gamma-1S" target="_blank" rel="noopener">код на GitHub</a>';
  tb.innerHTML = th;

  /* ---- 5. Table helpers ---- */
  function tableWrap(inner) {
    return '<div class="scroll"><table class="big">' + inner + '</table></div>';
  }

  /* ---- 5a. tbl-fits ---- */
  function buildFits(el) {
    var rows = ["k24g1t", "dtF", "dtP"];
    var h = "<thead><tr><th>Шаблоны</th><th>A(Ra-228)/паспорт</th><th>Δ к полной, %</th><th>A(Th-228)/паспорт</th><th>Δ к полной, %</th><th>χ²_ref/ν</th></tr></thead>";
    var b = "<tbody>";
    for (var i = 0; i < rows.length; i++) {
      var p = rows[i], f = D.fits[p];
      if (!f) continue;
      var lab = T.labels[p] || p;
      var raCell = fmt(f.ra, 3) + " ± " + fmt(f.ra * f.dra, 3);
      var thCell = fmt(f.th, 3) + " ± " + fmt(f.th * f.dth, 3);
      var chi = fmt(f.chi2, 2);
      var dRa, dTh;
      if (p === "dtF") { dRa = "\u2014"; dTh = "\u2014"; }
      else {
        var fd = D.fits["dtF"];
        dRa = signed((f.ra / fd.ra - 1) * 100, 1);
        dTh = signed((f.th / fd.th - 1) * 100, 1);
      }
      b += "<tr><td>" + lab + "</td><td class=\"num\">" + raCell + "</td><td class=\"num\">" + dRa + "</td><td class=\"num\">" + thCell + "</td><td class=\"num\">" + dTh + "</td><td class=\"num\">" + chi + "</td></tr>";
    }
    b += "</tbody>";
    el.innerHTML = tableWrap(h + b);
  }

  /* ---- 5b. tbl-nodes ---- */
  function buildNodes(el) {
    var rows = ["dtF", "dtH", "dtE1", "dtE2", "dtE3", "dtE4", "dtP"];
    var h = "<thead><tr><th>Вариант (от полной модели)</th><th>A(Ra-228)/паспорт</th><th>Δ, %</th><th>A(Th-228)/паспорт</th><th>Δ, %</th><th>χ²_ref/ν</th><th>пики, %</th><th>50–260 кэВ, %</th></tr></thead>";
    var b = "<tbody>";
    for (var i = 0; i < rows.length; i++) {
      var p = rows[i], f = D.fits[p];
      if (!f) continue;
      var lab = T.labels[p] || p;
      var raCell = fmt(f.ra, 3), thCell = fmt(f.th, 3), chi = fmt(f.chi2, 2);
      var dRa, dTh, pk, bd;
      if (p === "dtF") { dRa = "\u2014"; dTh = "\u2014"; pk = "\u2014"; bd = "\u2014"; }
      else {
        var fd = D.fits["dtF"];
        dRa = signed((f.ra / fd.ra - 1) * 100, 1);
        dTh = signed((f.th / fd.th - 1) * 100, 1);
        // peaks
        var pkVals = [];
        for (var c = 0; c < D.chains.length; c++) {
          var pd = D.peaks[D.chains[c]][p];
          if (pd && pd.peaks) { for (var k = 0; k < pd.peaks.length; k++) { if (pd.peaks[k].r != null) pkVals.push((pd.peaks[k].r - 1) * 100); } }
        }
        pk = pkVals.length ? signed(Math.min.apply(null, pkVals), 1) + "\u2026" + signed(Math.max.apply(null, pkVals), 1) : "\u2014";
        // bands lo 50 or 150
        var bdVals = [];
        for (var c2 = 0; c2 < D.chains.length; c2++) {
          var pd2 = D.peaks[D.chains[c2]][p];
          if (pd2 && pd2.bands) { for (var k2 = 0; k2 < pd2.bands.length; k2++) { var bn = pd2.bands[k2]; if ((bn.lo === 50 || bn.lo === 150) && bn.r != null) bdVals.push((bn.r - 1) * 100); } }
        }
        bd = bdVals.length ? signed(Math.min.apply(null, bdVals), 1) + "\u2026" + signed(Math.max.apply(null, bdVals), 1) : "\u2014";
      }
      b += "<tr><td>" + lab + "</td><td class=\"num\">" + raCell + "</td><td class=\"num\">" + dRa + "</td><td class=\"num\">" + thCell + "</td><td class=\"num\">" + dTh + "</td><td class=\"num\">" + chi + "</td><td class=\"num\">" + pk + "</td><td class=\"num\">" + bd + "</td></tr>";
    }
    b += "</tbody>";
    el.innerHTML = tableWrap(h + b);
  }

  /* ---- 5c. tbl-peaks ---- */
  function buildPeaks(el) {
    var html = "";
    for (var c = 0; c < D.chains.length; c++) {
      var chain = D.chains[c];
      var cap = T.chain_labels[chain] || chain;
      html += '<p class="tbl-cap">' + cap + '</p>';
      var rows = D.variants.concat([D.control]);
      // get peak/band structure from first available prefix
      var refP = null;
      for (var ri = 0; ri < rows.length; ri++) { if (D.peaks[chain][rows[ri]]) { refP = rows[ri]; break; } }
      if (!refP) continue;
      var refD = D.peaks[chain][refP];
      var peaks = refD.peaks || [], bands = refD.bands || [];
      var h = "<thead><tr><th>Вариант</th>";
      for (var pi = 0; pi < peaks.length; pi++) h += "<th>" + fmt(peaks[pi].E, 1) + " кэВ</th>";
      for (var bi = 0; bi < bands.length; bi++) h += "<th>" + bands[bi].lo + "–" + bands[bi].hi + " кэВ</th>";
      h += "</tr></thead>";
      var b = "<tbody>";
      for (var r = 0; r < rows.length; r++) {
        var p = rows[r], ld = D.peaks[chain][p];
        if (!ld) continue;
        var lab = T.labels[p] || p;
        b += "<tr><td>" + lab + "</td>";
        var ldPeaks = ld.peaks || [];
        for (var pi2 = 0; pi2 < peaks.length; pi2++) {
          var pv = (pi2 < ldPeaks.length) ? ldPeaks[pi2] : null;
          if (pv) b += "<td class=\"num\">" + fmt(pv.r, 3) + " ± " + fmt(pv.s, 3) + "</td>";
          else b += "<td class=\"num\">\u2014</td>";
        }
        var ldBands = ld.bands || [];
        for (var bi2 = 0; bi2 < bands.length; bi2++) {
          var bv = (bi2 < ldBands.length) ? ldBands[bi2] : null;
          if (bv && bv.r != null) b += "<td class=\"num\">" + fmt(bv.r, 3) + "</td>";
          else b += "<td class=\"num\">\u2014</td>";
        }
        b += "</tr>";
      }
      b += "</tbody>";
      html += tableWrap(h + b);
    }
    el.innerHTML = html;
  }

  /* ---- 6. Charts ---- */
  function buildCharts(el) {
    var COL = {
      dtF: "#16140f", dtP: "#d1495b", dtH: "#edae49", dtE1: "#66a182",
      dtE2: "#00798c", dtE3: "#8d5a97", dtE4: "#30638e", k24g1t: "#9a9a9a"
    };
    // read --ink for dtF in dark theme
    var cs = getComputedStyle(document.documentElement);
    var inkVar = cs.getPropertyValue("--ink").trim();
    if (inkVar) COL.dtF = inkVar;
    var inkColor = inkVar || "#16140f";
    var bgGrey = "#9a9a9a";

    var prefixes = ["dtF"].concat(D.variants);
    var M = {l: 62, r: 14, t: 12, b: 34};
    var bodyFont = cs.getPropertyValue("font-family").trim() || "sans-serif";
    var ruleSoft = cs.getPropertyValue("--rule-soft").trim() || "#ccc";
    var sunk = cs.getPropertyValue("--sunk").trim() || "#ddd";
    var dim = cs.getPropertyValue("--dim").trim() || "#888";
    var accentSoft = cs.getPropertyValue("--accent-soft").trim() || "rgba(100,150,200,0.3)";

    // Build DOM
    var ctrl = document.createElement("div"); ctrl.className = "ctrl";
    var seg = document.createElement("div"); seg.className = "seg"; seg.id = "rangeSeg";
    var ranges = [["30-3000","весь спектр"],["30-300","30–300 кэВ"],["300-1000","300–1000 кэВ"],["1000-2800","1–2,8 МэВ"]];
    for (var ri = 0; ri < ranges.length; ri++) {
      var btn = document.createElement("button");
      btn.className = "btn"; btn.setAttribute("data-range", ranges[ri][0]);
      btn.textContent = ranges[ri][1]; btn.setAttribute("aria-pressed", "false");
      seg.appendChild(btn);
    }
    ctrl.appendChild(seg);

    var chips = document.createElement("div"); chips.className = "chips"; chips.id = "varChips";
    // fixed entries
    var f1 = document.createElement("span"); f1.className = "chip fixed";
    f1.innerHTML = '<span class="sw" style="background:' + inkColor + '"></span><span>измерение</span>';
    chips.appendChild(f1);
    // checkbox entries
    for (var ci = 0; ci < prefixes.length; ci++) {
      var pfx = prefixes[ci];
      var lab = T.labels[pfx] || pfx;
      var lbl = document.createElement("label"); lbl.className = "chip";
      var inp = document.createElement("input"); inp.type = "checkbox"; inp.setAttribute("data-p", pfx);
      inp.checked = (pfx === "dtF" || pfx === "dtP");
      var sw = document.createElement("span"); sw.className = "sw"; sw.style.background = COL[pfx] || "#888";
      var txt = document.createElement("span"); txt.textContent = lab;
      lbl.appendChild(inp); lbl.appendChild(sw); lbl.appendChild(txt);
      chips.appendChild(lbl);
    }
    ctrl.appendChild(chips);
    el.appendChild(ctrl);

    var box1 = document.createElement("div"); box1.className = "chartbox";
    var cv1 = document.createElement("canvas"); cv1.id = "cvSpec"; cv1.setAttribute("height", "300");
    box1.appendChild(cv1); el.appendChild(box1);
    var box2 = document.createElement("div"); box2.className = "chartbox";
    var cv2 = document.createElement("canvas"); cv2.id = "cvRatio"; cv2.setAttribute("height", "260");
    box2.appendChild(cv2); el.appendChild(box2);

    var cursorLine = document.createElement("div"); cursorLine.className = "cursor-line"; cursorLine.id = "cursorLine";
    cursorLine.textContent = "наведи курсор на график";
    el.appendChild(cursorLine);

    // State
    var S = {zoom: null, drag: null, dragging: false, cursorE: null, on: {dtF: true, dtP: true, dtH: false, dtE1: false, dtE2: false, dtE3: false, dtE4: false}};

    // Data reference
    var meas = D.meas;
    if (!meas) { window.G1SD_UI = {state: S, redraw: function(){}}; return; }
    var eArr = meas.e, counts = meas.counts, bgArr = meas.bg, models = meas.model;
    var fitLo = meas.lo_keV || 50, fitHi = meas.hi_keV || 260;
    var xLoDefault = 30, xHiDefault = 3000;

    function curXLo() { return S.zoom ? S.zoom.xLo : xLoDefault; }
    function curXHi() { return S.zoom ? S.zoom.xHi : xHiDefault; }

    function eFromX(x, w) {
      var lo = curXLo(), hi = curXHi();
      return lo + (x - M.l) / (w - M.l - M.r) * (hi - lo);
    }

    function cssW(cv) { return cv.parentElement ? cv.parentElement.clientWidth : (cv.width / (window.devicePixelRatio || 1)); }
    function cssH(cv) { return cv.id === "cvSpec" ? 300 : 260; }

    function setupCanvas(cv) {
      var dpr = window.devicePixelRatio || 1;
      var w = cssW(cv) || 800, h = cssH(cv);
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(h * dpr);
      cv.style.width = "100%";
      cv.style.height = h + "px";
      var ctx = cv.getContext("2d");
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      return {ctx: ctx, w: w, h: h};
    }

    function niceStep(range) {
      var rough = range / 7;
      var mag = Math.pow(10, Math.floor(Math.log(rough) / Math.LN10));
      var norm = rough / mag;
      var nice;
      if (norm < 1.5) nice = 1; else if (norm < 3) nice = 2; else if (norm < 7) nice = 5; else nice = 10;
      return nice * mag;
    }

    function xTicks(lo, hi) {
      var step = niceStep(hi - lo);
      var start = Math.ceil(lo / step) * step;
      var ticks = [];
      for (var v = start; v <= hi + step * 0.01; v += step) {
        if (v >= lo - step * 0.01 && v <= hi + step * 0.01) ticks.push(v);
      }
      return ticks;
    }

    function drawAll() {
      // --- Top canvas: cvSpec (log y) ---
      var s1 = setupCanvas(cv1);
      var ctx1 = s1.ctx, W1 = s1.w, H1 = s1.h;
      ctx1.clearRect(0, 0, W1, H1);
      var lo = curXLo(), hi = curXHi();
      var plotW = W1 - M.l - M.r, plotH = H1 - M.t - M.b;
      function xP(E) { return M.l + (E - lo) / (hi - lo) * plotW; }

      // Determine visible indices
      var iLo = 0, iHi = eArr.length - 1;
      for (var vi = 0; vi < eArr.length; vi++) { if (eArr[vi] >= lo) { iLo = vi; break; } }
      for (var vj = eArr.length - 1; vj >= 0; vj--) { if (eArr[vj] <= hi) { iHi = vj; break; } }

      // Y range (positive values)
      var yMin = Infinity, yMax = 0;
      for (var yi = iLo; yi <= iHi; yi++) {
        if (counts[yi] > 0) { if (counts[yi] < yMin) yMin = counts[yi]; if (counts[yi] > yMax) yMax = counts[yi]; }
      }
      // include models
      var modelPrefixes = [];
      for (var mp = 0; mp < prefixes.length; mp++) { if (S.on[prefixes[mp]]) modelPrefixes.push(prefixes[mp]); }
      for (var mi = 0; mi < modelPrefixes.length; mi++) {
        var md = models[modelPrefixes[mi]];
        if (!md) continue;
        for (var yi2 = iLo; yi2 <= iHi; yi2++) {
          if (md[yi2] > 0) { if (md[yi2] < yMin) yMin = md[yi2]; if (md[yi2] > yMax) yMax = md[yi2]; }
        }
      }
      if (yMin === Infinity) { yMin = 1; yMax = 10; }
      var lin = !!S.lin;
      yMin = lin ? 0 : Math.max(0.1, yMin * 0.7);
      yMax = yMax * (lin ? 1.05 : 1.4);
      var logMin = Math.log(yMin) / Math.LN10, logMax = Math.log(yMax) / Math.LN10;
      function yP(v) {
        if (v <= 0 && !lin) return -1;
        if (lin) return H1 - M.b - (v - yMin) / (yMax - yMin) * plotH;
        var lv = Math.log(v) / Math.LN10;
        return H1 - M.b - (lv - logMin) / (logMax - logMin) * plotH;
      }

      // Grid and y-axis (powers of 10)
      ctx1.font = "11px " + bodyFont;
      ctx1.fillStyle = inkColor; ctx1.strokeStyle = ruleSoft;
      var expMin = Math.ceil(logMin), expMax = Math.floor(logMax), yTicks = [];
      if (lin) { var ys = niceStep(yMax - yMin); for (var tv = 0; tv <= yMax; tv += ys) yTicks.push(tv); }
      else { for (var te = expMin; te <= expMax; te++) yTicks.push(Math.pow(10, te)); }
      for (var ex = 0; ex < yTicks.length; ex++) {
        var val = yTicks[ex];
        var yy = yP(val);
        if (yy < M.t - 2 || yy > H1 - M.b + 2) continue;
        ctx1.beginPath(); ctx1.moveTo(M.l, yy); ctx1.lineTo(W1 - M.r, yy); ctx1.stroke();
        var label = lin ? thin(Math.round(val)) : "10" + sup(Math.round(Math.log(val) / Math.LN10));
        ctx1.textAlign = "right"; ctx1.textBaseline = "middle";
        ctx1.fillText(label, M.l - 4, yy);
      }
      // y-axis label (rotated)
      ctx1.save();
      ctx1.translate(12, M.t + plotH / 2);
      ctx1.rotate(-Math.PI / 2);
      ctx1.textAlign = "center"; ctx1.textBaseline = "middle";
      ctx1.fillText("отсчёты в канале", 0, 0);
      ctx1.restore();

      // Shade outside fit window
      ctx1.save();
      ctx1.globalAlpha = 0.5; ctx1.fillStyle = sunk;
      var fitXL = xP(fitLo), fitXR = xP(fitHi);
      if (fitXL > M.l) ctx1.fillRect(M.l, M.t, Math.min(fitXL, W1 - M.r) - M.l, plotH);
      if (fitXR < W1 - M.r) ctx1.fillRect(Math.max(fitXR, M.l), M.t, W1 - M.r - Math.max(fitXR, M.l), plotH);
      ctx1.restore();
      // label "вне окна подгонки" if wider than 60px
      var outsideW = (fitXL > M.l ? fitXL - M.l : 0) + (fitXR < W1 - M.r ? W1 - M.r - fitXR : 0);
      if (outsideW > 60) {
        ctx1.save(); ctx1.globalAlpha = 0.7; ctx1.fillStyle = dim; ctx1.font = "11px " + bodyFont;
        ctx1.textAlign = "center"; ctx1.fillText("вне окна подгонки", M.l + (fitXL - M.l) / 2, M.t + 14);
        ctx1.restore();
      }

      // Model lines
      for (var mi2 = 0; mi2 < modelPrefixes.length; mi2++) {
        var mp2 = modelPrefixes[mi2];
        var md2 = models[mp2];
        if (!md2) continue;
        ctx1.beginPath();
        ctx1.strokeStyle = COL[mp2] || "#888";
        ctx1.lineWidth = (mp2 === "dtF") ? 2 : 1.25;
        var mStarted = false;
        for (var mi3 = iLo; mi3 <= iHi; mi3++) {
          if (md2[mi3] > 0) {
            var mx = xP(eArr[mi3]), my = yP(md2[mi3]);
            if (!mStarted) { ctx1.moveTo(mx, my); mStarted = true; } else ctx1.lineTo(mx, my);
          } else { mStarted = false; }
        }
        ctx1.stroke();
      }

      // Measured counts as dots
      ctx1.fillStyle = inkColor;
      for (var di = iLo; di <= iHi; di++) {
        if (counts[di] > 0) {
          var dx = xP(eArr[di]), dy = yP(counts[di]);
          ctx1.fillRect(dx - 1.6, dy - 1.6, 3.2, 3.2);
        }
      }

      // X ticks (top canvas - no axis label)
      var xT = xTicks(lo, hi);
      ctx1.fillStyle = inkColor; ctx1.textAlign = "center"; ctx1.textBaseline = "top";
      ctx1.strokeStyle = ruleSoft;
      for (var xi = 0; xi < xT.length; xi++) {
        var xx = xP(xT[xi]);
        ctx1.beginPath(); ctx1.moveTo(xx, H1 - M.b); ctx1.lineTo(xx, H1 - M.b + 4); ctx1.stroke();
        ctx1.fillText(String(Math.round(xT[xi])), xx, H1 - M.b + 6);
      }

      // Drag rectangle
      if (S.drag) {
        ctx1.save(); ctx1.globalAlpha = 0.35; ctx1.fillStyle = accentSoft;
        var dx0 = Math.min(S.drag.x0, S.drag.x1), dx1 = Math.max(S.drag.x0, S.drag.x1);
        ctx1.fillRect(dx0, M.t, dx1 - dx0, plotH);
        ctx1.restore();
      }
      // Cursor line
      if (S.cursorE !== null) {
        var cx1 = xP(S.cursorE);
        if (cx1 >= M.l && cx1 <= W1 - M.r) {
          ctx1.beginPath(); ctx1.strokeStyle = dim; ctx1.lineWidth = 0.5;
          ctx1.moveTo(cx1, M.t); ctx1.lineTo(cx1, H1 - M.b); ctx1.stroke();
        }
      }

      // --- Bottom canvas: cvRatio (linear y) ---
      var s2 = setupCanvas(cv2);
      var ctx2 = s2.ctx, W2 = s2.w, H2 = s2.h;
      ctx2.clearRect(0, 0, W2, H2);
      var plotW2 = W2 - M.l - M.r, plotH2 = H2 - M.t - M.b;
      function xP2(E) { return M.l + (E - lo) / (hi - lo) * plotW2; }

      // Compute ratios (smoothed 5-channel)
      var ratioPrefixes = [];
      for (var rp = 0; rp < prefixes.length; rp++) { if (S.on[prefixes[rp]]) ratioPrefixes.push(prefixes[rp]); }

      // Y range from visible finite ratios inside fit window
      var rMin = Infinity, rMax = -Infinity;
      for (var rpi = 0; rpi < ratioPrefixes.length; rpi++) {
        var rmd = models[ratioPrefixes[rpi]];
        if (!rmd) continue;
        for (var rci = iLo; rci <= iHi; rci++) {
          if (eArr[rci] < fitLo || eArr[rci] > fitHi) continue;
          var sumM = 0, sumC = 0;
          for (var k = -2; k <= 2; k++) {
            var idx = rci + k;
            if (idx >= 0 && idx < eArr.length) { sumM += rmd[idx] || 0; sumC += counts[idx] || 0; }
          }
          if (sumC > 0 && sumM > 0) {
            var ratio = sumM / sumC;
            if (ratio < rMin) rMin = ratio;
            if (ratio > rMax) rMax = ratio;
          }
        }
      }
      if (rMin === Infinity) { rMin = 0.95; rMax = 1.05; }
      rMin = Math.max(0.7, Math.min(rMin, 1.3));
      rMax = Math.min(1.3, Math.max(rMax, 0.7));
      if (rMax - rMin < 0.1) { var mid = (rMin + rMax) / 2; rMin = mid - 0.05; rMax = mid + 0.05; }
      var pad2 = (rMax - rMin) * 0.1;
      rMin -= pad2; rMax += pad2;
      function yP2(v) { return H2 - M.b - (v - rMin) / (rMax - rMin) * plotH2; }

      // Y grid and ticks
      ctx2.font = "11px " + bodyFont;
      var yStep2 = niceStep(rMax - rMin);
      var yStart2 = Math.ceil(rMin / yStep2) * yStep2;
      ctx2.fillStyle = inkColor; ctx2.strokeStyle = ruleSoft;
      ctx2.textAlign = "right"; ctx2.textBaseline = "middle";
      for (var yv = yStart2; yv <= rMax + yStep2 * 0.01; yv += yStep2) {
        if (yv < rMin - yStep2 * 0.01) continue;
        var yy2 = yP2(yv);
        if (yy2 < M.t || yy2 > H2 - M.b) continue;
        ctx2.beginPath(); ctx2.moveTo(M.l, yy2); ctx2.lineTo(W2 - M.r, yy2); ctx2.stroke();
        ctx2.fillText(fmt(yv, 2), M.l - 4, yy2);
      }
      // y-axis label
      ctx2.save();
      ctx2.translate(12, M.t + plotH2 / 2);
      ctx2.rotate(-Math.PI / 2);
      ctx2.textAlign = "center"; ctx2.fillText("модель / измерение", 0, 0);
      ctx2.restore();

      // Horizontal line at 1.0 (dashed)
      var yOne = yP2(1.0);
      ctx2.save(); ctx2.setLineDash([5, 4]); ctx2.strokeStyle = inkColor; ctx2.lineWidth = 1;
      ctx2.beginPath(); ctx2.moveTo(M.l, yOne); ctx2.lineTo(W2 - M.r, yOne); ctx2.stroke();
      ctx2.restore();

      // Uncertainty band (grey, alpha 0.25)
      ctx2.save(); ctx2.globalAlpha = 0.25; ctx2.fillStyle = bgGrey;
      var bandPath = [];
      for (var bci = iLo; bci <= iHi; bci++) {
        if (eArr[bci] < fitLo || eArr[bci] > fitHi) continue;
        var sumC2 = 0;
        for (var k2 = -2; k2 <= 2; k2++) {
          var idx2 = bci + k2;
          if (idx2 >= 0 && idx2 < eArr.length) sumC2 += counts[idx2] || 0;
        }
        if (sumC2 > 0) {
          var unc = 1 / Math.sqrt(sumC2);
          var bx2 = xP2(eArr[bci]);
          bandPath.push({x: bx2, yLo: yP2(1 - unc), yHi: yP2(1 + unc)});
        }
      }
      if (bandPath.length > 1) {
        ctx2.beginPath();
        ctx2.moveTo(bandPath[0].x, bandPath[0].yHi);
        for (var bpi = 1; bpi < bandPath.length; bpi++) ctx2.lineTo(bandPath[bpi].x, bandPath[bpi].yHi);
        for (var bpi2 = bandPath.length - 1; bpi2 >= 0; bpi2--) ctx2.lineTo(bandPath[bpi2].x, bandPath[bpi2].yLo);
        ctx2.closePath(); ctx2.fill();
      }
      ctx2.restore();

      // Ratio lines
      for (var rpi2 = 0; rpi2 < ratioPrefixes.length; rpi2++) {
        var rmp = ratioPrefixes[rpi2];
        var rmd2 = models[rmp];
        if (!rmd2) continue;
        ctx2.beginPath();
        ctx2.strokeStyle = COL[rmp] || "#888";
        ctx2.lineWidth = (rmp === "dtF") ? 2 : 1.25;
        var rStarted = false;
        for (var rci2 = iLo; rci2 <= iHi; rci2++) {
          if (eArr[rci2] < fitLo || eArr[rci2] > fitHi) { rStarted = false; continue; }
          var sumM2 = 0, sumC3 = 0;
          for (var k3 = -2; k3 <= 2; k3++) {
            var idx3 = rci2 + k3;
            if (idx3 >= 0 && idx3 < eArr.length) { sumM2 += rmd2[idx3] || 0; sumC3 += counts[idx3] || 0; }
          }
          if (sumC3 > 0 && sumM2 > 0) {
            var ratio2 = sumM2 / sumC3;
            var rx = xP2(eArr[rci2]), ry = yP2(ratio2);
            if (!rStarted) { ctx2.moveTo(rx, ry); rStarted = true; } else ctx2.lineTo(rx, ry);
          } else { rStarted = false; }
        }
        ctx2.stroke();
      }

      // X ticks and axis label
      var xT2 = xTicks(lo, hi);
      ctx2.fillStyle = inkColor; ctx2.textAlign = "center"; ctx2.textBaseline = "top";
      ctx2.strokeStyle = ruleSoft;
      for (var xi2 = 0; xi2 < xT2.length; xi2++) {
        var xx2 = xP2(xT2[xi2]);
        ctx2.beginPath(); ctx2.moveTo(xx2, H2 - M.b); ctx2.lineTo(xx2, H2 - M.b + 4); ctx2.stroke();
        ctx2.fillText(String(Math.round(xT2[xi2])), xx2, H2 - M.b + 6);
      }
      ctx2.textAlign = "center"; ctx2.fillText("энергия, кэВ", M.l + plotW2 / 2, H2 - 10);

      // Drag rectangle
      if (S.drag) {
        ctx2.save(); ctx2.globalAlpha = 0.35; ctx2.fillStyle = accentSoft;
        var dxx0 = Math.min(S.drag.x0, S.drag.x1), dxx1 = Math.max(S.drag.x0, S.drag.x1);
        ctx2.fillRect(dxx0, M.t, dxx1 - dxx0, plotH2);
        ctx2.restore();
      }
      // Cursor line
      if (S.cursorE !== null) {
        var cx2 = xP2(S.cursorE);
        if (cx2 >= M.l && cx2 <= W2 - M.r) {
          ctx2.beginPath(); ctx2.strokeStyle = dim; ctx2.lineWidth = 0.5;
          ctx2.moveTo(cx2, M.t); ctx2.lineTo(cx2, H2 - M.b); ctx2.stroke();
        }
      }

      // Update aria-pressed on range buttons
      var btns = seg.querySelectorAll("button");
      for (var bni = 0; bni < btns.length; bni++) {
        var br = btns[bni].getAttribute("data-range");
        var pressed = false;
        if (S.zoom === null && br === "30-3000") pressed = true;
        else if (S.zoom !== null) {
          var parts2 = br.split("-");
          if (Math.abs(S.zoom.xLo - parseFloat(parts2[0])) < 0.01 && Math.abs(S.zoom.xHi - parseFloat(parts2[1])) < 0.01) pressed = true;
        }
        btns[bni].setAttribute("aria-pressed", pressed ? "true" : "false");
      }
    }

    // --- Interaction ---
    function onMouseMove(ev, cv) {
      var rect = cv.getBoundingClientRect();
      var x = ev.clientX - rect.left;
      var w = cssW(cv);
      if (S.dragging && S.drag) {
        var plotMin = M.l, plotMax = w - M.r;
        S.drag.x1 = Math.max(plotMin, Math.min(plotMax, x));
        drawAll();
        return;
      }
      if (x < M.l || x > w - M.r) {
        S.cursorE = null;
        cursorLine.textContent = "наведи курсор на график";
        drawAll();
        return;
      }
      S.cursorE = eFromX(x, w);
      // Find nearest channel
      var iNearest = 0, bestD = Infinity;
      for (var ni = 0; ni < eArr.length; ni++) {
        var d = Math.abs(eArr[ni] - S.cursorE);
        if (d < bestD) { bestD = d; iNearest = ni; }
      }
      var txt = "канал " + iNearest + " \u00B7 " + fmt(eArr[iNearest], 1) + " кэВ \u00B7 измерение " + Math.round(counts[iNearest]) + " \u00B7 фон " + fmt(bgArr[iNearest], 1);
      // checked prefixes
      for (var cp = 0; cp < prefixes.length; cp++) {
        var cpfx = prefixes[cp];
        if (!S.on[cpfx]) continue;
        var cm = models[cpfx];
        if (cm && cm[iNearest] != null) {
          txt += " \u00B7 " + (T.labels[cpfx] || cpfx) + " " + fmt(cm[iNearest], 1);
        }
      }
      cursorLine.textContent = txt;
      drawAll();
    }

    function onMouseLeave(cv) {
      if (S.dragging) return;
      S.cursorE = null;
      cursorLine.textContent = "наведи курсор на график";
      drawAll();
    }

    function onMouseDown(ev, cv) {
      var rect = cv.getBoundingClientRect();
      var x = ev.clientX - rect.left;
      var w = cssW(cv);
      if (x >= M.l && x <= w - M.r) {
        ev.preventDefault();
        S.dragging = true;
        S.drag = {x0: x, x1: x, w: w};
      }
    }

    function onMouseUp(ev) {
      if (!S.dragging) return;
      S.dragging = false;
      if (S.drag) {
        var dx = Math.abs(S.drag.x1 - S.drag.x0);
        if (dx < 6) { S.drag = null; drawAll(); return; }
        var w = S.drag.w;
        var e0 = eFromX(Math.min(S.drag.x0, S.drag.x1), w);
        var e1 = eFromX(Math.max(S.drag.x0, S.drag.x1), w);
        S.zoom = {xLo: Math.max(0, e0), xHi: e1};
        S.drag = null;
        drawAll();
      }
    }

    // Attach listeners to both canvases
    [cv1, cv2].forEach(function (cv) {
      cv.addEventListener("mousemove", function (ev) { onMouseMove(ev, cv); });
      cv.addEventListener("mouseleave", function () { onMouseLeave(cv); });
      cv.addEventListener("mousedown", function (ev) { onMouseDown(ev, cv); });
      cv.addEventListener("dblclick", function () { S.zoom = null; drawAll(); });
    });
    document.addEventListener("mouseup", onMouseUp);

    // Range buttons
    var rBtns = seg.querySelectorAll("button");
    for (var rbi = 0; rbi < rBtns.length; rbi++) {
      (function (btn) {
        btn.addEventListener("click", function () {
          var r = btn.getAttribute("data-range");
          if (r === "30-3000") { S.zoom = null; }
          else {
            var pp = r.split("-");
            S.zoom = {xLo: parseFloat(pp[0]), xHi: parseFloat(pp[1])};
          }
          drawAll();
        });
      })(rBtns[rbi]);
    }

    // Checkboxes
    var cbs = chips.querySelectorAll("input[type=checkbox]");
    for (var cbi = 0; cbi < cbs.length; cbi++) {
      (function (cb) {
        cb.addEventListener("change", function () {
          var p = cb.getAttribute("data-p");
          S.on[p] = cb.checked;
          drawAll();
        });
      })(cbs[cbi]);
    }

    // Resize
    window.addEventListener("resize", drawAll);

    // Initial draw
    drawAll();

    // Expose
    window.G1SD_UI = {state: S, redraw: drawAll};
  }

  /* ---- 2. Sections ---- */
  var main = document.getElementById("main");
  var builders = {"tbl-fits": buildFits, "tbl-nodes": buildNodes, "tbl-peaks": buildPeaks, "charts": buildCharts};
  for (var si = 0; si < T.sections.length; si++) {
    var sec = T.sections[si];
    var secEl = document.createElement("section");
    secEl.className = "sec"; secEl.id = sec.id;
    var h2 = document.createElement("h2"); h2.className = "sec-h"; h2.textContent = sec.h;
    secEl.appendChild(h2);
    var body = document.createElement("div"); body.className = "sec-body";
    body.innerHTML = sec.html || "";
    secEl.appendChild(body);
    var widget = document.createElement("div"); widget.className = "widget";
    secEl.appendChild(widget);
    main.appendChild(secEl);
    if (sec.widget && builders[sec.widget]) {
      builders[sec.widget](widget);
    }
  }

  /* ---- 4. Value spans ---- */
  var dataVs = document.querySelectorAll("[data-v]");
  for (var vi2 = 0; vi2 < dataVs.length; vi2++) {
    var vel = dataVs[vi2];
    var ftype = vel.getAttribute("data-fmt");
    var vpath = vel.getAttribute("data-v");

    if (ftype === "relrange") {
      var paths = vpath.split("|");
      var rvals = [];
      for (var pi3 = 0; pi3 < paths.length; pi3++) {
        var arr = get(paths[pi3]);
        if (arr && arr.length) {
          for (var ki = 0; ki < arr.length; ki++) {
            if (arr[ki] && arr[ki].r != null) rvals.push((arr[ki].r - 1) * 100);
          }
        }
      }
      if (rvals.length) {
        vel.textContent = signed(Math.min.apply(null, rvals), 1) + "\u2026" + signed(Math.max.apply(null, rvals), 1);
      } else { vel.textContent = "\u2014"; }
      continue;
    }

    if (ftype === "drel") {
      var dv = get(vpath);
      var refP = vel.getAttribute("data-ref");
      var rv = refP ? get(refP) : undefined;
      if (dv != null && rv != null && rv !== 0) {
        vel.textContent = signed((dv / rv - 1) * 100, 1);
      } else {
        vel.textContent = "\u2014";
        console.warn("data-v path not resolved: " + vpath);
      }
      continue;
    }

    var val = get(vpath);
    if (val === undefined || val === null) {
      vel.textContent = "\u2014";
      console.warn("data-v path not resolved: " + vpath);
      continue;
    }
    switch (ftype) {
      case "r3": vel.textContent = fmt(val, 3); break;
      case "r2": vel.textContent = fmt(val, 2); break;
      case "int": vel.textContent = thin(val); break;
      case "str":
        var sv = String(val);
        if (/^\d{4}-\d{2}-\d{2}$/.test(sv)) {
          var dp = sv.split("-");
          sv = dp[2] + "." + dp[1] + "." + dp[0];
        }
        vel.textContent = sv; break;
      case "pct2": vel.textContent = fmt(val * 100, 2); break;
      default: vel.textContent = String(val);
    }
  }
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", init);
} else {
  init();
}

})();
