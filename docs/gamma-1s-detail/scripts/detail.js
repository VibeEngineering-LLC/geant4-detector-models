(function () {
"use strict";

var D = window.G1SD;
var T = window.G1SD_TEXT;
var L = T.labels;
var CHAINS = D.chains;
var VARIANTS = D.variants;
var CONTROL = D.control;
var REF = D.ref;

/* ---------- 3. Number formatting ---------- */

function fmt(x, d) {
  if (x === null || x === undefined || (typeof x === "number" && !isFinite(x))) return "\u2014";
  var s = x.toFixed(d);
  s = s.replace(".", ",");
  s = s.replace(/^-/, "\u2212");
  return s;
}

function signed(x, d) {
  if (x === null || x === undefined || (typeof x === "number" && !isFinite(x))) return "\u2014";
  var rounded = x.toFixed(d);
  var num = parseFloat(rounded);
  if (num === 0) rounded = (0).toFixed(d); // no "−0,0"
  var s = rounded;
  s = s.replace(".", ",");
  if (num > 0) {
    if (s.charAt(0) !== "+") s = "+" + s;
  } else {
    s = s.replace(/^-/, "\u2212");
  }
  return s;
}

function rel(r) {
  return signed((r - 1) * 100, 1);
}

function thin(n) {
  n = Math.round(n);
  var neg = n < 0;
  var s = String(Math.abs(n));
  var parts = [];
  while (s.length > 3) {
    parts.unshift(s.slice(-3));
    s = s.slice(0, -3);
  }
  parts.unshift(s);
  return (neg ? "\u2212" : "") + parts.join("\u202F");
}

/* ---------- Helpers ---------- */

function el(tag, cls, txt) {
  var e = document.createElement(tag);
  if (cls) e.className = cls;
  if (txt !== undefined && txt !== null) e.textContent = txt;
  return e;
}

function get(path) {
  var parts = path.split(".");
  var cur = window.G1SD;
  for (var i = 0; i < parts.length; i++) {
    if (cur === null || cur === undefined) return undefined;
    cur = cur[parts[i]];
  }
  return cur;
}

/* ---------- 1. Header and toolbar ---------- */

function buildHeader() {
  document.title = T.title;
  var hdr = document.getElementById("hdr");
  if (hdr) {
    hdr.innerHTML = '<p class="eyebrow">' + T.eyebrow + "</p><h1>" + T.h1 + "</h1>" + T.stand;
  }
  var tb = document.getElementById("toolbar");
  if (tb) {
    tb.innerHTML = "";
    var secs = T.sections;
    for (var i = 0; i < secs.length; i++) {
      var a = el("a", "btn", secs[i].h);
      a.href = "#" + secs[i].id;
      tb.appendChild(a);
    }
    var sp = el("span", "sp");
    tb.appendChild(sp);
    var gh = el("a", "btn", "\u043a\u043e\u0434 \u043d\u0430 GitHub");
    gh.href = "https://github.com/VibeEngineering-LLC/geant4-detector-models/tree/main/detectors/Gamma-1S";
    gh.target = "_blank";
    gh.rel = "noopener";
    tb.appendChild(gh);
  }
}

/* ---------- 2. Sections ---------- */

function buildSections() {
  var main = document.getElementById("main");
  if (!main) return;
  main.innerHTML = "";
  var secs = T.sections;
  for (var i = 0; i < secs.length; i++) {
    var s = secs[i];
    var sec = el("section", "sec");
    sec.id = s.id;
    var h = el("h2", "sec-h", s.h);
    sec.appendChild(h);
    var body = el("div", "sec-body");
    body.innerHTML = s.html || "";
    sec.appendChild(body);
    var wdiv = el("div", "widget");
    sec.appendChild(wdiv);
    main.appendChild(sec);
    if (s.widget) {
      if (s.widget === "tbl-fits") buildTblFits(wdiv);
      else if (s.widget === "tbl-nodes") buildTblNodes(wdiv);
      else if (s.widget === "tbl-peaks") buildTblPeaks(wdiv);
      else if (s.widget === "charts") buildCharts(wdiv);
    }
  }
}

/* ---------- 4. Value spans ---------- */

function resolveFmt(node) {
  var f = node.getAttribute("data-fmt");
  var vpath = node.getAttribute("data-v");
  if (!f || !vpath) return null;

  if (f === "relrange") {
    var paths = vpath.split("|");
    var vals = [];
    for (var i = 0; i < paths.length; i++) {
      var arr = get(paths[i]);
      if (arr && typeof arr === "object") {
        if (arr.length !== undefined) {
          for (var j = 0; j < arr.length; j++) {
            if (arr[j] && typeof arr[j].r === "number") vals.push(arr[j].r);
          }
        }
      }
    }
    if (vals.length === 0) return "\u2014";
    var lo = Infinity, hi = -Infinity;
    for (var k = 0; k < vals.length; k++) {
      var pc = (vals[k] - 1) * 100;
      if (pc < lo) lo = pc;
      if (pc > hi) hi = pc;
    }
    return signed(lo, 1) + "\u2026" + signed(hi, 1);
  }

  var v = get(vpath);
  if (v === undefined || v === null) {
    console.warn("data-v path not resolved: " + vpath);
    return "\u2014";
  }

  if (f === "r3") return fmt(v, 3);
  if (f === "r2") return fmt(v, 2);
  if (f === "int") return thin(Math.round(v));
  if (f === "str") return String(v).replace(/^(\d{4})-(\d{2})-(\d{2})$/, "$3.$2.$1");
  if (f === "pct2") return fmt(v * 100, 2);
  if (f === "drel") {
    var refpath = node.getAttribute("data-ref");
    var rv = refpath ? get(refpath) : undefined;
    if (rv === undefined || rv === null || rv === 0) return "\u2014";
    return signed((v / rv - 1) * 100, 1);
  }
  return fmt(v, 2);
}

function fillValueSpans() {
  var nodes = document.querySelectorAll("[data-v]");
  for (var i = 0; i < nodes.length; i++) {
    var txt = resolveFmt(nodes[i]);
    if (txt !== null) nodes[i].textContent = txt;
  }
}

/* ---------- 5. Tables ---------- */

function deltaPct(x, ref) {
  if (ref === undefined || ref === null || ref === 0) return "\u2014";
  return signed((x / ref - 1) * 100, 1);
}

function makeTable(headers) {
  var wrap = el("div", "scroll");
  var tbl = el("table", "big");
  var thead = el("thead");
  var hr = el("tr");
  for (var i = 0; i < headers.length; i++) {
    var th = el("th", null, headers[i]);
    hr.appendChild(th);
  }
  thead.appendChild(hr);
  tbl.appendChild(thead);
  var tbody = el("tbody");
  tbl.appendChild(tbody);
  wrap.appendChild(tbl);
  return {wrap: wrap, tbody: tbody};
}

function addRow(tbody, cells, numericFrom) {
  var tr = el("tr");
  for (var i = 0; i < cells.length; i++) {
    var td = el("td", null, cells[i]);
    if (i >= numericFrom) td.className = "num";
    tr.appendChild(td);
  }
  tbody.appendChild(tr);
}

function buildTblFits(w) {
  var rows = [CONTROL, REF, "dtP"];
  var fFits = D.fits;
  var t = makeTable(["\u0428\u0430\u0431\u043b\u043e\u043d\u044b", "A(Ra-228)/\u043f\u0430\u0441\u043f\u043e\u0440\u0442", "\u0394 \u043a \u043f\u043e\u043b\u043d\u043e\u0439, %", "A(Th-228)/\u043f\u0430\u0441\u043f\u043e\u0440\u0442", "\u0394 \u043a \u043f\u043e\u043b\u043d\u043e\u0439, %", "\u03c7\u00b2_ref/\u03bd"]);
  for (var i = 0; i < rows.length; i++) {
    var p = rows[i];
    var f = fFits[p];
    if (!f) continue;
    var isRef = (p === REF);
    var raStr = fmt(f.ra, 3) + " \u00b1 " + fmt(f.ra * f.dra, 3);
    var thStr = fmt(f.th, 3) + " \u00b1 " + fmt(f.th * f.dth, 3);
    var drA = isRef ? "\u2014" : deltaPct(f.ra, fFits[REF].ra);
    var drT = isRef ? "\u2014" : deltaPct(f.th, fFits[REF].th);
    addRow(t.tbody, [L[p] || p, raStr, drA, thStr, drT, fmt(f.chi2, 2)], 1);
  }
  w.appendChild(t.wrap);
}

function collectPeakR(prefix) {
  var vals = [];
  for (var c = 0; c < CHAINS.length; c++) {
    var pk = D.peaks[CHAINS[c]][prefix];
    if (pk && pk.peaks) {
      for (var i = 0; i < pk.peaks.length; i++) {
        if (typeof pk.peaks[i].r === "number") vals.push(pk.peaks[i].r);
      }
    }
  }
  return vals;
}

function collectBandR(prefix) {
  var vals = [];
  for (var c = 0; c < CHAINS.length; c++) {
    var pk = D.peaks[CHAINS[c]][prefix];
    if (pk && pk.bands) {
      for (var i = 0; i < pk.bands.length; i++) {
        var b = pk.bands[i];
        if (b.lo >= 50 && b.hi <= 260 && typeof b.r === "number") vals.push(b.r);
      }
    }
  }
  return vals;
}

function rangeStr(vals) {
  if (!vals || vals.length === 0) return "\u2014";
  var lo = Infinity, hi = -Infinity;
  for (var i = 0; i < vals.length; i++) {
    var pc = (vals[i] - 1) * 100;
    if (pc < lo) lo = pc;
    if (pc > hi) hi = pc;
  }
  return signed(lo, 1) + "\u2026" + signed(hi, 1);
}

function buildTblNodes(w) {
  var rows = [REF, "dtH", "dtE1", "dtE2", "dtE3", "dtE4", "dtP"];
  var fFits = D.fits;
  var t = makeTable(["\u0412\u0430\u0440\u0438\u0430\u043d\u0442 (\u043e\u0442 \u043f\u043e\u043b\u043d\u043e\u0439 \u043c\u043e\u0434\u0435\u043b\u0438)", "A(Ra-228)/\u043f\u0430\u0441\u043f\u043e\u0440\u0442", "\u0394, %", "A(Th-228)/\u043f\u0430\u0441\u043f\u043e\u0440\u0442", "\u0394, %", "\u03c7\u00b2_ref/\u03bd", "\u043f\u0438\u043a\u0438, %", "50\u2013260 \u043a\u044d\u0412, %"]);
  for (var i = 0; i < rows.length; i++) {
    var p = rows[i];
    var f = fFits[p];
    if (!f) continue;
    var isRef = (p === REF);
    var drA = isRef ? "\u2014" : deltaPct(f.ra, fFits[REF].ra);
    var drT = isRef ? "\u2014" : deltaPct(f.th, fFits[REF].th);
    var pkR = isRef ? "\u2014" : rangeStr(collectPeakR(p));
    var bdR = isRef ? "\u2014" : rangeStr(collectBandR(p));
    addRow(t.tbody, [L[p] || p, fmt(f.ra, 3), drA, fmt(f.th, 3), drT, fmt(f.chi2, 2), pkR, bdR], 1);
  }
  w.appendChild(t.wrap);
}

function buildTblPeaks(w) {
  var prefixes = VARIANTS.concat([CONTROL]);
  for (var c = 0; c < CHAINS.length; c++) {
    var chain = CHAINS[c];
    var cap = el("p", "tbl-cap", T.chain_labels[chain] || chain);
    w.appendChild(cap);

    var pkData = D.peaks[chain];
    /* collect all peaks and bands from the first available prefix to determine columns */
    var peaksList = null;
    var bandsList = null;
    for (var pi = 0; pi < prefixes.length; pi++) {
      var pd = pkData[prefixes[pi]];
      if (pd) {
        if (!peaksList && pd.peaks) peaksList = pd.peaks;
        if (!bandsList && pd.bands) bandsList = pd.bands;
        if (peaksList && bandsList) break;
      }
    }
    if (!peaksList) peaksList = [];
    if (!bandsList) bandsList = [];

    var headers = ["\u0412\u0430\u0440\u0438\u0430\u043d\u0442"];
    for (var pi2 = 0; pi2 < peaksList.length; pi2++) {
      headers.push(fmt(peaksList[pi2].E, 1) + " \u043a\u044d\u0412");
    }
    for (var bi = 0; bi < bandsList.length; bi++) {
      headers.push(bandsList[bi].lo + "\u2013" + bandsList[bi].hi + " \u043a\u044d\u0412");
    }

    var t = makeTable(headers);
    for (var ri = 0; ri < prefixes.length; ri++) {
      var p = prefixes[ri];
      var pd = pkData[p];
      var cells = [L[p] || p];
      for (var pi3 = 0; pi3 < peaksList.length; pi3++) {
        if (pd && pd.peaks && pd.peaks[pi3]) {
          cells.push(fmt(pd.peaks[pi3].r, 3) + " \u00b1 " + fmt(pd.peaks[pi3].s, 3));
        } else {
          cells.push("\u2014");
        }
      }
      for (var bi2 = 0; bi2 < bandsList.length; bi2++) {
        if (pd && pd.bands && pd.bands[bi2]) {
          cells.push(fmt(pd.bands[bi2].r, 3));
        } else {
          cells.push("\u2014");
        }
      }
      addRow(t.tbody, cells, 1);
    }
    w.appendChild(t.wrap);
  }
}

/* ---------- 6. Charts ---------- */

var COL = {
  dtF: "#16140f",
  dtP: "#d1495b",
  dtH: "#edae49",
  dtE1: "#66a182",
  dtE2: "#00798c",
  dtE3: "#8d5a97",
  dtE4: "#30638e",
  k24g1t: "#9a9a9a"
};

var S; /* chart state */
var cvSpec, cvRatio, cursorLineEl;
var M = {l: 62, r: 14, t: 12, b: 34};
var SUPER = ["\u2070","\u00b9","\u00b2","\u00b3","\u2074","\u2075","\u2076","\u2077","\u00b8","\u2079"];

function getInk() {
  var v = getComputedStyle(document.documentElement).getPropertyValue("--ink");
  return v ? v.trim() : "#16140f";
}

function getRuleSoft() {
  var v = getComputedStyle(document.documentElement).getPropertyValue("--rule-soft");
  return v ? v.trim() : "#cccccc";
}

function getAccentSoft() {
  var v = getComputedStyle(document.documentElement).getPropertyValue("--accent-soft");
  return v ? v.trim() : "#d1495b33";
}

function getBodyFont() {
  var v = getComputedStyle(document.body).getPropertyValue("font-family");
  return v ? v.trim() : "sans-serif";
}

function niceStep(range, target) {
  var raw = range / target;
  var mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10));
  var frac = raw / mag;
  var nice;
  if (frac <= 1.5) nice = 1;
  else if (frac <= 3.5) nice = 2;
  else if (frac <= 7.5) nice = 5;
  else nice = 10;
  return nice * mag;
}

function supPow(n) {
  /* n is integer exponent, return "10" + superscript */
  if (n < 0) return "10" + "\u207b" + supDigits(-n);
  return "10" + supDigits(n);
}

function supDigits(n) {
  var s = "";
  n = Math.abs(Math.round(n));
  if (n === 0) return SUPER[0];
  var digits = [];
  while (n > 0) { digits.unshift(n % 10); n = Math.floor(n / 10); }
  for (var i = 0; i < digits.length; i++) s += SUPER[digits[i]];
  return s;
}

function expFmt(v) {
  /* toExponential(3) with comma and minus sign */
  var s = v.toExponential(3);
  s = s.replace(".", ",");
  s = s.replace(/e([+-])/, function (m, sign) {
    if (sign === "-") return "e\u2212";
    return "e+";
  });
  return s;
}

function getRange() {
  if (S.zoom) return {xLo: S.zoom.xLo, xHi: S.zoom.xHi};
  return {xLo: 30, xHi: 3000};
}

function xToE(x, w, xLo, xHi) {
  return xLo + (x - M.l) / (w - M.l - M.r) * (xHi - xLo);
}

function eToX(e, w, xLo, xHi) {
  return M.l + (e - xLo) / (xHi - xLo) * (w - M.l - M.r);
}

function setupCanvas(cv, cssW, cssH) {
  var dpr = window.devicePixelRatio || 1;
  cv.width = cssW * dpr;
  cv.height = cssH * dpr;
  cv.style.width = cssW + "px";
  cv.style.height = cssH + "px";
  var ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return ctx;
}

function drawAll() {
  if (!cvSpec || !cvRatio) return;
  var xLo, xHi;
  if (S.zoom) { xLo = S.zoom.xLo; xHi = S.zoom.xHi; }
  else { xLo = 30; xHi = 3000; }

  var chain = S.chain;
  var spec = D.spectra[chain];
  if (!spec) return;
  var eArr = spec.e;
  var yObj = spec.y;
  var ratioObj = spec.ratio;

  var ink = getInk();
  var ruleSoft = getRuleSoft();
  var accentSoft = getAccentSoft();
  var font = "11px " + getBodyFont();

  /* --- cvSpec (log y) --- */
  var w1 = cvSpec.parentElement ? cvSpec.parentElement.clientWidth : 800;
  var h1 = 300;
  var ctx1 = setupCanvas(cvSpec, w1, h1);
  ctx1.clearRect(0, 0, w1, h1);

  /* collect visible bins */
  var visIdx = [];
  for (var i = 0; i < eArr.length; i++) {
    if (eArr[i] >= xLo && eArr[i] <= xHi) visIdx.push(i);
  }

  /* determine y range from positive values of drawn series */
  var yMin = Infinity, yMax = -Infinity;
  function updateRange(arr) {
    for (var k = 0; k < visIdx.length; k++) {
      var v = arr[visIdx[k]];
      if (v !== null && v !== undefined && v > 0) {
        if (v < yMin) yMin = v;
        if (v > yMax) yMax = v;
      }
    }
  }
  updateRange(yObj[REF]);
  var allPrefixes = [REF].concat(VARIANTS).concat([CONTROL]);
  for (var ap = 0; ap < allPrefixes.length; ap++) {
    if (S.on[allPrefixes[ap]]) updateRange(yObj[allPrefixes[ap]]);
  }
  if (yMin === Infinity || yMax === -Infinity) { yMin = 1e-7; yMax = 1e-3; }
  yMin *= 0.7;
  yMax *= 1.4;
  var logMin = Math.log10(yMin);
  var logMax = Math.log10(yMax);

  function specY(v) {
    return M.t + (logMax - Math.log10(v)) / (logMax - logMin) * (h1 - M.t - M.b);
  }
  function specX(e) {
    return eToX(e, w1, xLo, xHi);
  }

  /* grid + y ticks (powers of 10) */
  ctx1.font = font;
  ctx1.textAlign = "right";
  ctx1.textBaseline = "middle";
  var tickLo = Math.ceil(logMin);
  var tickHi = Math.floor(logMax);
  for (var tk = tickLo; tk <= tickHi; tk++) {
    var ty = specY(Math.pow(10, tk));
    if (ty < M.t || ty > h1 - M.b) continue;
    ctx1.strokeStyle = ruleSoft;
    ctx1.lineWidth = 0.5;
    ctx1.beginPath();
    ctx1.moveTo(M.l, ty);
    ctx1.lineTo(w1 - M.r, ty);
    ctx1.stroke();
    ctx1.fillStyle = ink;
    ctx1.fillText(supPow(tk), M.l - 6, ty);
  }

  /* y axis label */
  ctx1.save();
  ctx1.translate(14, (M.t + h1 - M.b) / 2);
  ctx1.rotate(-Math.PI / 2);
  ctx1.textAlign = "center";
  ctx1.fillStyle = ink;
  ctx1.fillText("\u043e\u0442\u0441\u0447\u0451\u0442\u044b / \u0440\u0430\u0441\u043f\u0430\u0434 / \u043a\u044d\u0412", 0, 0);
  ctx1.restore();

  /* draw series as step lines */
  function drawStep(ctx, arr, color, lw) {
    ctx.strokeStyle = color;
    ctx.lineWidth = lw;
    ctx.beginPath();
    var started = false;
    for (var j = 0; j < visIdx.length; j++) {
      var idx = visIdx[j];
      var val = arr[idx];
      if (val === null || val === undefined || val <= 0) {
        started = false;
        continue;
      }
      var px = specX(eArr[idx]);
      var py = specY(val);
      if (!started) {
        ctx.moveTo(px, py);
        started = true;
      } else {
        /* step: horizontal then vertical */
        var prevIdx = visIdx[j - 1];
        var prevPx = specX(eArr[prevIdx]);
        ctx.lineTo(px, py); /* this is the vertical at current x going to current y */
        /* actually step: from prev point horizontal to current x, then vertical */
        /* Let me redo: I'll draw lineTo(prevPx, py) first then lineTo(px, py) */
      }
    }
    ctx.stroke();
  }

  /* Better step drawing */
  function drawStepLine(ctx, arr, color, lw) {
    ctx.strokeStyle = color;
    ctx.lineWidth = lw;
    ctx.beginPath();
    var lastX = null, lastY = null;
    for (var j = 0; j < visIdx.length; j++) {
      var idx = visIdx[j];
      var val = arr[idx];
      if (val === null || val === undefined || val <= 0) {
        lastX = null;
        continue;
      }
      var px = specX(eArr[idx]);
      var py = specY(val);
      if (lastX === null) {
        ctx.moveTo(px, py);
      } else {
        /* step: horizontal from lastX to px at current py, then we're at (px, py) */
        ctx.lineTo(px, lastY);
        ctx.lineTo(px, py);
      }
      lastX = px;
      lastY = py;
    }
    ctx.stroke();
  }

  /* dtF always */
  if (yObj[REF]) drawStepLine(ctx1, yObj[REF], COL.dtF, 2);
  /* checked variants + control */
  for (var ci = 0; ci < allPrefixes.length; ci++) {
    var cp = allPrefixes[ci];
    if (cp === REF) continue;
    if (S.on[cp] && yObj[cp]) {
      drawStepLine(ctx1, yObj[cp], COL[cp] || "#666", 1.25);
    }
  }

  /* --- cvRatio (linear y) --- */
  var w2 = cvRatio.parentElement ? cvRatio.parentElement.clientWidth : 800;
  var h2 = 260;
  var ctx2 = setupCanvas(cvRatio, w2, h2);
  ctx2.clearRect(0, 0, w2, h2);

  /* collect r±s for y range */
  var rLo = Infinity, rHi = -Infinity;
  for (var cp2 = 0; cp2 < allPrefixes.length; cp2++) {
    var pfx = allPrefixes[cp2];
    if (pfx === REF) continue;
    if (!S.on[pfx]) continue;
    var rArr = ratioObj && ratioObj[pfx] ? ratioObj[pfx].r : null;
    var sArr = ratioObj && ratioObj[pfx] ? ratioObj[pfx].s : null;
    if (!rArr) continue;
    for (var vi = 0; vi < visIdx.length; vi++) {
      var rv = rArr[visIdx[vi]];
      var sv = sArr ? sArr[visIdx[vi]] : 0;
      if (rv === null || rv === undefined) continue;
      var lo2 = rv - (sv || 0);
      var hi2 = rv + (sv || 0);
      if (isFinite(lo2) && lo2 < rLo) rLo = lo2;
      if (isFinite(hi2) && hi2 > rHi) rHi = hi2;
    }
  }
  if (rLo === Infinity) { rLo = 0.97; rHi = 1.03; }
  /* clamp to [0.8, 1.2] */
  if (rLo < 0.8) rLo = 0.8;
  if (rHi > 1.2) rHi = 1.2;
  /* ensure at least [0.97, 1.03] */
  if (rLo > 0.97) rLo = 0.97;
  if (rHi < 1.03) rHi = 1.03;
  /* pad 10% */
  var pad = (rHi - rLo) * 0.1;
  rLo -= pad;
  rHi += pad;

  function ratioY(v) {
    return M.t + (rHi - v) / (rHi - rLo) * (h2 - M.t - M.b);
  }
  function ratioX(e) {
    return eToX(e, w2, xLo, xHi);
  }

  /* grid + y ticks */
  var rStep = niceStep(rHi - rLo, 5);
  var rTickStart = Math.ceil(rLo / rStep) * rStep;
  ctx2.font = font;
  ctx2.textAlign = "right";
  ctx2.textBaseline = "middle";
  for (var rt = rTickStart; rt <= rHi + rStep * 0.01; rt += rStep) {
    if (rt < rLo - 0.001) continue;
    var ry = ratioY(rt);
    if (ry < M.t || ry > h2 - M.b) continue;
    ctx2.strokeStyle = ruleSoft;
    ctx2.lineWidth = 0.5;
    ctx2.beginPath();
    ctx2.moveTo(M.l, ry);
    ctx2.lineTo(w2 - M.r, ry);
    ctx2.stroke();
    ctx2.fillStyle = ink;
    ctx2.fillText(fmt(rt, 2), M.l - 6, ry);
  }

  /* horizontal line at 1.0 */
  var y1 = ratioY(1.0);
  ctx2.strokeStyle = ink;
  ctx2.lineWidth = 1;
  ctx2.setLineDash([4, 3]);
  ctx2.beginPath();
  ctx2.moveTo(M.l, y1);
  ctx2.lineTo(w2 - M.r, y1);
  ctx2.stroke();
  ctx2.setLineDash([]);

  /* y axis label */
  ctx2.save();
  ctx2.translate(14, (M.t + h2 - M.b) / 2);
  ctx2.rotate(-Math.PI / 2);
  ctx2.textAlign = "center";
  ctx2.fillStyle = ink;
  ctx2.fillText("\u043a \u043f\u043e\u043b\u043d\u043e\u0439 \u043c\u043e\u0434\u0435\u043b\u0438", 0, 0);
  ctx2.restore();

  /* draw ratio series with bands */
  for (var cp3 = 0; cp3 < allPrefixes.length; cp3++) {
    var pfx3 = allPrefixes[cp3];
    if (pfx3 === REF) continue;
    if (!S.on[pfx3]) continue;
    var rArr2 = ratioObj && ratioObj[pfx3] ? ratioObj[pfx3].r : null;
    var sArr2 = ratioObj && ratioObj[pfx3] ? ratioObj[pfx3].s : null;
    if (!rArr2) continue;
    var col = COL[pfx3] || "#666";

    /* band */
    ctx2.globalAlpha = 0.18;
    ctx2.fillStyle = col;
    ctx2.beginPath();
    var bandStarted = false;
    for (var bj = 0; bj < visIdx.length; bj++) {
      var bidx = visIdx[bj];
      var brv = rArr2[bidx];
      var bsv = sArr2 ? sArr2[bidx] : 0;
      if (brv === null || brv === undefined || !isFinite(brv)) {
        bandStarted = false;
        continue;
      }
      var bx = ratioX(eArr[bidx]);
      var btop = ratioY(brv + (bsv || 0));
      if (!bandStarted) {
        ctx2.moveTo(bx, btop);
        bandStarted = true;
      } else {
        ctx2.lineTo(bx, btop);
      }
    }
    /* bottom in reverse */
    for (var bj2 = visIdx.length - 1; bj2 >= 0; bj2--) {
      var bidx2 = visIdx[bj2];
      var brv2 = rArr2[bidx2];
      var bsv2 = sArr2 ? sArr2[bidx2] : 0;
      if (brv2 === null || brv2 === undefined || !isFinite(brv2)) continue;
      var bx2 = ratioX(eArr[bidx2]);
      var bbot = ratioY(brv2 - (bsv2 || 0));
      ctx2.lineTo(bx2, bbot);
    }
    ctx2.closePath();
    ctx2.fill();
    ctx2.globalAlpha = 1;

    /* line */
    ctx2.strokeStyle = col;
    ctx2.lineWidth = 1.25;
    ctx2.beginPath();
    var lStarted = false;
    for (var lj = 0; lj < visIdx.length; lj++) {
      var lidx = visIdx[lj];
      var lrv = rArr2[lidx];
      if (lrv === null || lrv === undefined || !isFinite(lrv)) {
        lStarted = false;
        continue;
      }
      var lpx = ratioX(eArr[lidx]);
      var lpy = ratioY(lrv);
      if (!lStarted) { ctx2.moveTo(lpx, lpy); lStarted = true; }
      else ctx2.lineTo(lpx, lpy);
    }
    ctx2.stroke();
  }

  /* --- X axis ticks (shared, drawn on both) --- */
  function drawXTicks(ctx, w, h, label) {
    var step = niceStep(xHi - xLo, 7);
    var tStart = Math.ceil(xLo / step) * step;
    ctx.font = font;
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (var xt = tStart; xt <= xHi + step * 0.01; xt += step) {
      if (xt < xLo - step * 0.01) continue;
      var xx = eToX(xt, w, xLo, xHi);
      if (xx < M.l - 1 || xx > w - M.r + 1) continue;
      ctx.strokeStyle = ruleSoft;
      ctx.lineWidth = 0.5;
      ctx.beginPath();
      ctx.moveTo(xx, h - M.b);
      ctx.lineTo(xx, h - M.b + 4);
      ctx.stroke();
      ctx.fillStyle = ink;
      ctx.fillText(String(Math.round(xt)), xx, h - M.b + 6);
    }
    if (label) {
      ctx.textAlign = "center";
      ctx.fillText(label, (M.l + w - M.r) / 2, h - 8);
    }
  }

  drawXTicks(ctx1, w1, h1, "");
  drawXTicks(ctx2, w2, h2, "\u044d\u043d\u0435\u0440\u0433\u0438\u044f, \u043a\u044d\u0412");

  /* --- Drag rectangle --- */
  if (S.drag) {
    var dx0 = Math.min(S.drag.x0, S.drag.x1);
    var dx1 = Math.max(S.drag.x0, S.drag.x1);
    ctx1.globalAlpha = 0.35;
    ctx1.fillStyle = accentSoft;
    ctx1.fillRect(dx0, M.t, dx1 - dx0, h1 - M.t - M.b);
    ctx1.globalAlpha = 1;
    ctx2.globalAlpha = 0.35;
    ctx2.fillStyle = accentSoft;
    ctx2.fillRect(dx0, M.t, dx1 - dx0, h2 - M.t - M.b);
    ctx2.globalAlpha = 1;
  }

  /* --- Cursor line --- */
  if (S.cursorE !== null && S.cursorE !== undefined) {
    var cx1 = eToX(S.cursorE, w1, xLo, xHi);
    if (cx1 >= M.l && cx1 <= w1 - M.r) {
      ctx1.strokeStyle = ink;
      ctx1.lineWidth = 0.5;
      ctx1.beginPath();
      ctx1.moveTo(cx1, M.t);
      ctx1.lineTo(cx1, h1 - M.b);
      ctx1.stroke();
    }
    var cx2 = eToX(S.cursorE, w2, xLo, xHi);
    if (cx2 >= M.l && cx2 <= w2 - M.r) {
      ctx2.strokeStyle = ink;
      ctx2.lineWidth = 0.5;
      ctx2.beginPath();
      ctx2.moveTo(cx2, M.t);
      ctx2.lineTo(cx2, h2 - M.b);
      ctx2.stroke();
    }
  }

  /* update range button aria-pressed */
  updateRangeButtons();
}

function updateRangeButtons() {
  var seg = document.getElementById("rangeSeg");
  if (!seg) return;
  var btns = seg.querySelectorAll("[data-range]");
  var xLo, xHi;
  if (S.zoom) { xLo = S.zoom.xLo; xHi = S.zoom.xHi; }
  else { xLo = 30; xHi = 3000; }
  for (var i = 0; i < btns.length; i++) {
    var rng = btns[i].getAttribute("data-range");
    if (rng === "30-3000") {
      btns[i].setAttribute("aria-pressed", String(!S.zoom));
    } else {
      var parts = rng.split("-");
      var rlo = parseFloat(parts[0]);
      var rhi = parseFloat(parts[1]);
      var match = S.zoom && Math.abs(S.zoom.xLo - rlo) < 0.001 && Math.abs(S.zoom.xHi - rhi) < 0.001;
      btns[i].setAttribute("aria-pressed", String(match));
    }
  }
}

function buildCharts(w) {
  S = {
    chain: "Th228chain",
    zoom: null,
    drag: null,
    dragging: false,
    cursorE: null,
    on: {dtP: true, dtE2: true, dtE4: true}
  };

  /* ctrl div */
  var ctrl = el("div", "ctrl");

  /* chain seg */
  var chainSeg = el("div", "seg");
  chainSeg.id = "chainSeg";
  var chainBtns = [
    {chain: "Th228chain", text: "Th-228 \u0438 \u0434\u043e\u0447\u0435\u0440\u043d\u0438\u0435"},
    {chain: "Ra228chain", text: "Ra-228 + Ac-228"}
  ];
  for (var cb = 0; cb < chainBtns.length; cb++) {
    var btn = el("button", "btn", chainBtns[cb].text);
    btn.setAttribute("data-chain", chainBtns[cb].chain);
    btn.setAttribute("aria-pressed", String(S.chain === chainBtns[cb].chain));
    btn.addEventListener("click", function () {
      S.chain = this.getAttribute("data-chain");
      var allBtns = document.getElementById("chainSeg").querySelectorAll("[data-chain]");
      for (var j = 0; j < allBtns.length; j++) {
        allBtns[j].setAttribute("aria-pressed", String(allBtns[j].getAttribute("data-chain") === S.chain));
      }
      drawAll();
    });
    chainSeg.appendChild(btn);
  }
  ctrl.appendChild(chainSeg);

  /* range seg */
  var rangeSeg = el("div", "seg");
  rangeSeg.id = "rangeSeg";
  var rangeDefs = [
    {range: "30-3000", text: "\u0432\u0435\u0441\u044c"},
    {range: "30-300", text: "30\u2013300 \u043a\u044d\u0412"},
    {range: "300-1000", text: "300\u20131000 \u043a\u044d\u0412"},
    {range: "1000-2800", text: "1\u20132,8 \u041c\u044d\u0412"}
  ];
  for (var rb = 0; rb < rangeDefs.length; rb++) {
    var rbtn = el("button", "btn", rangeDefs[rb].text);
    rbtn.setAttribute("data-range", rangeDefs[rb].range);
    rbtn.setAttribute("aria-pressed", String(rb === 0));
    rbtn.addEventListener("click", function () {
      var rng = this.getAttribute("data-range");
      if (rng === "30-3000") {
        S.zoom = null;
      } else {
        var parts = rng.split("-");
        S.zoom = {xLo: parseFloat(parts[0]), xHi: parseFloat(parts[1])};
      }
      drawAll();
    });
    rangeSeg.appendChild(rbtn);
  }
  ctrl.appendChild(rangeSeg);

  /* chips */
  var chips = el("div", "chips");
  chips.id = "varChips";
  var chipPrefixes = VARIANTS.concat([CONTROL]);
  for (var chp = 0; chp < chipPrefixes.length; chp++) {
    var pfx = chipPrefixes[chp];
    var label = el("label", "chip");
    var inp = document.createElement("input");
    inp.type = "checkbox";
    inp.setAttribute("data-p", pfx);
    inp.checked = !!S.on[pfx];
    inp.addEventListener("change", function () {
      S.on[this.getAttribute("data-p")] = this.checked;
      drawAll();
    });
    label.appendChild(inp);
    var sw = el("span", "sw");
    sw.style.background = COL[pfx] || "#666";
    label.appendChild(sw);
    label.appendChild(document.createTextNode(L[pfx] || pfx));
    chips.appendChild(label);
  }
  ctrl.appendChild(chips);
  w.appendChild(ctrl);

  /* canvas boxes */
  var box1 = el("div", "chartbox");
  cvSpec = document.createElement("canvas");
  cvSpec.id = "cvSpec";
  cvSpec.setAttribute("height", "300");
  box1.appendChild(cvSpec);
  w.appendChild(box1);

  var box2 = el("div", "chartbox");
  cvRatio = document.createElement("canvas");
  cvRatio.id = "cvRatio";
  cvRatio.setAttribute("height", "260");
  box2.appendChild(cvRatio);
  w.appendChild(box2);

  /* cursor line */
  cursorLineEl = el("div", "cursor-line", "\u043d\u0430\u0432\u0435\u0434\u0438 \u043a\u0443\u0440\u0441\u043e\u0440 \u043d\u0430 \u0433\u0440\u0430\u0444\u0438\u043a");
  cursorLineEl.id = "cursorLine";
  w.appendChild(cursorLineEl);

  /* interactions */
  attachCanvasEvents(cvSpec);
  attachCanvasEvents(cvRatio);

  /* resize */
  window.addEventListener("resize", drawAll);

  /* initial draw */
  drawAll();
}

function attachCanvasEvents(cv) {
  cv.addEventListener("mousemove", function (ev) {
    var rect = cv.getBoundingClientRect();
    var x = ev.clientX - rect.left;
    var w = cv.clientWidth;
    if (S.dragging && S.drag) {
      var lo = M.l, hi = w - M.r;
      if (x < lo) x = lo;
      if (x > hi) x = hi;
      S.drag.x1 = x;
      drawAll();
      return;
    }
    /* cursor */
    if (x < M.l || x > w - M.r) {
      S.cursorE = null;
      if (cursorLineEl) cursorLineEl.textContent = "\u043d\u0430\u0432\u0435\u0434\u0438 \u043a\u0443\u0440\u0441\u043e\u0440 \u043d\u0430 \u0433\u0440\u0430\u0444\u0438\u043a";
      drawAll();
      return;
    }
    var xLo, xHi;
    if (S.zoom) { xLo = S.zoom.xLo; xHi = S.zoom.xHi; } else { xLo = 30; xHi = 3000; }
    var eVal = xToE(x, w, xLo, xHi);
    S.cursorE = eVal;

    /* find nearest bin */
    var spec = D.spectra[S.chain];
    if (spec && spec.e) {
      var eArr = spec.e;
      var bestI = 0, bestD = Infinity;
      for (var i = 0; i < eArr.length; i++) {
        var d = Math.abs(eArr[i] - eVal);
        if (d < bestD) { bestD = d; bestI = i; }
      }
      var parts = ["\u0431\u0438\u043d " + fmt(eArr[bestI], 1) + " \u043a\u044d\u0412"];
      var yObj = spec.y;
      var ratioObj = spec.ratio;
      /* dtF */
      if (yObj[REF] && yObj[REF][bestI] !== null && yObj[REF][bestI] !== undefined) {
        parts.push((L[REF] || REF) + " " + expFmt(yObj[REF][bestI]));
      }
      /* checked prefixes */
      var allP = VARIANTS.concat([CONTROL]);
      for (var pi = 0; pi < allP.length; pi++) {
        if (!S.on[allP[pi]]) continue;
        if (yObj[allP[pi]] && yObj[allP[pi]][bestI] !== null && yObj[allP[pi]][bestI] !== undefined) {
          var yv = yObj[allP[pi]][bestI];
          var seg2 = (L[allP[pi]] || allP[pi]) + " " + expFmt(yv);
          var rv = ratioObj && ratioObj[allP[pi]] && ratioObj[allP[pi]].r ? ratioObj[allP[pi]].r[bestI] : null;
          if (rv !== null && rv !== undefined) seg2 += " (" + fmt(rv, 3) + ")";
          parts.push(seg2);
        }
      }
      if (cursorLineEl) cursorLineEl.textContent = parts.join(" \u00b7 ");
    }
    drawAll();
  });

  cv.addEventListener("mouseleave", function () {
    if (S.dragging) return;
    S.cursorE = null;
    if (cursorLineEl) cursorLineEl.textContent = "\u043d\u0430\u0432\u0435\u0434\u0438 \u043a\u0443\u0440\u0441\u043e\u0440 \u043d\u0430 \u0433\u0440\u0430\u0444\u0438\u043a";
    drawAll();
  });

  cv.addEventListener("mousedown", function (ev) {
    var rect = cv.getBoundingClientRect();
    var x = ev.clientX - rect.left;
    var w = cv.clientWidth;
    if (x < M.l || x > w - M.r) return;
    ev.preventDefault();
    S.dragging = true;
    S.drag = {x0: x, x1: x, w: w};
  });

  document.addEventListener("mouseup", function (ev) {
    if (!S.dragging) return;
    S.dragging = false;
    if (S.drag) {
      var dx0 = S.drag.x0;
      var dx1 = S.drag.x1;
      if (Math.abs(dx1 - dx0) < 6) {
        S.drag = null;
        drawAll();
        return;
      }
      var w = S.drag.w;
      var xLo, xHi;
      if (S.zoom) { xLo = S.zoom.xLo; xHi = S.zoom.xHi; } else { xLo = 30; xHi = 3000; }
      var eA = xToE(Math.min(dx0, dx1), w, xLo, xHi);
      var eB = xToE(Math.max(dx0, dx1), w, xLo, xHi);
      S.zoom = {xLo: Math.max(0, eA), xHi: eB};
      S.drag = null;
      drawAll();
    }
  });

  cv.addEventListener("dblclick", function () {
    S.zoom = null;
    drawAll();
  });
}

/* ---------- Init ---------- */

function init() {
  buildHeader();
  buildSections();
  fillValueSpans();
  window.G1SD_UI = {state: S, redraw: drawAll};
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", init);
} else {
  init();
}

})();
