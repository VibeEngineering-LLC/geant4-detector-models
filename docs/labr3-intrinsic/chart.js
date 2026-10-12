/* chart.js — log-y spectrum chart (ES2019, no libraries, no modules) */
(function () {
  "use strict";

  // ---- globals / DOM ----
  var BG = window.BG;
  var canvas = document.getElementById("cv");
  var tip = document.getElementById("tip");
  var status = document.getElementById("status");
  var rangeBox = document.getElementById("ctl-range");
  var ctx = canvas.getContext("2d");

  // ---- state (only xlo/xhi define the view; sel is transient) ----
  var xlo = BG.ranges.all[0];
  var xhi = BG.ranges.all[1];
  var sel = null; // {x0,x1} in CSS px while dragging

  // ---- geometry (CSS px) ----
  var W = 0, H = 0, dpr = 1;
  var M = { l: 56, r: 12, t: 10, b: 36 };
  var plotW = 0, plotH = 0;

  // ---- cached colours ----
  var C = { grid: "#ccc", ink: "#222", photo: "#2d5fa8", sum: "#d21f1f" };

  // ---- data ----
  var E = BG.E, y = BG.y, model = BG.model, comps = BG.comps || [];
  var n = E.length;

  // ---- i18n ----
  var I18N = BG.i18n || {};
  var LOCALE = (BG.lang === "ru") ? "ru-RU" : "en-US";
  var fmtInt = new Intl.NumberFormat(LOCALE);
  var fmtNum = new Intl.NumberFormat(LOCALE, { maximumFractionDigits: 4 });

  // ---- helpers ----
  function clamp(v, a, b) { return v < a ? a : (v > b ? b : v); }

  function cssVar(name, fallback) {
    var v = getComputedStyle(document.documentElement).getPropertyValue(name);
    v = v ? v.trim() : "";
    return v || fallback;
  }

  function readColors() {
    C.grid = cssVar("--grid", "#cccccc");
    C.ink = cssVar("--ink", "#222222");
    C.photo = cssVar("--ch-photo", "#2d5fa8");
    C.sum = cssVar("--sum-line", "#d21f1f");
  }

  // binary search: index of the energy nearest to target
  function nearestIndex(target) {
    var lo = 0, hi = n - 1;
    if (target <= E[0]) return 0;
    if (target >= E[hi]) return hi;
    while (lo < hi) {
      var mid = (lo + hi) >> 1;
      if (E[mid] < target) lo = mid + 1; else hi = mid;
    }
    if (lo > 0 && (target - E[lo - 1]) <= (E[lo] - target)) lo--;
    return lo;
  }

  // ---- coordinate transforms (CSS px) ----
  function xToPx(e) { return M.l + (e - xlo) / (xhi - xlo) * plotW; }
  function pxToX(px) { return xlo + (px - M.l) / plotW * (xhi - xlo); }

  // ---- view range (decades) ----
  function viewRange() {
    var i0 = nearestIndex(xlo), i1 = nearestIndex(xhi);
    if (i1 < i0) { var t = i0; i0 = i1; i1 = t; }
    var mn = Infinity, mx = -Infinity;
    for (var i = i0; i <= i1; i++) {
      var v = y[i];
      if (v > 0) { if (v < mn) mn = v; if (v > mx) mx = v; }
    }
    if (!isFinite(mn)) mn = 1e-5;
    if (!isFinite(mx)) mx = 1;
    var dmin = Math.max(Math.floor(Math.log10(mn)), -5);
    var dmax = Math.ceil(Math.log10(mx));
    if (dmax <= dmin) dmax = dmin + 1;
    return { dmin: dmin, dmax: dmax };
  }

  function yToPx(v, dmin, dmax) {
    if (v <= 0) v = 1e-300;
    var t = (Math.log10(v) - dmin) / (dmax - dmin);
    return M.t + (1 - t) * plotH;
  }

  // ---- nice ticks ----
  function niceStep(range, target) {
    var raw = range / target;
    var mag = Math.pow(10, Math.floor(Math.log10(raw)));
    var norm = raw / mag;
    var step;
    if (norm < 1.5) step = 1; else if (norm < 3) step = 2; else if (norm < 7) step = 5; else step = 10;
    return step * mag;
  }

  function xTicks() {
    var range = xhi - xlo;
    var step = niceStep(range, 6);
    var start = Math.ceil(xlo / step) * step;
    var out = [];
    for (var v = start; v <= xhi + step * 1e-6; v += step) {
      if (v >= xlo - step * 1e-6) out.push(v);
    }
    return out;
  }

  function fmtTick(v) {
    if (Math.abs(v) >= 1000) return fmtInt.format(Math.round(v));
    return fmtNum.format(v);
  }

  // ---- drawing ----
  function draw() {
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    var vr = viewRange();
    var dmin = vr.dmin, dmax = vr.dmax;

    // grid + Y decade labels (superscript digits)
    ctx.font = "12px system-ui, sans-serif";
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (var d = dmin; d <= dmax; d++) {
      var py = yToPx(Math.pow(10, d), dmin, dmax);
      ctx.strokeStyle = C.grid;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(M.l, py);
      ctx.lineTo(M.l + plotW, py);
      ctx.stroke();
      ctx.fillStyle = C.ink;
      ctx.fillText("10" + sup(d), M.l - 6, py);
    }

    // X ticks + labels
    var xt = xTicks();
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (var k = 0; k < xt.length; k++) {
      var px = xToPx(xt[k]);
      if (px < M.l - 0.5 || px > M.l + plotW + 0.5) continue;
      ctx.strokeStyle = C.grid;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(px, M.t);
      ctx.lineTo(px, M.t + plotH);
      ctx.stroke();
      ctx.fillStyle = C.ink;
      ctx.fillText(fmtTick(xt[k]), px, M.t + plotH + 6);
    }

    // axes
    ctx.strokeStyle = C.ink;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(M.l, M.t);
    ctx.lineTo(M.l, M.t + plotH);
    ctx.lineTo(M.l + plotW, M.t + plotH);
    ctx.stroke();

    // clip to plot area for the curves
    ctx.save();
    ctx.beginPath();
    ctx.rect(M.l, M.t, plotW, plotH);
    ctx.clip();

    // measured curve (thin line)
    ctx.strokeStyle = C.photo;
    ctx.lineWidth = 1;
    ctx.beginPath();
    var started = false;
    for (var i = 0; i < n; i++) {
      var ex = E[i];
      if (ex < xlo || ex > xhi) continue;
      var cx = xToPx(ex), cy = yToPx(y[i], dmin, dmax);
      if (!started) { ctx.moveTo(cx, cy); started = true; } else ctx.lineTo(cx, cy);
    }
    ctx.stroke();

    // model curve (red), only where model > 0
    ctx.strokeStyle = C.sum;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    var mstarted = false;
    for (var j = 0; j < n; j++) {
      var ej = E[j];
      if (ej < xlo || ej > xhi) { mstarted = false; continue; }
      if (model[j] > 0) {
        var mx2 = xToPx(ej), my2 = yToPx(model[j], dmin, dmax);
        if (!mstarted) { ctx.moveTo(mx2, my2); mstarted = true; } else ctx.lineTo(mx2, my2);
      } else {
        mstarted = false;
      }
    }
    ctx.stroke();

    // template/external component curves (thin coloured lines), BG.curves = [{y:[...], color:"#rrggbb"}]
    var cv_ = BG.curves || [];
    for (var q = 0; q < cv_.length; q++) {
      ctx.strokeStyle = cv_[q].color; ctx.lineWidth = 1; ctx.beginPath();
      var qs = false;
      for (var w = 0; w < n; w++) {
        if (E[w] < xlo || E[w] > xhi || !(cv_[q].y[w] > 0)) { qs = false; continue; }
        var qx = xToPx(E[w]), qy = yToPx(cv_[q].y[w], dmin, dmax);
        if (!qs) { ctx.moveTo(qx, qy); qs = true; } else ctx.lineTo(qx, qy);
      }
      ctx.stroke();
    }

    // component marker lines (dashed) + rotated labels
    for (var c = 0; c < comps.length; c++) {
      var ce = comps[c].E;
      if (ce < xlo || ce > xhi) continue;
      var cxp = xToPx(ce);
      ctx.save();
      ctx.strokeStyle = C.ink;
      ctx.setLineDash([4, 4]);
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(cxp, M.t);
      ctx.lineTo(cxp, M.t + plotH);
      ctx.stroke();
      ctx.setLineDash([]);
      // rotated label near the top
      var label = String(comps[c].name || "").slice(0, 14);
      ctx.translate(cxp + 4, M.t + 12);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillStyle = C.ink;
      ctx.font = "11px system-ui, sans-serif";
      ctx.fillText(label, 0, 0);
      ctx.restore();
    }

    // selection band (translucent vertical band, full plot height)
    if (sel) {
      var a = Math.min(sel.x0, sel.x1), b = Math.max(sel.x0, sel.x1);
      ctx.fillStyle = "rgba(45,95,168,0.18)";
      ctx.fillRect(a, M.t, b - a, plotH);
      ctx.strokeStyle = C.photo;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(a, M.t); ctx.lineTo(a, M.t + plotH);
      ctx.moveTo(b, M.t); ctx.lineTo(b, M.t + plotH);
      ctx.stroke();
    }

    // cursor crosshair
    if (cursorIdx != null) {
      var cex = E[cursorIdx];
      if (cex >= xlo && cex <= xhi) {
        var cpx = xToPx(cex);
        ctx.save();
        ctx.strokeStyle = C.ink;
        ctx.setLineDash([2, 3]);
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(cpx, M.t);
        ctx.lineTo(cpx, M.t + plotH);
        ctx.stroke();
        ctx.restore();
      }
    }

    ctx.restore(); // unclip
  }

  // superscript digits
  function sup(d) {
    var map = { "0": "\u2070", "1": "\u00B9", "2": "\u00B2", "3": "\u00B3", "4": "\u2074",
      "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079", "-": "\u207B" };
    var s = String(d), out = "";
    for (var i = 0; i < s.length; i++) out += map[s[i]] || s[i];
    return out;
  }

  // ---- cursor / tooltip ----
  var cursorIdx = null;

  function showTip(idx, clientX, clientY) {
    var ch = idx;
    var en = E[idx];
    var yv = y[idx];
    var mv = model[idx];
    var line = I18N.cursor + " " + fmtInt.format(ch) +
      " \u00B7 " + I18N.energy + " " + fmtNum.format(en) + (window.BG.lang === "ru" ? " кэВ" : " keV") +
      " \u00B7 " + I18N.rate + " " + fmtNum.format(yv) +
      " \u00B7 " + I18N.model + " " + (mv > 0 ? fmtNum.format(mv) : "\u2014");
    tip.textContent = line;
    tip.hidden = false;
    var box = canvas.parentElement;
    var br = box.getBoundingClientRect();
    var tx = clientX - br.left + 14;
    var ty = clientY - br.top + 14;
    var tw = tip.offsetWidth, th = tip.offsetHeight;
    if (tx + tw > br.width) tx = clientX - br.left - tw - 14;
    if (ty + th > br.height) ty = clientY - br.top - th - 14;
    if (tx < 0) tx = 0;
    if (ty < 0) ty = 0;
    tip.style.left = tx + "px";
    tip.style.top = ty + "px";
  }

  function hideTip() {
    tip.hidden = true;
    cursorIdx = null;
    draw();
  }

  // ---- range buttons ----
  function markButtons(key) {
    var btns = rangeBox.querySelectorAll(".btn[data-range]");
    for (var i = 0; i < btns.length; i++) {
      btns[i].setAttribute("aria-pressed", btns[i].getAttribute("data-range") === key ? "true" : "false");
    }
  }

  function setRange(key) {
    var r = BG.ranges[key];
    if (!r) return;
    xlo = r[0]; xhi = r[1];
    markButtons(key);
    draw();
  }

  // ---- navigation ----
  function zoomAround(factor) {
    var c = (xlo + xhi) / 2;
    var half = (xhi - xlo) / 2 * factor;
    xlo = c - half; xhi = c + half;
    markButtons(null);
    draw();
  }

  function reset() {
    xlo = BG.ranges.all[0];
    xhi = BG.ranges.all[1];
    markButtons("all");
    draw();
  }

  function localX(clientX) {
    var r = canvas.getBoundingClientRect();
    return clientX - r.left;
  }

  // ---- events ----
  canvas.addEventListener("pointerdown", function (e) {
    if (e.button !== 0) return;
    var lx = localX(e.clientX);
    sel = { x0: lx, x1: lx };
    canvas.setPointerCapture && canvas.setPointerCapture(e.pointerId);
    draw();
  });

  canvas.addEventListener("pointermove", function (e) {
    if (sel) {
      sel.x1 = localX(e.clientX);
      draw();
      return;
    }
    var lx = localX(e.clientX);
    if (lx < M.l || lx > M.l + plotW) { hideTip(); return; }
    var en = pxToX(lx);
    var idx = nearestIndex(en);
    cursorIdx = idx;
    showTip(idx, e.clientX, e.clientY);
    draw();
  });

  canvas.addEventListener("pointerleave", function () {
    if (!sel) hideTip();
  });

  // mouseup on document (so a drag that leaves the canvas still commits)
  document.addEventListener("mouseup", function (e) {
    if (!sel) return;
    var a = Math.min(sel.x0, sel.x1), b = Math.max(sel.x0, sel.x1);
    sel = null;
    if (b - a >= 6) {
      var e0 = pxToX(a), e1 = pxToX(b);
      if (e1 > e0) { xlo = e0; xhi = e1; markButtons(null); }
    }
    draw();
  });

  canvas.addEventListener("dblclick", function () { reset(); });

  canvas.addEventListener("keydown", function (e) {
    if (e.key === "Escape") { reset(); }
    else if (e.key === "+" || e.key === "=") { zoomAround(1 / 1.5); }
    else if (e.key === "-" || e.key === "_") { zoomAround(1.5); }
  });

  rangeBox.addEventListener("click", function (e) {
    var btn = e.target.closest(".btn[data-range]");
    if (!btn) return;
    setRange(btn.getAttribute("data-range"));
  });

  // ---- sizing / theme ----
  function resize() {
    dpr = window.devicePixelRatio || 1;
    var r = canvas.getBoundingClientRect();
    W = Math.max(1, Math.round(r.width));
    H = Math.max(1, Math.round(r.height));
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    plotW = Math.max(1, W - M.l - M.r);
    plotH = Math.max(1, H - M.t - M.b);
    readColors();
    draw();
  }

  window.addEventListener("resize", resize);
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  if (mq.addEventListener) mq.addEventListener("change", resize);
  else if (mq.addListener) mq.addListener(resize);

  // ---- init ----
  status.textContent = I18N.tip_hint || "";
  markButtons("all");
  resize();
})();
