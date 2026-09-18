(function () { "use strict";

const D = window.GAGG_MODEL;
if (!D) {
  document.getElementById("status").textContent = "Нет данных: файл gagg-data.js не загрузился.";
  return;
}

const E = D.E_keV;
const LT = D.meta.live_time_s;
const COLORS = ["--ch-photo", "--ch-compt_esc1", "--ch-pair_esc1", "--ch-brems_esc", "--ch-external", "--ch-xray_esc", "--ch-pair_esc2", "--ch-compt_escN", "--ch-pair_full", "--ch-compt_full", "--ch-other"];
const RANGES = { all: [0, 3000], low: [0, 600], alpha: [200, 500], high: [500, 3000] };
const M = { l: 78, r: 14, t: 36, b: 34 };
const LOW_FWHM_KEV = 238;
const MARKERS = [ { E: 349.8, t: "измерено RC-103G (экземпляр 1): 349,8" }, { E: 356.5, t: "максимум модели 356,5" } ];

const state = { scale: "log", smear: "on", unit: "rate", range: "all", zoom: null,
  on: {}, showTotal: true, cursor: null, drag: { active: false, kind: "", x0: 0, x1: 0 },
  meas: "on" };

const MEAS = window.GAGG_MEAS || null;
const showMeas = () => MEAS && state.meas === "on";

let k = 0;
const COMPS = D.components.map(c => ({
  id: c.id,
  label: c.label,
  act: c.activity_mBq_per_kg,
  weight: c.weight,
  rate: c.rate,
  rate_unsmeared: c.rate_unsmeared,
  inSum: c.weight > 0,
  color: c.weight > 0 ? COLORS[(k++) % COLORS.length] : "--acc-gray"
}));

COMPS.forEach(c => {
  state.on[c.id] = c.inSum;
});

const cv = document.getElementById("cv-main");
const tip = document.getElementById("tip-main");
const res = document.getElementById("cv-res");
const tipRes = document.getElementById("tip-res");

function currentRange() { return state.zoom ? [state.zoom.lo, state.zoom.hi] : RANGES[state.range]; }
function css(name, fallback) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback; }
function fmt(v, digits) {
  if (typeof v !== "number" || isNaN(v)) return "—";
  return Number(v).toLocaleString("ru-RU", { minimumFractionDigits: digits, maximumFractionDigits: digits }).replace(/-/g, "−");
}
function niceTicks(lo, hi, count) {
  const raw = (hi - lo) / count, pow = Math.floor(Math.log10(raw));
  let step = 10 * Math.pow(10, pow);
  for (const m of [1, 2, 5, 10]) { if (m * Math.pow(10, pow) >= raw) { step = m * Math.pow(10, pow); break; } }
  const t = [];
  for (let k = Math.ceil(lo / step); k <= Math.floor(hi / step); k++) t.push(k * step);
  return t;
}
function fitCanvas(cv) {
  const dpr = devicePixelRatio || 1, r = cv.getBoundingClientRect();
  const w = Math.max(200, Math.floor(r.width)), h = Math.max(120, Math.floor(r.height));
  cv.width = Math.floor(w * dpr); cv.height = Math.floor(h * dpr);
  const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, h);
  return { g, w, h };
}

function sci(v, sig) {
  if (v === 0) return "0";
  const exp = Math.floor(Math.log10(Math.abs(v)));
  const mantissa = v / Math.pow(10, exp);
  const mStr = fmt(mantissa, sig - 1);
  const expStr = exp.toString().replace(/-/g, "⁻").replace(/0/g, "⁰").replace(/1/g, "¹").replace(/2/g, "²").replace(/3/g, "³").replace(/4/g, "⁴").replace(/5/g, "⁵").replace(/6/g, "⁶").replace(/7/g, "⁷").replace(/8/g, "⁸").replace(/9/g, "⁹");
  return exp === 0 ? mStr : `${mStr}·10${expStr}`;
}

function mult() { return state.unit === "counts" ? LT : 1; }

function arr(c) { return state.smear === "on" ? c.rate : c.rate_unsmeared; }
function totalArr() { return state.smear === "on" ? D.total : D.total_unsmeared; }

function visibleIdx() {
  const [lo, hi] = currentRange();
  const idx = [];
  for (let i = 0; i < E.length; i++) if (E[i] >= lo && E[i] <= hi) idx.push(i);
  return idx;
}

function unitLabel() { return state.unit === "rate" ? "отсч./(с·кэВ)" : "отсч./кэВ за 668 023,5 с"; }

function drawChart() {
  const { g, w, h } = fitCanvas(cv);
  const x0 = M.l, x1 = w - M.r, y0 = M.t, y1 = h - M.b;
  const [xlo, xhi] = currentRange();
  const X = e => x0 + (e - xlo) / (xhi - xlo) * (x1 - x0);
  const idx = visibleIdx();
  if (idx.length < 2) return;

  const stackComps = COMPS.filter(c => c.inSum && state.on[c.id]).sort((a, b) => {
    let sa = 0, sb = 0;
    for (let i = 0; i < 3000; i++) { sa += arr(a)[i]; sb += arr(b)[i]; }
    return sb - sa;
  });

  const cum = [];
  if (stackComps.length > 0) {
    for (let k = 0; k < stackComps.length; k++) {
      cum[k] = [];
      for (let i = 0; i < 3000; i++) {
        cum[k][i] = (k === 0 ? 0 : cum[k - 1][i]) + arr(stackComps[k])[i] * mult();   // правка автора: стопка в единицах графика
      }
    }
  }

  let hiV = 0;
  if (state.showTotal) {
    for (let i of idx) hiV = Math.max(hiV, totalArr()[i] * mult());
  }
  for (let k = 0; k < stackComps.length; k++) {
    for (let i of idx) hiV = Math.max(hiV, cum[k][i]);
  }
  if (showMeas()) {
    for (let i of idx) hiV = Math.max(hiV, MEAS.rate[i] * mult());
  }
  for (let c of COMPS) {
    if (!c.inSum && state.on[c.id]) {
      for (let i of idx) hiV = Math.max(hiV, arr(c)[i] * mult());
    }
  }

  if (hiV <= 0) hiV = 1;

  let lo, hi, map;
  if (state.scale === "log") {
    hi = hiV * 2;
    lo = hi * 1e-6;
    map = v => y0 + (1 - (Math.log10(Math.max(v, lo)) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo))) * (y1 - y0);
  } else {
    lo = 0;
    hi = hiV * 1.08;
    map = v => y0 + (hi - v) / (hi - lo) * (y1 - y0);
  }

  if (xlo < LOW_FWHM_KEV) {
    g.fillStyle = css("--rule-soft", "#c9c4b4");
    g.globalAlpha = 0.18;
    g.fillRect(X(xlo), y0, X(Math.min(LOW_FWHM_KEV, xhi)) - X(xlo), y1 - y0);
    g.globalAlpha = 1;
  }

  // правка автора: в лог-режиме деления — степени десяти
  const yTicks = [];
  if (state.scale === "log") {
    for (let p = Math.ceil(Math.log10(lo)); p <= Math.floor(Math.log10(hi)); p++) yTicks.push(Math.pow(10, p));
  } else yTicks.push(...niceTicks(0, hi, 5));

  g.strokeStyle = css("--grid", "#ccc");
  g.lineWidth = 1;
  for (let t of yTicks) {
    const y = map(t);
    g.beginPath();
    g.moveTo(x0, y);
    g.lineTo(x1, y);
    g.stroke();
  }

  const xTicks = niceTicks(xlo, xhi, w < 520 ? 4 : 7);

  for (let t of xTicks) {
    const x = X(t);
    g.beginPath();
    g.moveTo(x, y0);
    g.lineTo(x, y1);
    g.stroke();
  }

  if (stackComps.length > 0) {
    for (let k = 0; k < stackComps.length; k++) {
      const color = stackComps[k].color;
      g.fillStyle = css(color);
      g.globalAlpha = 0.75;
      g.beginPath();
      let first = true;
      for (let i of idx) {
        if (first) { g.moveTo(X(E[i]), map(cum[k][i])); first = false; }
        else { g.lineTo(X(E[i]), map(cum[k][i])); }
      }
      if (k > 0) {
        for (let i = idx.length - 1; i >= 0; i--) {
          g.lineTo(X(E[idx[i]]), map(cum[k - 1][idx[i]]));
        }
      } else {
        const floor = state.scale === "log" ? lo : 0;
        for (let i = idx.length - 1; i >= 0; i--) {
          g.lineTo(X(E[idx[i]]), map(floor));
        }
      }
      g.closePath();
      g.fill();
      g.globalAlpha = 1;
      g.strokeStyle = css(color);
      g.lineWidth = 0.6;
      g.beginPath();
      first = true;
      for (let i of idx) {
        if (first) { g.moveTo(X(E[i]), map(cum[k][i])); first = false; }
        else { g.lineTo(X(E[i]), map(cum[k][i])); }
      }
      g.stroke();
    }
  }

  for (let c of COMPS) {
    if (!c.inSum && state.on[c.id]) {
      g.strokeStyle = css(c.color);
      g.lineWidth = 1.6;
      g.setLineDash([6, 4]);
      g.beginPath();
      let first = true;
      for (let i of idx) {
        const y = map(arr(c)[i] * mult());
        if (first) { g.moveTo(X(E[i]), y); first = false; }
        else { g.lineTo(X(E[i]), y); }
      }
      g.stroke();
      g.setLineDash([]);
    }
  }

  if (state.showTotal) {
    g.strokeStyle = css("--sum-line");
    g.lineWidth = 1.6;
    g.beginPath();
    let first = true;
    for (let i of idx) {
      const y = map(totalArr()[i] * mult());
      if (first) { g.moveTo(X(E[i]), y); first = false; }
      else { g.lineTo(X(E[i]), y); }
    }
    g.stroke();
  }

  if (showMeas()) {
    g.strokeStyle = css("--ink");
    g.lineWidth = 1;
    g.globalAlpha = 0.9;
    g.beginPath();
    let first = true;
    for (let i of idx) {
      const y = map(MEAS.rate[i] * mult());
      if (first) { g.moveTo(X(E[i] - 0.5), y); first = false; }
      else { g.lineTo(X(E[i] - 0.5), y); }
      g.lineTo(X(E[i] + 0.5), y);
    }
    g.stroke();
    g.globalAlpha = 1;
  }

  for (let m of MARKERS) {
    if (m.E >= xlo && m.E <= xhi) {
      const x = X(m.E);
      g.strokeStyle = css("--rule-soft", "#c9c4b4");
      g.lineWidth = 1;
      g.setLineDash([4, 3]);
      g.beginPath();
      g.moveTo(x, y0);
      g.lineTo(x, y1);
      g.stroke();
      g.setLineDash([]);
    }
  }

  const markerLabels = [];
  for (let m of MARKERS) {
    if (m.E >= xlo && m.E <= xhi) {
      const x = X(m.E);
      let row = 0;
      let lastEnd = -100;
      for (let i = 0; i < markerLabels.length; i++) {
        if (markerLabels[i].x + markerLabels[i].w < x - 6) continue;
        if (markerLabels[i].row === 0) { row = 1; lastEnd = markerLabels[i].x + markerLabels[i].w; }
        else if (markerLabels[i].row === 1 && x - 6 > lastEnd) { row = 1; lastEnd = markerLabels[i].x + markerLabels[i].w; }
        else { row = -1; break; }
      }
      if (row !== -1) {
        const text = m.t;
        g.textAlign = "center";
        g.textBaseline = "bottom";
        g.fillStyle = css("--dim");
        g.font = "11px monospace";
        const w = g.measureText(text).width;
        markerLabels.push({ x, y: y0 - 5 - row * 13, text, w, row });
      }
    }
  }

  for (let m of markerLabels) {
    g.fillText(m.text, m.x, m.y);
  }

  g.textAlign = "right";
  g.textBaseline = "middle";
  g.fillStyle = css("--dim");
  g.font = "11px monospace";
  for (let t of yTicks) {
    const y = map(t);
    g.fillText(sci(t, state.scale === "log" ? 1 : 2), x0 - 6, y);
  }

  g.textAlign = "center";
  g.textBaseline = "top";
  for (let t of xTicks) {
    const x = X(t);
    g.fillText(fmt(t, 0), x, y1 + 6);
  }

  g.textAlign = "right";
  g.textBaseline = "bottom";
  g.fillText("кэВ-экв.", x1, h - 2);

  g.textAlign = "left";
  g.textBaseline = "alphabetic";
  g.fillText(unitLabel(), x0, 14);

  g.strokeStyle = css("--rule");
  g.lineWidth = 1.2;
  g.beginPath();
  g.moveTo(x0, y0);
  g.lineTo(x0, y1);
  g.lineTo(x1, y1);
  g.stroke();

  if (state.cursor !== null && !state.drag.active) {
    const x = X(E[state.cursor]);
    g.strokeStyle = css("--rule");
    g.lineWidth = 1;
    g.setLineDash([4, 3]);
    g.globalAlpha = 0.7;
    g.beginPath();
    g.moveTo(x, y0);
    g.lineTo(x, y1);
    g.stroke();
    g.setLineDash([]);
    g.globalAlpha = 1;
  }

  if (state.drag.active && state.drag.kind === "main") {
    const xa = Math.max(x0, Math.min(x1, state.drag.x0));
    const xb = Math.max(x0, Math.min(x1, state.drag.x1));
    g.fillStyle = css("--accent-soft");
    g.globalAlpha = 0.22;
    g.fillRect(xa, y0, xb - xa, y1 - y0);
    g.globalAlpha = 1;
    g.strokeStyle = css("--rule");
    g.lineWidth = 1.5;
    g.setLineDash([4, 4]);
    g.strokeRect(xa, y0, xb - xa, y1 - y0);
    g.setLineDash([]);
  }
}

function drawRes() {
  if (!showMeas()) {
    const { g, w, h } = fitCanvas(res);
    const x0 = M.l, x1 = w - M.r, y0 = M.t, y1 = h - M.b;
    const [xlo, xhi] = currentRange();
    const X = e => x0 + (e - xlo) / (xhi - xlo) * (x1 - x0);
    g.textAlign = "center";
    g.textBaseline = "middle";
    g.fillStyle = css("--dim");
    g.font = "13px monospace";
    g.fillText("Измерение скрыто", (x0 + x1) / 2, (y0 + y1) / 2);
    return;
  }

  const { g, w, h } = fitCanvas(res);
  const x0 = M.l, x1 = w - M.r, y0 = M.t, y1 = h - M.b;
  const [xlo, xhi] = currentRange();
  const X = e => x0 + (e - xlo) / (xhi - xlo) * (x1 - x0);
  const idx = visibleIdx();

  let A = 0;
  for (let i of idx) {
    A = Math.max(A, Math.abs(MEAS.rate[i] - totalArr()[i]) * mult());
    A = Math.max(A, MEAS.rate_sigma[i] * mult());
  }
  if (A === 0) A = 1;
  A *= 1.15;

  const map = v => y0 + (A - v) / (2 * A) * (y1 - y0);

  if (xlo < LOW_FWHM_KEV) {
    g.fillStyle = css("--rule-soft", "#c9c4b4");
    g.globalAlpha = 0.18;
    g.fillRect(X(xlo), y0, X(Math.min(LOW_FWHM_KEV, xhi)) - X(xlo), y1 - y0);
    g.globalAlpha = 1;
  }

  const yTicks = niceTicks(-A, A, 4);

  g.strokeStyle = css("--grid", "#ccc");
  g.lineWidth = 1;
  for (let t of yTicks) {
    const y = map(t);
    g.beginPath();
    g.moveTo(x0, y);
    g.lineTo(x1, y);
    g.stroke();
  }

  g.fillStyle = css("--grid");
  g.globalAlpha = 0.9;
  g.beginPath();
  let first = true;
  for (let i of idx) {
    const y = map(MEAS.rate_sigma[i] * mult());
    if (first) { g.moveTo(X(E[i]), y); first = false; }
    else { g.lineTo(X(E[i]), y); }
  }
  for (let i = idx.length - 1; i >= 0; i--) {
    const y = map(-MEAS.rate_sigma[idx[i]] * mult());
    g.lineTo(X(E[idx[i]]), y);
  }
  g.closePath();
  g.fill();
  g.globalAlpha = 1;

  g.strokeStyle = css("--rule-soft");
  g.lineWidth = 1;
  g.beginPath();
  g.moveTo(x0, map(0));
  g.lineTo(x1, map(0));
  g.stroke();

  g.strokeStyle = css("--acc-blue");
  g.lineWidth = 1;
  g.beginPath();
  first = true;   // правка автора: повторное объявление let first
  for (let i of idx) {
    const y = map((MEAS.rate[i] - totalArr()[i]) * mult());
    if (first) { g.moveTo(X(E[i] - 0.5), y); first = false; }
    else { g.lineTo(X(E[i] - 0.5), y); }
    g.lineTo(X(E[i] + 0.5), y);
  }
  g.stroke();

  for (let m of MARKERS) {
    if (m.E >= xlo && m.E <= xhi) {
      const x = X(m.E);
      g.strokeStyle = css("--rule-soft", "#c9c4b4");
      g.lineWidth = 1;
      g.setLineDash([4, 3]);
      g.beginPath();
      g.moveTo(x, y0);
      g.lineTo(x, y1);
      g.stroke();
      g.setLineDash([]);
    }
  }

  g.textAlign = "right";
  g.textBaseline = "middle";
  g.fillStyle = css("--dim");
  g.font = "11px monospace";
  for (let t of yTicks) {
    const y = map(t);
    g.fillText(sci(t, 2), x0 - 6, y);
  }

  g.textAlign = "center";
  g.textBaseline = "top";
  const xTicks = niceTicks(xlo, xhi, w < 520 ? 4 : 7);
  for (let t of xTicks) {
    const x = X(t);
    g.fillText(fmt(t, 0), x, y1 + 6);
  }

  g.textAlign = "right";
  g.textBaseline = "bottom";
  g.fillText("кэВ-экв.", x1, h - 2);

  g.textAlign = "left";
  g.textBaseline = "alphabetic";
  g.fillText(unitLabel() + " · остаток = измерение − сумма модели", x0, 14);

  g.strokeStyle = css("--rule");
  g.lineWidth = 1.2;
  g.beginPath();
  g.moveTo(x0, y0);
  g.lineTo(x0, y1);
  g.lineTo(x1, y1);
  g.stroke();

  if (state.cursor !== null && !state.drag.active) {
    const x = X(E[state.cursor]);
    g.strokeStyle = css("--rule");
    g.lineWidth = 1;
    g.setLineDash([4, 3]);
    g.globalAlpha = 0.7;
    g.beginPath();
    g.moveTo(x, y0);
    g.lineTo(x, y1);
    g.stroke();
    g.setLineDash([]);
    g.globalAlpha = 1;
  }

  if (state.drag.active && state.drag.kind === "res") {
    const xa = Math.max(x0, Math.min(x1, state.drag.x0));
    const xb = Math.max(x0, Math.min(x1, state.drag.x1));
    g.fillStyle = css("--accent-soft");
    g.globalAlpha = 0.22;
    g.fillRect(xa, y0, xb - xa, y1 - y0);
    g.globalAlpha = 1;
    g.strokeStyle = css("--rule");
    g.lineWidth = 1.5;
    g.setLineDash([4, 4]);
    g.strokeRect(xa, y0, xb - xa, y1 - y0);
    g.setLineDash([]);
  }
}

function drawAll() {
  drawChart();
  drawRes();
}

// правка автора: без innerHTML, по строке на ряд; вне суммы — внутри цикла (было обращение к c вне цикла)
function tipLines(el, lines) {
  el.replaceChildren();
  lines.forEach((t, j) => { const n = document.createElement(j ? "div" : "b"); n.textContent = t; el.appendChild(n); });
}
function fillTip(i) {
  const L = [`E ${fmt(E[i], 1)} кэВ-экв.`];
  if (state.showTotal) L.push(`Сумма: ${sci(totalArr()[i] * mult(), 4)}`);
  if (showMeas()) L.push(`Измерение: ${sci(MEAS.rate[i] * mult(), 4)} ± ${sci(MEAS.rate_sigma[i] * mult(), 2)}`);
  for (let c of COMPS) if (c.inSum && state.on[c.id]) L.push(`${c.label}: ${sci(arr(c)[i] * mult(), 4)}`);
  for (let c of COMPS) if (!c.inSum && state.on[c.id]) L.push(`${c.label} (вне суммы): ${sci(arr(c)[i] * mult(), 4)}`);
  tipLines(tip, L);
}

function fillTipRes(i) {
  tipLines(tipRes, [`E ${fmt(E[i], 1)} кэВ-экв.`,
    `Остаток: ${sci((MEAS.rate[i] - totalArr()[i]) * mult(), 4)}`,
    `Измерение: ${sci(MEAS.rate[i] * mult(), 4)}`,
    `Сумма модели: ${sci(totalArr()[i] * mult(), 4)}`,
    `σ измерения: ${sci(MEAS.rate_sigma[i] * mult(), 2)}`]);
}

function pointermove(ev) {
  const r = cv.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;
  if (state.drag.active && state.drag.kind === "main") {
    state.drag.x1 = Math.max(M.l, Math.min(r.width - M.r, x));
    drawChart();
    return;
  }
  if (x < M.l || x > r.width - M.r) { hide(); return; }
  const [xlo, xhi] = currentRange();
  const e = xlo + (x - M.l) / (r.width - M.r - M.l) * (xhi - xlo);
  const idx = visibleIdx();
  let nearest = null;
  let minDist = Infinity;
  for (let i of idx) {
    const dist = Math.abs(E[i] - e);
    if (dist < minDist) { minDist = dist; nearest = i; }
  }
  state.cursor = nearest;
  drawAll();   // правка автора: курсор общий — перерисовать оба графика
  fillTip(nearest);
  tip.hidden = false;
  tip.style.left = x + "px";
  tip.style.top = Math.max(0, y) + "px";
}

function pointerleave() {
  if (!state.drag.active) hide();
}

function mousedown(ev) {
  const r = cv.getBoundingClientRect(), x = ev.clientX - r.left;
  if (x < M.l || x > r.width - M.r) return;
  ev.preventDefault();
  state.drag = { active: true, kind: "main", x0: x, x1: x };
  tip.hidden = true;
  drawChart();
}

function dblclick() {
  state.zoom = null;
  state.range = "all";
  syncRangeButtons();
  render();
}

function pointermoveRes(ev) {
  const r = res.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;
  if (state.drag.active && state.drag.kind === "res") {
    state.drag.x1 = Math.max(M.l, Math.min(r.width - M.r, x));
    drawRes();
    return;
  }
  if (x < M.l || x > r.width - M.r) { hide(); return; }
  const [xlo, xhi] = currentRange();
  const e = xlo + (x - M.l) / (r.width - M.r - M.l) * (xhi - xlo);
  const idx = visibleIdx();
  let nearest = null;
  let minDist = Infinity;
  for (let i of idx) {
    const dist = Math.abs(E[i] - e);
    if (dist < minDist) { minDist = dist; nearest = i; }
  }
  state.cursor = nearest;
  drawAll();
  if (!showMeas()) return;
  fillTipRes(nearest);
  tipRes.hidden = false;
  tipRes.style.left = x + "px";
  tipRes.style.top = Math.max(0, y) + "px";
}

function pointerleaveRes() {
  if (!state.drag.active) hide();
}

function mousedownRes(ev) {
  const r = res.getBoundingClientRect(), x = ev.clientX - r.left;
  if (x < M.l || x > r.width - M.r) return;
  ev.preventDefault();
  state.drag = { active: true, kind: "res", x0: x, x1: x };
  tipRes.hidden = true;
  drawRes();
}

function hide() {
  state.cursor = null;
  tip.hidden = true;
  tipRes.hidden = true;
  drawAll();
}

function keydown(ev) {
  if (state.cursor === null) {
    const idx = visibleIdx();
    state.cursor = idx[Math.floor(idx.length / 2)];
  }
  let i = state.cursor;
  if (ev.key === "ArrowLeft") {
    const idx = visibleIdx();
    const j = idx.indexOf(i);
    if (j > 0) {
      i = idx[j - 1];
      if (ev.shiftKey) {
        for (let k = 0; k < 10 && j - 1 - k >= 0; k++) {
          i = idx[j - 1 - k];
        }
      }
    }
  } else if (ev.key === "ArrowRight") {
    const idx = visibleIdx();
    const j = idx.indexOf(i);
    if (j < idx.length - 1) {
      i = idx[j + 1];
      if (ev.shiftKey) {
        for (let k = 0; k < 10 && j + 1 + k < idx.length; k++) {
          i = idx[j + 1 + k];
        }
      }
    }
  } else return;
  state.cursor = i;
  ev.preventDefault();
  drawAll();
  fillTip(i);
  const r = cv.getBoundingClientRect();
  const x = M.l + (E[i] - currentRange()[0]) / (currentRange()[1] - currentRange()[0]) * (r.width - M.l - M.r);
  tip.hidden = false;
  tip.style.left = x + "px";
  tip.style.top = "0px";
}

function bindSeg(id, prop, onSet) {
  document.querySelectorAll(`#${id} button[data-v]`).forEach(btn => {
    btn.addEventListener("click", () => {
      state[prop] = btn.dataset.v;
      if (onSet) onSet();
      document.querySelectorAll(`#${id} button`).forEach(b =>
        b.setAttribute("aria-pressed", String(b === btn)));
      render();
    });
  });
}

function syncRangeButtons() {
  document.querySelectorAll("#ctl-range button[data-v]").forEach(b =>
    b.setAttribute("aria-pressed", String(!state.zoom && b.dataset.v === state.range)));
}

function buildLegend() {
  const legend = document.getElementById("legend-main");
  legend.replaceChildren();

  const totalChip = document.createElement("label");
  totalChip.className = "chip";
  const totalInput = document.createElement("input");
  totalInput.type = "checkbox";
  totalInput.checked = state.showTotal;
  totalInput.addEventListener("change", () => {
    state.showTotal = totalInput.checked;
    render();
  });
  const totalSw = document.createElement("span");
  totalSw.className = "sw sum";
  const totalNm = document.createElement("span");
  totalNm.className = "nm";
  totalNm.textContent = "Сумма";
  totalChip.appendChild(totalInput);
  totalChip.appendChild(totalSw);
  totalChip.appendChild(totalNm);
  legend.appendChild(totalChip);

  for (let c of COMPS) {
    const chip = document.createElement("label");
    chip.className = "chip";
    const input = document.createElement("input");
    input.type = "checkbox";
    input.checked = state.on[c.id];
    input.addEventListener("change", () => {
      state.on[c.id] = input.checked;
      render();
    });
    const sw = document.createElement("span");
    sw.className = "sw";
    sw.style.backgroundColor = css(c.color, "#888");
    const nm = document.createElement("span");
    nm.className = "nm";
    nm.textContent = c.label + (c.weight === 0 ? " · вне суммы" : "");
    chip.appendChild(input);
    chip.appendChild(sw);
    chip.appendChild(nm);
    legend.appendChild(chip);
  }
}

function buildTable() {
  const box = document.getElementById("comp-table");   // правка автора: строки — в <table class="big">, не прямо в div
  const table = document.createElement("table");
  table.className = "big";

  const [lo, hi] = currentRange();
  const idx = visibleIdx();

  const rangeText = document.getElementById("table-range");
  rangeText.textContent = `Диапазон: ${fmt(lo, 1)}–${fmt(hi, 1)} кэВ-экв., размытие: ${state.smear === "on" ? "да" : "нет"}`;

  const thead = document.createElement("thead");
  const tr = document.createElement("tr");
  const ths = ["Нуклид / подцепочка", "Активность, мБк/кг", "Вес", "Скорость, с⁻¹", "Отсчёты за 668 023,5 с"];
  for (let t of ths) {
    const th = document.createElement("th");
    th.textContent = t;
    tr.appendChild(th);
  }
  thead.appendChild(tr);
  table.appendChild(thead);

  const tbody = document.createElement("tbody");

  let totalRate = 0;
  let totalCounts = 0;

  for (let c of COMPS) {
    const row = document.createElement("tr");
    if (!c.inSum) row.className = "off";

    const firstCell = document.createElement("td");
    const sw = document.createElement("span");
    sw.className = "sw";
    sw.style.backgroundColor = css(c.color, "#888");
    const label = document.createElement("span");
    label.textContent = c.label + (c.weight === 0 ? " (вне суммы)" : "");
    firstCell.appendChild(sw);
    firstCell.appendChild(label);
    row.appendChild(firstCell);

    const actCell = document.createElement("td");
    actCell.className = "n";
    actCell.textContent = fmt(c.act, 2);
    row.appendChild(actCell);

    const weightCell = document.createElement("td");
    weightCell.className = "n";
    weightCell.textContent = fmt(c.weight, 0);
    row.appendChild(weightCell);

    let rateSum = 0;
    for (let i of idx) rateSum += arr(c)[i];
    const rateCell = document.createElement("td");
    rateCell.className = "n";
    rateCell.textContent = sci(rateSum, 5);
    row.appendChild(rateCell);

    const countCell = document.createElement("td");
    countCell.className = "n";
    countCell.textContent = fmt(rateSum * LT, 1);
    row.appendChild(countCell);

    tbody.appendChild(row);

  }
  // правка автора: строка «Сумма» — из total модели, не из строк таблицы
  for (let i of idx) totalRate += totalArr()[i];
  totalCounts = totalRate * LT;

  const sumRow = document.createElement("tr");
  sumRow.className = "sum";
  const firstSumCell = document.createElement("td");
  firstSumCell.textContent = "Сумма (вес 1)";
  sumRow.appendChild(firstSumCell);

  const actSumCell = document.createElement("td");
  actSumCell.className = "n";
  actSumCell.textContent = "—";
  sumRow.appendChild(actSumCell);

  const weightSumCell = document.createElement("td");
  weightSumCell.className = "n";
  weightSumCell.textContent = "—";
  sumRow.appendChild(weightSumCell);

  const rateSumCell = document.createElement("td");
  rateSumCell.className = "n";
  rateSumCell.textContent = sci(totalRate, 5);
  sumRow.appendChild(rateSumCell);

  const countSumCell = document.createElement("td");
  countSumCell.className = "n";
  countSumCell.textContent = fmt(totalCounts, 1);
  sumRow.appendChild(countSumCell);

  tbody.appendChild(sumRow);

  if (MEAS) {   // правка автора: строки измерения в таблице не зависят от переключателя
    const measRow = document.createElement("tr");
    measRow.className = "meas";
    const firstMeasCell = document.createElement("td");
    firstMeasCell.textContent = `Измерение, ${MEAS.meta.serial}`;
    measRow.appendChild(firstMeasCell);

    const actMeasCell = document.createElement("td");
    actMeasCell.className = "n";
    actMeasCell.textContent = "—";
    measRow.appendChild(actMeasCell);

    const weightMeasCell = document.createElement("td");
    weightMeasCell.className = "n";
    weightMeasCell.textContent = "—";
    measRow.appendChild(weightMeasCell);

    let measRateSum = 0;
    for (let i of idx) measRateSum += MEAS.rate[i];
    const rateMeasCell = document.createElement("td");
    rateMeasCell.className = "n";
    rateMeasCell.textContent = sci(measRateSum, 5);
    measRow.appendChild(rateMeasCell);

    let measCountSum = 0;
    for (let i of idx) measCountSum += MEAS.counts[i];
    const countMeasCell = document.createElement("td");
    countMeasCell.className = "n";
    countMeasCell.textContent = fmt(measCountSum, 1);
    measRow.appendChild(countMeasCell);

    tbody.appendChild(measRow);

    const residRow = document.createElement("tr");
    residRow.className = "meas";
    const firstResidCell = document.createElement("td");
    firstResidCell.textContent = "Остаток = измерение − сумма";
    residRow.appendChild(firstResidCell);

    const actResidCell = document.createElement("td");
    actResidCell.className = "n";
    actResidCell.textContent = "—";
    residRow.appendChild(actResidCell);

    const weightResidCell = document.createElement("td");
    weightResidCell.className = "n";
    weightResidCell.textContent = "—";
    residRow.appendChild(weightResidCell);

    let residRateSum = 0;
    for (let i of idx) residRateSum += MEAS.rate[i] - totalArr()[i];
    const rateResidCell = document.createElement("td");
    rateResidCell.className = "n";
    rateResidCell.textContent = sci(residRateSum, 5);
    residRow.appendChild(rateResidCell);

    const countResidCell = document.createElement("td");
    countResidCell.className = "n";
    countResidCell.textContent = fmt(residRateSum * LT, 1);
    residRow.appendChild(countResidCell);

    tbody.appendChild(residRow);
  }

  table.appendChild(tbody);
  box.replaceChildren(table);
}

function render() {
  drawAll();
  buildTable();
}

function wire(kind) {
  const canvas = document.getElementById(`cv-${kind}`);
  const tipEl = document.getElementById(`tip-${kind}`);
  canvas.addEventListener("pointermove", kind === "main" ? pointermove : pointermoveRes);
  canvas.addEventListener("pointerleave", kind === "main" ? pointerleave : pointerleaveRes);
  canvas.addEventListener("mousedown", kind === "main" ? mousedown : mousedownRes);
  canvas.addEventListener("dblclick", dblclick);
  canvas.addEventListener("keydown", keydown);
  canvas.addEventListener("blur", hide);
}

document.addEventListener("mouseup", () => {
  if (!state.drag.active) return;
  state.drag.active = false;
  const r = document.getElementById(`cv-${state.drag.kind}`).getBoundingClientRect();
  if (Math.abs(state.drag.x1 - state.drag.x0) < 6) { drawAll(); return; }
  const [xlo, xhi] = currentRange();
  const toE = px => xlo + (px - M.l) / (r.width - M.r - M.l) * (xhi - xlo);
  const lo = Math.max(0, toE(Math.min(state.drag.x0, state.drag.x1)));
  const hi = Math.min(3000, toE(Math.max(state.drag.x0, state.drag.x1)));
  if (hi - lo < 1) { drawAll(); return; }
  state.zoom = { lo, hi };
  syncRangeButtons();
  render();
});

bindSeg("ctl-scale", "scale");
bindSeg("ctl-smear", "smear");
bindSeg("ctl-unit", "unit");
bindSeg("ctl-range", "range", () => { state.zoom = null; });
bindSeg("ctl-meas", "meas");

buildLegend();

wire("main");
wire("res");

const hash = location.hash;
if (hash && /zoom=([\d.]+)-([\d.]+)/.test(hash)) {
  const match = hash.match(/zoom=([\d.]+)-([\d.]+)/);
  const lo = parseFloat(match[1]);
  const hi = parseFloat(match[2]);
  if (hi - lo >= 1) state.zoom = { lo: Math.max(0, lo), hi: Math.min(3000, hi) };
}

syncRangeButtons();

let lastWidth = 0;
let frame = 0;

const ro = new ResizeObserver(entries => {
  for (let entry of entries) {
    const width = entry.contentRect.width;
    if (width === lastWidth) continue;
    lastWidth = width;
    if (frame) return;
    frame = requestAnimationFrame(() => {
      render();
      frame = 0;
    });
  }
});

ro.observe(document.querySelector(".app"));

window.GAGG_APP = { state, render, currentRange };

render();

})();
