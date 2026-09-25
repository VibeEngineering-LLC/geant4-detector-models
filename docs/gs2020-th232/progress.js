(function(){ "use strict";
const D = window.GS2020_PROGRESS;
if (!D) {
  document.getElementById("meta").textContent = "Нет данных: data.js не загрузился.";
  return;
}

const state = { scale:"log", range:"all", zoom:null, cursor:{raw:null, net:null}, drag:{active:false, kind:null, x0:0, x1:0} };
const RANGES = { all:[0,3000], low:[0,800], high:[1300,2800] };
function currentRange(){ return state.zoom ? [state.zoom.lo, state.zoom.hi] : RANGES[state.range]; }

function fmt(v, d) {
  if (v === null || v === undefined || isNaN(v)) return "—";
  return v.toFixed(d).replace(".", ",");
}
function css(name, fallback) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
}

const M = { l:64, r:16, t:14, b:36 };
const COLORS = {
  ink: css("--ink"),
  faint: css("--faint"),
  grid: css("--grid"),
  ruleSoft: css("--rule-soft"),
  accBlue: css("--acc-blue"),
  accOrange: css("--acc-orange"),
  sumLine: css("--sum-line")
};

function fitCanvas(cv) {
  const dpr = devicePixelRatio || 1;
  const r = cv.getBoundingClientRect();
  const w = Math.max(200, Math.floor(r.width));
  const h = Math.max(120, Math.floor(r.height));
  cv.width = Math.floor(w * dpr);
  cv.height = Math.floor(h * dpr);
  const g = cv.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);
  return { g, w, h };
}

function visibleIdx(kind) {
  const [xlo, xhi] = currentRange();
  const data = D.spectrum;
  const idxs = [];
  for (let i = 0; i < data.E_keV.length; i++) {
    if (data.E_keV[i] >= xlo && data.E_keV[i] <= xhi) idxs.push(i);
  }
  return idxs;
}

function nearestIdx(E) {
  const data = D.spectrum.E_keV;
  let min = Infinity, idx = -1;
  for (let i = 0; i < data.length; i++) {
    const d = Math.abs(data[i] - E);
    if (d < min) { min = d; idx = i; }
  }
  return idx;
}

function niceTicks(lo, hi, n) {
  const span = hi - lo;
  const step = span / n;
  const logStep = Math.log10(step);
  const pow = Math.floor(logStep);
  const base = Math.pow(10, pow);
  let tickStep = 1;
  if (step <= base * 2) tickStep = 2;
  else if (step <= base * 5) tickStep = 5;
  else tickStep = 10;
  const start = Math.ceil(lo / (base * tickStep)) * base * tickStep;
  const ticks = [];
  for (let t = start; t <= hi; t += base * tickStep) {
    if (t >= lo && t <= hi) ticks.push(t);
  }
  return ticks;
}

function logTicks(lo, hi) {
  if (lo <= 0) lo = 1;
  const logLo = Math.log10(lo);
  const logHi = Math.log10(hi);
  const n = Math.max(2, Math.floor(logHi - logLo));
  const ticks = [];
  for (let i = 0; i <= n; i++) {
    const t = Math.pow(10, logLo + i * (logHi - logLo) / n);
    if (t >= lo && t <= hi) ticks.push(t);
  }
  return ticks;
}

function drawEff() {
  const cv = document.getElementById("cv-eff");
  const { g, w, h } = fitCanvas(cv);
  const data = D.eff;
  const xlo = 150, xhi = 3000;
  let yMin = Infinity, yMax = -Infinity;
  for (let i = 0; i < data.length; i++) {
    const p = data[i].eff_peak * 100;
    const t = data[i].eff_total * 100;
    if (p < yMin) yMin = p;
    if (t < yMin) yMin = t;
    if (p > yMax) yMax = p;
    if (t > yMax) yMax = t;
  }
  const yLo = 0.5 * yMin, yHi = 2 * yMax;
  const xScale = (w - M.l - M.r) / Math.log10(xhi / xlo);
  const yScale = (h - M.t - M.b) / Math.log10(yHi / yLo);

  g.strokeStyle = COLORS.grid;
  g.lineWidth = 0.5;
  g.beginPath();
  for (let i = 0; i < data.length; i++) {
    const x = M.l + xScale * Math.log10(data[i].E_keV / xlo);
    const yPeak = h - M.b - yScale * Math.log10((data[i].eff_peak * 100) / yLo);
    const yTotal = h - M.b - yScale * Math.log10((data[i].eff_total * 100) / yLo);
    g.moveTo(x, yPeak);
    g.lineTo(x, yPeak + 4);
    g.stroke();
    g.beginPath();
    g.arc(x, yPeak, 4, 0, Math.PI * 2);
    g.fillStyle = COLORS.accBlue;
    g.fill();
    g.beginPath();
    g.rect(x - 3.5, yTotal - 3.5, 7, 7);
    g.fillStyle = COLORS.accOrange;
    g.fill();
  }

  g.strokeStyle = COLORS.ink;
  g.lineWidth = 1;
  g.font = "10px sans-serif";
  g.textAlign = "center";
  g.textBaseline = "top";

  const xTicks = logTicks(xlo, xhi);
  for (let i = 0; i < xTicks.length; i++) {
    const x = M.l + xScale * Math.log10(xTicks[i] / xlo);
    g.beginPath();
    g.moveTo(x, M.t);
    g.lineTo(x, h - M.b);
    g.stroke();
    g.fillText(xTicks[i].toFixed(0), x, h - M.b + 2);
  }

  const yTicks = logTicks(yLo, yHi);
  for (let i = 0; i < yTicks.length; i++) {
    const y = h - M.b - yScale * Math.log10(yTicks[i] / yLo);
    g.beginPath();
    g.moveTo(M.l, y);
    g.lineTo(w - M.r, y);
    g.stroke();
    g.fillText(yTicks[i].toFixed(1), M.l - 5, y - 4);
  }

  g.fillStyle = COLORS.ink;
  g.fillText("E, кэВ", (w - M.l - M.r) / 2 + M.l, h - 2);
  g.fillText("эффективность, %", M.l - 30, (h - M.t - M.b) / 2 + M.t);
  g.fillText("пик (узкое окно)", w - M.r - 50, M.t + 10);
  g.fillText("полная", w - M.r - 50, M.t + 25);
  g.fillStyle = COLORS.accBlue;
  g.beginPath();
  g.arc(w - M.r - 60, M.t + 14, 4, 0, Math.PI * 2);
  g.fill();
  g.fillStyle = COLORS.accOrange;
  g.fillRect(w - M.r - 63.5, M.t + 20, 7, 7);
}

function drawChart(kind) {
  const cv = document.getElementById("cv-" + kind);
  const { g, w, h } = fitCanvas(cv);
  const data = D.spectrum;
  const [xlo, xhi] = currentRange();
  const idxs = visibleIdx(kind);
  if (idxs.length === 0) return;

  let yMin = Infinity, yMax = -Infinity;
  for (let i = 0; i < idxs.length; i++) {
    const idx = idxs[i];
    if (kind === "raw") {
      yMin = 0;
      if (data.counts[idx] > yMax) yMax = data.counts[idx];
      if (data.bg[idx] > yMax) yMax = data.bg[idx];
    } else {
      if (data.net[idx] < yMin) yMin = data.net[idx];
      if (data.net[idx] > yMax) yMax = data.net[idx];
    }
  }

  const xScale = (w - M.l - M.r) / (xhi - xlo);
  const yScale = (h - M.t - M.b) / (yMax - yMin);
  const useLog = kind === "raw" && state.scale === "log";
  if (useLog) yMin = Math.max(1, yMax * 1e-4);
  const Y = v => useLog
    ? h - M.b - (h - M.t - M.b) * Math.log10(Math.max(v, yMin) / yMin) / Math.log10(yMax / yMin)
    : h - M.b - yScale * (v - yMin);

  g.strokeStyle = COLORS.grid;
  g.lineWidth = 0.5;
  g.beginPath();
  for (let i = 0; i < idxs.length; i++) {
    const idx = idxs[i];
    const x = M.l + xScale * (data.E_keV[idx] - xlo);
    if (i === 0) g.moveTo(x, h - M.b);
    else g.lineTo(x, h - M.b);
  }
  g.stroke();

  g.strokeStyle = COLORS.ink;
  g.lineWidth = 1.2;
  g.beginPath();
  for (let i = 0; i < idxs.length; i++) {
    const idx = idxs[i];
    const x = M.l + xScale * (data.E_keV[idx] - xlo);
    const y = Y(data[kind === "raw" ? "counts" : "net"][idx]);
    if (i === 0) g.moveTo(x, y);
    else g.lineTo(x, y);
  }
  g.stroke();

  if (kind === "raw") {
    g.strokeStyle = COLORS.faint;
    g.beginPath();
    for (let i = 0; i < idxs.length; i++) {
      const idx = idxs[i];
      const x = M.l + xScale * (data.E_keV[idx] - xlo);
      const y = Y(data.bg[idx]);
      if (i === 0) g.moveTo(x, y);
      else g.lineTo(x, y);
    }
    g.stroke();
  }

  if (kind === "net") {
    g.strokeStyle = COLORS.sumLine;
    g.lineWidth = 1.2;
    g.setLineDash([4, 3]);
    const yZero = Y(0);
    g.beginPath();
    g.moveTo(M.l, yZero);
    g.lineTo(w - M.r, yZero);
    g.stroke();
    g.setLineDash([]);
  }

  if (state.drag.active && state.drag.kind === kind) {
    const x0 = Math.max(M.l, Math.min(w - M.r, state.drag.x0));
    const x1 = Math.max(M.l, Math.min(w - M.r, state.drag.x1));
    g.fillStyle = COLORS.ruleSoft;
    g.globalAlpha = 0.25;
    g.fillRect(x0, M.t, x1 - x0, h - M.t - M.b);
    g.globalAlpha = 1;
  }

  if (state.cursor[kind] !== null) {
    const x = M.l + xScale * (data.E_keV[state.cursor[kind]] - xlo);
    g.strokeStyle = COLORS.ink;
    g.lineWidth = 1;
    g.beginPath();
    g.moveTo(x, M.t);
    g.lineTo(x, h - M.b);
    g.stroke();
  }

  g.strokeStyle = COLORS.ink;
  g.lineWidth = 1;
  g.font = "10px sans-serif";
  g.textAlign = "center";
  g.textBaseline = "top";

  const xTicks = niceTicks(xlo, xhi, 5);
  for (let i = 0; i < xTicks.length; i++) {
    const x = M.l + xScale * (xTicks[i] - xlo);
    g.beginPath();
    g.moveTo(x, M.t);
    g.lineTo(x, h - M.b);
    g.stroke();
    g.fillText(xTicks[i].toFixed(0), x, h - M.b + 2);
  }

  if (kind === "raw" && state.scale === "log") {
    const yTicks = logTicks(yMin, yMax);
    for (let i = 0; i < yTicks.length; i++) {
      const y = Y(yTicks[i]);
      g.beginPath();
      g.moveTo(M.l, y);
      g.lineTo(w - M.r, y);
      g.stroke();
      g.fillText(yTicks[i].toFixed(0), M.l - 5, y - 4);
    }
  } else {
    const yTicks = niceTicks(yMin, yMax, 5);
    for (let i = 0; i < yTicks.length; i++) {
      const y = Y(yTicks[i]);
      g.beginPath();
      g.moveTo(M.l, y);
      g.lineTo(w - M.r, y);
      g.stroke();
      g.fillText(yTicks[i].toFixed(0), M.l - 5, y - 4);
    }
  }

  g.fillStyle = COLORS.ink;
  g.fillText("E, кэВ", (w - M.l - M.r) / 2 + M.l, h - 2);
}

function show(idx, x, y) {
  const data = D.spectrum;
  const E = data.E_keV[idx];
  const rawTip = document.getElementById("tip-raw");
  const netTip = document.getElementById("tip-net");
  if (state.cursor.raw !== null) {
    rawTip.textContent = `E ${E.toFixed(1)} кэВ · измерение ${data.counts[idx]} · фон·k ${data.bg[idx]} (на кэВ)`;
    rawTip.style.left = x + "px";
    rawTip.style.top = y + "px";
    rawTip.hidden = false;
  }
  if (state.cursor.net !== null) {
    netTip.textContent = `E ${E.toFixed(1)} кэВ · нетто ${data.net[idx]} на кэВ`;
    netTip.style.left = x + "px";
    netTip.style.top = y + "px";
    netTip.hidden = false;
  }
}

function hide() {
  document.getElementById("tip-raw").hidden = true;
  document.getElementById("tip-net").hidden = true;
}

function attachCursor(kind) {
  const cv = document.getElementById("cv-" + kind);
  cv.addEventListener("pointermove", (ev) => {
    const r = cv.getBoundingClientRect();
    const x = ev.clientX - r.left;
    const y = ev.clientY - r.top;
    if (state.drag.active && state.drag.kind === kind) { state.drag.x1 = x; drawChart(kind); return; }
    if (x < M.l || x > r.width - M.r) { hide(); return; }
    const [xlo, xhi] = currentRange();
    const E = xlo + (x - M.l) / (r.width - M.r - M.l) * (xhi - xlo);
    const idx = nearestIdx(E);
    state.cursor[kind] = idx;
    show(idx, x, y);
    render();
  });
  cv.addEventListener("pointerleave", () => {
    hide();
    state.cursor[kind] = null;
    render();
  });
  cv.addEventListener("mousedown", (ev) => {
    const r = cv.getBoundingClientRect();
    const x = ev.clientX - r.left;
    if (x < M.l || x > r.width - M.r) return;
    ev.preventDefault();
    state.drag = { active: true, kind: kind, x0: x, x1: x };
    hide();
    drawChart(kind);
  });
  cv.addEventListener("dblclick", () => {
    state.zoom = null;
    state.range = "all";
    syncRangeButtons();
    render();
  });
}

function syncRangeButtons() {
  const rangeBtns = document.querySelectorAll("#ctl-range button");
  rangeBtns.forEach(btn => {
    btn.setAttribute("aria-pressed", btn.dataset.v === state.range ? "true" : "false");
  });
  const scaleBtns = document.querySelectorAll("#ctl-scale button");
  scaleBtns.forEach(btn => {
    btn.setAttribute("aria-pressed", btn.dataset.v === state.scale ? "true" : "false");
  });
}

function bindSeg(id, prop, onSet) {
  const container = document.getElementById(id);
  container.addEventListener("click", (ev) => {
    if (!ev.target.matches("button")) return;
    const val = ev.target.dataset.v;
    state[prop] = val;
    if (onSet) onSet();
    syncRangeButtons();
    render();
  });
}

function render() {
  drawEff();
  drawChart("raw");
  drawChart("net");
  (window.GS2020_VIEW.listeners || []).forEach(fn => fn());
}
// общий диапазон энергии для прочих графиков страницы (#CHART-1 п.5: синхронный зум)
window.GS2020_VIEW = { listeners: [], range: currentRange,
  setZoom(lo, hi) { state.zoom = { lo, hi }; syncRangeButtons(); render(); },
  reset() { state.zoom = null; state.range = "all"; syncRangeButtons(); render(); } };

document.getElementById("meta").textContent =
  `Выгрузка: ${D.generated}. Живое время спектра ${fmt(D.live_s, 0)} с, фона ${fmt(D.live_bg_s, 0)} с, k фона = ${fmt(D.k_bg, 5)} (фиксирован). Паспорт КИ: ${fmt(D.passport.bq_per_kg, 0)} Бк/кг ± ${fmt(D.passport.unc_pct, 0)} %, ${fmt(D.passport.mass_g, 0)} г → ${fmt(D.activity_passport_bq, 1)} Бк.`;

const stages = document.getElementById("stages");
D.stages.forEach(stage => {
  const tr = document.createElement("tr");
  const statusDict = { done:"готово", rejected:"отбраковано", work:"в работе", wait:"ждёт" };
  tr.innerHTML = `<td>${stage.n}</td><td>${stage.name}</td><td><span class="st st-${stage.status}">${statusDict[stage.status]}</span></td><td>${stage.note}</td>`;
  stages.appendChild(tr);
});

const eff = document.getElementById("eff");
D.eff.forEach(row => {
  const tr = document.createElement("tr");
  tr.innerHTML = `<td>${row.nuclide}</td><td class="num">${fmt(row.E_keV, 2)}</td><td class="num">${fmt(row.eff_peak * 100, 4)} ± ${fmt(row.d_eff_peak * 100, 4)}</td><td class="num">${fmt(row.eff_total * 100, 3)}</td><td class="num">${row.n_events.toLocaleString("ru-RU")}</td><td class="num">${row.geomnav || "—"} </td>`;
  eff.appendChild(tr);
});

attachCursor("raw");
attachCursor("net");

bindSeg("ctl-range", "range", () => { state.zoom = null; });
bindSeg("ctl-scale", "scale");

document.addEventListener("mouseup", () => {
  if (!state.drag.active) return;
  const kind = state.drag.kind;
  state.drag.active = false;
  if (Math.abs(state.drag.x1 - state.drag.x0) < 6) { drawChart(kind); return; }
  const r = document.getElementById("cv-" + kind).getBoundingClientRect();
  const [xlo, xhi] = currentRange();
  const toE = (px) => xlo + (px - M.l) / (r.width - M.r - M.l) * (xhi - xlo);
  const lo = Math.max(0, toE(Math.min(state.drag.x0, state.drag.x1)));
  const hi = Math.min(3000, toE(Math.max(state.drag.x0, state.drag.x1)));
  if (hi - lo < 1) { drawChart(kind); return; }
  state.zoom = { lo, hi };
  syncRangeButtons();
  render();
});

new ResizeObserver(() => {
  requestAnimationFrame(render);
}).observe(document.querySelector(".app"));

render();

})();
