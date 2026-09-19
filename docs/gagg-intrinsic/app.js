(function () { "use strict";

// разделы «Исходный состав» и «Что и как считалось» — отдельные фрагменты, вставляются перед «Оговорками»
[["sec-composition", "sec-composition.html"], ["sec-method", "sec-method.html"]].forEach(([id, file]) => {
  const box = document.getElementById(id);
  if (!box) return;
  fetch(file).then(r => { if (!r.ok) throw new Error(String(r.status)); return r.text(); })
    .then(t => box.insertAdjacentHTML("beforeend", t))
    .catch(() => { box.textContent = "раздел не загрузился"; });
});

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
// максимум размытой суммы в окне 300–400 кэВ — из данных страницы (model_sum_B2.json), не набран руками
const PK = (() => { let b = E.findIndex(e => e >= 300); for (let i = b; i < E.length && E[i] <= 400; i++) if (D.total[i] > D.total[b]) b = i; return E[b]; })();
const MARKERS = [ { E: 349.8, t: "измерено RC-103G (экземпляр 1): 349,8", short: "измерено 349,8" },
                  { E: PK, t: `максимум модели ${fmt(PK, 1)}`, short: `модель ${fmt(PK, 1)}` } ];

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

// Основание нормировки активности каждого шаблона (Belli = Belli et al., J. Phys. G 53 (2026) 045101, табл. 1)
const SRC = {
  T1: "расчёт: природный Gd, доля Gd-152 0,20 % (IUPAC), T½ 3,408·10²¹ с",
  T2: "Belli: сумма Gd-152 + Sm-147 = 1555 мБк/кг минус расчёт Gd-152 (примесь конкретного кристалла)",
  T3: "верхний предел по измерению RC-103G; у Belli 11 026 мБк/кг",
  T4: "Belli, табл. 1",
  T5: "= Th-228 (допущение; Belli для Ra-228: ≤ 4 мБк/кг)",
  T6: "Belli, табл. 1",
  T7: "Belli, табл. 1",
  T8: "Belli, табл. 1",
  T9: "Belli, табл. 1",
  T10: "Belli, табл. 1",
  T11: "Belli, табл. 1 (верхний предел)",
  T12: "Belli, табл. 1"
};

// Тип распада и что он даёт в спектре; положения — свет α варианта B2, кэВ-экв. (results/birks_decay/lines_table.csv)
const PHYS = {
  T1: ["α", "α 2146 кэВ поглощается целиком, свет гасится (Биркс) → пик ≈ 349,5"],
  T2: ["α", "α 2248 кэВ, гашение → пик ≈ 367,5; вместе с Gd-152 дают главный пик ~357"],
  T3: ["β⁻ + γ", "β⁻ на уровень 597 кэВ ¹⁷⁶Hf, каскад γ 307/202/88 поглощается в том же импульсе → широкий горб 150–850, узких пиков нет"],
  T4: ["α", "α 4011 кэВ → ≈ 726,5; ветвь на уровень 64 кэВ — гашёная α плюс разрядка уровня (не гасится)"],
  T5: ["β⁻ + γ", "β⁻ Ra-228 и Ac-228 с γ-каскадом в одном импульсе → континуум от нуля, вклад в пик ~12"],
  T6: ["α, β⁻ + γ", "α Th-228, Ra-224, Rn-220, Po-216 → ≈ 1072,5 / 1141,5 / 1305,5 / 1443,5; β Pb-212, Bi-212, Tl-208 — континуум; Bi-212→Po-212 слипаются в один импульс (90,5 %), остальные Po-212 → ≈ 2054,5"],
  T7: ["α, β⁻", "α U-238 4198 кэВ → ≈ 769,5; β Th-234 и Pa-234m — континуум от нуля, главный вклад в пик ~12"],
  T8: ["α", "α 4774 кэВ → ≈ 907,5; ветвь на уровень 53 кэВ — гашёная α плюс разрядка уровня"],
  T9: ["α", "α 4687 кэВ → ≈ 886,5; ветвь на уровень 68 кэВ — гашёная α плюс разрядка уровня"],
  T10: ["α, β⁻ + γ", "α Ra-226, Rn-222, Po-218, Po-214 → ≈ 910,5 / 1089,5 / 1226,5 / 1711,5 (Po-214 не слипается, τ 236 мкс); β Pb-214, Bi-214 — континуум"],
  T11: ["α, β⁻", "α 4397 кэВ (57 %) на уровень 205 кэВ: гашёная α ≈ 816,5 плюс 205 кэВ разрядки; β Th-231 — мягкий континуум"],
  T12: ["α, β⁻ + γ", "α Th-227, Ra-223, Rn-219, Bi-211, Po-215 → 1100–1620 (Rn-219 ≈ 1455,5, Po-215 ≈ 1621,5); β Ac-227, Pb-211, Tl-207 — континуум"]
};

const cv = document.getElementById("cv-main");
const tip = document.getElementById("tip-main");

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

function visRate(c) { let s = 0; for (const i of visibleIdx()) s += arr(c)[i]; return s; }   // вклад компоненты в видимом диапазоне

function unitLabel() { return state.unit === "rate" ? "отсч./(с·кэВ)" : "отсч./кэВ за 668 023,5 с"; }

function drawChart() {
  const { g, w, h } = fitCanvas(cv);
  const x0 = M.l, x1 = w - M.r, y0 = M.t, y1 = h - M.b;
  const [xlo, xhi] = currentRange();
  const X = e => x0 + (e - xlo) / (xhi - xlo) * (x1 - x0);
  const idx = visibleIdx();
  if (idx.length < 2) return;

  // стопка: снизу МЕНЬШИЙ вклад в видимом диапазоне, сверху крупнейший (Lu-176 не закрывает малые компоненты в лог-шкале)
  const stackComps = COMPS.filter(c => c.inSum && state.on[c.id]).sort((a, b) => visRate(a) - visRate(b));

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

  // Подписи маркеров — ВНУТРИ графика: левая линия подписывается слева, правая — справа (не пересекаются).
  g.font = "11px monospace"; g.textBaseline = "top"; g.fillStyle = css("--dim");
  MARKERS.forEach((m, k) => {
    if (m.E < xlo || m.E > xhi) return;
    const x = X(m.E), w = g.measureText(m.short).width;
    let left = (k === 0), yy = y0 + 4;
    if (left && x - 4 - w < x0 + 2) { left = false; yy = y0 + 18; }
    if (!left && x + 4 + w > x1 - 2) { left = true; yy = y0 + 18; }
    g.textAlign = left ? "right" : "left";
    g.fillText(m.short, left ? x - 4 : x + 4, yy);
  });

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

function drawAll() {
  drawChart();
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

function hide() {
  state.cursor = null;
  tip.hidden = true;
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

  const total = c => c.rate.reduce((s, v) => s + v, 0);
  for (let c of COMPS.slice().sort((a, b) => (b.inSum - a.inSum) || (total(b) - total(a)))) {   // легенда — по вкладу во всём диапазоне
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

// образец линии графика для строк «Сумма» и «Измерение»: рисунок + название + тип линии словами
function lineSample(kind, name, type) {
  const NS = "http://www.w3.org/2000/svg", wrap = document.createElement("span");
  const svg = document.createElementNS(NS, "svg"), p = document.createElementNS(NS, "polyline");
  svg.setAttribute("width", "26"); svg.setAttribute("height", "12"); svg.setAttribute("class", "ln");
  p.setAttribute("fill", "none");
  p.setAttribute("stroke", kind === "sum" ? css("--sum-line") : css("--ink"));
  p.setAttribute("stroke-width", kind === "sum" ? "2" : "1");
  p.setAttribute("points", kind === "sum" ? "0,6 26,6" : "0,9 4,9 4,4 8,4 8,8 12,8 12,3 16,3 16,7 20,7 20,5 26,5");
  svg.appendChild(p);
  const t = document.createElement("span"); t.className = "lt"; t.textContent = ` — ${type}`;
  wrap.append(svg, document.createTextNode(name), t);
  return wrap;
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

  const byContribution = COMPS.slice().sort((a, b) => (b.inSum - a.inSum) || (visRate(b) - visRate(a)));   // таблица — по вкладу, крупнейший сверху
  for (let c of byContribution) {
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
    if (PHYS[c.id]) {   // тип распада — метка у названия, механизм — строкой ниже
      const t = document.createElement("span"); t.className = "dtype"; t.textContent = PHYS[c.id][0];
      const m = document.createElement("span"); m.className = "srcnote mech"; m.textContent = PHYS[c.id][1];
      firstCell.append(t, m);
    }
    if (SRC[c.id]) {   // основание нормировки активности — второй строкой под названием
      const src = document.createElement("span");
      src.className = "srcnote";
      src.textContent = SRC[c.id];
      firstCell.appendChild(src);
    }
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
  firstSumCell.appendChild(lineSample("sum", "Сумма (вес 1)", "красная сплошная линия"));
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
    firstMeasCell.appendChild(lineSample("meas", `Измерение, ${MEAS.meta.serial}`, "чёрная ступенчатая линия"));
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
  canvas.addEventListener("pointermove", pointermove);
  canvas.addEventListener("pointerleave", pointerleave);
  canvas.addEventListener("mousedown", mousedown);
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
