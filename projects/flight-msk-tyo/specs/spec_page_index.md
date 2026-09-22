Write ONE complete file `index.html` (HTML + inline CSS + inline JavaScript, no external libraries, no network). Output ONLY the file content (no markdown fences, no prose). Text of the page in Russian. It loads `data.js` with `<script src="data.js"></script>` BEFORE the inline script; `data.js` defines `window.SPECTRA = {E0, dE, unit, series:[{id,label,group,y:[...],sig?:[...]}], lines:[{E,reaction,net,sig,signif,origin}], unlisted:[{E,net,sig,signif,origin}]}`; the energy of point i is `E0 + i*dE` keV; `group` is `total`, `component` or `origin`; `y` is counts per hour per keV (may contain zeros).

## Page content
* `<title>Спектр на эшелоне</title>`, viewport meta. Header: h1 `Спектр на эшелоне: рейс DAD–SGN`; paragraph `A320/A321, 20.09.2026, высота до 9800 м, столик у окна, эконом. Расчёт Geant4 11.4.2 + PARMA/EXPACS; прибор АтомНано 16 (CsI, 16,2 см³), разрешение 41,6 кэВ на 662 кэВ. Единицы: отсчёты в час на кэВ.`
* Canvas `#cv` (width 100 %, height 420 px, `devicePixelRatio` aware), toolbar buttons: presets `20–100`, `100–511`, `511–1500`, `1,5–3 МэВ`, `3–10 МэВ`, `весь` (data attributes `data-lo`, `data-hi`; `весь` = 20..10000) and a toggle `лог / лин` (log Y by default). The preset that is active is highlighted (class `active`); with a free zoom NO preset is highlighted.
* One checkbox PER SERIES (from `series`, in the order of the data), grouped under three captions `Итог`, `Компоненты`, `Происхождение`; default checked: the `total` series only. Each series has its own colour taken from this fixed palette by its index in `series`: `['#222222','#d62728','#1f77b4','#2ca02c','#ff7f0e','#9467bd','#8c564b','#e377c2','#17becf','#bcbd22','#7f7f7f','#393b79','#e6550d','#31a354','#756bb1','#636363','#b15928','#fb9a99','#6a3d9a','#33a02c','#a6cee3']`; the colour swatch is shown next to the label. `total` is drawn 2 px thick, others 1.2 px.
* Readout line under the chart (`#readout`): on pointer move over the plot area `E = <энергия> кэВ`, then for every ticked series `<label>: <значение> ± <σ>` (σ only for `total`), values with 3 significant digits.
* Two tables: `Линии, найденные в спектре` (columns E, кэВ · Реакция · Чистая площадь ± σ, отсч./ч · Значимость · Происхождение) and `Пики вне перечня` (E · площадь ± σ · значимость · происхождение). Clicking a row zooms to `E ± 8 %`. Build the tables ONCE at load (not inside `draw`). Escape text with a helper `esc()` before `innerHTML`.
* `<details><summary>Как читать</summary>` 3 short sentences: colours per series, drag to select and double click to reset, meaning of `сорт / процесс / материал`.
* Colour tokens as CSS variables on `:root`, redefined for dark mode with `@media (prefers-color-scheme: dark)`; body background explicit; no horizontal scroll at 375 px width; 16 px side gutter; table container `overflow-x: auto`.

## Colours on the canvas (IMPORTANT)
Canvas cannot use `var(--x)`. Read the theme colours once per draw: `const css = getComputedStyle(document.documentElement); const C = (n) => css.getPropertyValue(n).trim();` and use `C('--grid')`, `C('--axis')`, `C('--text')`.

## Drawing (reference skeleton — follow it)
```js
function ticksLin(lo, hi, n) {                    // «красивые» деления
  const raw = (hi - lo) / n, p = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / p;
  const step = (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p, out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * step; v += step) out.push(v);
  return out;
}
function draw() {
  const w = cv.clientWidth, h = 420, dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);               // не накапливать scale
  const [xlo, xhi] = currentRange(), pw = w - M.l - M.r, ph = h - M.t - M.b;
  const X = (E) => M.l + (E - xlo) / (xhi - xlo) * pw;
  const vis = series.filter((s) => shown[s.id]);
  // диапазон Y — только по ВИДИМЫМ рядам и ВИДИМОМУ диапазону энергий
  let ymax = 0, ymin = Infinity;
  for (const s of vis) for (let j = idx(xlo); j <= idx(xhi); j++) { const v = s.y[j]; if (v > 0) { ymax = Math.max(ymax, v); ymin = Math.min(ymin, v); } }
```
Continue in the same spirit: define `idx(E) = clamp(Math.round((E - E0)/dE), 0, n-1)`; if `ymax <= 0` draw only the axes and return. Log mode: `lo = max(ymin, ymax*1e-6)`, `Y(v) = M.t + ph*(1 - (log10(max(v, lo)) - log10(lo))/(log10(ymax*1.3) - log10(lo)))` with decade tick marks (`10^k`) and labels; linear mode: `Y(v) = M.t + ph*(1 - v/(ymax*1.05))` with `ticksLin`. X ticks with `ticksLin(xlo, xhi, 8)`, grid lines in `C('--grid')`, labels with `C('--text')`, axis titles `E, кэВ` and `отсч./ч/кэВ`. Series polylines: for every pixel column keep the MAXIMUM of the y values falling into it when more than one point falls into the same column (so peaks survive zoom-out), else the plain polyline. Draw the selection rectangle (semi-transparent) while dragging.

## Navigation (MANDATORY; reference behaviour)
State: `zoom = null` or `{lo, hi}`; `drag = {active:false, x0:0, x1:0}`; margins `M = {l: 64, r: 16, t: 12, b: 34}`; `currentRange()` = zoom if set, else the range of the active preset.
```js
cv.addEventListener("mousedown", (ev) => {
  const r = cv.getBoundingClientRect(), x = ev.clientX - r.left;
  if (x < M.l || x > r.width - M.r) return;
  ev.preventDefault(); drag = { active: true, x0: x, x1: x }; draw();
});
cv.addEventListener("pointermove", (ev) => {
  const r = cv.getBoundingClientRect(), x = ev.clientX - r.left;
  if (drag.active) { drag.x1 = Math.max(M.l, Math.min(r.width - M.r, x)); draw(); return; }
  updateReadout(x);
});
cv.addEventListener("dblclick", () => { zoom = null; setPreset(0, 10000); });
document.addEventListener("mouseup", () => {          // один слушатель на document
  if (!drag.active) return; drag.active = false;
  if (Math.abs(drag.x1 - drag.x0) < 6) { draw(); return; }       // короткий клик — не зум
  const r = cv.getBoundingClientRect(), [xlo, xhi] = currentRange();
  const toE = (px) => xlo + (px - M.l) / (r.width - M.r - M.l) * (xhi - xlo);
  const lo = Math.max(20, toE(Math.min(drag.x0, drag.x1))), hi = Math.min(10000, toE(Math.max(drag.x0, drag.x1)));
  if (hi - lo < 1) { draw(); return; }
  zoom = { lo, hi }; markPreset(null); draw();
});
```
`setPreset(lo, hi)` (used by preset buttons: they read `data-lo`/`data-hi`) clears the zoom, stores the preset range, highlights the matching button and redraws (`весь` uses lo 20, hi 10000 — make the dblclick call consistent with this). Redraw on `ResizeObserver` of the canvas's parent only when the width changed. If `!window.SPECTRA` show a message instead of the chart. No `localStorage`, no network.
