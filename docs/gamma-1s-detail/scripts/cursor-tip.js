/* Подсказка у курсора: энергия и отсчёты ближайшего канала (поля графика — как в detail.js) */
(function () {
  "use strict";
  var tip = document.createElement("div");
  tip.style.cssText = "position:fixed;z-index:20;pointer-events:none;padding:2px 7px;font:600 12px system-ui,sans-serif;" +
    "background:var(--paper,#fff);color:var(--ink,#000);border:2px solid var(--ink,#000);white-space:nowrap;display:none";
  document.body.appendChild(tip);
  function onMove(ev) {
    var r = ev.target.getBoundingClientRect(), x = ev.clientX - r.left, U = window.G1SD_UI, m = window.G1SD && window.G1SD.meas;
    if (!U || !m || x < 62 || x > r.width - 14) { tip.style.display = "none"; return; }
    var z = U.state.zoom, lo = z ? z.xLo : 30, hi = z ? z.xHi : 3000, k = 0, best = 1e9;
    var E = lo + (x - 62) / (r.width - 76) * (hi - lo);
    for (var i = 0; i < m.e.length; i++) if (Math.abs(m.e[i] - E) < best) { best = Math.abs(m.e[i] - E); k = i; }
    tip.textContent = E.toFixed(1).replace(".", ",") + " кэВ · канал " + k + " · измерение " + Math.round(m.counts[k]);
    tip.style.display = "block"; tip.style.left = (ev.clientX + 14) + "px"; tip.style.top = (ev.clientY - 30) + "px";
  }
  document.addEventListener("mousemove", function (ev) {
    if (ev.target.id === "cvSpec" || ev.target.id === "cvRatio") onMove(ev); else tip.style.display = "none";
  });
})();
