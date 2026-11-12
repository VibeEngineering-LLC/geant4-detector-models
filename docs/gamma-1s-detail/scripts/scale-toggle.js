/* Лин/лог вертикальной оси спектра: состояние G1SD_UI.state.lin, перерисовка G1SD_UI.redraw */
(function () {
  "use strict";
  function addSeg() {
    var rs = document.getElementById("rangeSeg");
    if (!rs || document.getElementById("scaleSeg")) return;
    var ys = document.createElement("div"); ys.className = "seg"; ys.id = "scaleSeg";
    ys.innerHTML = '<button type="button" class="btn" data-scale="log" aria-pressed="true">шкала: логарифмическая</button>' +
      '<button type="button" class="btn" data-scale="lin" aria-pressed="false">шкала: линейная</button>';
    rs.parentNode.insertBefore(ys, rs.nextSibling);
  }
  document.addEventListener("click", function (ev) {
    var b = ev.target.closest ? ev.target.closest("[data-scale]") : null;
    if (!b || !window.G1SD_UI) return;
    window.G1SD_UI.state.lin = b.getAttribute("data-scale") === "lin";
    var all = b.parentNode.querySelectorAll("[data-scale]");
    for (var i = 0; i < all.length; i++) all[i].setAttribute("aria-pressed", all[i] === b ? "true" : "false");
    window.G1SD_UI.redraw();
  });
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", addSeg); else addSeg();
})();
