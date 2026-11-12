/* Увеличение рисунка по клику: затемнённый фон, рисунок по высоте окна, закрытие кликом или Esc */
(function () {
  "use strict";
  var box = null;
  function close() { if (box) { box.remove(); box = null; } }
  document.addEventListener("click", function (ev) {
    var img = ev.target.closest ? ev.target.closest(".figpair img") : null;
    if (box) { close(); return; }
    if (!img) return;
    box = document.createElement("div");
    box.style.cssText = "position:fixed;inset:0;background:rgba(0,0,0,.8);display:flex;align-items:center;justify-content:center;z-index:50;cursor:zoom-out";
    var big = document.createElement("img"); big.src = img.src; big.alt = img.alt;
    big.style.cssText = "max-height:96vh;max-width:96vw;background:#fff";
    box.appendChild(big); document.body.appendChild(box);
  });
  document.addEventListener("keydown", function (ev) { if (ev.key === "Escape") close(); });
})();
