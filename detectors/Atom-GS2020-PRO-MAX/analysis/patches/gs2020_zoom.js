function zEfromX(x, rectWidth) {
  var e = D.spectrum.e_of_ch;
  var xLo = ST.zoom ? ST.zoom.xLo : 0;
  var xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];
  var m = { l: 62, r: 14 };
  return xLo + ((x - m.l) / (rectWidth - m.r - m.l)) * (xHi - xLo);
}

function wireZoom(cvId) {
  var cv = document.getElementById(cvId);
  if (!cv) return;
  var dragging = false;
  var dragState = null;
  var m = { l: 62, r: 14 };

  cv.addEventListener("mousedown", function (ev) {
    var r = cv.getBoundingClientRect();
    var x = ev.clientX - r.left;
    if (x < m.l || x > r.width - m.r) return;
    ev.preventDefault();
    dragging = true;
    dragState = { x0: x, x1: x, w: r.width };
  });

  document.addEventListener("mousemove", function (ev) {
    if (!dragging) return;
    var r = cv.getBoundingClientRect();
    var x = ev.clientX - r.left;
    if (x < m.l) x = m.l;
    if (x > r.width - m.r) x = r.width - m.r;
    dragState.x1 = x;
  });

  document.addEventListener("mouseup", function () {
    if (!dragging) return;
    dragging = false;
    var x0 = dragState.x0, x1 = dragState.x1, w = dragState.w;
    dragState = null;
    if (Math.abs(x1 - x0) < 6) return;
    ST.zoom = { xLo: Math.max(0, zEfromX(Math.min(x0, x1), w)),
                xHi: zEfromX(Math.max(x0, x1), w) };
    redraw();
  });

  cv.addEventListener("dblclick", function () {
    ST.zoom = null;
    redraw();
  });
}
