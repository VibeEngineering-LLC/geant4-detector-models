function cmpCol(kind) {
  if (kind === "ref") return pal().ink;
  if (kind === "m1") return "#0f5aa8";
  if (kind === "m2") return "#c8541c";
  if (kind === "area") return "#c8541c";
  if (kind === "ster") return "#7a3b12";
  return pal().ink;
}

function cmpItems() {
  if (!D.cmp || !D.cmp.cs) return [];
  return D.cmp.cs.map(function (it) {
    return { lab: it.lab, A: it.A, dA: it.dA == null ? 0 : it.dA, col: cmpCol(it.kind) };
  });
}

function cmpItemsK() {
  if (!D.cmp || !D.cmp.k) return [];
  return D.cmp.k.map(function (it) {
    return { lab: it.lab, A: it.A, dA: it.dA == null ? 0 : it.dA, col: cmpCol(it.kind) };
  });
}

function fillSummaries() {
  var el = document.getElementById(PFX + "sumM1");
  if (!el) return;
  var m1 = D.method1;
  var pn = m1.per_nuclide;
  var cs = pn.CS137;
  var k = pn.K40;
  var html = "";
  html += cell("активность Cs-137 (± с поправкой Бирге)", cnt(cs.A_Bq) + " Бк <em>± " + cnt(cs.dA_Bq) + " Бк</em> <em>· " + cnt(cs.per_kg) + " Бк/кг</em>", true);
  html += cell("активность K-40, предварительно (± с поправкой Бирге)", cnt(k.A_Bq) + " Бк <em>± " + cnt(k.dA_Bq) + " Бк</em> <em>· " + cnt(k.per_kg) + " Бк/кг</em>", true);
  html += cell("к Бета-1С (СпектраЛайн; его значение пересчитано на дату нашего измерения 09.10.2026)", "Cs-137 " + num(cs.ref_ratio, 3) + " (" + signedPct(cs.ref_ratio) + "); K-40 " + num(k.ref_ratio, 3) + " (" + signedPct(k.ref_ratio) + ")");
  html += cell("χ²/ν", num(m1.chi2_ndof, 2));
  html += cell("число степеней свободы ν", cnt(m1.ndof));
  html += cell("окно подгонки", num(m1.E_fit_lo, 0) + "–" + num(m1.E_fit_hi, 0) + " кэВ");
  el.innerHTML = html;
  fillM2Summary();   // #GS-74: метод 2 посчитан (berry_m2_block.js)
  buildM2Nuc();
}

function buildM1() {
  var tbl = document.getElementById(PFX + "tblM1");
  if (!tbl) return;
  var pn = D.method1.per_nuclide;
  var html = "<thead><tr><th>нуклид</th><th class='num'>амплитуда, Бк</th><th class='num'>к опорному (Бета-1С)</th><th class='num'>доля в спектре</th><th>пояснение</th></tr></thead><tbody>";
  D.nuclides.forEach(function (nu) {
    var e = pn[nu.key];
    if (!e) return;
    var amp;
    if (e.nuisance) {
      amp = "вспомогательный параметр подгонки, активность не определяется";
    } else if (e.A_Bq === 0) {
      amp = "0 Бк (амплитуда подгонки равна нулю)";
    } else if (e.dA_Bq == null) {
      amp = cnt(e.A_Bq) + " Бк";
    } else {
      amp = cnt(e.A_Bq) + " ± " + cnt(e.dA_Bq) + " Бк";
    }
    var ref = e.ref_ratio != null ? num(e.ref_ratio, 3) : "—";
    var share = num(100 * e.share, 1) + " %";
    var note = esc(nu.note || "");
    html += "<tr><td><span class='sw' style='background:" + nu.color + "'></span>" + esc(nu.label_ru) + "</td>";
    html += "<td class='num'>" + amp + "</td>";
    html += "<td class='num'>" + ref + "</td>";
    html += "<td class='num'>" + share + "</td>";
    html += "<td>" + note + "</td></tr>";
  });
  html += "</tbody>";
  tbl.innerHTML = html;
}

function fillCmpTable() {
  var el = document.getElementById(PFX + "cmpTable");
  if (!el) return;
  var massKg = D.cmp.mass_kg;
  var html = "";
  function group(items, label) {
    var refA = null;
    for (var i = 0; i < items.length; i++) {
      if (items[i].kind === "ref") { refA = items[i].A; break; }
    }
    var h = "<div class='cmp-row'><span class='cmp-lab'>" + label + "</span><span class='cmp-val'></span><span class='cmp-note'></span></div>";
    for (var j = 0; j < items.length; j++) {
      var it = items[j];
      var cls = it.kind === "ref" ? "cmp-pass" : (it.kind === "m1" ? "cmp-m1" : "cmp-m2");
      var val = cnt(it.A) + " Бк";
      if (it.dA != null && it.dA !== 0) val += " <em>± " + cnt(it.dA) + " Бк</em>";
      val += " <em>· " + cnt(it.A / massKg) + " Бк/кг</em>";
      var note = esc(it.note || "");
      if ((it.kind === "m1" || it.kind === "m2") && refA != null) {
        note += "; к Бета-1С " + num(it.A / refA, 3) + " (" + signedPct(it.A / refA) + ")";
      }
      h += "<div class='cmp-row " + cls + "'><span class='cmp-lab'>" + esc(it.lab) + "</span><span class='cmp-val big-num'>" + val + "</span><span class='cmp-note'>" + note + "</span></div>";
    }
    return h;
  }
  html += group(D.cmp.cs, "Cs-137");
  html += group(D.cmp.k, "K-40");
  el.innerHTML = html;
}

function drawCmpOne(cvId, items) {
  var cv = document.getElementById(PFX + cvId);
  if (!cv) return;
  var p = pal();
  var f = fit(cv);
  var g = f.g, W = f.w, H = f.h;
  var m = { l: 26, r: 20, t: 18, b: 34 };
  var lo = Infinity, hi = -Infinity;
  items.forEach(function (it) {
    lo = Math.min(lo, it.A - it.dA);
    hi = Math.max(hi, it.A + it.dA);
  });
  var pad = (hi - lo) * 0.15 + 1;
  lo -= pad; hi += pad; lo = Math.max(0, lo);

  g.strokeStyle = p.rule; g.lineWidth = 2;
  g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);

  var range = hi - lo;
  var stp = Math.pow(10, Math.floor(Math.log10(range / 4)));
  var s = Math.max(stp, Math.ceil(range / 5 / stp) * stp);
  g.strokeStyle = p.grid; g.beginPath();
  g.fillStyle = p.faint; g.font = "11px system-ui, sans-serif";
  g.textAlign = "center"; g.textBaseline = "top";
  for (var v = Math.ceil(lo / s) * s; v <= hi; v += s) {
    var x = mapX(v, lo, hi, m.l, W - m.r);
    g.moveTo(x, m.t); g.lineTo(x, H - m.b);
    g.fillText(cnt(v), x, H - m.b + 4);
  }
  g.stroke();
  g.textAlign = "center"; g.textBaseline = "bottom";
  g.fillText("активность, Бк", (m.l + W - m.r) / 2, H - 2);

  var innerH = H - m.b - m.t;
  var rowH = innerH / items.length;
  for (var j = 0; j < items.length; j++) {
    var it = items[j];
    var yc = m.t + rowH * (j + 0.5);
    var xl = mapX(it.A - it.dA, lo, hi, m.l, W - m.r);
    var xr = mapX(it.A + it.dA, lo, hi, m.l, W - m.r);
    var xm = mapX(it.A, lo, hi, m.l, W - m.r);
    g.fillStyle = it.col; g.globalAlpha = 0.28;
    g.fillRect(xl, yc - rowH * 0.28, Math.max(xr - xl, 1), rowH * 0.56);
    g.globalAlpha = 1;
    g.strokeStyle = it.col; g.lineWidth = 3;
    g.beginPath();
    g.moveTo(xm, yc - rowH * 0.36); g.lineTo(xm, yc + rowH * 0.36);
    g.stroke();
    g.fillStyle = it.col;
    g.textAlign = "left"; g.textBaseline = "middle";
    g.font = "bold 13px system-ui, sans-serif";
    g.fillText(it.lab, m.l + 6, yc - rowH * 0.28);
    g.fillStyle = p.ink;
    g.font = "12px ui-monospace, Menlo, monospace";
    g.fillText(cnt(it.A) + " ± " + cnt(it.dA) + " Бк",
               Math.min(xr + 8, W - m.r - 130), yc);
  }
}

function drawCmp() {
  drawCmpOne("cvCmp", cmpItems());
  drawCmpOne("cvCmpK", cmpItemsK());
}
