function fillM2Summary() {
    var el = document.getElementById(PFX + "sumM2");
    if (!el) return;
    var m2 = D.method2;
    var cs = m2.per_nuclide.CS137;
    var k = m2.per_nuclide.K40;
    var html = "";
    html += cell("активность Cs-137 (± с поправкой Бирге)", cnt(cs.A_Bq) + " Бк <em>± " + cnt(cs.dA_Bq) + " Бк</em> <em>· " + cnt(cs.per_kg) + " Бк/кг</em>", true);
    html += cell("активность K-40, предварительно (± с поправкой Бирге)", cnt(k.A_Bq) + " Бк <em>± " + cnt(k.dA_Bq) + " Бк</em> <em>· " + cnt(k.per_kg) + " Бк/кг</em>", true);
    html += cell("к Бета-1С (СпектраЛайн; его значение пересчитано на дату нашего измерения 09.10.2026)", "Cs-137 " + num(cs.ref_ratio, 3) + " (" + signedPct(cs.ref_ratio) + "); K-40 " + num(k.ref_ratio, 3) + " (" + signedPct(k.ref_ratio) + ")");
    html += cell("χ²/ν", num(m2.chi2_ndof, 2));
    html += cell("число степеней свободы ν", cnt(m2.ndof));
    html += cell("окно подгонки", num(m2.E_fit_lo, 0) + "–" + num(m2.E_fit_hi, 0) + " кэВ");
    html += cell("линий в библиотеке", cnt(m2.n_lines) + " (в таблице " + cnt(m2.n_lines_shown) + ": ниже 25 кэВ не показаны " + cnt(m2.n_below_25) + "; рентгеновских линий в библиотеке " + cnt(m2.n_xray_energies) + ", в таблице " + cnt(m2.n_xray_shown) + ")");
    html += cell("узлов сетки откликов", cnt(m2.n_nodes));
    el.innerHTML = html;
}

function buildM2Nuc() {
    var tbl = document.getElementById(PFX + "tblM2N");
    if (!tbl) return;
    var html = "<thead><tr><th>нуклид</th><th class='num'>амплитуда, Бк</th><th class='num'>удельная, Бк/кг</th><th class='num'>к опорному (Бета-1С)</th><th class='num'>доля в спектре</th><th class='num'>из неё тормозное излучение β/e⁻</th><th>пояснение</th></tr></thead><tbody>";
    var keys = ["CS137", "K40", "SRY90"];
    for (var i = 0; i < keys.length; i++) {
        var key = keys[i];
        var nuclide = null;
        for (var j = 0; j < D.nuclides.length; j++) {
            if (D.nuclides[j].key === key) {
                nuclide = D.nuclides[j];
                break;
            }
        }
        var e = D.method2.per_nuclide[key];
        if (!nuclide || !e) continue;
        var ampHtml;
        if (e.nuisance) {
            ampHtml = "вспомогательный параметр подгонки, активность не определяется";
        } else if (e.dA_Bq === null) {
            ampHtml = cnt(e.A_Bq) + " Бк";
        } else {
            ampHtml = cnt(e.A_Bq) + " ± " + cnt(e.dA_Bq) + " Бк";
        }
        var specHtml = e.per_kg == null ? "—" : cnt(e.per_kg);
        var refHtml = e.ref_ratio == null ? "—" : num(e.ref_ratio, 3) + " (" + signedPct(e.ref_ratio) + ")";
        var shareHtml = num(100 * e.share, 1) + " %";
        var betaHtml = e.beta_frac == null ? "—" : num(100 * e.beta_frac, 1) + " %";
        var noteHtml = esc(e.note_ru || "");
        html += "<tr><td><span class='sw' style='background:" + nuclide.color + "'></span>" + esc(nuclide.label_ru) + "</td>";
        html += "<td class='num'>" + ampHtml + "</td>";
        html += "<td class='num'>" + specHtml + "</td>";
        html += "<td class='num'>" + refHtml + "</td>";
        html += "<td class='num'>" + shareHtml + "</td>";
        html += "<td class='num'>" + betaHtml + "</td>";
        html += "<td>" + noteHtml + "</td></tr>";
    }
    html += "</tbody>";
    tbl.innerHTML = html;
}
