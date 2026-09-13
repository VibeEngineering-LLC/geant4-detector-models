(function () {
  "use strict";
  if (!window.AMTICSEU || !window.G1SA) {
    console.error("Не загружены данные AMTICSEU или G1SA");
    return;
  }

  var D = window.AMTICSEU;
  // Помощники основного файла (генерация обращалась к ним напрямую, без
  // привязки — ReferenceError на каждом числе; поймано при отрисовке 11.09).
  var G = window.G1SA, num = G.num, cnt = G.cnt, esc = G.esc,
      fit = G.fit, pal = G.pal, mapX = G.mapX;

  function labelRu(key) {
    var n = D.nuclides;
    for (var i = 0; i < n.length; i++) {
      if (n[i].key === key) return n[i].label_ru;
    }
    return key;
  }

  function colorOf(key) {
    var n = D.nuclides;
    for (var i = 0; i < n.length; i++) {
      if (n[i].key === key) return n[i].color;
    }
    return "#000";
  }

  function sw(color) {
    return "<span class='sw' style='background:" + color + "'></span>";
  }

  function mainNet(key) {
    var lines = D.netarea.lines;
    for (var i = 0; i < lines.length; i++) {
      if (lines[i].nuclide === key && lines[i].main) return lines[i];
    }
    return null;
  }

  function interp(arr, x) {
    if (!arr || arr.length < 2) return 0;
    for (var i = 0; i < arr.length - 1; i++) {
      var a = arr[i], b = arr[i + 1];
      if (a[0] <= x && x <= b[0]) {
        var t = (x - a[0]) / (b[0] - a[0]);
        return a[1] + t * (b[1] - a[1]);
      }
    }
    return arr[0][1];
  }

  function fillValues() {
    var els = document.querySelectorAll("[data-v]");
    for (var i = 0; i < els.length; i++) {
      var el = els[i], key = el.getAttribute("data-v");
      var val = "?";
      try {
        switch (key) {
          case "n_refs":
            val = D.meta.cal_sample.refs ? cnt(D.meta.cal_sample.refs.length) : "?";
            break;
          case "n_fwhm_points":
            val = D.fwhm_cal.points ? cnt(D.fwhm_cal.points.length) : "?";
            break;
          case "m1_chi2":
            val = num(D.method1.chi2_ndof, 1);
            break;
          // Единая метрика (дисперсия A1): ею сравнимы между собой разные критерии, поэтому
          // в тексте страницы согласие модели названо именно ею, а не χ² рабочего критерия.
          case "m1_chi2_ref":
            val = num(D.method1.chi2_ref_ndof, 1);
            break;
          case "m1_e1_tv":
            val = D.method1_e1 ? num(D.method1_e1.tv, 4) : "?";
            break;
          case "m2_chi2":
            // D-020 п.3: единая метрика (дисперсия измерения), сопоставимая с путём 1
            val = num(D.method2.chi2_ref_ndof, 1);
            break;
          case "m2_zone":
            val = num(D.method2.zone_50_70_ratio, 3);
            break;
          case "m2_nodes":
            val = cnt(D.method2.n_nodes);
            break;
          case "zone_ratio":
            val = D.method1.zone_50_70 ? num(D.method1.zone_50_70.ratio, 3) : "?";
            break;
          case "decays":
            var ds = D.meta.model.decays;
            if (ds && ds.length > 0) {
              var allSame = true;
              for (var j = 1; j < ds.length; j++) {
                if (ds[j] !== ds[0]) {
                  allSame = false;
                  break;
                }
              }
              val = allSame ? cnt(ds[0]) : ds.map(cnt).join(" / ");
            } else val = "?";
            break;
          case "tail_T":
            val = D.meta.tail_T ? num(D.meta.tail_T, 2) : "?";
            break;
          case "net_win":
            val = D.netarea.definition.win_fwhm ? num(D.netarea.definition.win_fwhm, 2) : "?";
            break;
          case "net_side":
            var side = D.netarea.definition.side_fwhm;
            if (side && side.length >= 2) {
              val = num(side[0], 1) + "…" + num(side[1], 1);
            } else val = "?";
            break;
          case "blend_k":
            val = D.netarea.definition.blend_fwhm ? num(D.netarea.definition.blend_fwhm, 1) : "?";
            break;
          case "am_fw_marinelli":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Am241" && g[j].record === "маринелли") {
                val = num(g[j].fwhm_ratio, 3);
                break;
              }
            }
            break;
          case "am_fw_point":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Am241" && g[j].record === "точечный 5 см") {
                val = num(g[j].fwhm_ratio, 3);
                break;
              }
            }
            break;
          case "am_fw_petri":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Am241" && g[j].record === "Петри-60") {
                val = num(g[j].fwhm_ratio, 3);
                break;
              }
            }
            break;
          case "am_point_ratio":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Am241" && g[j].record === "точечный 5 см") {
                val = num(g[j].ratio, 3);
                break;
              }
            }
            break;
          case "am_petri_ratio":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Am241" && g[j].record === "Петри-60") {
                val = num(g[j].ratio, 3);
                break;
              }
            }
            break;
          case "am_fw_curve_marinelli":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Am241" && g[j].record === "маринелли") {
                val = num(g[j].fwhm_curve_keV, 2);
                break;
              }
            }
            break;
          case "am_marinelli_ratio":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Am241" && g[j].record === "маринелли") {
                val = num(g[j].ratio, 3);
                break;
              }
            }
            break;
          case "cs_point_ratio":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Cs137chain" && g[j].record === "точечный 5 см") {
                val = num(g[j].ratio, 3);
                break;
              }
            }
            break;
          case "cs_petri_ratio":
            var g = D.cross_geometry;
            for (var j = 0; j < g.length; j++) {
              if (g[j].nuclide === "Cs137chain" && g[j].record === "Петри-60") {
                val = num(g[j].ratio, 3);
                break;
              }
            }
            break;
          case "geo_spread":
            var maxRatio = 0;
            for (var j = 0; j < D.cross_geometry.length; j++) {
              var g = D.cross_geometry[j];
              if (g.nuclide === "Am241" && g.record === "маринелли") continue;
              if (g.ratio) {
                var diff = Math.abs(g.ratio - 1);
                if (diff > maxRatio) maxRatio = diff;
              }
            }
            val = num(100 * maxRatio, 0);
            break;
          case "bg_rms":
            val = D.meta.cal_bg.rms_keV ? num(D.meta.cal_bg.rms_keV, 2) : "?";
            break;
          case "n_bg_nat":
            val = D.meta.cal_bg.anchors ? cnt(D.meta.cal_bg.anchors.length - 1) : "?";
            break;
          case "am_shift_fwhm":
            var refs = D.meta.cal_sample.refs;
            for (var j = 0; j < refs.length; j++) {
              if (Math.abs(refs[j].E_true - 59.541) < 0.01) {
                val = num(Math.abs(refs[j].shift_fwhm), 2);
                break;
              }
            }
            break;
          case "sqrt59":
            val = D.fwhm_cal.sqrt_law ? num(interp(D.fwhm_cal.sqrt_law, 59.541), 2) : "?";
            break;
          case "meas59":
            var pts = D.fwhm_cal.points;
            for (var j = 0; j < pts.length; j++) {
              if (Math.abs(pts[j].E_keV - 59.541) < 0.01) {
                val = num(pts[j].fwhm_keV, 2);
                break;
              }
            }
            break;
        }
      } catch (e) {
        console.warn("Ошибка при заполнении значения", key, e);
      }
      el.textContent = val;
    }
  }

  function renderTwoPaths() {
    var table = document.getElementById("tblTwoPaths");
    if (!table) return;
    table.innerHTML = "";
    var thead = document.createElement("thead");
    var headRow = document.createElement("tr");
    // Две оценки пути 1 публикуются рядом (D-020): A2 — рабочий критерий, E1 — вторая мера.
    // Расхождение между ними — факт о модели, а не повод показать одно удобное число.
    ["нуклид", "паспорт, Бк", "путь 1, критерий A2", "путь 1, вторая мера E1",
     "путь 2: к паспорту (по линии)", "расхождение A2 и пути 2, %", "примечание"].forEach(function (h) {
      var th = document.createElement("th");
      th.textContent = h;
      headRow.appendChild(th);
    });
    thead.appendChild(headRow);
    table.appendChild(thead);
    var tbody = document.createElement("tbody");

    for (var i = 0; i < D.nuclides.length; i++) {
      var key = D.nuclides[i].key;
      var main = mainNet(key);
      if (!main) continue;

      var row = document.createElement("tr");
      var cell1 = document.createElement("td");
      cell1.innerHTML = sw(colorOf(key)) + " " + esc(labelRu(key));
      row.appendChild(cell1);

      var p = D.passport[key];
      if (p) {
        var cell2 = document.createElement("td");
        cell2.textContent = cnt(p.A_Bq) + " ± " + cnt(p.dA_Bq);
        row.appendChild(cell2);
      } else {
        var cell2 = document.createElement("td");
        cell2.textContent = "?";
        row.appendChild(cell2);
      }

      var m1 = D.method1.groups[key];
      if (m1 && m1.A_over_passport !== null) {
        var cell3 = document.createElement("td");
        cell3.textContent = num(m1.A_over_passport, 3);
        row.appendChild(cell3);
      } else {
        var cell3 = document.createElement("td");
        cell3.textContent = "?";
        row.appendChild(cell3);
      }

      // Вторая оценка того же пути 1 — критерий E1 (минимум невязки формы).
      var e1g = D.method1_e1 && D.method1_e1.groups ? D.method1_e1.groups[key] : null;
      var cellE1 = document.createElement("td");
      cellE1.textContent = (e1g && e1g.A_over_passport !== null && e1g.A_over_passport !== undefined)
        ? num(e1g.A_over_passport, 3) : "?";
      row.appendChild(cellE1);

      // Путь 2 — нетто-площадь главной линии (D.netarea), а не путь 1 с фоном
      // канал в канал: прежняя привязка давала расхождение 0,0 при 0,738 против 1,475.
      var p2 = (main.ratio !== null && main.ratio !== undefined) ? main.ratio : null;
      var cell4 = document.createElement("td");
      cell4.textContent = p2 === null ? "?" : num(p2, 3) + " ± " + num(main.d_ratio_stat, 3) +
        " (по " + num(main.E_keV, 1) + ")";
      row.appendChild(cell4);

      if (m1 && m1.A_over_passport !== null && p2 !== null) {
        var diff = Math.abs(p2 - m1.A_over_passport) / m1.A_over_passport * 100;
        var cell5 = document.createElement("td");
        cell5.textContent = num(diff, 1);
        row.appendChild(cell5);

        var cell6 = document.createElement("td");
        if (main.blends && main.blends.length > 0) {
          var list = [];
          for (var j = 0; j < main.blends.length; j++) {
            var b = main.blends[j];
            if (b.nuclide) {
              list.push(num(b.E_keV, 1) + " (" + labelRu(b.nuclide) + ")");
            } else {
              list.push(num(b.E_keV, 1) + " (сумма)");
            }
          }
          cell6.textContent = "бленд: " + list.join(", ");
        } else if (diff <= 5) {
          cell6.textContent = "сходятся";
        } else {
          cell6.textContent = "-";
        }
        row.appendChild(cell6);
      } else {
        var cell5 = document.createElement("td");
        cell5.textContent = "?";
        row.appendChild(cell5);

        var cell6 = document.createElement("td");
        cell6.textContent = "?";
        row.appendChild(cell6);
      }

      tbody.appendChild(row);
    }

    table.appendChild(tbody);
  }

  // Метод 2 (γ-линии × прямой отклик) рядом с путями 1 и 2 — только числа выгрузки.
  function renderMethod2() {
    var table = document.getElementById("tblMethod2");
    if (!table || !D.method2) return;
    // D-020 п.3 (13.09.2026): метод 2 тем же критерием A2; вторая мера E1 — рядом, как у пути 1.
    var rows = ["<thead><tr><th>нуклид</th><th>метод 2: активность, Бк</th>" +
      "<th>метод 2: к паспорту</th><th>метод 2, вторая мера E1: к паспорту</th>" +
      "<th>путь 1: к паспорту</th><th>путь 2: к паспорту</th></tr></thead><tbody>"];
    for (var i = 0; i < D.nuclides.length; i++) {
      var key = D.nuclides[i].key, g2 = D.method2.groups[key], g1 = D.method1.groups[key];
      if (!g2) continue;
      var g2e = D.method2_e1 ? D.method2_e1.groups[key] : null;
      var main = mainNet(key);
      rows.push("<tr><td>" + sw(colorOf(key)) + " " + esc(labelRu(key)) + "</td><td>" +
        cnt(g2.A_Bq) + " ± " + cnt(g2.dA_Bq) + "</td><td>" + num(g2.A_over_passport, 3) +
        "</td><td>" + (g2e ? num(g2e.A_over_passport, 3) : "?") +
        "</td><td>" + (g1 ? num(g1.A_over_passport, 3) : "?") + "</td><td>" +
        (main && main.ratio != null ? num(main.ratio, 3) : "?") + "</td></tr>");
    }
    table.innerHTML = rows.join("") + "</tbody>";
  }

  function renderBands() {
    var table = document.getElementById("tblBands");
    if (!table) return;
    table.innerHTML = "";
    var thead = document.createElement("thead");
    var headRow = document.createElement("tr");
    var th1 = document.createElement("th");
    th1.textContent = "вклады полос в χ²";
    th1.colSpan = 3;
    headRow.appendChild(th1);
    thead.appendChild(headRow);
    table.appendChild(thead);

    var tbody = document.createElement("tbody");

    if (D.method1.bands) {
      for (var i = 0; i < D.method1.bands.length; i++) {
        var band = D.method1.bands[i];
        var row = document.createElement("tr");
        var cell1 = document.createElement("td");
        cell1.textContent = num(band.lo, 1) + "–" + num(band.hi, 1) + " кэВ";
        row.appendChild(cell1);

        var cell2 = document.createElement("td");
        cell2.textContent = num(band.chi2, 0);
        row.appendChild(cell2);

        var cell3 = document.createElement("td");
        cell3.textContent = num(band.share_pct, 1) + "%";
        row.appendChild(cell3);

        tbody.appendChild(row);
      }
    }

    if (D.method1.zone_50_70) {
      var row = document.createElement("tr");
      var cell1 = document.createElement("td");
      cell1.textContent = "зона 50–70 кэВ: модель / измерение";
      row.appendChild(cell1);

      var cell2 = document.createElement("td");
      cell2.textContent = num(D.method1.zone_50_70.ratio, 3);
      row.appendChild(cell2);

      var cell3 = document.createElement("td");
      cell3.textContent = "максимумы: изм. " + num(D.method1.zone_50_70.peak_meas_keV, 1) +
        " / мод. " + num(D.method1.zone_50_70.peak_model_keV, 1) + " кэВ";
      row.appendChild(cell3);

      tbody.appendChild(row);
    }

    if (D.method1.chi2_ndof !== undefined) {
      // Двух χ² здесь не избежать, и подписать их обязательно: у критерия A2 в знаменателе
      // стоит ещё и дисперсия шаблонов, поэтому его χ²/ν меньше — но это шире знаменатель,
      // а не лучше согласие. Сопоставимая между критериями величина — вторая.
      var row = document.createElement("tr");
      var cell1 = document.createElement("td");
      cell1.textContent = "χ²/ν: по дисперсии критерия A2 / единая метрика; √(χ²_ref/ν)";
      row.appendChild(cell1);

      var cell2 = document.createElement("td");
      cell2.textContent = num(D.method1.chi2_ndof, 1) +
        (D.method1.chi2_ref_ndof !== undefined ? " / " + num(D.method1.chi2_ref_ndof, 1) : "");
      row.appendChild(cell2);

      var cell3 = document.createElement("td");
      cell3.textContent = num(D.method1.birge, 2) +
        (D.method1_e1 && D.method1_e1.tv !== undefined
          ? " (E1: невязка формы " + num(D.method1_e1.tv, 4) + ")" : "");
      row.appendChild(cell3);

      tbody.appendChild(row);
    }

    if (D.method1_bg_by_channel.chi2_ndof !== undefined) {
      var row = document.createElement("tr");
      var cell1 = document.createElement("td");
      cell1.textContent = "фон канал в канал (как в отчёте 11.09): зона 50–70, χ²/ν";
      row.appendChild(cell1);

      var cell2 = document.createElement("td");
      cell2.textContent = num(D.method1_bg_by_channel.zone_ratio, 3);
      row.appendChild(cell2);

      var cell3 = document.createElement("td");
      cell3.textContent = num(D.method1_bg_by_channel.chi2_ndof, 1);
      row.appendChild(cell3);

      tbody.appendChild(row);
    }

    table.appendChild(tbody);
  }

  function renderNet() {
    var sumM2 = document.getElementById("sumM2");
    if (sumM2) {
      sumM2.innerHTML = "";
      var div = document.createElement("div");
      div.innerHTML = "<span class='lab'>окно</span><span class='val'>" + num(D.netarea.definition.win_fwhm, 2) + " ПШПВ</span>";
      sumM2.appendChild(div);

      var side = D.netarea.definition.side_fwhm;
      if (side && side.length >= 2) {
        div = document.createElement("div");
        div.innerHTML = "<span class='lab'>боковые полосы</span><span class='val'>" + num(side[0], 1) + "…" + num(side[1], 1) + "</span>";
        sumM2.appendChild(div);
      }

      div = document.createElement("div");
      div.innerHTML = "<span class='lab'>сетка модели</span><span class='val'>" + num(D.netarea.definition.grid_keV, 1) + " кэВ</span>";
      sumM2.appendChild(div);

      div = document.createElement("div");
      div.innerHTML = "<span class='lab'>порог бленда</span><span class='val'>" + num(D.netarea.definition.blend_fwhm, 1) + " ПШПВ</span>";
      sumM2.appendChild(div);
    }

    var tblM2 = document.getElementById("tblM2");
    if (!tblM2) return;
    tblM2.innerHTML = "";
    var thead = document.createElement("thead");
    var headRow = document.createElement("tr");
    var th1 = document.createElement("th");
    th1.textContent = "нуклид";
    headRow.appendChild(th1);
    var th2 = document.createElement("th");
    th2.textContent = "E, кэВ";
    headRow.appendChild(th2);
    var th3 = document.createElement("th");
    th3.textContent = "пик прибора, кэВ";
    headRow.appendChild(th3);
    var th4 = document.createElement("th");
    th4.textContent = "ПШПВ прибора / по кривой, кэВ";
    headRow.appendChild(th4);
    var th5 = document.createElement("th");
    th5.textContent = "отношение ширин";
    headRow.appendChild(th5);
    var th6 = document.createElement("th");
    th6.textContent = "площадь ±";
    headRow.appendChild(th6);
    var th7 = document.createElement("th");
    th7.textContent = "эфф. измер., %";
    headRow.appendChild(th7);
    var th8 = document.createElement("th");
    th8.textContent = "эфф. модели, %";
    headRow.appendChild(th8);
    var th9 = document.createElement("th");
    th9.textContent = "изм./модель ±";
    headRow.appendChild(th9);
    var th10 = document.createElement("th");
    th10.textContent = "путь 1";
    headRow.appendChild(th10);
    var th11 = document.createElement("th");
    th11.textContent = "бленды";
    headRow.appendChild(th11);
    thead.appendChild(headRow);
    tblM2.appendChild(thead);

    var tbody = document.createElement("tbody");

    var groups = {};
    for (var i = 0; i < D.netarea.lines.length; i++) {
      var line = D.netarea.lines[i];
      if (!groups[line.nuclide]) groups[line.nuclide] = [];
      groups[line.nuclide].push(line);
    }

    for (var i = 0; i < D.nuclides.length; i++) {
      var key = D.nuclides[i].key;
      if (!groups[key]) continue;

      var lines = groups[key];
      for (var j = 0; j < lines.length; j++) {
        var line = lines[j];
        var row = document.createElement("tr");

        if (j === 0) {
          var cell1 = document.createElement("td");
          cell1.innerHTML = sw(colorOf(key)) + " " + esc(labelRu(key));
          row.appendChild(cell1);
        } else {
          var cell1 = document.createElement("td");
          cell1.textContent = "";
          row.appendChild(cell1);
        }

        var cell2 = document.createElement("td");
        if (line.main) {
          var b = document.createElement("b");
          b.textContent = num(line.E_keV, 3);
          cell2.appendChild(b);
        } else {
          cell2.textContent = num(line.E_keV, 3);
        }
        row.appendChild(cell2);

        if (line.device) {
          var cell3 = document.createElement("td");
          cell3.textContent = num(line.device.E_keV, 1);
          row.appendChild(cell3);

          var cell4 = document.createElement("td");
          cell4.textContent = num(line.device.fwhm_keV, 2) + " / " + num(line.fwhm_curve_keV, 2);
          row.appendChild(cell4);

          var cell5 = document.createElement("td");
          cell5.textContent = num(line.fwhm_ratio, 3);
          if (line.fwhm_ratio < 0.9 || line.fwhm_ratio > 1.1) {
            row.className = "row-dirty";
          }
          row.appendChild(cell5);

          var cell6 = document.createElement("td");
          cell6.textContent = cnt(line.device.area) + " ± " + cnt(line.device.d_area);
          row.appendChild(cell6);

          var cell7 = document.createElement("td");
          cell7.textContent = num(line.eff_meas_pct, 4);
          row.appendChild(cell7);

          var cell8 = document.createElement("td");
          cell8.textContent = num(line.eff_model_pct, 4);
          row.appendChild(cell8);

          var cell9 = document.createElement("td");
          cell9.textContent = num(line.ratio, 3) + " ± " + num(line.d_ratio_stat, 3);
          row.appendChild(cell9);

          if (line.main) {
            // B-14: здесь стояло area/d_area; столбец «путь 1» — отношение пути 1 к паспорту.
            var cell10 = document.createElement("td");
            var g1 = D.method1.groups[key];
            cell10.textContent = g1 ? num(g1.A_over_passport, 3) : "?";
            row.appendChild(cell10);
          } else {
            var cell10 = document.createElement("td");
            cell10.textContent = "";
            row.appendChild(cell10);
          }

          if (line.blends && line.blends.length > 0) {
            var list = [];
            for (var k = 0; k < line.blends.length; k++) {
              var b = line.blends[k];
              if (b.nuclide) {
                list.push(num(b.E_keV, 1) + " (" + labelRu(b.nuclide) + ")");
              } else {
                list.push(num(b.E_keV, 1) + " (сумма)");
              }
            }
            var cell11 = document.createElement("td");
            cell11.textContent = list.join(", ");
            row.appendChild(cell11);
          } else {
            var cell11 = document.createElement("td");
            cell11.textContent = "-";
            row.appendChild(cell11);
          }

        } else {
          var cell3 = document.createElement("td");
          cell3.textContent = "нет пика в таблице прибора";
          row.appendChild(cell3);

          var cell4 = document.createElement("td");
          cell4.textContent = "";
          row.appendChild(cell4);

          var cell5 = document.createElement("td");
          cell5.textContent = "";
          row.appendChild(cell5);

          var cell6 = document.createElement("td");
          cell6.textContent = "";
          row.appendChild(cell6);

          var cell7 = document.createElement("td");
          cell7.textContent = "";
          row.appendChild(cell7);

          var cell8 = document.createElement("td");
          cell8.textContent = "";
          row.appendChild(cell8);

          var cell9 = document.createElement("td");
          cell9.textContent = "";
          row.appendChild(cell9);

          var cell10 = document.createElement("td");
          cell10.textContent = "";
          row.appendChild(cell10);

          var cell11 = document.createElement("td");
          cell11.textContent = "";
          row.appendChild(cell11);
        }

        tbody.appendChild(row);
      }
    }

    tblM2.appendChild(tbody);
  }

  function renderAmCheck() {
    var table = document.getElementById("tblAmCheck");
    if (!table) return;
    table.innerHTML = "";
    var thead = document.createElement("thead");
    var headRow = document.createElement("tr");
    var th1 = document.createElement("th");
    th1.textContent = "запись";
    headRow.appendChild(th1);
    var th2 = document.createElement("th");
    th2.textContent = "нуклид";
    headRow.appendChild(th2);
    var th3 = document.createElement("th");
    th3.textContent = "E, кэВ";
    headRow.appendChild(th3);
    var th4 = document.createElement("th");
    th4.textContent = "пик прибора, кэВ";
    headRow.appendChild(th4);
    var th5 = document.createElement("th");
    th5.textContent = "ПШПВ прибора";
    headRow.appendChild(th5);
    var th6 = document.createElement("th");
    th6.textContent = "ПШПВ по кривой";
    headRow.appendChild(th6);
    var th7 = document.createElement("th");
    th7.textContent = "отношение ширин";
    headRow.appendChild(th7);
    var th8 = document.createElement("th");
    th8.textContent = "A, Бк";
    headRow.appendChild(th8);
    var th9 = document.createElement("th");
    th9.textContent = "эфф. измер., %";
    headRow.appendChild(th9);
    var th10 = document.createElement("th");
    th10.textContent = "эфф. модели, %";
    headRow.appendChild(th10);
    var th11 = document.createElement("th");
    th11.textContent = "изм./модель";
    headRow.appendChild(th11);
    thead.appendChild(headRow);
    table.appendChild(thead);

    var tbody = document.createElement("tbody");

    for (var i = 0; i < D.cross_geometry.length; i++) {
      var g = D.cross_geometry[i];
      var row = document.createElement("tr");

      var cell1 = document.createElement("td");
      cell1.textContent = esc(g.record);
      row.appendChild(cell1);

      var cell2 = document.createElement("td");
      cell2.innerHTML = sw(colorOf(g.nuclide)) + " " + esc(labelRu(g.nuclide));
      row.appendChild(cell2);

      var cell3 = document.createElement("td");
      cell3.textContent = num(g.E_keV, 1);
      row.appendChild(cell3);

      if (g.device) {
        var cell4 = document.createElement("td");
        cell4.textContent = num(g.device.E_keV, 2);
        row.appendChild(cell4);

        var cell5 = document.createElement("td");
        cell5.textContent = num(g.device.fwhm_keV, 2);
        row.appendChild(cell5);

        var cell6 = document.createElement("td");
        cell6.textContent = num(g.fwhm_curve_keV, 2);
        row.appendChild(cell6);

        var cell7 = document.createElement("td");
        cell7.textContent = num(g.fwhm_ratio, 3);
        if (g.fwhm_ratio < 0.9 || g.fwhm_ratio > 1.1) {
          row.className = "row-dirty";
        }
        row.appendChild(cell7);

        var cell8 = document.createElement("td");
        cell8.textContent = cnt(g.A_Bq);
        row.appendChild(cell8);

        var cell9 = document.createElement("td");
        cell9.textContent = num(g.eff_meas_pct, 4);
        row.appendChild(cell9);

        var cell10 = document.createElement("td");
        cell10.textContent = num(g.eff_model_pct, 4);
        row.appendChild(cell10);

        var cell11 = document.createElement("td");
        cell11.textContent = num(g.ratio, 3);
        row.appendChild(cell11);
      } else {
        var cell4 = document.createElement("td");
        cell4.textContent = "";
        row.appendChild(cell4);

        var cell5 = document.createElement("td");
        cell5.textContent = "";
        row.appendChild(cell5);

        var cell6 = document.createElement("td");
        cell6.textContent = "";
        row.appendChild(cell6);

        var cell7 = document.createElement("td");
        cell7.textContent = "";
        row.appendChild(cell7);

        var cell8 = document.createElement("td");
        cell8.textContent = "";
        row.appendChild(cell8);

        var cell9 = document.createElement("td");
        cell9.textContent = "";
        row.appendChild(cell9);

        var cell10 = document.createElement("td");
        cell10.textContent = "";
        row.appendChild(cell10);

        var cell11 = document.createElement("td");
        cell11.textContent = "";
        row.appendChild(cell11);
      }

      tbody.appendChild(row);
    }

    table.appendChild(tbody);
  }

  function coefsHtml(arr) {
    if (!arr || !arr.length) return "";
    var html = [];
    for (var i = 0; i < arr.length; i++) {
      var c = arr[i];
      var s;
      if (Math.abs(c) >= 0.01 && Math.abs(c) <= 10000) {
        s = num(c, 6);
      } else {
        s = c.toExponential(4).replace(".", ",");
      }
      html.push("<span class='mono'>c" + i + " = " + s + "</span>");
    }
    return html.join("<br>");
  }

  function renderCalTables() {
    var tblCal = document.getElementById("tblCal");
    if (!tblCal) return;
    tblCal.innerHTML = "";
    var tbody = document.createElement("tbody");

    var row1 = document.createElement("tr");
    var cell1 = document.createElement("td");
    cell1.textContent = "живое время, с";
    row1.appendChild(cell1);
    var cell2 = document.createElement("td");
    cell2.textContent = num(D.meta.live_s, 2);
    row1.appendChild(cell2);
    var cell3 = document.createElement("td");
    cell3.textContent = num(D.meta.bg_live_s, 2);
    row1.appendChild(cell3);
    tbody.appendChild(row1);

    var row2 = document.createElement("tr");
    var cell1 = document.createElement("td");
    cell1.textContent = "реальное время, с";
    row2.appendChild(cell1);
    var cell2 = document.createElement("td");
    cell2.textContent = num(D.meta.real_s, 2);
    row2.appendChild(cell2);
    var cell3 = document.createElement("td");
    cell3.textContent = num(D.meta.bg_real_s, 2);
    row2.appendChild(cell3);   // B-17: дважды добавлялся cell2
    tbody.appendChild(row2);

    var row3 = document.createElement("tr");
    var cell1 = document.createElement("td");
    cell1.textContent = "мёртвое время, %";
    row3.appendChild(cell1);
    var cell2 = document.createElement("td");
    cell2.textContent = num((D.meta.real_s - D.meta.live_s) / D.meta.real_s * 100, 3);
    row3.appendChild(cell2);
    var cell3 = document.createElement("td");
    cell3.textContent = num((D.meta.bg_real_s - D.meta.bg_live_s) / D.meta.bg_real_s * 100, 3);
    row3.appendChild(cell3);
    tbody.appendChild(row3);

    var row4 = document.createElement("tr");
    var cell1 = document.createElement("td");
    cell1.textContent = "шкала файла";
    row4.appendChild(cell1);
    var cell2 = document.createElement("td");
    cell2.innerHTML = coefsHtml(D.meta.cal_sample.coefs_file);
    row4.appendChild(cell2);
    var cell3 = document.createElement("td");
    // B-15: здесь стояли коэффициенты рефита под подписью «шкала файла».
    cell3.innerHTML = coefsHtml(D.meta.cal_bg.coefs_file) + "<br>рефит по якорям: " +
      coefsHtml(D.meta.cal_bg.coefs_refit) + "<br>СКО " + num(D.meta.cal_bg.rms_keV, 2) + " кэВ";
    row4.appendChild(cell3);
    tbody.appendChild(row4);

    var row5 = document.createElement("tr");
    var cell1 = document.createElement("td");
    cell1.textContent = "масштаб фона по времени";
    row5.appendChild(cell1);
    var cell5v = document.createElement("td");   // B-18: значение не выводилось
    cell5v.colSpan = 2;
    cell5v.textContent = D.meta.bg_scale_time != null ? num(D.meta.bg_scale_time, 5) : "?";
    row5.appendChild(cell5v);
    tbody.appendChild(row5);

    tblCal.appendChild(tbody);

    var tblRefs = document.getElementById("tblRefs");
    if (!tblRefs) return;
    tblRefs.innerHTML = "";
    var thead = document.createElement("thead");
    var headRow = document.createElement("tr");
    var th1 = document.createElement("th");
    th1.textContent = "реперы шкалы пробы";
    th1.colSpan = 5;
    headRow.appendChild(th1);
    thead.appendChild(headRow);
    tblRefs.appendChild(thead);

    var tbody = document.createElement("tbody");

    if (D.meta.cal_sample.refs) {
      for (var i = 0; i < D.meta.cal_sample.refs.length; i++) {
        var ref = D.meta.cal_sample.refs[i];
        var row = document.createElement("tr");
        var cell1 = document.createElement("td");
        cell1.textContent = num(ref.E_true, 3);
        row.appendChild(cell1);
        var cell2 = document.createElement("td");
        cell2.textContent = num(ref.ch, 2);
        row.appendChild(cell2);
        var cell3 = document.createElement("td");
        cell3.textContent = num(ref.E_file, 2);
        row.appendChild(cell3);
        var cell4 = document.createElement("td");
        cell4.textContent = num(ref.shift_keV, 2);
        row.appendChild(cell4);
        var cell5 = document.createElement("td");
        cell5.textContent = num(ref.shift_fwhm, 2);
        row.appendChild(cell5);
        tbody.appendChild(row);
      }
    }

    // B-16: якоря фона — своим блоком с заголовками столбцов, не в thead и не в tbody реперов пробы.
    tblRefs.appendChild(tbody);
    tbody = document.createElement("tbody");
    var hr2 = document.createElement("tr");
    ["якоря шкалы фона: E, кэВ", "канал", "нетто-счёт", "невязка, кэВ"].forEach(function (h) {
      var th = document.createElement("th");
      th.textContent = h;
      hr2.appendChild(th);
    });
    tbody.appendChild(hr2);

    if (D.meta.cal_bg.anchors) {
      for (var i = 0; i < D.meta.cal_bg.anchors.length; i++) {
        var anchor = D.meta.cal_bg.anchors[i];
        var row = document.createElement("tr");
        var cell1 = document.createElement("td");
        cell1.textContent = num(anchor.E_keV, 2);
        row.appendChild(cell1);
        var cell2 = document.createElement("td");
        cell2.textContent = num(anchor.ch, 2);
        row.appendChild(cell2);
        var cell3 = document.createElement("td");
        cell3.textContent = cnt(anchor.netsum);
        row.appendChild(cell3);
        var cell4 = document.createElement("td");
        cell4.textContent = num(anchor.resid_keV, 2);
        row.appendChild(cell4);
        tbody.appendChild(row);
      }
    }

    tblRefs.appendChild(tbody);
  }

  function renderFwhm() {
    var tblFwhm = document.getElementById("tblFwhm");
    if (!tblFwhm) return;
    tblFwhm.innerHTML = "";
    var thead = document.createElement("thead");
    var headRow = document.createElement("tr");
    var th1 = document.createElement("th");
    th1.textContent = "E, кэВ";
    headRow.appendChild(th1);
    var th2 = document.createElement("th");
    th2.textContent = "источник";
    headRow.appendChild(th2);
    var th3 = document.createElement("th");
    th3.textContent = "ПШПВ измер., кэВ";
    headRow.appendChild(th3);
    var th4 = document.createElement("th");
    th4.textContent = "закон √E, кэВ";
    headRow.appendChild(th4);
    var th5 = document.createElement("th");
    th5.textContent = "отклонение √E, %";
    headRow.appendChild(th5);
    thead.appendChild(headRow);
    tblFwhm.appendChild(thead);

    var tbody = document.createElement("tbody");

    if (D.fwhm_cal.points) {
      for (var i = 0; i < D.fwhm_cal.points.length; i++) {
        var p = D.fwhm_cal.points[i];
        var row = document.createElement("tr");
        var cell1 = document.createElement("td");
        cell1.textContent = num(p.E_keV, 1);
        row.appendChild(cell1);

        var cell2 = document.createElement("td");
        cell2.textContent = esc(p.source);
        row.appendChild(cell2);

        var cell3 = document.createElement("td");
        cell3.textContent = num(p.fwhm_keV, 2) + " ± " + num(p.d_fwhm_keV, 2);
        row.appendChild(cell3);

        var cell4 = document.createElement("td");
        cell4.textContent = num(interp(D.fwhm_cal.sqrt_law, p.E_keV), 2);
        row.appendChild(cell4);

        var cell5 = document.createElement("td");
        var sqrtE = interp(D.fwhm_cal.sqrt_law, p.E_keV);
        if (sqrtE > 0) {
          var diff = (sqrtE - p.fwhm_keV) / p.fwhm_keV * 100;
          cell5.textContent = (diff >= 0 ? "+" : "") + num(diff, 1);
        } else {
          cell5.textContent = "?";
        }
        row.appendChild(cell5);

        tbody.appendChild(row);
      }
    }

    tblFwhm.appendChild(tbody);

    var tblFwhmDev = document.getElementById("tblFwhmDev");
    if (!tblFwhmDev) return;
    tblFwhmDev.innerHTML = "";
    var thead = document.createElement("thead");
    var headRow = document.createElement("tr");
    var th1 = document.createElement("th");
    th1.textContent = "пики этой записи (таблица прибора)";
    th1.colSpan = 4;
    headRow.appendChild(th1);
    thead.appendChild(headRow);
    tblFwhmDev.appendChild(thead);

    var tbody = document.createElement("tbody");

    if (D.fwhm_cal.device_peaks) {
      for (var i = 0; i < D.fwhm_cal.device_peaks.length; i++) {
        var p = D.fwhm_cal.device_peaks[i];
        var row = document.createElement("tr");
        var cell1 = document.createElement("td");
        cell1.textContent = num(p.E_keV, 1);
        row.appendChild(cell1);

        var cell2 = document.createElement("td");
        cell2.textContent = num(p.fwhm_keV, 2);
        row.appendChild(cell2);

        var cell3 = document.createElement("td");
        cell3.textContent = num(p.fwhm_curve_keV, 2);
        row.appendChild(cell3);

        var cell4 = document.createElement("td");
        cell4.textContent = num(p.ratio, 3);
        if (p.ratio < 0.9 || p.ratio > 1.1) {
          row.className = "row-dirty";
        }
        row.appendChild(cell4);

        tbody.appendChild(row);
      }
    }

    tblFwhmDev.appendChild(tbody);

    var cv = document.getElementById("cvFwhm");
    if (!cv) return;
    var fit = window.G1SA.fit(cv);
    var g = fit.g, w = fit.w, h = fit.h;

    var max = 0;
    if (D.fwhm_cal.points) {
      for (var i = 0; i < D.fwhm_cal.points.length; i++) {
        if (D.fwhm_cal.points[i].fwhm_keV > max) max = D.fwhm_cal.points[i].fwhm_keV;
      }
    }
    if (D.fwhm_cal.curve) {
      for (var i = 0; i < D.fwhm_cal.curve.length; i++) {
        if (D.fwhm_cal.curve[i][1] > max) max = D.fwhm_cal.curve[i][1];
      }
    }
    if (D.fwhm_cal.sqrt_law) {
      for (var i = 0; i < D.fwhm_cal.sqrt_law.length; i++) {
        if (D.fwhm_cal.sqrt_law[i][1] > max) max = D.fwhm_cal.sqrt_law[i][1];
      }
    }
    if (D.fwhm_cal.device_peaks) {
      for (var i = 0; i < D.fwhm_cal.device_peaks.length; i++) {
        if (D.fwhm_cal.device_peaks[i].fwhm_keV > max) max = D.fwhm_cal.device_peaks[i].fwhm_keV;
      }
    }

    g.clearRect(0, 0, w, h);
    var m = {l: 62, r: 16, t: 14, b: 34};
    var xScale = (w - m.l - m.r) / 3000;
    var yScale = (h - m.t - m.b) / (max * 1.15);

    g.strokeStyle = window.G1SA.pal().grid;
    g.lineWidth = 1;

    for (var x = 0; x <= 3000; x += 500) {
      var px = m.l + x * xScale;
      g.beginPath();
      g.moveTo(px, m.t);
      g.lineTo(px, h - m.b);
      g.stroke();
    }

    for (var y = 0; y <= max * 1.15; y += 20) {
      var py = h - m.b - y * yScale;
      g.beginPath();
      g.moveTo(m.l, py);
      g.lineTo(w - m.r, py);
      g.stroke();
    }

    g.fillStyle = window.G1SA.pal().faint;
    g.font = "11px system-ui";
    g.textAlign = "center";

    for (var x = 0; x <= 3000; x += 500) {
      var px = m.l + x * xScale;
      g.fillText(x, px, h - m.b + 12);
    }

    g.textAlign = "right";
    g.textBaseline = "middle";

    for (var y = 0; y <= max * 1.15; y += 20) {
      var py = h - m.b - y * yScale;
      g.fillText(num(y, 0), m.l - 4, py);
    }

    g.textAlign = "center";
    g.textBaseline = "top";
    g.font = "bold 11px system-ui";

    g.fillText("энергия, кэВ", w / 2, h - 4);
    g.save();
    g.translate(4, h / 2);
    g.rotate(-Math.PI / 2);
    g.fillText("ПШПВ, кэВ", 0, 0);
    g.restore();

    g.strokeStyle = window.G1SA.pal().rule;
    g.lineWidth = 2;
    g.beginPath();
    g.moveTo(m.l, m.t);
    g.lineTo(m.l, h - m.b);
    g.lineTo(w - m.r, h - m.b);
    g.stroke();

    if (D.fwhm_cal.curve) {
      g.strokeStyle = "#0f5aa8";
      g.lineWidth = 2;
      g.beginPath();
      for (var i = 0; i < D.fwhm_cal.curve.length; i++) {
        var p = D.fwhm_cal.curve[i];
        var x = m.l + p[0] * xScale;
        var y = h - m.b - p[1] * yScale;
        if (i === 0) g.moveTo(x, y);
        else g.lineTo(x, y);
      }
      g.stroke();
    }

    if (D.fwhm_cal.sqrt_law) {
      g.strokeStyle = window.G1SA.pal().faint;
      g.lineWidth = 1.5;
      g.setLineDash([5, 4]);
      g.beginPath();
      for (var i = 0; i < D.fwhm_cal.sqrt_law.length; i++) {
        var p = D.fwhm_cal.sqrt_law[i];
        var x = m.l + p[0] * xScale;
        var y = h - m.b - p[1] * yScale;
        if (i === 0) g.moveTo(x, y);
        else g.lineTo(x, y);
      }
      g.stroke();
      g.setLineDash([]);
    }

    if (D.fwhm_cal.points) {
      g.fillStyle = "#c8541c";
      for (var i = 0; i < D.fwhm_cal.points.length; i++) {
        var p = D.fwhm_cal.points[i];
        var x = m.l + p.E_keV * xScale;
        var y = h - m.b - p.fwhm_keV * yScale;
        g.beginPath();
        g.arc(x, y, 4, 0, Math.PI * 2);
        g.fill();
      }
    }

    if (D.fwhm_cal.device_peaks) {
      g.strokeStyle = window.G1SA.pal().ink;
      g.fillStyle = "#fff";
      g.lineWidth = 1;
      for (var i = 0; i < D.fwhm_cal.device_peaks.length; i++) {
        var p = D.fwhm_cal.device_peaks[i];
        var x = m.l + p.E_keV * xScale;
        var y = h - m.b - p.fwhm_keV * yScale;
        g.beginPath();
        g.arc(x, y, 4, 0, Math.PI * 2);
        g.stroke();
      }
    }

    g.fillStyle = window.G1SA.pal().ink;
    g.font = "600 11px system-ui, sans-serif";
    g.textAlign = "left";

    var legendY = m.t + 4;
    g.fillText("измеренные точки комплекта", m.l + 4, legendY);
    legendY += 12;
    g.fillText("кривая (интерполяция в логарифмах)", m.l + 4, legendY);
    legendY += 12;
    g.fillText("прежний закон √E", m.l + 4, legendY);
    legendY += 12;
    g.fillText("пики этой записи (прибор)", m.l + 4, legendY);
  }

  function showNet() {
    renderNet();
  }

  function showCmp() {
    renderAmCheck();
  }

  function showCal() {
    renderCalTables();
    renderFwhm();
  }

  function redraw() {
    if (document.getElementById("viewCal") && document.getElementById("viewCal").style.display !== "none") {
      renderFwhm();
    }
  }

  window.G1SA_RERUN = {
    showNet: showNet,
    showCmp: showCmp,
    showCal: showCal,
    redraw: redraw
  };

  fillValues();
  renderTwoPaths();
  renderMethod2();
  renderBands();

  window.addEventListener("resize", redraw);
})();
