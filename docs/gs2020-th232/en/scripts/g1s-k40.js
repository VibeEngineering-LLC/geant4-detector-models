(function () {
    "use strict";
    var D = window.GS_K40;
    var PFX = "k40-";
    if (!D) {
        console.error("window.GS_K40 is not defined");
        return;
    }
    function num(x, d) {
        if (d === undefined)
            d = 1;
        return (Number(x)).toFixed(d).replace(/^(-?)(\d{5,})/, function (m, s, i) { return s + i.replace(/\B(?=(\d{3})+(?!\d))/g, "\u202F"); }).replace(/^-/, "\u2212");
    }
    function cnt(x) {
        var s = String(Math.round(Number(x))), out = "", neg = s[0] === "-";
        if (neg)
            s = s.slice(1);
        if (s.length > 4) while (s.length > 3) {
            out = " " + s.slice(-3) + out;
            s = s.slice(0, -3);
        }
        return (neg ? "−" : "") + s + out;
    }
    function signedPct(ratio) {
        var s = 100 * (ratio - 1);
        return (s < 0 ? "−" : "+") + num(Math.abs(s), 1) + " %";
    }
    function esc(s) {
        return String(s === undefined || s === null ? "" : s)
            .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
    }
    var ST = {
        on: {},
        log: true,
        sum: true,
        fwhmLaw: "cs",
        lib: "fixed",
        m2sort: "contrib",
        cursorE: null,
    };
    D.nuclides.forEach(function (n) { ST.on[n.key] = true; });
    function SRC() { return ST.fwhmLaw === "cs" ? D.cs : D; }
    function SPEC() { return ST.fwhmLaw === "cs" ? D.cs.spectrum : D.spectrum; }
    function M1() { return SRC().method1; }
    function M2() {
        return ST.lib === "full" ? SRC().method2_full : SRC().method2;
    }
    function STACK1() { return SPEC().stack; }
    function STACK2() {
        return ST.lib === "full" ? SPEC().stack2_full : SPEC().stack2;
    }
    function STACK2CHAN() {
        return ST.lib === "full" ? SPEC().stack2_chan_full : SPEC().stack2_chan;
    }
    function TRUSTED() { return SPEC().trusted || D.spectrum.trusted || null; }
    function TRUST_MASK(key) {
        var t = TRUSTED();
        var v = t && t[key];
        return (v && v.length === D.spectrum.e_of_ch.length) ? v : null;
    }
    function fit(cv) {
        var dpr = window.devicePixelRatio || 1;
        var r = cv.getBoundingClientRect();
        var w = Math.max(200, Math.floor(r.width));
        var h = Math.max(120, Math.floor(r.height || cv.height || 260));
        cv.width = Math.floor(w * dpr);
        cv.height = Math.floor(h * dpr);
        var g = cv.getContext("2d");
        g.setTransform(dpr, 0, 0, dpr, 0, 0);
        g.clearRect(0, 0, w, h);
        return { g: g, w: w, h: h };
    }
    function css(name, fallback) {
        var v = getComputedStyle(document.documentElement).getPropertyValue(name);
        return (v || "").trim() || fallback;
    }
    function pal() {
        return {
            ink: css("--ink", "#16140f"),
            rule: css("--rule", "#16140f"),
            faint: css("--faint", "#6a6558"),
            grid: css("--grid", "#dcd7c8"),
            paper: css("--paper", "#f5f2ea"),
            sum: css("--sum-line", "#d21f1f"),
        };
    }
    function mapX(v, lo, hi, x0, x1) {
        return x0 + (v - lo) / (hi - lo) * (x1 - x0);
    }
    function makeY(logY, lo, hi) {
        if (logY) {
            lo = Math.max(0.5, lo);
            hi = Math.max(lo * 10, hi);
            var l0 = Math.log10(lo), l1 = Math.log10(hi);
            return {
                lo: lo, hi: hi, log: true,
                map: function (v, y0, y1) {
                    var t = (Math.log10(Math.max(v, lo)) - l0) / (l1 - l0);
                    return y0 + (1 - t) * (y1 - y0);
                }
            };
        }
        return {
            lo: 0, hi: hi, log: false,
            map: function (v, y0, y1) { return y0 + (1 - v / hi) * (y1 - y0); }
        };
    }
    function supNum(n) {
        var d = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵",
            "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
        return String(n).split("").map(function (c) { return d[c] || c; }).join("");
    }
    function stackTotal(stk, i) {
        var acc = 0;
        for (var j = 0; j < D.nuclides.length; j++) {
            var k = D.nuclides[j].key;
            if (!ST.on[k] || !stk[k])
                continue;
            acc += stk[k][i];
        }
        return acc;
    }
    function drawSpectrum(cvId, stk, extra) {
        var cv = document.getElementById(PFX + cvId);
        if (!cv)
            return null;
        var p = pal();
        var f = fit(cv);
        var g = f.g, W = f.w, H = f.h;
        var m = { l: 62, r: 14, t: 12, b: 34 };
        var e = D.spectrum.e_of_ch;
        var yy = D.spectrum.counts;
        var xLo = ST.zoom ? ST.zoom.xLo : 25, xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];
        var vMax = 1;
        for (var i0 = 0; i0 < e.length; i0++) {
            if (e[i0] < 30 || e[i0] > 3500)
                continue;
            var v0 = Math.max(yy[i0], stackTotal(stk, i0));
            if (v0 > vMax)
                vMax = v0;
        }
        var Y = makeY(ST.log, ST.log ? 0.5 : 0, vMax * (ST.log ? 2.0 : 1.1));
        g.strokeStyle = p.grid;
        g.lineWidth = 1;
        g.beginPath();
        var xTicks = [250, 500, 750, 1000, 1250, 1500, 1750, 2000,
            2250, 2500, 2750, 3000, 3250, 3500];
        for (var xi = 0; xi < xTicks.length; xi++) {
            if (xTicks[xi] > xHi)
                break;
            var xx = mapX(xTicks[xi], xLo, xHi, m.l, W - m.r);
            g.moveTo(xx, m.t);
            g.lineTo(xx, H - m.b);
        }
        var yTicks = [];
        if (ST.log) {
            for (var d = 0; d <= Math.ceil(Math.log10(Y.hi)); d++)
                yTicks.push(Math.pow(10, d));
        }
        else {
            var stp = Math.pow(10, Math.floor(Math.log10(Y.hi / 4)));
            var s0 = Math.ceil(Y.hi / 4 / stp) * stp;
            for (var v = s0; v < Y.hi; v += s0)
                yTicks.push(v);
        }
        for (var yi = 0; yi < yTicks.length; yi++) {
            var yg = Y.map(yTicks[yi], m.t, H - m.b);
            g.moveTo(m.l, yg);
            g.lineTo(W - m.r, yg);
        }
        g.stroke();
        g.fillStyle = p.faint;
        g.font = "11px system-ui, sans-serif";
        g.textAlign = "center";
        g.textBaseline = "top";
        for (var xj = 0; xj < xTicks.length; xj++) {
            if (xTicks[xj] > xHi)
                break;
            g.fillText(String(xTicks[xj]), mapX(xTicks[xj], xLo, xHi, m.l, W - m.r), H - m.b + 4);
        }
        g.textAlign = "right";
        g.textBaseline = "middle";
        for (var yj = 0; yj < yTicks.length; yj++) {
            var lbl = ST.log ? "10" + supNum(Math.round(Math.log10(yTicks[yj])))
                : cnt(yTicks[yj]);
            g.fillText(lbl, m.l - 4, Y.map(yTicks[yj], m.t, H - m.b));
        }
        g.textAlign = "center";
        g.textBaseline = "bottom";
        g.fillText("energy, keV", (m.l + W - m.r) / 2, H - 2);
        g.strokeStyle = p.rule;
        g.lineWidth = 2;
        g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);
        g.save();
        g.beginPath();
        g.rect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);
        g.clip();
        var order = [];
        for (var ni = 0; ni < D.nuclides.length; ni++) {
            var nk = D.nuclides[ni].key;
            if (!ST.on[nk] || !stk[nk])
                continue;
            var s = 0;
            for (var si = 0; si < e.length; si++) {
                if (e[si] < xLo || e[si] > xHi)
                    continue;
                s += stk[nk][si];
            }
            order.push({ nuc: D.nuclides[ni], area: s });
        }
        order.sort(function (a, b) { return b.area - a.area; });
        function edgeLo(i) {
            return mapX(i > 0 ? (e[i - 1] + e[i]) / 2 : e[0], xLo, xHi, m.l, W - m.r);
        }
        function edgeHi(i) {
            return mapX(i + 1 < e.length ? (e[i] + e[i + 1]) / 2 : e[e.length - 1], xLo, xHi, m.l, W - m.r);
        }
        function runs(msk, want) {
            var out = [], i2 = 0, n2 = e.length;
            while (i2 < n2) {
                var ok = msk ? !!msk[i2] : true;
                var j2 = i2;
                while (j2 + 1 < n2 && (msk ? !!msk[j2 + 1] : true) === ok)
                    j2++;
                if (ok === want)
                    out.push([i2, j2]);
                i2 = j2 + 1;
            }
            return out;
        }
        function segPath(vec, a, b) {
            g.beginPath();
            for (var q = a; q <= b; q++) {
                var x = mapX(e[q], xLo, xHi, m.l, W - m.r);
                var y = Y.map(vec[q] > Y.lo ? vec[q] : Y.lo, m.t, H - m.b);
                if (q === a)
                    g.moveTo(edgeLo(a), y);
                g.lineTo(x, y);
                if (q === b)
                    g.lineTo(edgeHi(b), y);
            }
        }
        var yBase = Y.map(Y.lo, m.t, H - m.b);
        function fillRuns(vec, rr, style) {
            g.fillStyle = style;
            for (var ri = 0; ri < rr.length; ri++) {
                segPath(vec, rr[ri][0], rr[ri][1]);
                g.lineTo(edgeHi(rr[ri][1]), yBase);
                g.lineTo(edgeLo(rr[ri][0]), yBase);
                g.closePath();
                g.fill();
            }
        }
        for (var oi = 0; oi < order.length; oi++) {
            var kf = order[oi].nuc.key;
            var vec = stk[kf];
            var mf = TRUST_MASK(kf);
            fillRuns(vec, [[0, e.length - 1]], order[oi].nuc.color);
            strokeLayer(oi);
        }
        function strokeLayer(oj) {
            var ks = order[oj].nuc.key;
            var vs = stk[ks];
            var msk = TRUST_MASK(ks);
            g.strokeStyle = order[oj].nuc.color;
            var solid = runs(msk, true);
            g.lineWidth = 1.2;
            g.setLineDash([]);
            for (var sj = 0; sj < solid.length; sj++) {
                segPath(vs, solid[sj][0], solid[sj][1]);
                g.stroke();
            }
            var noisy = runs(msk, false);
            g.lineWidth = 0.8;
            g.setLineDash([2, 3]);
            for (var nj = 0; nj < noisy.length; nj++) {
                segPath(vs, noisy[nj][0], noisy[nj][1]);
                g.stroke();
            }
            g.setLineDash([]);
        }
        function trace(getV, color, width, dash) {
            g.strokeStyle = color;
            g.lineWidth = width;
            if (dash)
                g.setLineDash(dash);
            g.beginPath();
            var st = false;
            for (var k = 0; k < e.length; k++) {
                if (e[k] < xLo || e[k] > xHi)
                    continue;
                var v = getV(k);
                var x = mapX(e[k], xLo, xHi, m.l, W - m.r);
                var y = Y.map(v > Y.lo ? v : Y.lo, m.t, H - m.b);
                if (!st) {
                    g.moveTo(x, y);
                    st = true;
                }
                else
                    g.lineTo(x, y);
            }
            g.stroke();
            if (dash)
                g.setLineDash([]);
        }
        if (ST.sum) {
            trace(function (i) { return stackTotal(stk, i); }, p.sum, 1.8, [6, 3]);
        }
        (extra || []).forEach(function (ex) {
            if (!ex.arr)
                return;
            trace(function (i) { return ex.arr[i]; }, ex.color, ex.width || 1.4, ex.dash);
        });
        trace(function (i) { return yy[i]; }, p.ink, 0.8);
        if (ST.cursorE !== null) {
            var xC = mapX(ST.cursorE, xLo, xHi, m.l, W - m.r);
            g.strokeStyle = p.rule;
            g.lineWidth = 1;
            g.setLineDash([4, 4]);
            g.beginPath();
            g.moveTo(xC, m.t);
            g.lineTo(xC, H - m.b);
            g.stroke();
            g.setLineDash([]);
        }
        if (ST.drag) {
            var xa3 = mapX(Math.min(ST.drag.e0, ST.drag.e1), xLo, xHi, m.l, W - m.r);
            var xb3 = mapX(Math.max(ST.drag.e0, ST.drag.e1), xLo, xHi, m.l, W - m.r);
            g.fillStyle = "rgba(246,211,28,.22)";
            g.fillRect(xa3, m.t, xb3 - xa3, H - m.b - m.t);
            g.strokeStyle = "#16140f";
            g.lineWidth = 1.5;
            g.setLineDash([4, 4]);
            g.strokeRect(xa3, m.t, xb3 - xa3, H - m.b - m.t);
            g.setLineDash([]);
        }
        g.restore();
        g.strokeStyle = p.rule;
        g.lineWidth = 2;
        g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);
        return { m: m, W: W, H: H, xLo: xLo, xHi: xHi, Y: Y };
    }
    function nucOf(key) {
        for (var i = 0; i < D.nuclides.length; i++)
            if (D.nuclides[i].key === key)
                return D.nuclides[i];
        return null;
    }
    function noiseFrac(key) {
        var nf = SPEC().noise_frac || D.spectrum.noise_frac;
        if (!nf || !(key in nf))
            return null;
        return nf[key];
    }
    function guardHint(nf) {
        return "the MC template statistics are sufficient for " + num(100 * (1 - nf), 1) + " %"
            + " of the component contribution; the remainder corresponds to regions below "
            + (D.spectrum.n_eff_min || 0) + " MC counts per channel, where the component is shown "
            + "as a dashed line without shading (the nuclide fraction is determined by template noise)";
    }
    function buildLegend(elId) {
        var el = document.getElementById(PFX + elId);
        if (!el || el.dataset.built)
            return;
        el.dataset.built = "1";
        var html = "";
        var stkL = (elId === 'legendM2') ? STACK2() : STACK1();
        var arrL = D.nuclides.map(function (n, i) {
            var v = stkL && stkL[n.key], s = 0;
            if (v)
                for (var q = 0; q < v.length; q++)
                    s += v[q];
            return { n: n, s: s, i: i };
        });
        arrL.sort(function (a, b) { return (b.s - a.s) || (a.i - b.i); });
        arrL.map(function (o) { return o.n; }).forEach(function (nuc) {
            if (nuc.key === "SECOND")
                return;
            var nf = noiseFrac(nuc.key);
            var sw = nuc.color;
            var hint = (nf === null) ? "" : " title='" + esc(guardHint(nf)) + "'";
            html += "<label class='chip' data-nuc='" + nuc.key + "'" + hint + ">"
                + "<input type='checkbox' " + (ST.on[nuc.key] ? "checked" : "") + ">"
                + "<span class='sw' style='background:" + sw + "'></span>"
                + "<span class='nm'>" + esc(nuc.label_ru) + "</span>"
                + "</label>";
        });
        html += "<label class='chip'>"
            + "<input type='checkbox' class='c-sum' " + (ST.sum ? "checked" : "") + ">"
            + "<span class='sw sum'></span>"
            + "<span class='nm'>sum</span></label>";
        html += "<label class='chip toggle'>"
            + "<input type='checkbox' class='c-log' " + (ST.log ? "checked" : "") + ">"
            + "<span class='sw log'></span><span class='nm'>logarithmic scale</span></label>";
        el.innerHTML = html;
        el.querySelectorAll("label.chip[data-nuc] input").forEach(function (inp) {
            var key = inp.parentNode.getAttribute("data-nuc");
            inp.addEventListener("change", function (ev) {
                ST.on[key] = ev.target.checked;
                syncLegends();
                redraw();
            });
        });
        el.querySelectorAll(".c-sum").forEach(function (inp) {
            inp.addEventListener("change", function (ev) {
                ST.sum = ev.target.checked;
                syncLegends();
                redraw();
            });
        });
        el.querySelectorAll(".c-log").forEach(function (inp) {
            inp.addEventListener("change", function (ev) {
                ST.log = ev.target.checked;
                syncLegends();
                redraw();
            });
        });
    }
    function syncLegends() {
        ["legendM1", "legendM2"].forEach(function (id) {
            var el = document.getElementById(PFX + id);
            if (!el || !el.dataset.built)
                return;
            el.querySelectorAll("label.chip[data-nuc]").forEach(function (lab) {
                var k = lab.getAttribute("data-nuc");
                lab.querySelector("input").checked = !!ST.on[k];
                var nf = noiseFrac(k), nuc = nucOf(k);
                if (nf === null || !nuc)
                    return;
                lab.querySelector(".sw").style.background = nuc.color;
                lab.title = guardHint(nf);
            });
            el.querySelectorAll(".c-sum").forEach(function (i) { i.checked = ST.sum; });
            el.querySelectorAll(".c-log").forEach(function (i) { i.checked = ST.log; });
        });
    }
    function cell(lab, val, big, hint) {
        return "<div" + (hint ? " title='" + hint + "'" : "") + "><span class='lab'>"
            + lab + "</span><span class='val"
            + (big ? " big-num" : "") + "'>" + val + "</span></div>";
    }
    var CONT_LAB = "multiplier for the background scaled to live time (continuum correction)";
    var CONT_HINT = "coefficient for the second, non-nuclide term in the fit "
        + "(background scaled to live time); a value markedly greater than one indicates a correction "
        + "for the continuum, not a multiple of the actual background; details are given under "
        + "“how it is calculated”.";
    function fillSummaries() {
        var m1 = M1(), m2 = M2(), pass = D.passport;
        var s1 = document.getElementById(PFX + "sumM1");
        if (s1) {
            s1.innerHTML =
                cell("K-40 activity (± with Birge ratio correction)", cnt(m1.A_Bq) + " Bq <em>± "
                    + cnt(m1.dA_Bq) + " Bq</em>", true)
                    + cell("ratio to the activity expected from the sample mass", num(m1.A_Bq / pass.A_Bq, 3) + " ("
                        + signedPct(m1.A_Bq / pass.A_Bq) + ")")
                    + cell("χ²/ν", num(m1.chi2_ndof, 2))
                    + cell("channels in the fit", cnt(m1.ndof))
                    + cell("fit window", num(m1.E_fit_lo, 0) + "–" + num(m1.E_fit_hi, 0)
                        + " keV")
                    + cell(CONT_LAB, num(m1.bg_amplitude, 2), false, CONT_HINT);
        }
        var s2 = document.getElementById(PFX + "sumM2");
        if (s2) {
            s2.innerHTML =
                cell("K-40 activity (± with Birge ratio correction)", cnt(m2.A_Bq) + " Bq <em>± "
                    + cnt(m2.dA_Bq) + " Bq</em>", true)
                    + cell("ratio to the activity expected from the sample mass", num(m2.A_Bq / pass.A_Bq, 3) + " ("
                        + signedPct(m2.A_Bq / pass.A_Bq) + ")")
                    + cell("χ²/ν", num(m2.chi2_ndof, 2))
                    + cell("lines in the model", cnt(m2.n_lines) + " (1 γ-ray line and " + cnt(m2.n_lines - 1) + " Ar X-ray lines; no sum peaks)")
                    + cell("channels in the fit", cnt(m2.n_channels_fit))
                    + cell(CONT_LAB, num(m2.bg_amplitude, 2), false, CONT_HINT);
        }
    }
    function buildM1() {
        var tbl = document.getElementById(PFX + "tblM1");
        if (!tbl)
            return;
        var m1 = M1();
        var stk = STACK1();
        var head = "<thead><tr><th>nuclide</th><th class='num'>amplitude, Bq</th>"
            + "<th class='num'>ratio to the activity expected from the sample mass</th><th class='num'>fraction in the spectrum</th>"
            + "<th>note</th></tr></thead>";
        var body = "<tbody>";
        var grand = 0;
        D.nuclides.forEach(function (nuc) {
            var arr = stk[nuc.key];
            if (!arr)
                return;
            for (var i = 0; i < arr.length; i++)
                grand += arr[i];
        });
        D.nuclides.forEach(function (nuc) {
            var arr = stk[nuc.key];
            if (!arr)
                return;
            var sum = 0;
            for (var i = 0; i < arr.length; i++)
                sum += arr[i];
            var amp, damp, tag;
            var v = m1.per_nuclide[nuc.key];
            if (!v)
                return;
            amp = v.A_Bq;
            damp = v.dA_Bq;
            tag = "";
            body += "<tr>"
                + "<td><span class='sw' style='background:" + nuc.color + "'></span>"
                + esc(nuc.label_ru) + "</td>"
                + "<td class='num'>" + cnt(amp) + " ± " + cnt(damp) + "</td>"
                + "<td class='num'>" + num(amp / D.passport.A_Bq, 3) + "</td>"
                + "<td class='num'>" + num(100 * sum / Math.max(grand, 1e-9), 1)
                + " %</td>"
                + "<td>" + esc(nuc.note) + tag + "</td></tr>";
        });
        tbl.innerHTML = head + body + "</tbody>";
    }
    function buildM2Chan() {
        var tbl = document.getElementById(PFX + "tblM2Chan");
        if (!tbl || !D.channels)
            return;
        var sc = STACK2CHAN();
        if (!sc) {
            tbl.innerHTML = "";
            return;
        }
        var rows = D.channels.map(function (ch) {
            var arr = sc[ch.key] || [];
            var sum = 0;
            for (var i = 0; i < arr.length; i++)
                sum += arr[i];
            return { ch: ch, sum: sum };
        });
        var total = rows.reduce(function (a, r) { return a + r.sum; }, 0);
        rows.sort(function (a, b) { return b.sum - a.sum; });
        var head = "<thead><tr><th>channel</th><th class='num'>fraction</th>"
            + "<th class='num'>count</th></tr></thead>";
        var body = "<tbody>";
        rows.forEach(function (r) {
            var pct = total > 0 ? 100 * r.sum / total : 0;
            body += "<tr>"
                + "<td><span class='sw' style='background:" + r.ch.color + "'></span>"
                + esc(r.ch.label_ru) + "</td>"
                + "<td class='num'>" + num(pct, 1) + " %</td>"
                + "<td class='num'>" + cnt(r.sum) + "</td></tr>";
        });
        tbl.innerHTML = head + body + "</tbody>";
    }
    function buildM2() {
        var tbl = document.getElementById(PFX + "tblM2");
        if (!tbl)
            return;
        var nucCol = {}, nucLab = {}, nucIdx = {};
        D.nuclides.forEach(function (n, i) {
            nucCol[n.key] = n.color;
            nucLab[n.key] = n.label_ru;
            nucIdx[n.key] = i;
        });
        function predicted(r) { return r.predicted_net; }
        var sorters = {
            contrib: function (a, b) { return predicted(b) - predicted(a); },
            nuclide: function (a, b) {
                var d = (nucIdx[a.nuclide] || 0) - (nucIdx[b.nuclide] || 0);
                return d !== 0 ? d : a.E_keV - b.E_keV;
            },
            energy: function (a, b) { return a.E_keV - b.E_keV; },
        };
        var rows = M2().lines.slice().sort(sorters[ST.m2sort] || sorters.contrib);
        var head = "<thead><tr>"
            + "<th>line</th><th>nuclide</th>"
            + "<th class='num'>emission probability I<sub>γ</sub> per decay of the nuclide</th>"
            + "<th class='num'>ε<sub>FEP</sub></th>"
            + "<th class='num'>model counts</th><th>note</th></tr></thead>";
        var body = "<tbody>";
        for (var i = 0; i < rows.length; i++) {
            var r = rows[i];
            var lineTxt, iTxt;
            if (r.kind === "sum") {
                lineTxt = num(r.E1_keV, 1) + "+" + num(r.E2_keV, 1) + " = "
                    + num(r.E_keV, 1) + " keV";
                iTxt = num(r.I1_pct, 2) + " % × " + num(r.I2_pct, 2) + " %";
            }
            else if (r.kind === "xray") {
                lineTxt = "K-series, centroid " + num(r.E_keV, 1) + " keV";
                iTxt = num(r.I_gamma_pct, 1) + " % per decay of the nuclide";
            }
            else {
                lineTxt = num(r.E_keV, 1) + " keV";
                iTxt = num(r.I_gamma_pct, r.I_gamma_pct < 0.1 ? 4 : 2) + " %";
                if (typeof r.branch === "number" && r.branch < 0.999)
                    iTxt += " <em>× " + num(100 * r.branch, 2) + " % (branching)</em>";
            }
            var tag = r.kind === "sum" ? " · sum peak"
                : (r.kind === "xray" ? " · X-ray" : "");
            var epsTxt = (typeof r.eps_peak === "number")
                ? r.eps_peak.toExponential(3).replace(/^-/, "\u2212").replace(/e([+-])(\d+)$/, function (m, s, p) { return "\u00B710" + (s === "-" ? "\u207B" : "") + p.replace(/\d/g, function (c) { return "\u2070\u00B9\u00B2\u00B3\u2074\u2075\u2076\u2077\u2078\u2079"[c]; }); }) : "—";
            body += "<tr>"
                + "<td><span class='sw' style='background:"
                + (nucCol[r.nuclide] || "#999") + "'></span>" + lineTxt + "</td>"
                + "<td>" + esc(nucLab[r.nuclide] || r.nuclide) + tag + "</td>"
                + "<td class='num'>" + iTxt + "</td>"
                + "<td class='num'>" + epsTxt + "</td>"
                + "<td class='num'>" + cnt(predicted(r)) + "</td>"
                + "<td>" + esc(r.note || "") + "</td></tr>";
        }
        tbl.innerHTML = head + body + "</tbody>";
    }
    function cursorText(elId, stk, tipId) {
        var el = document.getElementById(PFX + elId);
        if (!el)
            return;
        if (ST.cursorE === null) {
            el.textContent = "hover over a spectrum channel";
            return;
        }
        var e = D.spectrum.e_of_ch;
        var i = 0, best = Infinity;
        for (var k = 0; k < e.length; k++) {
            var dd = Math.abs(e[k] - ST.cursorE);
            if (dd < best) {
                best = dd;
                i = k;
            }
        }
        var meas = D.spectrum.counts[i];
        var contribs = [];
        D.nuclides.forEach(function (nuc) {
            if (!ST.on[nuc.key] || !stk[nuc.key])
                return;
            var v = stk[nuc.key][i];
            if (v > 0.5)
                contribs.push({ nuc: nuc, v: v });
        });
        contribs.sort(function (a, b) { return b.v - a.v; });
        var top = contribs.slice(0, 4).map(function (c) {
            return c.nuc.label_ru + " " + cnt(c.v);
        }).join(" · ");
        var txt = num(e[i], 0) + " keV — measured " + cnt(meas)
            + ", model " + cnt(stackTotal(stk, i));
        if (top)
            txt += " — " + top;
        el.textContent = txt;
        var tip = document.getElementById(PFX + tipId);
        if (tip)
            tip.textContent = num(e[i], 0) + " keV · " + cnt(meas);
    }
    function attachCursor(cvId, tipId, onMove) {
        var cv = document.getElementById(PFX + cvId);
        var tip = document.getElementById(PFX + tipId);
        if (!cv)
            return;
        cv.addEventListener("pointermove", function (ev) {
            var r = cv.getBoundingClientRect();
            var x = ev.clientX - r.left, y = ev.clientY - r.top;
            var m = { l: 62, r: 14 };
            var e = D.spectrum.e_of_ch;
            var xLo = ST.zoom ? ST.zoom.xLo : 25, xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];
            if (x < m.l || x > r.width - m.r)
                ST.cursorE = null;
            else
                ST.cursorE = xLo + ((x - m.l) / (r.width - m.r - m.l)) * (xHi - xLo);
            onMove();
            if (tip) {
                if (ST.cursorE === null)
                    tip.hidden = true;
                else {
                    tip.hidden = false;
                    tip.style.left = x + "px";
                    tip.style.top = Math.max(0, y) + "px";
                }
            }
        });
        cv.addEventListener("pointerleave", function () {
            ST.cursorE = null;
            onMove();
            if (tip)
                tip.hidden = true;
        });
    }
    function zEfromX(x, rectWidth) {
        var e = D.spectrum.e_of_ch;
        var xLo = ST.zoom ? ST.zoom.xLo : 25;
        var xHi = ST.zoom ? ST.zoom.xHi : e[e.length - 1];
        var m = { l: 62, r: 14 };
        return xLo + ((x - m.l) / (rectWidth - m.r - m.l)) * (xHi - xLo);
    }
    function wireZoom(cvId) {
        var cv = document.getElementById(PFX + cvId);
        if (!cv)
            return;
        var dragging = false;
        var dragState = null;
        var m = { l: 62, r: 14 };
        cv.addEventListener("mousedown", function (ev) {
            var r = cv.getBoundingClientRect();
            var x = ev.clientX - r.left;
            if (x < m.l || x > r.width - m.r)
                return;
            ev.preventDefault();
            dragging = true;
            dragState = { x0: x, x1: x, w: r.width };
            ST.drag = { e0: zEfromX(x, r.width), e1: zEfromX(x, r.width) };
            redraw();
        });
        document.addEventListener("mousemove", function (ev) {
            if (!dragging)
                return;
            var r = cv.getBoundingClientRect();
            var x = ev.clientX - r.left;
            if (x < m.l)
                x = m.l;
            if (x > r.width - m.r)
                x = r.width - m.r;
            dragState.x1 = x;
            ST.drag.e1 = zEfromX(x, r.width);
            redraw();
        });
        document.addEventListener("mouseup", function () {
            if (!dragging)
                return;
            dragging = false;
            var x0 = dragState.x0, x1 = dragState.x1, w = dragState.w;
            dragState = null;
            ST.drag = null;
            if (Math.abs(x1 - x0) < 6) {
                redraw();
                return;
            }
            ST.zoom = { xLo: Math.max(0, zEfromX(Math.min(x0, x1), w)),
                xHi: zEfromX(Math.max(x0, x1), w) };
            redraw();
        });
        cv.addEventListener("dblclick", function () {
            ST.zoom = null;
            redraw();
        });
    }
    function redraw() {
        var vis = document.querySelector("#src-k40 .stage:not([hidden])");
        if (!vis)
            return;
        if (vis.id === PFX + "viewM1") {
            drawSpectrum("cvM1", STACK1(), []);
            cursorText("cursorM1", STACK1(), "m1-tip");
        }
        else if (vis.id === PFX + "viewM2") {
            drawSpectrum("cvM2", STACK2(), []);
            cursorText("cursorM2", STACK2(), "m2-tip");
        }
        else if (vis.id === PFX + "viewCmp") {
            drawCmp();
            fillCmpTable();
        }
        else if (vis.id === PFX + "viewCal") {
            drawCal();
            drawFwhm();
        }
    }
    function refreshAll() {
        syncLegends();
        fillSummaries();
        buildM1();
        buildM2();
        buildM2Chan();
        fillCmpTable();
        redraw();
    }
    function openPop(id) {
        var tpl = document.getElementById(PFX + id);
        if (!tpl)
            return;
        var pop = document.getElementById(PFX + "pop");
        var scrim = document.getElementById(PFX + "scrim");
        pop.innerHTML = "";
        pop.appendChild(tpl.content.cloneNode(true));
        var close = document.createElement("button");
        close.type = "button";
        close.className = "pop-close";
        close.setAttribute("aria-label", "close");
        close.textContent = "×";
        close.addEventListener("click", closePop);
        pop.appendChild(close);
        pop.hidden = false;
        scrim.hidden = false;
    }
    function closePop() {
        document.getElementById(PFX + "pop").hidden = true;
        document.getElementById(PFX + "scrim").hidden = true;
    }
    function switchTab(name) {
        document.querySelectorAll("#k40-tabs .tab").forEach(function (t) {
            t.setAttribute("aria-selected", t.getAttribute("data-tab") === name ? "true" : "false");
        });
        var map = { m1: "viewM1", m2: "viewM2", cmp: "viewCmp", cal: "viewCal" };
        Object.keys(map).forEach(function (k) {
            var el = document.getElementById(PFX + map[k]);
            if (el)
                el.hidden = (k !== name);
        });
        if (name === "m1") {
            buildLegend("legendM1");
            fillSummaries();
            buildM1();
        }
        if (name === "m2") {
            buildLegend("legendM2");
            fillSummaries();
            buildM2();
            buildM2Chan();
        }
        if (name === "cal")
            buildCal();
        redraw();
    }
    function buildCal() {
        var tbl = document.getElementById(PFX + "tblCal");
        if (!tbl || tbl.dataset.built)
            return;
        tbl.dataset.built = "1";
        var m = D.meta;
        function coefsHtml(coefs) {
            return coefs.map(function (c, i) {
                var abs = Math.abs(c), s;
                if (abs === 0)
                    s = "0";
                else if (abs >= 0.01 && abs < 10000)
                    s = num(c, 6);
                else
                    s = c.toExponential(4).replace(/^-/, "\u2212").replace(/e([+-])(\d+)$/, function (m, s, p) { return "\u00B710" + (s === "-" ? "\u207B" : "") + p.replace(/\d/g, function (c) { return "\u2070\u00B9\u00B2\u00B3\u2074\u2075\u2076\u2077\u2078\u2079"[c]; }); });
                return "<span class='mono'>c" + i + " = " + s + "</span>";
            }).join("<br>");
        }
        var head = "<thead><tr><th>parameter</th>"
            + "<th>sample (KCl)</th><th>Marinelli + water background</th></tr></thead>";
        var body = "<tbody>"
            + "<tr><td>channels</td><td class='num'>" + m.cal_sample.n_channels
            + "</td><td class='num'>" + m.cal_bg.n_channels + "</td></tr>"
            + "<tr><td>live time, s</td><td class='num'>" + num(m.live_s, 2)
            + "</td><td class='num'>" + num(m.bg_live_s, 2) + "</td></tr>"
            + "<tr><td>real time, s</td><td class='num'>" + num(m.real_s, 2)
            + "</td><td class='num'>" + num(m.bg_real_s, 2) + "</td></tr>"
            + "<tr><td>dead time, %</td><td class='num'>"
            + num(100 * (m.real_s - m.live_s) / m.real_s, 3) + "</td><td class='num'>"
            + num(100 * (m.bg_real_s - m.bg_live_s) / m.bg_real_s, 3) + "</td></tr>"
            + "<tr><td>degree of the E(channel) polynomial</td><td class='num'>"
            + m.cal_sample.order + "</td><td class='num'>" + m.cal_bg.order + "</td></tr>"
            + "<tr><td>polynomial representing this spectrum’s own energy calibration (from the shapes of its peak groups)</td><td>" + coefsHtml(m.cal_sample.coefs)
            + "</td><td>" + coefsHtml(m.cal_bg.coefs) + "</td></tr>"
            + "<tr><td>background scaling factor (t_sample / t_bg)</td>"
            + "<td class='num' colspan='2'>" + num(m.bg_scale_time, 4) + "</td></tr>"
            + "</tbody>";
        tbl.innerHTML = head + body;
        buildFwhmTable();
    }
    var CAL = { smp: true, bg: true, diff: false, log: true, anch: true,
        zoom: null, drag: null, dragging: false };
    var CAL_wired = false;
    var CAL_MARGIN = { l: 62, r: 14, t: 12, b: 32 };
    var CAL_MARK_H = 7;
    var CAL_MARK_HIT_PX = 5;
    function drawCal() {
        var cv = document.getElementById(PFX + "cvCal");
        if (!cv)
            return;
        if (!CAL_wired)
            wireCal();
        var p = pal();
        var f = fit(cv);
        var g = f.g, W = f.w, H = f.h;
        var m = CAL_MARGIN;
        var e = D.spectrum.e_of_ch;
        var y = D.spectrum.counts;
        var b = D.spectrum.bg_counts;
        var xLo = CAL.zoom ? CAL.zoom.xLo : 0;
        var xHi = CAL.zoom ? CAL.zoom.xHi : e[e.length - 1];
        var vMax = 1, series = [];
        if (CAL.smp)
            series.push(y);
        if (CAL.bg)
            series.push(b);
        if (CAL.diff) {
            var diff = new Array(e.length);
            for (var i = 0; i < e.length; i++)
                diff[i] = y[i] - b[i];
            series.push(diff);
        }
        series.forEach(function (arr) {
            for (var i = 0; i < arr.length; i++) {
                if (e[i] < xLo || e[i] > xHi)
                    continue;
                var v = CAL.log ? Math.abs(arr[i]) : arr[i];
                if (v > vMax)
                    vMax = v;
            }
        });
        var Y = makeY(CAL.log, CAL.log ? 0.5 : 0, vMax * (CAL.log ? 2 : 1.1));
        g.strokeStyle = p.grid;
        g.lineWidth = 1;
        g.beginPath();
        var xTicks = [250, 500, 750, 1000, 1250, 1500, 1750, 2000,
            2250, 2500, 2750, 3000];
        for (var xi = 0; xi < xTicks.length; xi++) {
            if (xTicks[xi] > xHi)
                break;
            var xx = mapX(xTicks[xi], xLo, xHi, m.l, W - m.r);
            g.moveTo(xx, m.t);
            g.lineTo(xx, H - m.b);
        }
        var yTicks = [];
        if (CAL.log) {
            for (var d = 0; d <= Math.ceil(Math.log10(Y.hi)); d++)
                yTicks.push(Math.pow(10, d));
        }
        else {
            var stp = Math.pow(10, Math.floor(Math.log10(Y.hi / 4)));
            var s0 = Math.ceil(Y.hi / 4 / stp) * stp;
            for (var v = s0; v < Y.hi; v += s0)
                yTicks.push(v);
        }
        for (var yi = 0; yi < yTicks.length; yi++) {
            var yg = Y.map(yTicks[yi], m.t, H - m.b);
            g.moveTo(m.l, yg);
            g.lineTo(W - m.r, yg);
        }
        g.stroke();
        g.fillStyle = p.faint;
        g.font = "11px system-ui, sans-serif";
        g.textAlign = "center";
        g.textBaseline = "top";
        for (var xj = 0; xj < xTicks.length; xj++) {
            if (xTicks[xj] > xHi)
                break;
            g.fillText(String(xTicks[xj]), mapX(xTicks[xj], xLo, xHi, m.l, W - m.r), H - m.b + 4);
        }
        g.textAlign = "right";
        g.textBaseline = "middle";
        for (var yj = 0; yj < yTicks.length; yj++) {
            var label = CAL.log ? "10" + supNum(Math.round(Math.log10(yTicks[yj])))
                : cnt(yTicks[yj]);
            g.fillText(label, m.l - 4, Y.map(yTicks[yj], m.t, H - m.b));
        }
        g.textAlign = "center";
        g.textBaseline = "bottom";
        g.fillText("energy, keV", (m.l + W - m.r) / 2, H - 2);
        g.strokeStyle = p.rule;
        g.lineWidth = 2;
        g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);
        if (CAL.anch && D.reference_lines) {
            var nucCol = {};
            D.nuclides.forEach(function (n) { nucCol[n.label_ru] = n.color; });
            D.reference_lines.forEach(function (r) {
                var E = r[0], nuc = r[1];
                if (E < xLo || E > xHi)
                    return;
                var xa = mapX(E, xLo, xHi, m.l, W - m.r);
                var col = nucCol[nuc] || p.faint;
                g.strokeStyle = col;
                g.lineWidth = 1;
                g.setLineDash([3, 3]);
                g.beginPath();
                g.moveTo(xa, m.t + CAL_MARK_H + 2);
                g.lineTo(xa, H - m.b);
                g.stroke();
                g.setLineDash([]);
                g.fillStyle = col;
                g.beginPath();
                g.moveTo(xa - 4, m.t);
                g.lineTo(xa + 4, m.t);
                g.lineTo(xa, m.t + CAL_MARK_H);
                g.closePath();
                g.fill();
                g.strokeStyle = p.paper;
                g.lineWidth = 1;
                g.stroke();
            });
        }
        function drawTrace(arr, color, allowNeg) {
            g.strokeStyle = color;
            g.lineWidth = 1.2;
            g.beginPath();
            var started = false;
            for (var k = 0; k < e.length; k++) {
                if (e[k] < xLo || e[k] > xHi)
                    continue;
                var v = arr[k];
                if (!allowNeg && v < 0)
                    v = 0;
                var vv = CAL.log ? Math.max(v, Y.lo) : v;
                var xt = mapX(e[k], xLo, xHi, m.l, W - m.r);
                var yt = Y.map(vv, m.t, H - m.b);
                if (!started) {
                    g.moveTo(xt, yt);
                    started = true;
                }
                else
                    g.lineTo(xt, yt);
            }
            g.stroke();
        }
        if (CAL.bg)
            drawTrace(b, "#0f5aa8", false);
        if (CAL.diff)
            drawTrace((function () {
                var dd = new Array(e.length);
                for (var i = 0; i < e.length; i++)
                    dd[i] = y[i] - b[i];
                return dd;
            })(), "#c8541c", CAL.log ? false : true);
        if (CAL.smp)
            drawTrace(y, p.ink, false);
        if (CAL.drag) {
            var xa2 = Math.min(CAL.drag.x0, CAL.drag.x1);
            var xb2 = Math.max(CAL.drag.x0, CAL.drag.x1);
            g.fillStyle = "rgba(246,211,28,.22)";
            g.fillRect(xa2, m.t, xb2 - xa2, H - m.b - m.t);
            g.strokeStyle = "#16140f";
            g.lineWidth = 1.5;
            g.setLineDash([4, 4]);
            g.strokeRect(xa2, m.t, xb2 - xa2, H - m.b - m.t);
            g.setLineDash([]);
        }
    }
    function calEfromX(x, rectWidth) {
        var e = D.spectrum.e_of_ch;
        var xLo = CAL.zoom ? CAL.zoom.xLo : 0;
        var xHi = CAL.zoom ? CAL.zoom.xHi : e[e.length - 1];
        var m = CAL_MARGIN;
        return xLo + ((x - m.l) / (rectWidth - m.r - m.l)) * (xHi - xLo);
    }
    function calRefLineAt(x, y, rectWidth) {
        if (!CAL.anch || !D.reference_lines)
            return null;
        if (y < CAL_MARGIN.t - 2 || y > CAL_MARGIN.t + CAL_MARK_H + 3)
            return null;
        var e = D.spectrum.e_of_ch;
        var xLo = CAL.zoom ? CAL.zoom.xLo : 0;
        var xHi = CAL.zoom ? CAL.zoom.xHi : e[e.length - 1];
        var best = null, bestD = CAL_MARK_HIT_PX;
        D.reference_lines.forEach(function (r) {
            var E = r[0];
            if (E < xLo || E > xHi)
                return;
            var xa = mapX(E, xLo, xHi, CAL_MARGIN.l, rectWidth - CAL_MARGIN.r);
            var d = Math.abs(xa - x);
            if (d <= bestD) {
                bestD = d;
                best = r;
            }
        });
        return best;
    }
    function wireCal() {
        ["c-smp", "c-bg", "c-diff", "c-log", "c-anch"].forEach(function (id) {
            var el = document.getElementById(PFX + id);
            if (!el)
                return;
            var key = { "c-smp": "smp", "c-bg": "bg", "c-diff": "diff",
                "c-log": "log", "c-anch": "anch" }[id];
            el.addEventListener("change", function (ev) {
                CAL[key] = ev.target.checked;
                drawCal();
            });
        });
        var calReset = document.getElementById(PFX + "cal-reset");
        if (calReset)
            calReset.addEventListener("click", function () {
                CAL.zoom = null;
                drawCal();
            });
        var cv = document.getElementById(PFX + "cvCal");
        var ro = document.getElementById(PFX + "cal-ro");
        var tip = document.getElementById(PFX + "cal-tip");
        if (cv && ro) {
            cv.addEventListener("pointermove", function (ev) {
                var r = cv.getBoundingClientRect();
                var x = ev.clientX - r.left, y = ev.clientY - r.top;
                if (CAL.dragging) {
                    CAL.drag.x1 = Math.max(CAL_MARGIN.l, Math.min(r.width - CAL_MARGIN.r, x));
                    drawCal();
                    return;
                }
                if (x < CAL_MARGIN.l || x > r.width - CAL_MARGIN.r) {
                    ro.textContent = "";
                    if (tip)
                        tip.hidden = true;
                    return;
                }
                var ref = calRefLineAt(x, y, r.width);
                if (ref) {
                    var refTxt = ref[1] + " · " + num(ref[0], 1) + " keV";
                    ro.textContent = refTxt;
                    if (tip) {
                        tip.hidden = false;
                        tip.textContent = refTxt;
                        tip.style.left = x + "px";
                        tip.style.top = Math.max(0, y) + "px";
                    }
                    return;
                }
                var e = D.spectrum.e_of_ch;
                var E = calEfromX(x, r.width);
                var i = 0, best = Infinity;
                for (var k = 0; k < e.length; k++) {
                    var dd = Math.abs(e[k] - E);
                    if (dd < best) {
                        best = dd;
                        i = k;
                    }
                }
                var smp = D.spectrum.counts[i], bgv = D.spectrum.bg_counts[i];
                ro.textContent = num(e[i], 0) + " keV · sample " + cnt(smp)
                    + " · background " + cnt(bgv) + " · difference " + cnt(smp - bgv);
                if (tip) {
                    tip.hidden = false;
                    tip.textContent = num(e[i], 0) + " keV · " + cnt(smp);
                    tip.style.left = x + "px";
                    tip.style.top = Math.max(0, y) + "px";
                }
            });
            cv.addEventListener("pointerleave", function () {
                if (!CAL.dragging) {
                    ro.textContent = "";
                    if (tip)
                        tip.hidden = true;
                }
            });
            cv.addEventListener("mousedown", function (ev) {
                var r = cv.getBoundingClientRect();
                var x = ev.clientX - r.left;
                if (x < CAL_MARGIN.l || x > r.width - CAL_MARGIN.r)
                    return;
                ev.preventDefault();
                CAL.dragging = true;
                CAL.drag = { x0: x, x1: x };
                drawCal();
            });
            document.addEventListener("mouseup", function () {
                if (!CAL.dragging)
                    return;
                CAL.dragging = false;
                var r = cv.getBoundingClientRect();
                var x0 = CAL.drag.x0, x1 = CAL.drag.x1;
                CAL.drag = null;
                if (Math.abs(x1 - x0) < 6) {
                    drawCal();
                    return;
                }
                CAL.zoom = { xLo: Math.max(0, calEfromX(Math.min(x0, x1), r.width)),
                    xHi: calEfromX(Math.max(x0, x1), r.width) };
                drawCal();
            });
            cv.addEventListener("dblclick", function () {
                CAL.zoom = null;
                drawCal();
            });
        }
        CAL_wired = true;
    }
    function fwSrc(q) {
        var c = function (v, u) {
            return "<td class='num'>" + (v == null ? "—" : num(v, 2) + (u ? " ± " + num(u, 2) : "")) + "</td>";
        };
        return c(q.own_keV, q.own_unc_keV) + c(q.becqmoni_keV) + c(q.fwhm_sl_keV);
    }
    function buildFwhmTable() {
        var tbl = document.getElementById(PFX + "tblFwhm");
        if (!tbl || !D.fwhm_cal)
            return;
        var fw = D.fwhm_cal;
        var head = "<thead><tr><th>line, keV</th><th class='num'>centroid</th>"
            + "<th class='num'>FWHM, keV</th><th class='num'>resolution</th>"
            + "<th class='num'>own (measured in this spectrum)</th><th class='num'>BecqMoni</th><th class='num'>SpectraLine</th><th class='num'>power law k·E<sup>p</sup></th>"
            + "<th class='num'>deviation</th><th>status</th></tr></thead>";
        var body = "<tbody>";
        fw.points.forEach(function (q) {
            if (!q.used) {
                body += "<tr class='row-dirty'><td>" + num(q.E_nominal, 1) + "</td>"
                    + "<td class='num'>—</td><td class='num'>—</td><td class='num'>—</td>" + fwSrc(q)
                    + "<td class='num'>—</td><td class='num'>—</td>"
                    + "<td>not used: " + esc(q.reject) + "</td></tr>";
                return;
            }
            body += "<tr><td>" + num(q.E_nominal, 1) + "</td>"
                + "<td class='num'>" + num(q.E_centroid, 1) + "</td>"
                + "<td class='num'>" + num(q.fwhm_keV, 2) + " ± "
                + num(q.d_fwhm_keV, 2) + "</td>"
                + "<td class='num'>" + num(q.res_pct, 2) + " %</td>"
                + fwSrc(q)
                + "<td class='num'>" + num(q.fwhm_model_keV, 2) + "</td>"
                + "<td class='num'>" + (q.dev_pct >= 0 ? "+" : "−")
                + num(Math.abs(q.dev_pct), 1) + " %</td>"
                + "<td>in the fit: " + (q.own_used ? "own width" : "SpectraLine × multiplier") + "</td></tr>";
        });
        body += "<tr class='sum'><td>power law</td>"
            + "<td class='num' colspan='2'>FWHM = " + num(fw.k, 3) + "·E<sup>"
            + num(fw.p, 4) + "</sup></td>"
            + "<td class='num'>" + num(fw.res662_pct, 2) + " % at 662 keV</td>"
            + "<td class='num' colspan='3'>" + fw.n_used + " of " + fw.n_anchors + " points</td>"
            + "<td class='num'>" + num(fw.fwhm662_law, 1) + " keV</td>"
            + "<td class='num'>RMS " + num(fw.rms_dev_pct, 1) + " %</td>"
            + "<td>FWHM values used in the convolution: thorium reference spectrum (SpectraLine × multiplier)</td></tr>";
        tbl.innerHTML = head + body + "</tbody>";
    }
    function drawFwhm() {
        var cv = document.getElementById(PFX + "cvFwhm");
        if (!cv || !D.fwhm_cal)
            return;
        var fw = D.fwhm_cal;
        var p = pal();
        var f = fit(cv);
        var g = f.g, W = f.w, H = f.h;
        var m = { l: 62, r: 16, t: 14, b: 34 };
        var used = fw.points.filter(function (q) { return q.used; });
        if (!used.length)
            return;
        var xLo = 0, xHi = 2900;
        var vMax = 0;
        used.forEach(function (q) { vMax = Math.max(vMax, q.fwhm_keV); });
        vMax = Math.max(vMax, fw.k * Math.pow(xHi, fw.p)) * 1.15;
        g.strokeStyle = p.grid;
        g.lineWidth = 1;
        g.beginPath();
        var xTicks = [500, 1000, 1500, 2000, 2500];
        xTicks.forEach(function (t) {
            var x = mapX(t, xLo, xHi, m.l, W - m.r);
            g.moveTo(x, m.t);
            g.lineTo(x, H - m.b);
        });
        var yTicks = [25, 50, 75, 100, 125];
        yTicks.forEach(function (t) {
            if (t > vMax)
                return;
            var y = m.t + (1 - t / vMax) * (H - m.b - m.t);
            g.moveTo(m.l, y);
            g.lineTo(W - m.r, y);
        });
        g.stroke();
        g.fillStyle = p.faint;
        g.font = "11px system-ui, sans-serif";
        g.textAlign = "center";
        g.textBaseline = "top";
        xTicks.forEach(function (t) {
            g.fillText(String(t), mapX(t, xLo, xHi, m.l, W - m.r), H - m.b + 4);
        });
        g.textAlign = "right";
        g.textBaseline = "middle";
        yTicks.forEach(function (t) {
            if (t > vMax)
                return;
            g.fillText(String(t), m.l - 4, m.t + (1 - t / vMax) * (H - m.b - m.t));
        });
        g.textAlign = "center";
        g.textBaseline = "bottom";
        g.fillText("energy, keV", (m.l + W - m.r) / 2, H - 2);
        g.save();
        g.translate(12, (m.t + H - m.b) / 2);
        g.rotate(-Math.PI / 2);
        g.textAlign = "center";
        g.textBaseline = "top";
        g.fillText("FWHM, keV", 0, 0);
        g.restore();
        g.strokeStyle = p.rule;
        g.lineWidth = 2;
        g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);
        function yOf(v) { return m.t + (1 - v / vMax) * (H - m.b - m.t); }
        g.strokeStyle = "#0f5aa8";
        g.lineWidth = 2;
        g.beginPath();
        for (var E = 40; E <= xHi; E += 10) {
            var x = mapX(E, xLo, xHi, m.l, W - m.r);
            var y = yOf(fw.k * Math.pow(E, fw.p));
            if (E === 40)
                g.moveTo(x, y);
            else
                g.lineTo(x, y);
        }
        g.stroke();
        g.strokeStyle = p.faint;
        g.lineWidth = 1.5;
        g.setLineDash([5, 4]);
        g.beginPath();
        for (var E2 = 40; E2 <= xHi; E2 += 10) {
            var x2 = mapX(E2, xLo, xHi, m.l, W - m.r);
            var y2 = yOf(fw.fwhm662_cs * Math.sqrt(E2 / 661.657));
            if (E2 === 40)
                g.moveTo(x2, y2);
            else
                g.lineTo(x2, y2);
        }
        g.stroke();
        g.setLineDash([]);
        g.fillStyle = "#c8541c";
        g.strokeStyle = "#c8541c";
        g.lineWidth = 1.5;
        used.forEach(function (q) {
            var x = mapX(q.E_centroid, xLo, xHi, m.l, W - m.r);
            var y = yOf(q.fwhm_keV);
            g.beginPath();
            g.moveTo(x, yOf(q.fwhm_keV - q.d_fwhm_keV));
            g.lineTo(x, yOf(q.fwhm_keV + q.d_fwhm_keV));
            g.stroke();
            g.beginPath();
            g.arc(x, y, 4, 0, 2 * Math.PI);
            g.fill();
        });
        g.font = "600 11px system-ui, sans-serif";
        g.textAlign = "left";
        g.textBaseline = "top";
        g.fillStyle = "#c8541c";
        g.fillText("FWHM values used in the model convolution (from SpectraLine)", m.l + 10, m.t + 8);
        g.fillStyle = "#0f5aa8";
        g.fillText("FWHM(E) power law", m.l + 10, m.t + 24);
        g.fillStyle = p.faint;
        g.fillText("power law with exponent 0.5 through the 662 keV point", m.l + 10, m.t + 40);
    }
    function cmpItems() {
        var p = pal();
        return [
            { lab: "expected", A: D.passport.A_Bq, dA: D.passport.dA_Bq, col: p.ink },
            { lab: "Method 1", A: M1().A_Bq, dA: M1().dA_Bq, col: "#0f5aa8" },
            { lab: "Method 2", A: SRC().method2.A_Bq, dA: SRC().method2.dA_Bq, col: "#c8541c" }
        ];
    }
    function fillCmpTable() {
        var el = document.getElementById(PFX + "cmpTable");
        if (!el)
            return;
        var pass = D.passport, m1 = M1(), m2 = SRC().method2, m2f = SRC().method2_full;
        var modeTxt = ST.fwhmLaw === "cs" ? "background: Marinelli + water"
            : "background as acquired (without the beaker)";
        var massKg = pass.mass_g / 1000;
        function row(cls, lab, A, dA, note) {
            return "<div class='cmp-row " + cls + "'>"
                + "<span class='cmp-lab'>" + lab + "</span>"
                + "<span class='cmp-val big-num'>" + cnt(A) + " Bq <em>± "
                + cnt(dA) + " Bq</em> <em>· " + cnt(A / massKg) + " Bq/kg ± "
                + cnt(dA / massKg) + " Bq/kg</em></span>"
                + "<span class='cmp-note'>" + note + "</span></div>";
        }
        el.innerHTML =
            row("cmp-pass", "expected (from mass)", pass.A_Bq, pass.dA_Bq, "calculation from KCl mass " + num(pass.mass_g, 0) + " g: potassium fraction × K-40 abundance × ln2/T½; "
                + "uncertainty " + num(pass.unc_pct, 2) + " % — K-40 abundance, half-life, reagent purity")
                + row("cmp-m1", "Method 1: K-40 MC template, " + num(m1.E_fit_lo, 0) + "–"
                    + num(m1.E_fit_hi, 0) + " keV", m1.A_Bq, m1.dA_Bq, num(m1.A_Bq / pass.A_Bq, 3) + " of the expected value, "
                    + signedPct(m1.A_Bq / pass.A_Bq) + "; χ²/ν = " + num(m1.chi2_ndof, 2)
                    + " over " + cnt(m1.ndof) + " channels; " + modeTxt)
                + row("cmp-m2", "Method 2: FEP efficiency function + 1460.8 keV γ-ray line", m2.A_Bq, m2.dA_Bq, num(m2.A_Bq / pass.A_Bq, 3) + " of the expected value, "
                    + signedPct(m2.A_Bq / pass.A_Bq) + "; χ²/ν = " + num(m2.chi2_ndof, 2)
                    + " over " + cnt(m2.n_channels_fit) + " channels; " + modeTxt)
                + "<div class='cmp-row'><span class='cmp-lab'>discrepancy between methods</span>"
                + "<span class='cmp-val big-num'>" + signedPct(m1.A_Bq / m2.A_Bq)
                + "</span><span class='cmp-note'>Method 1 relative to Method 2</span></div>";
    }
    function drawCmp() {
        var cv = document.getElementById(PFX + "cvCmp");
        if (!cv)
            return;
        var p = pal();
        var f = fit(cv);
        var g = f.g, W = f.w, H = f.h;
        var m = { l: 26, r: 20, t: 18, b: 34 };
        var items = cmpItems();
        var lo = Infinity, hi = -Infinity;
        items.forEach(function (it) {
            lo = Math.min(lo, it.A - it.dA);
            hi = Math.max(hi, it.A + it.dA);
        });
        var pad = (hi - lo) * 0.15 + 1;
        lo -= pad;
        hi += pad;
        lo = Math.max(0, lo);
        g.strokeStyle = p.rule;
        g.lineWidth = 2;
        g.strokeRect(m.l, m.t, W - m.r - m.l, H - m.b - m.t);
        var range = hi - lo;
        var stp = Math.pow(10, Math.floor(Math.log10(range / 4)));
        var s = Math.max(stp, Math.ceil(range / 5 / stp) * stp);
        g.strokeStyle = p.grid;
        g.beginPath();
        g.fillStyle = p.faint;
        g.font = "11px system-ui, sans-serif";
        g.textAlign = "center";
        g.textBaseline = "top";
        for (var v = Math.ceil(lo / s) * s; v <= hi; v += s) {
            var x = mapX(v, lo, hi, m.l, W - m.r);
            g.moveTo(x, m.t);
            g.lineTo(x, H - m.b);
            g.fillText(cnt(v), x, H - m.b + 4);
        }
        g.stroke();
        g.textAlign = "center";
        g.textBaseline = "bottom";
        g.fillText("activity, Bq", (m.l + W - m.r) / 2, H - 2);
        var innerH = H - m.b - m.t;
        var rowH = innerH / items.length;
        for (var j = 0; j < items.length; j++) {
            var it = items[j];
            var yc = m.t + rowH * (j + 0.5);
            var xl = mapX(it.A - it.dA, lo, hi, m.l, W - m.r);
            var xr = mapX(it.A + it.dA, lo, hi, m.l, W - m.r);
            var xm = mapX(it.A, lo, hi, m.l, W - m.r);
            g.fillStyle = it.col;
            g.globalAlpha = 0.28;
            g.fillRect(xl, yc - rowH * 0.28, Math.max(xr - xl, 1), rowH * 0.56);
            g.globalAlpha = 1;
            g.strokeStyle = it.col;
            g.lineWidth = 3;
            g.beginPath();
            g.moveTo(xm, yc - rowH * 0.36);
            g.lineTo(xm, yc + rowH * 0.36);
            g.stroke();
            g.fillStyle = it.col;
            g.textAlign = "left";
            g.textBaseline = "middle";
            g.font = "bold 13px system-ui, sans-serif";
            g.fillText(it.lab, m.l + 6, yc - rowH * 0.28);
            g.fillStyle = p.ink;
            g.font = "12px ui-monospace, Menlo, monospace";
            g.fillText(cnt(it.A) + " ± " + cnt(it.dA) + " Bq", Math.min(xr + 8, W - m.r - 130), yc);
        }
    }
    function wire() {
        document.querySelectorAll("#k40-tabs .tab").forEach(function (t) {
            t.addEventListener("click", function () {
                switchTab(t.getAttribute("data-tab"));
            });
        });
        document.querySelectorAll("#k40-fwhmseg .btn").forEach(function (b) {
            b.addEventListener("click", function () {
                ST.fwhmLaw = b.getAttribute("data-fwhmlaw");
                document.querySelectorAll("#k40-fwhmseg .btn").forEach(function (o) {
                    o.setAttribute("aria-pressed", o.getAttribute("data-fwhmlaw") === ST.fwhmLaw ? "true" : "false");
                });
                refreshAll();
            });
        });
        document.querySelectorAll("#k40-libseg .btn").forEach(function (b) {
            b.addEventListener("click", function () {
                ST.lib = b.getAttribute("data-lib");
                document.querySelectorAll("#k40-libseg .btn").forEach(function (o) {
                    o.setAttribute("aria-pressed", o.getAttribute("data-lib") === ST.lib ? "true" : "false");
                });
                refreshAll();
            });
        });
        document.querySelectorAll("#k40-m2sortseg .btn").forEach(function (b) {
            b.addEventListener("click", function () {
                ST.m2sort = b.getAttribute("data-sort");
                document.querySelectorAll("#k40-m2sortseg .btn").forEach(function (o) {
                    o.setAttribute("aria-pressed", o.getAttribute("data-sort") === ST.m2sort ? "true" : "false");
                });
                buildM2();
            });
        });
        document.querySelectorAll("[data-k40pop]").forEach(function (b) {
            b.addEventListener("click", function () {
                openPop(b.getAttribute("data-k40pop"));
            });
        });
        var scrim = document.getElementById(PFX + "scrim");
        if (scrim)
            scrim.addEventListener("click", closePop);
        document.addEventListener("keydown", function (ev) {
            if (ev.key === "Escape")
                closePop();
        });
        window.addEventListener("resize", redraw);
        wireZoom("cvM1");
        wireZoom("cvM2");
        attachCursor("cvM1", "m1-tip", function () {
            cursorText("cursorM1", STACK1(), "m1-tip");
            drawSpectrum("cvM1", STACK1(), []);
        });
        attachCursor("cvM2", "m2-tip", function () {
            cursorText("cursorM2", STACK2(), "m2-tip");
            drawSpectrum("cvM2", STACK2(), []);
        });
    }
    document.addEventListener("DOMContentLoaded", function () {
        wire();
        switchTab("m1");
        var host = document.querySelector("#src-k40");
        if (host && window.MutationObserver)
            new MutationObserver(function () {
                if (!host.hidden)
                    redraw();
            }).observe(host, { attributes: true, attributeFilter: ["hidden"] });
    });
})();
