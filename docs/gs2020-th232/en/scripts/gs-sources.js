(function () {
    "use strict";
    const $ = (sel, ctx = document) => ctx.querySelector(sel);
    const $$ = (sel, ctx = document) => Array.from(ctx.querySelectorAll(sel));
    function fmt(x, d) {
        if (x == null || isNaN(x))
            return "—";
        let s = x.toFixed(d);
        let parts = s.split(".");
        let intPart = parts[0].length > 4 ? parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, "\u202F") : parts[0];
        return intPart + (parts[1] !== undefined ? "." + parts[1] : "");
    }
    function logTicks(yMin, yMax) {
        const sup = { "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
        const out = [];
        for (let p = Math.ceil(Math.log10(yMin)); p <= Math.floor(Math.log10(yMax)); p++) {
            out.push({ v: Math.pow(10, p), label: p <= 3 ? String(Math.pow(10, p)) : "10" + String(p).split("").map(c => sup[c]).join("") });
        }
        return out;
    }
    function escHtml(str) {
        if (!str)
            return "";
        return String(str).replace(/[&<>"']/g, function (m) {
            return ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;' })[m];
        });
    }
    function niceStep(range, targetTicks = 5) {
        if (range <= 0)
            return 1;
        let rough = range / targetTicks;
        let pow = Math.pow(10, Math.floor(Math.log10(rough)));
        let frac = rough / pow;
        let nice;
        if (frac <= 1.5)
            nice = 1;
        else if (frac <= 3.5)
            nice = 2;
        else if (frac <= 7.5)
            nice = 5;
        else
            nice = 10;
        return nice * pow;
    }
    function lowerBound(arr, val) {
        let lo = 0, hi = arr.length;
        while (lo < hi) {
            let mid = (lo + hi) >> 1;
            if (arr[mid] < val)
                lo = mid + 1;
            else
                hi = mid;
        }
        return lo;
    }
    function getCssVar(name, fallback) {
        let v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
        return v || fallback;
    }
    const GS_SRC = window.GS_SRC;
    if (!GS_SRC)
        return;
    const order = GS_SRC.order || [];
    const sources = GS_SRC.sources || {};
    const srcbar = $("#srcbar");
    if (!srcbar)
        return;
    const buttons = {};
    order.forEach(key => {
        const src = sources[key];
        if (!src)
            return;
        const btn = document.createElement("button");
        btn.type = "button";
        btn.className = "btn";
        btn.setAttribute("role", "tab");
        btn.setAttribute("data-src", key);
        let label = src.label || key;
        if (src.status === "soon") {
            label += " (to be added)";
            btn.disabled = true;
            btn.setAttribute("aria-disabled", "true");
            btn.title = "spectrum to be added";
        }
        btn.textContent = label;
        if (src.status !== "soon")
            btn.addEventListener("click", () => select(key));
        srcbar.appendChild(btn);
        buttons[key] = btn;
    });
    function select(key) {
        if (!buttons[key])
            return;
        order.forEach(k => {
            const b = buttons[k];
            if (b)
                b.setAttribute("aria-selected", k === key ? "true" : "false");
        });
        order.forEach(k => {
            const el = document.getElementById("src-" + k);
            if (el) {
                if (k === key) {
                    el.hidden = false;
                    if (k === "th") {
                        window.dispatchEvent(new Event("resize"));
                    }
                    else {
                        const block = blocks[k];
                        if (block)
                            block.draw();
                    }
                }
                else {
                    el.hidden = true;
                }
            }
        });
        history.replaceState(null, "", "#src=" + key);
    }
    function initSelection() {
        let hashKey = null;
        const hashMatch = location.hash.match(/^#src=(.+)$/);
        if (hashMatch) {
            const k = hashMatch[1];
            if (order.includes(k) && sources[k] && sources[k].status === "ready") {
                hashKey = k;
            }
        }
        const startKey = hashKey || order[0];
        select(startKey);
    }
    window.addEventListener("load", initSelection);
    const blocks = {};
    order.forEach(key => {
        const src = sources[key];
        if (!src || !src.spectrum)
            return;
        const container = document.getElementById("src-" + key);
        if (!container)
            return;
        const resTable = $('[data-role="res"]', container);
        if (resTable && src.results) {
            let html = `<thead><tr>
            <th>estimate</th>
            <th class="num">activity, Bq</th>
            <th class="num">ratio to expected</th>
            <th>note</th>
        </tr></thead><tbody>`;
            src.results.forEach(r => {
                let actStr = fmt(r.A, 1);
                if (r.dA != null)
                    actStr += " ± " + fmt(r.dA, 1);
                let ratioStr = r.ratio ? fmt(r.ratio, 4) : "—";
                html += `<tr>
                <td>${escHtml(r.lab)}</td>
                <td class="num">${actStr}</td>
                <td class="num">${ratioStr}</td>
                <td>${escHtml(r.note)}</td>
            </tr>`;
            });
            html += `</tbody>`;
            resTable.innerHTML = html;
        }
        const chartBox = $('[data-role="chart"]', container);
        if (!chartBox)
            return;
        const spec = src.spectrum;
        const ctrlsDiv = document.createElement("div");
        ctrlsDiv.className = "ctrls";
        const seriesConfig = [
            { id: "meas", label: "measurement", colorVar: "--ink", checked: true },
            { id: "bg", label: "Marinelli + water background (scaled)", color: "#0f5aa8", checked: true },
            { id: "model", label: "model + background", color: "#c8541c", checked: true },
            { id: "diff", label: "measurement − background", color: "#6a6558", checked: false }
        ];
        (spec.comps || []).forEach((c, k) => seriesConfig.push({ id: "comp" + k, label: c.label, color: c.color, checked: true }));
        const checkboxes = {};
        seriesConfig.forEach(cfg => {
            const lbl = document.createElement("label");
            const sw = document.createElement("span");
            sw.className = "sw";
            if (cfg.colorVar) {
                sw.style.background = getCssVar(cfg.colorVar, "#16140f");
            }
            else {
                sw.style.background = cfg.color;
            }
            const inp = document.createElement("input");
            inp.type = "checkbox";
            inp.checked = cfg.checked;
            inp.addEventListener("change", () => block.draw());
            lbl.appendChild(sw);
            lbl.appendChild(inp);
            lbl.appendChild(document.createTextNode(cfg.label));
            ctrlsDiv.appendChild(lbl);
            checkboxes[cfg.id] = inp;
        });
        const logLbl = document.createElement("label");
        const logInp = document.createElement("input");
        logInp.type = "checkbox";
        logInp.checked = true;
        logInp.addEventListener("change", () => block.draw());
        checkboxes.log = logInp;
        logLbl.appendChild(logInp);
        logLbl.appendChild(document.createTextNode("logarithmic scale"));
        ctrlsDiv.appendChild(logLbl);
        const resetBtn = document.createElement("button");
        resetBtn.type = "button";
        resetBtn.className = "btn";
        resetBtn.textContent = "full range";
        resetBtn.addEventListener("click", () => { block.zoom = null; block.draw(); });
        ctrlsDiv.appendChild(resetBtn);
        chartBox.appendChild(ctrlsDiv);
        const cvMain = document.createElement("canvas");
        cvMain.className = "src-spec";
        cvMain.setAttribute("role", "img");
        cvMain.setAttribute("aria-label", "spectrum: measurement, background, model; drag to zoom, double-click to reset");
        const cvResid = document.createElement("canvas");
        cvResid.className = "src-resid";
        cvResid.setAttribute("role", "img");
        cvResid.setAttribute("aria-label", "normalized residuals");
        chartBox.appendChild(cvMain);
        chartBox.appendChild(cvResid);
        const readout = document.createElement("div");
        readout.className = "readout";
        chartBox.appendChild(readout);
        const block = {
            key: key,
            spec: spec,
            cvMain: cvMain,
            cvResid: cvResid,
            readout: readout,
            checkboxes: checkboxes,
            zoom: null,
            drag: null,
            cursorE: null,
            draw() {
                if (container.hidden)
                    return;
                const w = this.cvMain.clientWidth;
                if (w === 0)
                    return;
                this.setupCanvas(this.cvMain);
                this.setupCanvas(this.cvResid);
                const ctxM = this.cvMain.getContext("2d");
                const ctxR = this.cvResid.getContext("2d");
                let xRange = this.zoom ? [this.zoom.lo, this.zoom.hi] : spec.view;
                const iStart = Math.max(0, lowerBound(spec.e, xRange[0]) - 1);
                const iEnd = Math.min(spec.e.length - 1, lowerBound(spec.e, xRange[1]));
                const ink = getCssVar("--ink", "#16140f");
                const dim = getCssVar("--dim", "#4c493f");
                const ruleSoft = getCssVar("--rule-soft", "#c9c4b4");
                const paper = getCssVar("--paper", "#f5f2ea");
                const padM = { l: 64, r: 14, t: 12, b: 34 };
                const padR = { l: 64, r: 14, t: 6, b: 30 };
                const plotW = w - padM.l - padM.r;
                const plotHm = this.cvMain.clientHeight - padM.t - padM.b;
                const plotHr = this.cvResid.clientHeight - padR.t - padR.b;
                [ctxM, ctxR].forEach(ctx => {
                    ctx.clearRect(0, 0, w, ctx.canvas.clientHeight);
                    ctx.fillStyle = paper;
                    ctx.fillRect(0, 0, w, ctx.canvas.clientHeight);
                });
                const xToPx = (e) => padM.l + (e - xRange[0]) / (xRange[1] - xRange[0]) * plotW;
                const pxToE = (px) => xRange[0] + (px - padM.l) / plotW * (xRange[1] - xRange[0]);
                ctxM.fillStyle = ruleSoft;
                ctxM.globalAlpha = 0.35;
                const clampPx = (px) => Math.max(padM.l, Math.min(padM.l + plotW, px));
                if (spec.lo > xRange[0]) {
                    const pxLo = clampPx(xToPx(spec.lo));
                    ctxM.fillRect(padM.l, padM.t, pxLo - padM.l, plotHm);
                }
                if (spec.hi < xRange[1]) {
                    const pxHi = clampPx(xToPx(spec.hi));
                    ctxM.fillRect(pxHi, padM.t, padM.l + plotW - pxHi, plotHm);
                }
                ctxM.globalAlpha = 1;
                let yMin = Infinity, yMax = -Infinity;
                const isLog = this.checkboxes.log.checked;
                const seriesData = [];
                const cb = this.checkboxes;
                let yPos = Infinity;
                if (cb.meas.checked)
                    seriesData.push({ color: ink, width: 1.2, calc: (i) => spec.meas[i] });
                if (cb.bg.checked)
                    seriesData.push({ color: "#0f5aa8", width: 1.2, calc: (i) => spec.bg[i] });
                if (cb.model.checked)
                    seriesData.push({ color: "#c8541c", width: 1.6, calc: (i) => spec.model[i] + spec.bg[i] });
                if (cb.diff.checked)
                    seriesData.push({ color: "#6a6558", width: 1.2, calc: (i) => spec.meas[i] - spec.bg[i] });
                (spec.comps || []).forEach((c, k) => { const b = cb["comp" + k]; if (b && b.checked)
                    seriesData.push({ color: c.color, width: 1.2, calc: (i) => c.data[i] }); });
                seriesData.forEach(s => {
                    for (let i = iStart; i <= iEnd; i++) {
                        const v = s.calc(i);
                        if (v < yMin)
                            yMin = v;
                        if (v > yMax)
                            yMax = v;
                        if (v > 0 && v < yPos)
                            yPos = v;
                    }
                });
                if (!isFinite(yMax)) {
                    yMin = 0;
                    yMax = 1;
                }
                if (isLog) {
                    yMin = Math.max(0.5, 0.7 * (isFinite(yPos) ? yPos : 1));
                    yMax = Math.max(1.5 * yMax, 10 * yMin);
                }
                else {
                    yMin = Math.min(0, yMin);
                    yMax = Math.max(yMax, yMin) + (yMax > yMin ? 0.08 * (yMax - yMin) : 1);
                }
                const yToPx = (v) => {
                    if (isLog) {
                        if (v <= yMin)
                            v = yMin;
                        const logMin = Math.log10(yMin);
                        const logMax = Math.log10(yMax);
                        const logV = Math.log10(v);
                        return padM.t + plotHm - (logV - logMin) / (logMax - logMin) * plotHm;
                    }
                    else {
                        return padM.t + plotHm - (v - yMin) / (yMax - yMin) * plotHm;
                    }
                };
                ctxM.save();
                ctxM.rect(padM.l, padM.t, plotW, plotHm);
                ctxM.clip();
                const xStep = niceStep(xRange[1] - xRange[0], plotW / 90);
                const xStartTick = Math.ceil(xRange[0] / xStep) * xStep;
                ctxM.strokeStyle = ruleSoft;
                ctxM.lineWidth = 1;
                ctxM.beginPath();
                for (let e = xStartTick; e <= xRange[1]; e += xStep) {
                    const px = xToPx(e);
                    ctxM.moveTo(px, padM.t);
                    ctxM.lineTo(px, padM.t + plotHm);
                }
                ctxM.stroke();
                if (isLog) {
                    logTicks(yMin, yMax).forEach(t => {
                        const py = yToPx(t.v);
                        ctxM.moveTo(padM.l, py);
                        ctxM.lineTo(padM.l + plotW, py);
                    });
                }
                else {
                    const yStep = niceStep(yMax - yMin, 5);
                    for (let v = Math.ceil(yMin / yStep) * yStep; v <= yMax; v += yStep) {
                        const py = yToPx(v);
                        ctxM.moveTo(padM.l, py);
                        ctxM.lineTo(padM.l + plotW, py);
                    }
                }
                ctxM.stroke();
                if (spec.anchors) {
                    ctxM.strokeStyle = dim;
                    ctxM.setLineDash([4, 4]);
                    spec.anchors.forEach(([eVal, label]) => {
                        if (eVal >= xRange[0] && eVal <= xRange[1]) {
                            const px = xToPx(eVal);
                            ctxM.beginPath();
                            ctxM.moveTo(px, padM.t);
                            ctxM.lineTo(px, padM.t + plotHm);
                            ctxM.stroke();
                            ctxM.fillStyle = dim;
                            ctxM.font = "12px ui-monospace, Consolas, monospace";
                            ctxM.textAlign = "center";
                            ctxM.fillText(label, px, padM.t + 12);
                        }
                    });
                    ctxM.setLineDash([]);
                }
                seriesData.forEach(s => {
                    ctxM.strokeStyle = s.color;
                    ctxM.lineWidth = s.width;
                    ctxM.beginPath();
                    let started = false;
                    for (let i = iStart; i <= iEnd; i++) {
                        const e = spec.e[i];
                        const v = s.calc ? s.calc(i) : s.data[i];
                        let drawV = v;
                        if (isLog && v <= yMin)
                            drawV = yMin;
                        const px = xToPx(e);
                        const py = yToPx(drawV);
                        if (!started) {
                            ctxM.moveTo(px, py);
                            started = true;
                        }
                        else {
                            ctxM.lineTo(px, py);
                        }
                    }
                    ctxM.stroke();
                });
                if (this.cursorE != null && this.cursorE >= xRange[0] && this.cursorE <= xRange[1]) {
                    const px = xToPx(this.cursorE);
                    ctxM.strokeStyle = dim;
                    ctxM.setLineDash([2, 2]);
                    ctxM.beginPath();
                    ctxM.moveTo(px, padM.t);
                    ctxM.lineTo(px, padM.t + plotHm);
                    ctxM.stroke();
                    ctxM.setLineDash([]);
                }
                if (this.drag && this.drag.active) {
                    const x0 = Math.max(padM.l, Math.min(padM.l + plotW, this.drag.x0));
                    const x1 = Math.max(padM.l, Math.min(padM.l + plotW, this.drag.x1));
                    ctxM.fillStyle = "rgba(246, 211, 28, .22)";
                    ctxM.fillRect(Math.min(x0, x1), padM.t, Math.abs(x1 - x0), plotHm);
                    ctxM.strokeStyle = "#16140f";
                    ctxM.setLineDash([4, 4]);
                    ctxM.strokeRect(Math.min(x0, x1), padM.t, Math.abs(x1 - x0), plotHm);
                    ctxM.setLineDash([]);
                }
                ctxM.restore();
                ctxM.strokeStyle = ink;
                ctxM.lineWidth = 1.5;
                ctxM.strokeRect(padM.l, padM.t, plotW, plotHm);
                ctxM.fillStyle = ink;
                ctxM.font = "12px ui-monospace, Consolas, monospace";
                ctxM.textAlign = "center";
                for (let e = xStartTick; e <= xRange[1]; e += xStep) {
                    const px = xToPx(e);
                    ctxM.fillText(fmt(e, 0), px, padM.t + plotHm + 16);
                }
                ctxM.fillText("energy, keV", padM.l + plotW / 2, padM.t + plotHm + 30);
                ctxM.textAlign = "right";
                if (isLog) {
                    logTicks(yMin, yMax).forEach(t => ctxM.fillText(t.label, padM.l - 6, yToPx(t.v) + 4));
                }
                else {
                    const yStep = niceStep(yMax - yMin, 5);
                    for (let v = Math.ceil(yMin / yStep) * yStep; v <= yMax; v += yStep) {
                        const py = yToPx(v);
                        ctxM.fillText(fmt(v, 0), padM.l - 6, py + 4);
                    }
                }
                ctxM.save();
                ctxM.translate(14, padM.t + plotHm / 2);
                ctxM.rotate(-Math.PI / 2);
                ctxM.textAlign = "center";
                ctxM.fillText("counts per channel", 0, 0);
                ctxM.restore();
                const padRl = padR.l;
                const plotWr = w - padR.l - padR.r;
                const plotHrVal = this.cvResid.clientHeight - padR.t - padR.b;
                ctxR.save();
                ctxR.rect(padRl, padR.t, plotWr, plotHrVal);
                ctxR.clip();
                ctxR.fillStyle = ruleSoft;
                ctxR.globalAlpha = 0.35;
                if (spec.lo > xRange[0]) {
                    const pxLo = clampPx(xToPx(spec.lo));
                    ctxR.fillRect(padRl, padR.t, pxLo - padRl, plotHrVal);
                }
                if (spec.hi < xRange[1]) {
                    const pxHi = clampPx(xToPx(spec.hi));
                    ctxR.fillRect(pxHi, padR.t, padRl + plotWr - pxHi, plotHrVal);
                }
                ctxR.globalAlpha = 1;
                ctxR.strokeStyle = ruleSoft;
                ctxR.lineWidth = 1;
                ctxR.beginPath();
                for (let e = xStartTick; e <= xRange[1]; e += xStep) {
                    const px = xToPx(e);
                    ctxR.moveTo(px, padR.t);
                    ctxR.lineTo(px, padR.t + plotHrVal);
                }
                ctxR.stroke();
                const yResidToPx = (v) => {
                    let cv = Math.max(-5, Math.min(5, v));
                    return padR.t + plotHrVal - (cv + 5) / 10 * plotHrVal;
                };
                ctxR.strokeStyle = dim;
                ctxR.beginPath();
                ctxR.moveTo(padRl, yResidToPx(0));
                ctxR.lineTo(padRl + plotWr, yResidToPx(0));
                ctxR.stroke();
                ctxR.strokeStyle = ruleSoft;
                ctxR.setLineDash([4, 4]);
                ctxR.beginPath();
                ctxR.moveTo(padRl, yResidToPx(2));
                ctxR.lineTo(padRl + plotWr, yResidToPx(2));
                ctxR.moveTo(padRl, yResidToPx(-2));
                ctxR.lineTo(padRl + plotWr, yResidToPx(-2));
                ctxR.stroke();
                ctxR.setLineDash([]);
                ctxR.fillStyle = ink;
                for (let i = iStart; i <= iEnd; i++) {
                    if (spec.resid[i] == null)
                        continue;
                    const px = xToPx(spec.e[i]);
                    const py = yResidToPx(spec.resid[i]);
                    ctxR.fillRect(px - 1, py - 1, 2, 2);
                }
                if (this.cursorE != null && this.cursorE >= xRange[0] && this.cursorE <= xRange[1]) {
                    const px = xToPx(this.cursorE);
                    ctxR.strokeStyle = dim;
                    ctxR.setLineDash([2, 2]);
                    ctxR.beginPath();
                    ctxR.moveTo(px, padR.t);
                    ctxR.lineTo(px, padR.t + plotHrVal);
                    ctxR.stroke();
                    ctxR.setLineDash([]);
                }
                if (this.drag && this.drag.active) {
                    const x0 = Math.max(padRl, Math.min(padRl + plotWr, this.drag.x0));
                    const x1 = Math.max(padRl, Math.min(padRl + plotWr, this.drag.x1));
                    ctxR.fillStyle = "rgba(246, 211, 28, .22)";
                    ctxR.fillRect(Math.min(x0, x1), padR.t, Math.abs(x1 - x0), plotHrVal);
                    ctxR.strokeStyle = "#16140f";
                    ctxR.setLineDash([4, 4]);
                    ctxR.strokeRect(Math.min(x0, x1), padR.t, Math.abs(x1 - x0), plotHrVal);
                    ctxR.setLineDash([]);
                }
                ctxR.restore();
                ctxR.strokeStyle = ink;
                ctxR.lineWidth = 1.5;
                ctxR.strokeRect(padRl, padR.t, plotWr, plotHrVal);
                ctxR.fillStyle = ink;
                ctxR.font = "12px ui-monospace, Consolas, monospace";
                ctxR.textAlign = "center";
                for (let e = xStartTick; e <= xRange[1]; e += xStep) {
                    const px = xToPx(e);
                    ctxR.fillText(fmt(e, 0), px, padR.t + plotHrVal + 16);
                }
                ctxR.fillText("energy, keV", padRl + plotWr / 2, padR.t + plotHrVal + 30);
                ctxR.textAlign = "right";
                [-5, -2, 0, 2, 5].forEach(v => {
                    const py = yResidToPx(v);
                    ctxR.fillText(String(v), padRl - 6, py + 4);
                });
                ctxR.save();
                ctxR.translate(14, padR.t + plotHrVal / 2);
                ctxR.rotate(-Math.PI / 2);
                ctxR.textAlign = "center";
                ctxR.fillText("(measured − model)/σ", 0, 0);
                ctxR.restore();
                this.updateReadout(pxToE);
            },
            setupCanvas(cv) {
                const dpr = window.devicePixelRatio || 1;
                const w = cv.clientWidth;
                const h = cv.clientHeight;
                if (w === 0 || h === 0)
                    return;
                cv.width = Math.round(w * dpr);
                cv.height = Math.round(h * dpr);
                const ctx = cv.getContext("2d");
                ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            },
            updateReadout(pxToE) {
                if (this.cursorE == null) {
                    this.readout.textContent = "hover over the spectrum; drag to zoom, double-click for the full range";
                    return;
                }
                const idx = lowerBound(spec.e, this.cursorE);
                let bestIdx = idx;
                if (idx > 0) {
                    const distCurr = Math.abs(spec.e[idx] - this.cursorE);
                    const distPrev = Math.abs(spec.e[idx - 1] - this.cursorE);
                    if (distPrev < distCurr)
                        bestIdx = idx - 1;
                }
                if (bestIdx >= spec.e.length)
                    bestIdx = spec.e.length - 1;
                if (bestIdx < 0)
                    bestIdx = 0;
                const chNum = spec.ch0 + bestIdx;
                const eVal = spec.e[bestIdx];
                const measVal = spec.meas[bestIdx];
                const bgVal = spec.bg[bestIdx];
                const modelVal = spec.model[bestIdx] + spec.bg[bestIdx];
                const residVal = spec.resid[bestIdx];
                this.readout.textContent =
                    `channel ${chNum} · ${fmt(eVal, 1)} keV · measured ${fmt(measVal, 1)} · background ${fmt(bgVal, 1)} · model + background ${fmt(modelVal, 1)} · residual ${fmt(residVal, 2)}`;
            }
        };
        blocks[key] = block;
        const handlePointerDown = (e) => {
            if (e.button !== 0)
                return;
            e.preventDefault();
            const rect = e.target.getBoundingClientRect();
            const x = e.clientX - rect.left;
            if (x < 64 || x > rect.width - 14)
                return;
            block.drag = { active: true, x0: x, x1: x, cv: e.target };
        };
        const handlePointerMove = (e) => {
            if (block.drag && block.drag.active && block.drag.cv === e.target) {
                const rect = e.target.getBoundingClientRect();
                const x = e.clientX - rect.left;
                block.drag.x1 = Math.max(64, Math.min(rect.width - 14, x));
                block.draw();
            }
            else {
                const rect = e.target.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const padL = 64;
                const plotW = e.target.clientWidth - 64 - 14;
                if (x >= padL && x <= padL + plotW) {
                    let xRange = block.zoom ? [block.zoom.lo, block.zoom.hi] : spec.view;
                    const pxToE = (px) => xRange[0] + (px - padL) / plotW * (xRange[1] - xRange[0]);
                    block.cursorE = pxToE(x);
                }
                else {
                    block.cursorE = null;
                }
                block.draw();
            }
        };
        const handlePointerLeave = () => {
            if (!block.drag || !block.drag.active) {
                block.cursorE = null;
                block.draw();
            }
        };
        [block.cvMain, block.cvResid].forEach(cv => {
            cv.addEventListener("mousedown", handlePointerDown);
            cv.addEventListener("pointermove", handlePointerMove);
            cv.addEventListener("pointerleave", handlePointerLeave);
            cv.addEventListener("dblclick", () => {
                block.zoom = null;
                block.draw();
            });
        });
        document.addEventListener("mouseup", (e) => {
            if (!block.drag || !block.drag.active)
                return;
            block.drag.active = false;
            const dx = Math.abs(block.drag.x1 - block.drag.x0);
            if (dx < 6) {
                block.draw();
                return;
            }
            let xRange = block.zoom ? [block.zoom.lo, block.zoom.hi] : spec.view;
            const plotW = block.drag.cv.clientWidth - 64 - 14;
            const pxToE = (px) => xRange[0] + (px - 64) / plotW * (xRange[1] - xRange[0]);
            let lo = Math.max(spec.e[0], Math.min(pxToE(block.drag.x0), pxToE(block.drag.x1)));
            let hi = Math.min(spec.e[spec.e.length - 1], Math.max(pxToE(block.drag.x0), pxToE(block.drag.x1)));
            if (hi - lo >= 1) {
                block.zoom = { lo, hi };
            }
            block.draw();
        });
        const resizeObs = new ResizeObserver(() => {
            requestAnimationFrame(() => block.draw());
        });
        resizeObs.observe(chartBox);
        window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
            block.draw();
        });
    });
})();
