function yRange(arrays, j0, j1) {
    let ymax = 0;
    let ymin = Infinity;
    for (let i = 0; i < arrays.length; i++) {
        const arr = arrays[i];
        for (let j = j0; j <= j1; j++) {
            const v = arr[j];
            if (v > 0) {
                ymax = Math.max(ymax, v);
                ymin = Math.min(ymin, v);
            }
        }
    }
    return { ymax, ymin };
}

function drawStackLayers(ctx, st, o) {
    const { j0, j1, E0, dE, X, Y, color } = o;
    for (let i = 0; i < st.length; i++) {
        const l = st[i];
        ctx.fillStyle = color(l.id);
        ctx.globalAlpha = 0.9;
        ctx.beginPath();
        // Верхняя граница
        for (let j = j0; j <= j1; j++) {
            if (j === j0) ctx.moveTo(X(E0 + j * dE), Y(l.upper[j]));
            else ctx.lineTo(X(E0 + j * dE), Y(l.upper[j]));
        }
        // Нижняя граница
        for (let j = j1; j >= j0; j--) {
            ctx.lineTo(X(E0 + j * dE), Y(l.lower[j]));
        }
        ctx.closePath();
        ctx.fill();
    }
    ctx.globalAlpha = 1;
}

function legendRows(series, st, o) {
    const { dE, hidden, colorOf, nameOf } = o;
    let sum = 0;
    for (let i = 0; i < st.length; i++) {
        sum += st[i].contrib;
    }
    if (sum === 0) sum = 1;

    const rows = [];
    for (let i = 0; i < st.length; i++) {
        const l = st[i];
        const s = series.find(s => s.id === l.id);
        if (!s) continue;
        const counts = l.contrib * dE;
        const share = 100 * l.contrib / sum;
        rows.push({
            id: l.id,
            name: nameOf(s),
            color: colorOf(l.id),
            counts,
            share,
            on: !hidden[l.id]
        });
    }

    rows.reverse();   // st по возрастанию вклада -> строки по убыванию
    return rows;
}

function legendHtml(rows, totalOn, totalColor) {
    let html = `<div class="series-group">`;
    html += `<label class="series-row"><span class="color-swatch" style="background:${totalColor}"></span>`;
    html += `<input type="checkbox" data-id="total" ${totalOn ? "checked" : ""}> `;
    html += `Полный спектр (сумма слоёв)</label></div>`;
    html += `<div class="series-group"><h3>Слои по убыванию вклада в показанном диапазоне</h3>`;
    for (let i = 0; i < rows.length; i++) {
        const r = rows[i];
        const name = r.name.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
        const counts = Math.round(r.counts);
        const share = r.share.toFixed(1).replace(".", ",");
        html += `<label class="series-row"><span class="color-swatch" style="background:${r.color}"></span>`;
        html += `<input type="checkbox" data-id="${r.id}" ${r.on ? "checked" : ""}> `;
        html += `${name} — ${counts} отсч./ч, ${share} %</label><br>`;
    }
    html += `</div>`;
    return html;
}

// Самотест
if (typeof require !== "undefined" && require.main === module) {
    const s = require("./stack.js");   // мутация: отключить reverse — тест порядка должен краснеть
    global.buildStack = s.buildStack;
    global.russianOrigin = s.russianOrigin;

    // Тест yRange
    let res = yRange([[0,2,5],[1,0,0]], 0, 2);
    if (res.ymax !== 5 || res.ymin !== 1) {
        console.log("SELFTEST FAIL: yRange");
        process.exit(1);
    }
    res = yRange([[0,0]], 0, 1);
    if (res.ymax !== 0 || res.ymin !== Infinity) {
        console.log("SELFTEST FAIL: yRange2");
        process.exit(1);
    }

    // Тест legendRows
    const st = buildStack([
        {id: "a", y: [1,1]},
        {id: "b", y: [5,5]},
        {id: "c", y: [2,2]}
    ], 0, 1);
    const series = [
        {id: "a", label: "a", group: "component"},
        {id: "b", label: "b", group: "component"},
        {id: "c", label: "c", group: "component"}
    ];
    const rows = legendRows(series, st, {
        dE: 1,
        hidden: {c: true},
        colorOf: id => "#" + id,
        nameOf: s => s.label.toUpperCase()
    });
    if (rows[0].id !== "b" || rows[1].id !== "c" || rows[2].id !== "a") {
        console.log("SELFTEST FAIL: legendRows order");
        process.exit(1);
    }
    if (Math.abs(rows[0].share - 62.5) > 1e-9) {
        console.log("SELFTEST FAIL: legendRows share");
        process.exit(1);
    }
    if (rows[1].on !== false) {
        console.log("SELFTEST FAIL: legendRows on");
        process.exit(1);
    }
    if (rows[0].name !== "B") {
        console.log("SELFTEST FAIL: legendRows name");
        process.exit(1);
    }

    // Тест legendHtml
    const html = legendHtml(rows, true, "#000");
    if (!html.includes('data-id="total" checked')) {
        console.log("SELFTEST FAIL: legendHtml totalOn");
        process.exit(1);
    }
    if (!html.includes("B — 10 отсч./ч, 62,5 %")) {
        console.log("SELFTEST FAIL: legendHtml content");
        process.exit(1);
    }

    // Тест drawStackLayers
    const ctx = {
        fillCount: 0,
        moveToCount: 0,
        lineToCount: 0,
        globalAlpha: 0.9,
        beginPath() {},
        closePath() {},
        fill() { this.fillCount++; },
        moveTo(x, y) { this.moveToCount++; },
        lineTo(x, y) { this.lineToCount++; }
    };
    const fakeSt = [
        {id: "a", contrib: 1, lower: new Float64Array([1,2,3]), upper: new Float64Array([1,2,3])},
        {id: "b", contrib: 2, lower: new Float64Array([4,5,6]), upper: new Float64Array([4,5,6])}
    ];
    drawStackLayers(ctx, fakeSt, {
        j0: 0,
        j1: 2,
        E0: 0,
        dE: 1,
        X: x => x,
        Y: y => y,
        color: id => "#000"
    });
    if (ctx.fillCount !== 2 || ctx.moveToCount !== 2 || ctx.lineToCount !== 10) {
        console.log("SELFTEST FAIL: drawStackLayers calls");
        process.exit(1);
    }
    if (ctx.globalAlpha !== 1) {
        console.log("SELFTEST FAIL: drawStackLayers globalAlpha");
        process.exit(1);
    }

    console.log("SELFTEST PASS");
    process.exit(0);
}

if (typeof module !== "undefined") module.exports = { yRange, drawStackLayers, legendRows, legendHtml };
