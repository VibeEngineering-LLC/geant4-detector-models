(function () { "use strict";

// --- Утилиты ---

const $ = (sel, ctx = document) => ctx.querySelector(sel);
const $$ = (sel, ctx = document) => Array.from(ctx.querySelectorAll(sel));

// Форматирование числа: 1234.56 -> 1 234,56
function fmt(x, d) {
    if (x == null || isNaN(x)) return "—";
    // toFixed возвращает строку с точкой
    let s = x.toFixed(d);
    let parts = s.split(".");
    // Целая часть: группировка по 3 через U+202F
    let intPart = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, "\u202F");
    return intPart + (parts[1] !== undefined ? "," + parts[1] : "");
}

// Деления лог-шкалы: степени 10 внутри [yMin, yMax]; подписи 1, 10, 100, 1000, 10⁴…
function logTicks(yMin, yMax) {
    const sup = { "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
    const out = [];
    for (let p = Math.ceil(Math.log10(yMin)); p <= Math.floor(Math.log10(yMax)); p++) {
        out.push({ v: Math.pow(10, p), label: p <= 3 ? String(Math.pow(10, p)) : "10" + String(p).split("").map(c => sup[c]).join("") });
    }
    return out;
}

// Экранирование HTML для безопасной вставки текста
function escHtml(str) {
    if (!str) return "";
    return String(str).replace(/[&<>"']/g, function(m) {
        return ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;' })[m];
    });
}

// Красивый шаг для осей (1, 2, 5 * 10^n)
function niceStep(range, targetTicks = 5) {
    if (range <= 0) return 1;
    let rough = range / targetTicks;
    let pow = Math.pow(10, Math.floor(Math.log10(rough)));
    let frac = rough / pow;
    let nice;
    if (frac <= 1.5) nice = 1;
    else if (frac <= 3.5) nice = 2;
    else if (frac <= 7.5) nice = 5;
    else nice = 10;
    return nice * pow;
}

// Двоичный поиск индекса, где e[i] >= val
function lowerBound(arr, val) {
    let lo = 0, hi = arr.length;
    while (lo < hi) {
        let mid = (lo + hi) >> 1;
        if (arr[mid] < val) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

// Получение CSS переменной с фолбэком
function getCssVar(name, fallback) {
    let v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return v || fallback;
}

// --- Основная логика ---

const GS_SRC = window.GS_SRC;
if (!GS_SRC) return; // Нет данных — выходим

const order = GS_SRC.order || [];
const sources = GS_SRC.sources || {};

// 1. Переключатель (Tabs)
const srcbar = $("#srcbar");
if (!srcbar) return;

const buttons = {};
order.forEach(key => {
    const src = sources[key];
    if (!src) return;
    
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "btn";
    btn.setAttribute("role", "tab");
    btn.setAttribute("data-src", key);
    
    let label = src.label || key;
    if (src.status === "soon") {
        label += " — скоро";
        btn.disabled = true;
        btn.setAttribute("aria-disabled", "true");
        btn.title = "спектр будет добавлен";
    }
    
    btn.textContent = label;
    if (src.status !== "soon") btn.addEventListener("click", () => select(key));   // правка: в генерации клик не был привязан
    srcbar.appendChild(btn);
    buttons[key] = btn;
});

// Функция выбора источника
function select(key) {
    if (!buttons[key]) return;

    // Обновляем кнопки
    order.forEach(k => {
        const b = buttons[k];
        if (b) b.setAttribute("aria-selected", k === key ? "true" : "false");
    });

    // Скрываем/показываем контейнеры
    order.forEach(k => {
        const el = document.getElementById("src-" + k);
        if (el) {
            if (k === key) {
                el.hidden = false;
                // Если это торий — триггерим ресайз для его приложения
                if (k === "th") {
                    window.dispatchEvent(new Event("resize"));
                } else {
                    // Для остальных — запускаем отрисовку, если есть блок
                    const block = blocks[k];
                    if (block) block.draw();
                }
            } else {
                el.hidden = true;
            }
        }
    });

    // Обновляем хэш
    history.replaceState(null, "", "#src=" + key);
}

// Инициализация выбора
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

// Ждем загрузки окна, чтобы приложение тория успело инициализироваться
window.addEventListener("load", initSelection);


// 2. Блок источника (Графики)
const blocks = {}; // Хранилище состояний блоков

order.forEach(key => {
    const src = sources[key];
    if (!src || !src.spectrum) return;
    
    const container = document.getElementById("src-" + key);
    if (!container) return;

    // 2.1 Таблица результатов
    const resTable = $('[data-role="res"]', container);
    if (resTable && src.results) {
        let html = `<thead><tr>
            <th>оценка</th>
            <th class="num">активность, Бк</th>
            <th class="num">к ожидаемой</th>
            <th>примечание</th>
        </tr></thead><tbody>`;
        
        src.results.forEach(r => {
            let actStr = fmt(r.A, 1);
            if (r.dA != null) actStr += " ± " + fmt(r.dA, 1);
            
            let ratioStr = r.ratio ? fmt(r.ratio, 4) : "—";
            
            html += `<tr>
                <td>${escHtml(r.lab)}</td>
                <td class="num">${actStr}</td>
                <td class="num">${ratioStr}</td>
                <td>${escHtml(r.note)}</td>
            </tr>`;
        });
        html += `</tbody>`;
        
        // Используем innerHTML только для статической структуры таблицы, данные экранированы
        resTable.innerHTML = html;
    }

    // 2.2 Панель управления и канвы
    const chartBox = $('[data-role="chart"]', container);
    if (!chartBox) return;

    const spec = src.spectrum;
    
    // Создаем структуру внутри chartBox
    const ctrlsDiv = document.createElement("div");
    ctrlsDiv.className = "ctrls";
    
    // Флажки серий
    const seriesConfig = [
        { id: "meas", label: "измерение", colorVar: "--ink", checked: true },
        { id: "bg", label: "фон воды (приведён)", color: "#0f5aa8", checked: true },
        { id: "model", label: "модель + фон", color: "#c8541c", checked: true },
        { id: "diff", label: "измерение − фон", color: "#6a6558", checked: false }
    ];

    const checkboxes = {};
    
    seriesConfig.forEach(cfg => {
        const lbl = document.createElement("label");
        const sw = document.createElement("span");
        sw.className = "sw";
        if (cfg.colorVar) {
            // Цвет из CSS переменной, но для фона span лучше взять конкретный hex или rgba
            // В ТЗ сказано: цвет --ink. Получим его значение.
            sw.style.background = getCssVar(cfg.colorVar, "#16140f");
        } else {
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

    // Флажок лог-шкалы
    const logLbl = document.createElement("label");
    const logInp = document.createElement("input");
    logInp.type = "checkbox";
    logInp.checked = true; // По умолчанию включён
    logInp.addEventListener("change", () => block.draw());
    checkboxes.log = logInp;   // правка: draw() читает checkboxes.log
    logLbl.appendChild(logInp);
    logLbl.appendChild(document.createTextNode("лог-шкала"));
    ctrlsDiv.appendChild(logLbl);

    // Кнопка сброса
    const resetBtn = document.createElement("button");
    resetBtn.type = "button";
    resetBtn.className = "btn";
    resetBtn.textContent = "весь диапазон";
    resetBtn.addEventListener("click", () => { block.zoom = null; block.draw(); });
    ctrlsDiv.appendChild(resetBtn);

    chartBox.appendChild(ctrlsDiv);

    // Канвы
    const cvMain = document.createElement("canvas");
    cvMain.className = "src-spec";
    cvMain.setAttribute("role", "img");
    cvMain.setAttribute("aria-label", "спектр: измерение, фон, модель; протяжка мышью — приближение, двойной клик — сброс");
    
    const cvResid = document.createElement("canvas");
    cvResid.className = "src-resid";
    cvResid.setAttribute("role", "img");
    cvResid.setAttribute("aria-label", "нормированные остатки");

    chartBox.appendChild(cvMain);
    chartBox.appendChild(cvResid);

    // Readout
    const readout = document.createElement("div");
    readout.className = "readout";
    chartBox.appendChild(readout);

    // --- Логика блока (Block State & Logic) ---
    
    const block = {
        key: key,
        spec: spec,
        cvMain: cvMain,
        cvResid: cvResid,
        readout: readout,
        checkboxes: checkboxes,
        zoom: null, // {lo, hi} или null
        drag: null, // {active, x0, x1, cv}
        cursorE: null,
        
        // Методы рисования
        draw() {
            if (container.hidden) return; // Не рисуем скрытый
            
            const w = this.cvMain.clientWidth;
            if (w === 0) return; // Контейнер скрыт или нулевой ширины

            this.setupCanvas(this.cvMain);
            this.setupCanvas(this.cvResid);

            const ctxM = this.cvMain.getContext("2d");
            const ctxR = this.cvResid.getContext("2d");

            // Определяем диапазон X
            let xRange = this.zoom ? [this.zoom.lo, this.zoom.hi] : spec.view;
            
            // Индексы видимых точек
            const iStart = Math.max(0, lowerBound(spec.e, xRange[0]) - 1);   // правка: точка за краем — линия до рамки
            const iEnd = Math.min(spec.e.length - 1, lowerBound(spec.e, xRange[1]));   // правка: включительно, без выхода за массив
            
            // Цвета и стили
            const ink = getCssVar("--ink", "#16140f");
            const dim = getCssVar("--dim", "#4c493f");
            const ruleSoft = getCssVar("--rule-soft", "#c9c4b4");
            const paper = getCssVar("--paper", "#f5f2ea");

            // Поля
            const padM = { l: 64, r: 14, t: 12, b: 34 };
            const padR = { l: 64, r: 14, t: 6, b: 30 };

            // Размеры областей построения
            const plotW = w - padM.l - padM.r;
            const plotHm = this.cvMain.clientHeight - padM.t - padM.b;
            const plotHr = this.cvResid.clientHeight - padR.t - padR.b;

            // Очистка и фон
            [ctxM, ctxR].forEach(ctx => {
                ctx.clearRect(0, 0, w, ctx.canvas.clientHeight);
                ctx.fillStyle = paper;
                ctx.fillRect(0, 0, w, ctx.canvas.clientHeight);
            });

            // Функция маппинга X -> px
            const xToPx = (e) => padM.l + (e - xRange[0]) / (xRange[1] - xRange[0]) * plotW;
            // Функция маппинга px -> E
            const pxToE = (px) => xRange[0] + (px - padM.l) / plotW * (xRange[1] - xRange[0]);

            // --- Основная канва ---
            
            // 1. Заливка вне диапазона подгонки [lo, hi]
            ctxM.fillStyle = ruleSoft;
            ctxM.globalAlpha = 0.35;
            const clampPx = (px) => Math.max(padM.l, Math.min(padM.l + plotW, px));   // правка: заливка не шире рамки при зуме
            if (spec.lo > xRange[0]) {
                const pxLo = clampPx(xToPx(spec.lo));
                ctxM.fillRect(padM.l, padM.t, pxLo - padM.l, plotHm);
            }
            if (spec.hi < xRange[1]) {
                const pxHi = clampPx(xToPx(spec.hi));
                ctxM.fillRect(pxHi, padM.t, padM.l + plotW - pxHi, plotHm);
            }
            ctxM.globalAlpha = 1;

            // 2. Подготовка данных Y
            let yMin = Infinity, yMax = -Infinity;
            const isLog = this.checkboxes.log.checked;
            
            // Собираем значения для определения диапазона Y
            const seriesData = [];
            
            // правка: прежние циклы сравнивали v > yMin при yMin = Infinity — нижняя граница не находилась никогда
            const cb = this.checkboxes; let yPos = Infinity;
            if (cb.meas.checked) seriesData.push({ color: ink, width: 1.2, calc: (i) => spec.meas[i] });
            if (cb.bg.checked) seriesData.push({ color: "#0f5aa8", width: 1.2, calc: (i) => spec.bg[i] });
            if (cb.model.checked) seriesData.push({ color: "#c8541c", width: 1.6, calc: (i) => spec.model[i] + spec.bg[i] });
            if (cb.diff.checked) seriesData.push({ color: "#6a6558", width: 1.2, calc: (i) => spec.meas[i] - spec.bg[i] });
            seriesData.forEach(s => { for (let i = iStart; i <= iEnd; i++) { const v = s.calc(i);
                if (v < yMin) yMin = v; if (v > yMax) yMax = v; if (v > 0 && v < yPos) yPos = v; } });
            if (!isFinite(yMax)) { yMin = 0; yMax = 1; }   // все серии выключены
            if (isLog) {
                yMin = Math.max(0.5, 0.7 * (isFinite(yPos) ? yPos : 1));
                yMax = Math.max(1.5 * yMax, 10 * yMin);
            } else {
                yMin = Math.min(0, yMin);
                yMax = Math.max(yMax, yMin) + (yMax > yMin ? 0.08 * (yMax - yMin) : 1);
            }

            // Функция маппинга Y -> px
            const yToPx = (v) => {
                if (isLog) {
                    if (v <= yMin) v = yMin;
                    const logMin = Math.log10(yMin);
                    const logMax = Math.log10(yMax);
                    const logV = Math.log10(v);
                    return padM.t + plotHm - (logV - logMin) / (logMax - logMin) * plotHm;
                } else {
                    return padM.t + plotHm - (v - yMin) / (yMax - yMin) * plotHm;
                }
            };

            // 3. Сетка и оси X
            ctxM.save();
            ctxM.rect(padM.l, padM.t, plotW, plotHm);
            ctxM.clip();

            // Сетка X
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

            // Сетка Y
            if (isLog) {
                logTicks(yMin, yMax).forEach(t => {   // правка: только деления внутри диапазона
                    const py = yToPx(t.v);
                    ctxM.moveTo(padM.l, py);
                    ctxM.lineTo(padM.l + plotW, py);
                });
            } else {
                const yStep = niceStep(yMax - yMin, 5);
                for (let v = Math.ceil(yMin / yStep) * yStep; v <= yMax; v += yStep) {
                    const py = yToPx(v);
                    ctxM.moveTo(padM.l, py);
                    ctxM.lineTo(padM.l + plotW, py);
                }
            }
            ctxM.stroke();

            // 4. Реперы (anchors)
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
                        
                        // Подпись
                        ctxM.fillStyle = dim;
                        ctxM.font = "12px ui-monospace, Consolas, monospace";
                        ctxM.textAlign = "center";
                        ctxM.fillText(label, px, padM.t + 12);
                    }
                });
                ctxM.setLineDash([]);
            }

            // 5. Рисование серий
            seriesData.forEach(s => {
                ctxM.strokeStyle = s.color;
                ctxM.lineWidth = s.width;
                ctxM.beginPath();
                let started = false;
                
                for (let i = iStart; i <= iEnd; i++) {
                    const e = spec.e[i];
                    const v = s.calc ? s.calc(i) : s.data[i];
                    
                    // Для лог-шкалы значения <= ymin рисуем на ymin
                    let drawV = v;
                    if (isLog && v <= yMin) drawV = yMin;
                    
                    const px = xToPx(e);
                    const py = yToPx(drawV);
                    
                    if (!started) {
                        ctxM.moveTo(px, py);
                        started = true;
                    } else {
                        ctxM.lineTo(px, py);
                    }
                }
                ctxM.stroke();
            });

            // 6. Курсор (вертикальная линия)
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

            // 7. Рамка выделения (drag)
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

            ctxM.restore(); // Снимаем clip

            // Рамка области построения (поверх всего)
            ctxM.strokeStyle = ink;
            ctxM.lineWidth = 1.5;
            ctxM.strokeRect(padM.l, padM.t, plotW, plotHm);

            // Подписи осей X и Y (вне clip)
            ctxM.fillStyle = ink;
            ctxM.font = "12px ui-monospace, Consolas, monospace";
            
            // Подписи X
            ctxM.textAlign = "center";
            for (let e = xStartTick; e <= xRange[1]; e += xStep) {
                const px = xToPx(e);
                ctxM.fillText(fmt(e, 0), px, padM.t + plotHm + 16);
            }
            
            // Подпись оси X (центр снизу)
            ctxM.fillText("энергия, кэВ", padM.l + plotW / 2, padM.t + plotHm + 30);

            // Подписи Y
            ctxM.textAlign = "right";
            if (isLog) {
                logTicks(yMin, yMax).forEach(t => ctxM.fillText(t.label, padM.l - 6, yToPx(t.v) + 4));
            } else {
                const yStep = niceStep(yMax - yMin, 5);
                for (let v = Math.ceil(yMin / yStep) * yStep; v <= yMax; v += yStep) {
                    const py = yToPx(v);
                    ctxM.fillText(fmt(v, 0), padM.l - 6, py + 4);
                }
            }

            // Подпись оси Y (повёрнута)
            ctxM.save();
            ctxM.translate(14, padM.t + plotHm / 2);
            ctxM.rotate(-Math.PI / 2);
            ctxM.textAlign = "center";
            ctxM.fillText("отсчётов на канал", 0, 0);
            ctxM.restore();


            // --- Канва остатков ---
            
            const padRl = padR.l;
            const plotWr = w - padR.l - padR.r;
            const plotHrVal = this.cvResid.clientHeight - padR.t - padR.b;

            ctxR.save();
            ctxR.rect(padRl, padR.t, plotWr, plotHrVal);
            ctxR.clip();

            // Заливка вне диапазона подгонки
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

            // Сетка X (такая же как на основной)
            ctxR.strokeStyle = ruleSoft;
            ctxR.lineWidth = 1;
            ctxR.beginPath();
            for (let e = xStartTick; e <= xRange[1]; e += xStep) {
                const px = xToPx(e);
                ctxR.moveTo(px, padR.t);
                ctxR.lineTo(px, padR.t + plotHrVal);
            }
            ctxR.stroke();

            // Оси Y: -5...+5. Горизонтали 0 (сплошная), ±2 (пунктир)
            const yResidToPx = (v) => {
                // clamp v to [-5, 5]
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

            // Точки остатков
            ctxR.fillStyle = ink;
            for (let i = iStart; i <= iEnd; i++) {
                if (spec.resid[i] == null) continue;
                const px = xToPx(spec.e[i]);
                const py = yResidToPx(spec.resid[i]);
                ctxR.fillRect(px - 1, py - 1, 2, 2); // Квадратик 2x2
            }

            // Курсор на канве остатков
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

            // Рамка выделения (drag) на канве остатков
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

            // Рамка области построения остатков
            ctxR.strokeStyle = ink;
            ctxR.lineWidth = 1.5;
            ctxR.strokeRect(padRl, padR.t, plotWr, plotHrVal);

            // Подписи осей X (повторяем)
            ctxR.fillStyle = ink;
            ctxR.font = "12px ui-monospace, Consolas, monospace";
            ctxR.textAlign = "center";
            for (let e = xStartTick; e <= xRange[1]; e += xStep) {
                const px = xToPx(e);
                ctxR.fillText(fmt(e, 0), px, padR.t + plotHrVal + 16);
            }
            
            // Подпись оси X
            ctxR.fillText("энергия, кэВ", padRl + plotWr / 2, padR.t + plotHrVal + 30);

            // Подписи Y остатков (-5, -2, 0, 2, 5)
            ctxR.textAlign = "right";
            [-5, -2, 0, 2, 5].forEach(v => {
                const py = yResidToPx(v);
                ctxR.fillText(String(v), padRl - 6, py + 4);
            });

            // Подпись оси Y остатков
            ctxR.save();
            ctxR.translate(14, padR.t + plotHrVal / 2);
            ctxR.rotate(-Math.PI / 2);
            ctxR.textAlign = "center";
            ctxR.fillText("(изм − модель)/σ", 0, 0);
            ctxR.restore();

            // --- Readout ---
            this.updateReadout(pxToE);
        },

        setupCanvas(cv) {
            const dpr = window.devicePixelRatio || 1;
            const w = cv.clientWidth;
            const h = cv.clientHeight;
            if (w === 0 || h === 0) return;
            
            cv.width = Math.round(w * dpr);
            cv.height = Math.round(h * dpr);
            const ctx = cv.getContext("2d");
            ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        },

        updateReadout(pxToE) {
            if (this.cursorE == null) {
                this.readout.textContent = "наведи курсор на спектр; протяжка — приближение, двойной клик — весь диапазон";
                return;
            }

            // Ближайший индекс к cursorE
            const idx = lowerBound(spec.e, this.cursorE);
            // Проверяем предыдущий и текущий, чтобы найти ближайший
            let bestIdx = idx;
            if (idx > 0) {
                const distCurr = Math.abs(spec.e[idx] - this.cursorE);
                const distPrev = Math.abs(spec.e[idx-1] - this.cursorE);
                if (distPrev < distCurr) bestIdx = idx - 1;
            }
            
            // Проверка границ
            if (bestIdx >= spec.e.length) bestIdx = spec.e.length - 1;
            if (bestIdx < 0) bestIdx = 0;

            const chNum = spec.ch0 + bestIdx;
            const eVal = spec.e[bestIdx];
            const measVal = spec.meas[bestIdx];
            const bgVal = spec.bg[bestIdx];
            const modelVal = spec.model[bestIdx] + spec.bg[bestIdx];
            const residVal = spec.resid[bestIdx];

            this.readout.textContent = 
                `канал ${chNum} · ${fmt(eVal, 1)} кэВ · изм. ${fmt(measVal, 1)} · фон ${fmt(bgVal, 1)} · модель + фон ${fmt(modelVal, 1)} · остаток ${fmt(residVal, 2)}`;
        }
    };

    blocks[key] = block;

    // --- Навигация и события ---

    const handlePointerDown = (e) => {
        if (e.button !== 0) return; // Только ЛКМ
        e.preventDefault();
        
        const rect = e.target.getBoundingClientRect();
        const x = e.clientX - rect.left;
        if (x < 64 || x > rect.width - 14) return;   // правка: протяжка только из области построения

        block.drag = { active: true, x0: x, x1: x, cv: e.target };
    };

    const handlePointerMove = (e) => {
        if (block.drag && block.drag.active && block.drag.cv === e.target) {
            const rect = e.target.getBoundingClientRect();
            const x = e.clientX - rect.left;
            block.drag.x1 = Math.max(64, Math.min(rect.width - 14, x));   // правка: зажим в область построения
            block.draw();
        } else {
            // Курсор
            const rect = e.target.getBoundingClientRect();
            const x = e.clientX - rect.left;
            
            // Проверяем, внутри ли области построения
            const padL = 64;
            const plotW = e.target.clientWidth - 64 - 14;
            
            if (x >= padL && x <= padL + plotW) {
                // Вычисляем энергию
                let xRange = block.zoom ? [block.zoom.lo, block.zoom.hi] : spec.view;
                const pxToE = (px) => xRange[0] + (px - padL) / plotW * (xRange[1] - xRange[0]);
                block.cursorE = pxToE(x);
            } else {
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

    // Слушатели на канвах
    [block.cvMain, block.cvResid].forEach(cv => {
        cv.addEventListener("mousedown", handlePointerDown);
        cv.addEventListener("pointermove", handlePointerMove);
        cv.addEventListener("pointerleave", handlePointerLeave);
        
        // Двойной клик — сброс зума
        cv.addEventListener("dblclick", () => {
            block.zoom = null;
            block.draw();
        });
    });

    // Глобальный mouseup для завершения drag
    document.addEventListener("mouseup", (e) => {
        if (!block.drag || !block.drag.active) return;
        
        block.drag.active = false;
        
        const dx = Math.abs(block.drag.x1 - block.drag.x0);
        if (dx < 6) {
            // Короткий клик — просто перерисовка (уже сделана в move/up, но на всякий случай)
            block.draw();
            return;
        }

        // Пересчет в энергии
        let xRange = block.zoom ? [block.zoom.lo, block.zoom.hi] : spec.view;
        const plotW = block.drag.cv.clientWidth - 64 - 14;   // правка: ширина канвы протяжки, не элемента под курсором
        const pxToE = (px) => xRange[0] + (px - 64) / plotW * (xRange[1] - xRange[0]);

        let lo = Math.max(spec.e[0], Math.min(pxToE(block.drag.x0), pxToE(block.drag.x1)));
        let hi = Math.min(spec.e[spec.e.length - 1], Math.max(pxToE(block.drag.x0), pxToE(block.drag.x1)));
        
        if (hi - lo >= 1) { // Минимум 1 кэВ
            block.zoom = { lo, hi };
        }
        
        block.draw();
    });

    // ResizeObserver
    const resizeObs = new ResizeObserver(() => {
        requestAnimationFrame(() => block.draw());
    });
    resizeObs.observe(chartBox);

    // Смена темы
    window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
        block.draw();
    });

});

})();
