// Итоговый раздел страницы GS2020 Th-232: таблицы и графики «нетто против моделей». Спека: scripts\specs\SPEC-result_js_gs2020.md
(function(){ "use strict";

    // Проверка наличия данных
    const D = window.GS2020_RESULT;
    if (!D) {
        const el = document.getElementById('res-lede');
        if (el) el.textContent = "Нет данных: result.js не загрузился.";
        return;
    }

    // Глобальный объект вида (опционально)
    const V = window.GS2020_VIEW;

    // --- Вспомогательные функции ---

    // Форматирование чисел: null/NaN -> "—", иначе toFixed с заменой точки на запятую
    function fmt(v, d) {
        if (v === null || v === undefined || Number.isNaN(v)) return "—";
        return v.toFixed(d).replace(".", ",");
    }

    // Получение CSS переменной
    function css(name) {
        return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    }

    // Цвета из CSS переменных
    const COLORS = {
        ink: css('--ink'),
        faint: css('--faint'),
        grid: css('--grid'),
        ruleSoft: css('--rule-soft'),
        m1: css('--acc-blue'),
        m2: css('--acc-orange'),
        one: css('--sum-line')
    };

    // Отступы графика
    const M = {l:64, r:16, t:14, b:36};

    // Инициализация Canvas с учетом devicePixelRatio
    function fitCanvas(cv) {
        const dpr = window.devicePixelRatio || 1;
        const r = cv.getBoundingClientRect();
        const w = Math.max(200, Math.floor(r.width));
        const h = Math.max(120, Math.floor(r.height));
        
        cv.width = Math.floor(w * dpr);
        cv.height = Math.floor(h * dpr);
        
        const g = cv.getContext("2d");
        g.setTransform(dpr, 0, 0, dpr, 0, 0);
        g.clearRect(0, 0, w, h);
        return {g, w, h};
    }

    // Текущий диапазон энергий
    function range() {
        if (V) return V.range();
        return [D.lo, D.hi];
    }

    // Генерация "красивых" делений шкалы
    function niceTicks(lo, hi, n) {
        const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / n)));
        const candidates = [1, 2, 5, 10];
        let bestStep = step;
        for (let c of candidates) {
            if (Math.ceil((hi - lo) / (step * c)) <= n) {
                bestStep = step * c;
                break;
            }
        }
        
        const ticks = [];
        let start = Math.ceil(lo / bestStep) * bestStep;
        for (let v = start; v <= hi + 1e-9; v += bestStep) {
            ticks.push(v);
        }
        return ticks;
    }

    // Логарифмические деления (степени 10)
    function logTicks(lo, hi) {
        const ticks = [];
        let p = Math.pow(10, Math.floor(Math.log10(lo)));
        while (p <= hi * 1.0001) {
            if (p >= lo * 0.9999) ticks.push(p);
            p *= 10;
        }
        return ticks;
    }

    // --- Заполнение таблиц и текста ---

    // Таблица результатов (#res-rows)
    const tbodyRows = document.getElementById('res-rows');
    if (tbodyRows && D.rows) {
        D.rows.forEach((row, idx) => {
            const tr = document.createElement('tr');
            if (idx === 0) tr.classList.add('main');
            
            // title (bold)
            const tdTitle = document.createElement('td');
            tdTitle.innerHTML = `<b>${row.title}</b>`;
            tr.appendChild(tdTitle);

            // what
            const tdWhat = document.createElement('td');
            tdWhat.textContent = row.what;
            tr.appendChild(tdWhat);

            // A ± dA
            const tdA = document.createElement('td');
            tdA.className = 'num';
            tdA.textContent = fmt(row.A, 1) + " ± " + fmt(row.dA, 1);
            tr.appendChild(tdA);

            // E1
            const tdE1 = document.createElement('td');
            tdE1.className = 'num';
            tdE1.textContent = fmt(row.E1, 1);
            tr.appendChild(tdE1);

            // ratio
            const tdRatio = document.createElement('td');
            tdRatio.className = 'num';
            tdRatio.textContent = fmt(row.ratio, 3);
            tr.appendChild(tdRatio);

            // chi2_ndof
            const tdChi = document.createElement('td');
            tdChi.className = 'num';
            tdChi.textContent = fmt(row.chi2_ndof, 2);
            tr.appendChild(tdChi);

            tbodyRows.appendChild(tr);
        });
    }

    // Таблица ссылок (#res-links)
    const tbodyLinks = document.getElementById('res-links');
    if (tbodyLinks && D.links) {
        D.links.forEach(link => {
            const tr = document.createElement('tr');
            
            // Method label
            const tdMethod = document.createElement('td');
            tdMethod.textContent = link.method === 'm1' ? "М1" : "М2";
            tr.appendChild(tdMethod);

            // Link text
            const tdLink = document.createElement('td');
            tdLink.textContent = link.link;
            tr.appendChild(tdLink);

            // A ± dA
            const tdA = document.createElement('td');
            tdA.className = 'num';
            let aStr = fmt(link.A, 1);
            if (link.dA !== null) aStr += " ± " + fmt(link.dA, 1);
            tdA.textContent = aStr;
            tr.appendChild(tdA);

            // chain_eq
            const tdChain = document.createElement('td');
            tdChain.className = 'num';
            tdChain.textContent = fmt(link.chain_eq, 1);
            tr.appendChild(tdChain);

            tbodyLinks.appendChild(tr);
        });
    }

    // Вывод текста (#res-lede)
    const ledeEl = document.getElementById('res-lede');
    if (ledeEl && D.rows.length >= 2) {
        const r0 = D.rows[0];
        const r1 = D.rows[1];
        const passportVal = fmt(D.passport_Bq, 1);
        const passportUnc = fmt(D.passport_Bq * D.passport_unc_pct / 100, 1);
        
        ledeEl.textContent = `Паспорт КИ: ${passportVal} ± ${passportUnc} Бк. ${r0.title} (фон без ослабления): ${fmt(r0.A, 1)} ± ${fmt(r0.dA, 1)} Бк, отношение к паспорту ${fmt(r0.ratio, 3)}. ${r1.title}: ${fmt(r1.A, 1)} Бк, отношение ${fmt(r1.ratio, 3)}.`;
    }

    // --- Состояние и управление масштабом ---

    const st = {
        scale: "log",
        cursor: null,
        drag: { active: false, kind: null, x0: 0, x1: 0 }
    };

    // Переключатель масштаба
    const ctlFitscale = document.getElementById('ctl-fitscale');
    if (ctlFitscale) {
        const btns = ctlFitscale.querySelectorAll('[data-v]');
        btns.forEach(btn => {
            const v = btn.getAttribute('data-v');
            // Инициализация aria-pressed
            btn.setAttribute('aria-pressed', v === st.scale ? 'true' : 'false');
            
            btn.addEventListener('click', () => {
                st.scale = v;
                btns.forEach(b => b.setAttribute('aria-pressed', b.getAttribute('data-v') === v ? 'true' : 'false'));
                draw();
            });
        });
    }

    // --- Отрисовка графиков ---

    const cvFit = document.getElementById('cv-fit');
    const cvRatio = document.getElementById('cv-ratio');
    const tipFit = document.getElementById('tip-fit');

    function draw() {
        if (cvFit) drawFit();
        if (cvRatio) drawRatio();
        
        // Отрисовка курсора и драга поверх всего, если нужно
        if (st.cursor !== null) {
            drawCursor(cvFit, st.cursor);
            drawCursor(cvRatio, st.cursor);
        }
        if (st.drag.active) {
            drawDragRect(st.drag.kind === 'fit' ? cvFit : cvRatio);
        }
    }

    // График "Fit"
    function drawFit() {
        const ctx = fitCanvas(cvFit);
        const {g, w, h} = ctx;
        const [xlo, xhi] = range();

        // Фильтрация видимых данных
        const visible = D.fit.E_keV.map((e, i) => ({e, net: D.fit.net[i], m1: D.fit.m1[i], m2: D.fit.m2[i]}))
                                  .filter(d => d.e >= xlo && d.e <= xhi);

        if (visible.length === 0) return;

        // Y-диапазон
        let yMax = -Infinity, yMinLin = Infinity;
        visible.forEach(d => {
            [d.net, d.m1, d.m2].forEach(v => {
                if (v !== null && v !== undefined) {
                    if (v > yMax) yMax = v;
                    if (v < yMinLin) yMinLin = v;
                }
            });
        });

        let yMin, isLog = (st.scale === 'log');
        
        if (isLog) {
            yMin = Math.max(1e-3, yMax * 1e-4);
        } else {
            yMin = Math.min(0, yMinLin);
        }

        // Функции маппинга координат
        const mapX = (E) => M.l + (E - xlo) * (w - M.l - M.r) / (xhi - xlo);
        
        const mapY = (v) => {
            let val = v;
            if (isLog && val < yMin) val = yMin; // clamp for log
            
            if (isLog) {
                return h - M.b - (h - M.t - M.b) * Math.log10(val / yMin) / Math.log10(yMax / yMin);
            } else {
                return h - M.b - (h - M.t - M.b) * (val - yMin) / (yMax - yMin);
            }
        };

        // Сетка и оси
        g.strokeStyle = COLORS.grid;
        g.lineWidth = 1;
        g.fillStyle = COLORS.ink;
        g.font = '10px sans-serif';
        g.textAlign = 'center';

        // X ticks
        const xTicks = niceTicks(xlo, xhi, 6);
        xTicks.forEach(t => {
            const x = mapX(t);
            if (x >= M.l && x <= w - M.r) {
                g.beginPath();
                g.moveTo(x, M.t);
                g.lineTo(x, h - M.b);
                g.stroke();
                g.fillText(Math.round(t), x, h - M.b + 14);
            }
        });

        // Y ticks
        g.textAlign = 'right';
        const yTicks = isLog ? logTicks(yMin, yMax) : niceTicks(yMin, yMax, 5);
        yTicks.forEach(t => {
            const y = mapY(t);
            if (y >= M.t && y <= h - M.b) {
                g.beginPath();
                g.moveTo(M.l, y);
                g.lineTo(w - M.r, y);
                g.stroke();
                // Форматирование метки Y
                let label = t >= 1 ? t.toFixed(0) : t.toPrecision(1);
                g.fillText(label, M.l - 4, y + 3);
            }
        });

        // Подписи осей
        g.fillStyle = COLORS.ink;
        g.textAlign = 'center';
        g.fillText("E, кэВ", w / 2, h - 4);
        
        g.save();
        g.translate(12, h / 2);
        g.rotate(-Math.PI / 2);
        g.fillText("отсчёты на кэВ", 0, 0);
        g.restore();

        // Рамка графика
        g.strokeStyle = COLORS.ink;
        g.strokeRect(M.l, M.t, w - M.l - M.r, h - M.t - M.b);

        // Серии данных
        const drawLine = (dataArr, color, width, dash) => {
            g.strokeStyle = color;
            g.lineWidth = width;
            g.setLineDash(dash || []);
            g.beginPath();
            let started = false;
            visible.forEach(d => {
                const val = dataArr.includes('net') ? d.net : (dataArr.includes('m1') ? d.m1 : d.m2);
                if (val !== null && val !== undefined) {
                    const x = mapX(d.e);
                    const y = mapY(val);
                    if (!started) { g.moveTo(x, y); started = true; }
                    else g.lineTo(x, y);
                }
            });
            g.stroke();
        };

        drawLine(['net'], COLORS.ink, 1);
        drawLine(['m1'], COLORS.m1, 1.4);
        drawLine(['m2'], COLORS.m2, 1.4, [5, 3]);

        // Легенда
        g.setLineDash([]);
        g.lineWidth = 1;
        const lx = w - M.r - 10;
        let ly = M.t + 10;
        g.textAlign = 'right';
        
        const legendItems = [
            { label: "нетто (измерение − k·фон)", color: COLORS.ink, dash: [] },
            { label: "модель М1", color: COLORS.m1, dash: [] },
            { label: "модель М2", color: COLORS.m2, dash: [5, 3] }
        ];

        legendItems.forEach(item => {
            g.strokeStyle = item.color;
            g.setLineDash(item.dash);
            g.beginPath();
            g.moveTo(lx - 40, ly);
            g.lineTo(lx, ly);
            g.stroke();
            g.fillStyle = COLORS.ink;
            g.fillText(item.label, lx - 45, ly + 3);
            ly += 12;
        });
    }

    // График "Ratio"
    function drawRatio() {
        const ctx = fitCanvas(cvRatio);
        const {g, w, h} = ctx;
        const [xlo, xhi] = range();

        // Фильтрация видимых данных
        const visible = D.ratio.E_keV.map((e, i) => ({
            e, 
            m1: D.ratio.m1[i], m1_sd: D.ratio.m1_sd[i],
            m2: D.ratio.m2[i], m2_sd: D.ratio.m2_sd[i]
        })).filter(d => d.e >= xlo && d.e <= xhi);

        if (visible.length === 0) return;

        // Y-диапазон
        let yMin = Infinity, yMax = -Infinity;
        visible.forEach(d => {
            [d.m1, d.m2].forEach((v, idx) => {
                const sd = idx === 0 ? d.m1_sd : d.m2_sd;
                if (v !== null && v !== undefined) {
                    const lo = v - (sd || 0);
                    const hi = v + (sd || 0);
                    if (lo < yMin) yMin = lo;
                    if (hi > yMax) yMax = hi;
                }
            });
        });

        // Padding 5% и включение 1
        const pad = (yMax - yMin) * 0.05;
        yMin -= pad;
        yMax += pad;
        if (yMin > 1) yMin = 1;
        if (yMax < 1) yMax = 1;

        // Clamp [0, 2.5]
        yMin = Math.max(0, yMin);
        yMax = Math.min(2.5, yMax);

        const mapX = (E) => M.l + (E - xlo) * (w - M.l - M.r) / (xhi - xlo);
        const mapY = (v) => h - M.b - (h - M.t - M.b) * (v - yMin) / (yMax - yMin);

        // Сетка и оси
        g.strokeStyle = COLORS.grid;
        g.lineWidth = 1;
        g.fillStyle = COLORS.ink;
        g.font = '10px sans-serif';
        g.textAlign = 'center';

        const xTicks = niceTicks(xlo, xhi, 6);
        xTicks.forEach(t => {
            const x = mapX(t);
            if (x >= M.l && x <= w - M.r) {
                g.beginPath();
                g.moveTo(x, M.t);
                g.lineTo(x, h - M.b);
                g.stroke();
                g.fillText(Math.round(t), x, h - M.b + 14);
            }
        });

        g.textAlign = 'right';
        const yTicks = niceTicks(yMin, yMax, 5);
        yTicks.forEach(t => {
            const y = mapY(t);
            if (y >= M.t && y <= h - M.b) {
                g.beginPath();
                g.moveTo(M.l, y);
                g.lineTo(w - M.r, y);
                g.stroke();
                g.fillText(t.toFixed(1), M.l - 4, y + 3);
            }
        });

        // Подписи осей
        g.textAlign = 'center';
        g.fillText("E, кэВ", w / 2, h - 4);
        g.save();
        g.translate(12, h / 2);
        g.rotate(-Math.PI / 2);
        g.fillText("нетто / модель", 0, 0);
        g.restore();

        // Рамка
        g.strokeStyle = COLORS.ink;
        g.strokeRect(M.l, M.t, w - M.l - M.r, h - M.t - M.b);

        // Линия y=1
        const y1 = mapY(1);
        if (y1 >= M.t && y1 <= h - M.b) {
            g.strokeStyle = COLORS.one;
            g.setLineDash([5, 5]);
            g.beginPath();
            g.moveTo(M.l, y1);
            g.lineTo(w - M.r, y1);
            g.stroke();
            g.setLineDash([]);
        }

        // Точки и ошибки
        const drawPoints = (valKey, sdKey, color, shiftX) => {
            g.fillStyle = color;
            g.strokeStyle = color;
            visible.forEach(d => {
                const v = d[valKey];
                const sd = d[sdKey];
                if (v !== null && v !== undefined) {
                    const x = mapX(d.e) + shiftX;
                    const y = mapY(v);
                    
                    // Ошибка
                    if (sd !== null && sd !== undefined) {
                        g.beginPath();
                        g.moveTo(x, mapY(v - sd));
                        g.lineTo(x, mapY(v + sd));
                        g.stroke();
                    }

                    // Точка (квадрат 4x4)
                    g.fillRect(x - 2, y - 2, 4, 4);
                }
            });
        };

        drawPoints('m1', 'm1_sd', COLORS.m1, 0);
        drawPoints('m2', 'm2_sd', COLORS.m2, 2);

        // Легенда
        g.textAlign = 'left';
        const lx = M.l + 5;
        let ly = M.t + 10;
        
        g.fillStyle = COLORS.m1;
        g.fillRect(lx, ly - 4, 8, 8);
        g.fillStyle = COLORS.ink;
        g.fillText("М1", lx + 12, ly + 3);

        ly += 12;
        g.fillStyle = COLORS.m2;
        g.fillRect(lx, ly - 4, 8, 8);
        g.fillStyle = COLORS.ink;
        g.fillText("М2", lx + 12, ly + 3);
    }

    // --- Взаимодействие: Курсор и Zoom ---

    function getEFromX(cv, x) {
        const r = cv.getBoundingClientRect();
        const w = r.width;
        const [xlo, xhi] = range();
        // Учитываем отступы M.l и M.r относительно ширины канваса (которая равна r.width в CSS пикселях)
        // Но fitCanvas возвращает w/h в CSS пикселях.
        // Координата x передается относительно элемента.
        const plotW = w - M.l - M.r;
        if (plotW <= 0) return null;
        const relX = x - M.l;
        if (relX < 0 || relX > plotW) return null;
        return xlo + relX * (xhi - xlo) / plotW;
    }

    function drawCursor(cv, E) {
        if (!cv) return;
        const r = cv.getBoundingClientRect();
        const w = r.width;
        const h = r.height;
        
        // Перерисовка не нужна полностью, просто линия поверх? 
        // В данной реализации мы перерисовываем всё через draw(), но линию курсора добавим отдельно.
        // Чтобы не мигало, лучше рисовать на том же контексте, но так как мы используем clearRect в fitCanvas,
        // нам нужно рисовать линию ПОСЛЕ очистки. 
        // Однако draw() вызывает fitCanvas который чистит. 
        // Поэтому логика курсора должна быть внутри drawFit/drawRatio или вызываться после них без очистки.
        // В текущей структуре draw() вызывает drawFit (который чистит), затем drawCursor.
        
        const ctx = cv.getContext("2d");
        const dpr = window.devicePixelRatio || 1;
        // Восстанавливаем трансформацию, так как fitCanvas её сбрасывает
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        
        const [xlo, xhi] = range();
        const plotW = w - M.l - M.r;
        const relX = (E - xlo) * plotW / (xhi - xlo);
        const x = M.l + relX;

        ctx.strokeStyle = COLORS.ink;
        ctx.lineWidth = 1;
        ctx.setLineDash([2, 2]);
        ctx.beginPath();
        ctx.moveTo(x, M.t);
        ctx.lineTo(x, h - M.b);
        ctx.stroke();
        ctx.setLineDash([]);

        // Обновление текста подсказки
        if (cv === cvFit && tipFit) {
            // Найти ближайший индекс в D.fit.E_keV
            let idx = 0;
            let minDiff = Infinity;
            D.fit.E_keV.forEach((e, i) => {
                const diff = Math.abs(e - E);
                if (diff < minDiff) { minDiff = diff; idx = i; }
            });
            
            const net = D.fit.net[idx];
            const m1 = D.fit.m1[idx];
            const m2 = D.fit.m2[idx];
            
            tipFit.textContent = `E ${fmt(E, 1)} кэВ · нетто ${fmt(net, 1)} · М1 ${fmt(m1, 1)} · М2 ${fmt(m2, 1)} на кэВ`;
            tipFit.style.display = 'block';
            
            // Позиционирование подсказки
            // Получаем координаты мыши из события pointermove, но здесь у нас только E.
            // Нужно хранить координаты мыши в st.cursor? Или передавать x,y.
            // Упростим: позиционируем относительно центра графика или просто показываем.
            // По ТЗ: "set left/top from pointer". Значит нужно знать координаты мыши.
            // Добавим mx, my в стейт курсора.
        } else if (cv === cvRatio && tipFit) {
             // Для ratio тоже используем тот же элемент подсказки? ТЗ говорит "show #tip-fit".
             // Найдем ближайший бин для ratio
            let idx = 0;
            let minDiff = Infinity;
            D.ratio.E_keV.forEach((e, i) => {
                const diff = Math.abs(e - E);
                if (diff < minDiff) { minDiff = diff; idx = i; }
            });

            const r1 = D.ratio.m1[idx];
            const r2 = D.ratio.m2[idx];
            const binE = D.ratio.E_keV[idx];

            tipFit.textContent = `E ${fmt(binE, 0)} кэВ · нетто/М1 ${fmt(r1, 3)} · нетто/М2 ${fmt(r2, 3)}`;
            tipFit.style.display = 'block';
        }
    }

    function drawDragRect(cv) {
        if (!cv || !st.drag.active) return;
        const ctx = cv.getContext("2d");
        const dpr = window.devicePixelRatio || 1;
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        
        const r = cv.getBoundingClientRect();
        const w = r.width;
        const h = r.height;

        // Ограничиваем прямоугольник областью графика
        let x1 = Math.max(M.l, Math.min(w - M.r, st.drag.x0));
        let x2 = Math.max(M.l, Math.min(w - M.r, st.drag.x1));
        
        const left = Math.min(x1, x2);
        const width = Math.abs(x2 - x1);

        ctx.fillStyle = COLORS.ruleSoft;
        ctx.globalAlpha = 0.25;
        ctx.fillRect(left, M.t, width, h - M.t - M.b);
        ctx.globalAlpha = 1.0;
    }

    // Обработчики событий для графиков
    const setupInteraction = (cv, kind) => {
        if (!cv) return;

        cv.addEventListener('pointermove', (e) => {
            const rect = cv.getBoundingClientRect();
            const x = e.clientX - rect.left;
            const y = e.clientY - rect.top;

            if (st.drag.active && st.drag.kind === kind) {
                st.drag.x1 = x;
                draw(); // Перерисовка для обновления прямоугольника
                return;
            }

            const E = getEFromX(cv, x);
            if (E === null) {
                st.cursor = null;
                if (tipFit) tipFit.style.display = 'none';
                draw();
            } else {
                // Сохраняем координаты мыши для позиционирования подсказки
                st.cursor = { E, x: e.clientX, y: e.clientY };
                
                // Позиционирование подсказки
                if (tipFit) {
                    tipFit.style.left = (e.clientX + 10) + 'px';
                    tipFit.style.top = (e.clientY - 30) + 'px';
                }
                
                draw();
            }
        });

        cv.addEventListener('pointerleave', () => {
            if (st.drag.active) return; // Не сбрасываем при драге, если мышь ушла за пределы (хотя обычно драг завершается mouseup)
            st.cursor = null;
            if (tipFit) tipFit.style.display = 'none';
            draw();
        });

        cv.addEventListener('mousedown', (e) => {
            if (e.button !== 0) return; // Только левая кнопка
            const rect = cv.getBoundingClientRect();
            const x = e.clientX - rect.left;
            
            // Проверка, что клик внутри области графика
            if (x >= M.l && x <= rect.width - M.r) {
                e.preventDefault();
                st.drag = { active: true, kind: kind, x0: x, x1: x };
            }
        });

        cv.addEventListener('dblclick', () => {
            if (V) V.reset();
            else draw();
        });
    };

    setupInteraction(cvFit, 'fit');
    setupInteraction(cvRatio, 'ratio');

    // Глобальный mouseup для завершения драга
    document.addEventListener('mouseup', (e) => {
        if (!st.drag.active) return;
        
        st.drag.active = false;
        
        const rect = cvFit.getBoundingClientRect(); // Используем fit как референс, или тот где был драг
        // Но нам нужно знать x1 относительно того же канваса. 
        // В pointermove мы обновляли st.drag.x1.
        
        // Проверяем разницу
        if (Math.abs(st.drag.x1 - st.drag.x0) < 6) {
            draw();
            return;
        }

        // Конвертация координат в E
        // Нужно знать, на каком канвасе был драг, чтобы правильно масштабировать? 
        // Диапазон X одинаковый для обоих.
        const cv = st.drag.kind === 'fit' ? cvFit : cvRatio;
        if (!cv) return;
        
        const r = cv.getBoundingClientRect();
        const w = r.width;
        const [xlo, xhi] = range();
        const plotW = w - M.l - M.r;
        
        const getE = (px) => {
            const relX = px - M.l;
            return xlo + relX * (xhi - xlo) / plotW;
        };

        const e1 = getE(st.drag.x0);
        const e2 = getE(st.drag.x1);
        
        let lo = Math.max(D.lo, Math.min(e1, e2));
        let hi = Math.min(D.hi, Math.max(e1, e2));

        if (V) {
            V.setZoom(lo, hi);
        } else {
            draw();
        }
    });

    // --- Инициализация и наблюдатели ---

    if (V) V.listeners.push(draw);

    const fitParent = document.getElementById("cv-fit").parentElement;
    if (fitParent) {
        new ResizeObserver(() => requestAnimationFrame(draw)).observe(fitParent);
    }

    // Первый вызов отрисовки
    draw();

})();
