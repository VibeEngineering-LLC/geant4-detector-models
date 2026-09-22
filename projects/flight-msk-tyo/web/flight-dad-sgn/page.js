    const LANG = window.PAGE_LANG === 'en' ? 'en' : 'ru';
    const STR = {
      ru: { errData: 'Ошибка: данные не загружены', axisE: 'E, кэВ', axisY: 'отсч./ч/кэВ', keV: 'кэВ',
        detectorPath: 'нейтроны, дошедшие до прибора: реакции в самом приборе (кристалл, корпус, плата)',
        full: 'полный', unidentified: 'неотождествлённый пик', near: 'рядом', notSignif: 'не значим',
        candidateDb: 'кандидат в базе', reaction: 'реакция', noMatch: 'в базе реакций рядом ничего нет — пик неотождествлён',
        noInList: 'нет в перечне; ближайший кандидат в базе —', shift: 'сдвиг', unconfirmed: 'не подтверждён',
        thReaction: 'Реакция', thArea: 'Площадь ± σ, отсч./ч', thSignif: 'Значимость', thOrigin: 'Происхождение',
        modeComp: 'Спектр разложен по типу частицы космического излучения, пришедшей снаружи. Цветные слои сложены друг на друга и в сумме дают полный спектр (чёрная линия); меньший вклад — внизу, больший — сверху.',
        modeOrig: 'Другой разрез тех же отсчётов: какая частица дошла до прибора, каким процессом и где она родилась. Это не дополнение к разложению «по частицам поля», а тот же спектр, разрезанный иначе, — смотрите одно из двух.' },
      en: { errData: 'Error: data not loaded', axisE: 'E, keV', axisY: 'counts/h/keV', keV: 'keV',
        detectorPath: 'neutrons reaching the instrument: reactions inside the instrument itself (crystal, housing, PCB)',
        full: 'total', unidentified: 'unidentified peak', near: 'nearby', notSignif: 'not significant',
        candidateDb: 'candidate in the database', reaction: 'reaction', noMatch: 'no candidate nearby in the reaction database — peak unidentified',
        noInList: 'not in the list; nearest candidate in the database —', shift: 'shift', unconfirmed: 'unconfirmed',
        thReaction: 'Reaction', thArea: 'Area ± σ, counts/h', thSignif: 'Significance', thOrigin: 'Origin',
        modeComp: 'The spectrum is broken down by the type of cosmic-ray particle arriving from outside. Colored layers are stacked and together give the full spectrum (black line); smaller contributions are at the bottom, larger ones at the top.',
        modeOrig: 'A different cut of the same counts: which particle reached the instrument, by which process, and where it originated. This is not additional to the "field particle" breakdown — it is the same spectrum cut differently; look at one or the other.' }
    }[LANG];
    if (!window.SPECTRA) {
      document.body.innerHTML = '<h1>' + STR.errData + '</h1>';
      throw new Error("SPECTRA not defined");
    }

    const { E0, dE, unit, series, lines, unlisted, candidates } = window.SPECTRA;
    const n = (series.find((s) => s.id === 'total') || series[0]).y.length;   // длина реального массива данных (формула по E0/dE давала n+1 — NaN на правом краю)

    function idx(E) {
      return Math.max(0, Math.min(n - 1, Math.round((E - E0) / dE)));
    }

    function ticksLin(lo, hi, n) {
      const raw = (hi - lo) / n, p = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / p;
      const step = (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p, out = [];
      for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * step; v += step) out.push(v);
      return out;
    }

    const M = { l: 64, r: 16, t: 12, b: 34 };
    let zoom = null;
    let drag = { active: false, x0: 0, x1: 0 };
    let logY = true;
    let hoverPx = null;
    let mode = 'comp', legendKey = '';   // разложение: comp — по частицам поля, orig — по происхождению
    const hiddenLayer = {};
    const PAL = ['#222222','#d62728','#1f77b4','#2ca02c','#ff7f0e','#9467bd','#8c564b','#e377c2','#17becf','#bcbd22','#7f7f7f','#393b79','#e6550d','#31a354','#756bb1','#636363','#b15928','#fb9a99','#6a3d9a','#33a02c','#a6cee3','#9e9e9e'];
    const colorOf = (id) => PAL[series.findIndex((q) => q.id === id) % PAL.length];
    const grpOf = () => (mode === 'comp' ? 'component' : 'origin');
    const layersNow = (all) => series.filter((s) => s.group === grpOf() && (all || !hiddenLayer[s.id])).map((s) => ({ id: s.id, y: s.y }));
    const nameOf = (s) => (s.label === 'neutron / primary / World' ? STR.detectorPath : (s.group === 'origin' ? russianOrigin(s.label, LANG) : s.label));

    const cv = document.getElementById('cv');
    const ctx = cv.getContext('2d');
    const readout = document.getElementById('readout');
    const linesBar = document.getElementById('lines-bar'), linesDetail = document.getElementById('lines-detail');
    const unlistedBar = document.getElementById('unlisted-bar'), unlistedDetail = document.getElementById('unlisted-detail');

    const shown = {};
    series.forEach((s, i) => {
      shown[s.id] = s.group === 'total';
    });

    function esc(s) {
      return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
    }

    function setPreset(lo, hi) {
      zoom = null;
      const btns = document.querySelectorAll('#toolbar button');
      for (let i = 0; i < btns.length; i++) {
        if (parseInt(btns[i].dataset.lo) === lo && parseInt(btns[i].dataset.hi) === hi) {
          btns[i].classList.add('active');
        } else {
          btns[i].classList.remove('active');
        }
      }
      draw();
    }

    function currentRange() {
      if (zoom) return [zoom.lo, zoom.hi];
      const btn = document.querySelector('#toolbar button.active');
      if (btn) {
        return [parseInt(btn.dataset.lo), parseInt(btn.dataset.hi)];
      }
      return [20, 10000];
    }

    function markPreset(id) {
      const btns = document.querySelectorAll('#toolbar button');
      for (let i = 0; i < btns.length; i++) {
        if (id && parseInt(btns[i].dataset.lo) === id.lo && parseInt(btns[i].dataset.hi) === id.hi) {
          btns[i].classList.add('active');
        } else {
          btns[i].classList.remove('active');
        }
      }
    }

    function draw() {
      const w = cv.clientWidth, h = 420, dpr = window.devicePixelRatio || 1;
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(h * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

      const [xlo, xhi] = currentRange();
      const pw = w - M.l - M.r, ph = h - M.t - M.b;
      const X = (E) => M.l + (E - xlo) / (xhi - xlo) * pw;

      const key = xlo + '|' + xhi + '|' + mode;
      if (key !== legendKey) { legendKey = key; buildLegend(); }
      const j0 = idx(xlo), j1 = idx(xhi);
      const st = buildStack(layersNow(false), j0, j1);   // по возрастанию вклада: меньший снизу
      const totS = series.find((s) => s.id === 'total');
      let { ymax, ymin } = yRange(st.map((l) => l.upper).concat(shown.total ? [totS.y] : []), j0, j1);

      const css = getComputedStyle(document.documentElement);
      const C = (n) => css.getPropertyValue(n).trim();

      ctx.clearRect(0, 0, w, h);

      if (ymax <= 0) {
        ctx.strokeStyle = C('--axis');
        ctx.beginPath();
        ctx.moveTo(M.l, h - M.b);
        ctx.lineTo(w - M.r, h - M.b);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(M.l, M.t);
        ctx.lineTo(M.l, h - M.b);
        ctx.stroke();
        return;
      }

      const Y = logY ? (v) => {
        const lo = Math.max(ymin, ymax * 1e-6);
        return M.t + ph * (1 - (Math.log10(Math.max(v, lo)) - Math.log10(lo)) / (Math.log10(ymax * 1.3) - Math.log10(lo)));
      } : (v) => {
        return M.t + ph * (1 - v / (ymax * 1.05));
      };

      // grid
      ctx.strokeStyle = C('--grid');
      ctx.lineWidth = 1;
      const xTicks = ticksLin(xlo, xhi, 8);
      for (const x of xTicks) {
        const px = X(x);
        ctx.beginPath();
        ctx.moveTo(px, M.t);
        ctx.lineTo(px, h - M.b);
        ctx.stroke();
      }

      let yTicks;
      if (logY) {
        const lo = Math.max(ymin, ymax * 1e-6);
        const hi = ymax * 1.3;
        const logLo = Math.log10(lo), logHi = Math.log10(hi);
        const step = Math.ceil((logHi - logLo) / 4);
        yTicks = [];
        for (let i = Math.ceil(logLo / step) * step; i <= logHi + 1e-9 * step; i += step) {
          yTicks.push(Math.pow(10, i));
        }
      } else {
        yTicks = ticksLin(0, ymax * 1.05, 6);
      }

      for (const y of yTicks) {
        const py = Y(y);
        ctx.beginPath();
        ctx.moveTo(M.l, py);
        ctx.lineTo(w - M.r, py);
        ctx.stroke();
      }

      // axis
      ctx.strokeStyle = C('--axis');
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(M.l, h - M.b);
      ctx.lineTo(w - M.r, h - M.b);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(M.l, M.t);
      ctx.lineTo(M.l, h - M.b);
      ctx.stroke();

      // labels
      ctx.fillStyle = C('--text');
      ctx.font = '12px sans-serif';
      ctx.textAlign = 'center';
      for (const x of xTicks) {
        const px = X(x);
        ctx.fillText(x.toFixed(0), px, h - M.b + 16);
      }
      ctx.textAlign = 'right';
      for (const y of yTicks) {
        const py = Y(y);
        ctx.fillText(y.toExponential(2), M.l - 4, py + 4);
      }

      // axis title
      ctx.textAlign = 'center';
      ctx.fillText(STR.axisE, w / 2, h - 4);
      ctx.save();
      ctx.translate(8, h / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(STR.axisY, 0, 0);
      ctx.restore();

      // series (цвет — colorOf()/PAL, строка 27; здесь была неиспользуемая копия того же массива)
      const drawSeries = (s, color, width) => {
        ctx.beginPath();
        let lastX = -1;
        let lastY = -1;
        for (let j = idx(xlo); j <= idx(xhi); j++) {
          const v = s.y[j];
          if (v <= 0) continue;
          const x = X(E0 + j * dE);
          const y = Y(v);
          if (lastX === -1) {
            ctx.moveTo(x, y);
          } else {
            ctx.lineTo(x, y);
          }
          lastX = x;
          lastY = y;
        }
        ctx.strokeStyle = color;
        ctx.lineWidth = width;
        ctx.stroke();
      };

      drawStackLayers(ctx, st, { j0, j1, E0, dE, X, Y: (v) => Y(Math.max(v, 1e-30)), color: colorOf });
      if (shown.total) drawSeries(totS, C('--ink'), 2);

      if (hoverPx !== null && !drag.active) {   // курсор: вертикальный пунктир и подпись энергии
        const Eh = xlo + (hoverPx - M.l) / pw * (xhi - xlo);
        ctx.save(); ctx.setLineDash([4, 4]); ctx.strokeStyle = C('--ink'); ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(hoverPx, M.t); ctx.lineTo(hoverPx, h - M.b); ctx.stroke(); ctx.restore();
        const jh = idx(Eh), txt = [`${Eh.toFixed(1)} ${STR.keV}`];
        if (shown.total && totS.y[jh] > 0) txt.push(`${STR.full}: ${totS.y[jh].toPrecision(3)} ${STR.axisY}`);
        const topL = st.slice().reverse().slice(0, 3);
        for (const l of topL) { const v = l.upper[jh] - l.lower[jh]; if (v > 0) txt.push(`${nameOf(series.find((q) => q.id === l.id))}: ${v.toPrecision(3)}`); }
        const nearest = (arr, key) => { let b = null, bd = Infinity; for (const x of arr) { const d = Math.abs(Eh - x[key]); if (d < bd && d <= Math.max(3, x[key] * 0.01)) { bd = d; b = x; } } return b; };
        const trunc = (s) => s.length > 70 ? s.slice(0, 68) + '…' : s;
        const nearLine = nearest(lines, 'E');
        if (nearLine) { txt.push(trunc(nearLine.reaction)); }
        else {
          const nu = nearest(unlisted, 'E');
          // Δ кандидата считать от энергии САМОГО ПИКА (nu.E), не от позиции курсора Eh —
          // иначе число расходится с таблицей «Пики вне перечня» при наведении мимо центра пика.
          const candRef = nu ? nu.E : Eh;
          const candNear = (() => { let b = null, bd = Infinity; for (const x of (candidates || [])) { const d = Math.abs(candRef - x.E); if (d < bd && d <= Math.max(3, x.E * 0.01)) { bd = d; b = x; } } return b; })();
          if (nu) txt.push(`${STR.unidentified}, ${nu.signif.toFixed(1)}σ` + (candNear ? ` (${STR.near}: ${trunc(candNear.reaction)}, Δ=${(candNear.E - candRef).toFixed(1)})` : ''));
          else if (candNear) txt.push(`${STR.candidateDb}: ${trunc(candNear.reaction)} (Δ=${(candNear.E - candRef).toFixed(1)} ${STR.keV}, ${STR.notSignif})`);
        }
        ctx.font = '12px sans-serif';
        const bwMax = w - M.l - M.r - 16;
        const bw = Math.min(bwMax, Math.max(...txt.map((t) => ctx.measureText(t).width)) + 12), bh = 16 * txt.length + 6;
        const bx = Math.max(M.l, Math.min(w - M.r - bw, hoverPx + 8 + bw > w - M.r ? hoverPx - 8 - bw : hoverPx + 8));
        ctx.fillStyle = C('--paper'); ctx.globalAlpha = 0.92; ctx.fillRect(bx, M.t + 4, bw, bh); ctx.globalAlpha = 1;
        ctx.strokeStyle = C('--ink'); ctx.strokeRect(bx, M.t + 4, bw, bh);
        ctx.fillStyle = C('--ink'); ctx.textAlign = 'left';
        txt.forEach((t, i) => ctx.fillText(t, bx + 6, M.t + 20 + 16 * i));
      }

      // selection
      if (drag.active) {
        ctx.fillStyle = 'rgba(0, 0, 255, 0.2)';
        const x0 = Math.min(drag.x0, drag.x1), x1 = Math.max(drag.x0, drag.x1);
        ctx.fillRect(x0, M.t, x1 - x0, ph);
      }
    }

    function updateReadout(px) {
      if (px < M.l || px > cv.clientWidth - M.r) return;
      const [xlo, xhi] = currentRange();
      const E = xlo + (px - M.l) / (cv.clientWidth - M.r - M.l) * (xhi - xlo);
      const items = [];
      for (const s of series) {
        if (!(s.id === 'total' ? shown.total : (s.group === grpOf() && !hiddenLayer[s.id]))) continue;
        const i = idx(E);
        if (i >= 0 && i < s.y.length) {
          const v = s.y[i];
          if (v > 0) {
            let item = `${esc(nameOf(s))}: ${v.toFixed(3)}`;
            if (s.id === 'total' && s.sig && i < s.sig.length) item += ` ± ${s.sig[i].toFixed(3)}`;
            items.push(item);
          }
        }
      }
      let best = null, bestD = Infinity;
      for (const l of lines) {
        const d = Math.abs(E - l.E);
        if (d < bestD && d <= Math.max(3, l.E * 0.01)) { bestD = d; best = l; }
      }
      if (best) items.push(`${STR.reaction}: ${esc(best.reaction)}`);
      readout.innerHTML = `<b>E = ${E.toFixed(1)} ${STR.keV}</b>` + (items.length ? `<ul>${items.map((t) => `<li>${t}</li>`).join('')}</ul>` : '');
    }

    function lineTabs(bar, detail, rows, hasReaction) {   // #CHART-1: клик по вкладке — и деталь, и зум графика
      bar.innerHTML = '';
      rows.forEach((r, i) => {
        const b = document.createElement('button');
        b.type = 'button'; b.className = 'btn'; b.setAttribute('role', 'tab');
        b.setAttribute('aria-selected', i === 0 ? 'true' : 'false');
        b.textContent = r.E.toFixed(1);
        b.onclick = () => {
          bar.querySelectorAll('button').forEach((x) => x.setAttribute('aria-selected', 'false'));
          b.setAttribute('aria-selected', 'true');
          detail.innerHTML = lineDetailHtml(r, hasReaction);
          zoom = { lo: Math.max(20, r.E * 0.92), hi: Math.min(10000, r.E * 1.08) };
          markPreset(null);
          draw();
        };
        bar.appendChild(b);
      });
      detail.innerHTML = rows.length ? lineDetailHtml(rows[0], hasReaction) : '';
    }

    function lineDetailHtml(r, hasReaction) {
      const reactionRow = hasReaction ? `<tr><th>${STR.thReaction}</th><td>${esc(r.reaction)}</td></tr>`
        : `<tr><th>${STR.thReaction}</th><td>${r.nearReaction
          ? `${STR.noInList} ${esc(r.nearReaction)} (${STR.shift} ${r.nearDE > 0 ? '+' : ''}${r.nearDE} ${STR.keV}, ${STR.unconfirmed})`
          : STR.noMatch}</td></tr>`;
      return `<table style="width:auto"><tbody>
        <tr><th>${STR.axisE}</th><td class="num">${r.E.toFixed(2)}</td></tr>
        ${reactionRow}
        <tr><th>${STR.thArea}</th><td class="num">${r.net.toFixed(3)} ± ${r.sig.toFixed(3)}</td></tr>
        <tr><th>${STR.thSignif}</th><td class="num">${r.signif.toFixed(2)}</td></tr>
        <tr><th>${STR.thOrigin}</th><td>${esc(humanOrigin(r.origin)).replace(/\n/g, '<br>')}</td></tr>
      </tbody></table>`;
    }

    function buildTables() {
      lineTabs(linesBar, linesDetail, lines, true);
      lineTabs(unlistedBar, unlistedDetail, unlisted, false);
    }

    document.querySelectorAll('#toolbar button').forEach(btn => {
      btn.addEventListener('click', () => {
        const lo = parseInt(btn.dataset.lo);
        const hi = parseInt(btn.dataset.hi);
        setPreset(lo, hi);
      });
    });

    document.getElementById('logy-toggle').addEventListener('change', (ev) => {
      logY = ev.target.checked;
      draw();
    });

    cv.addEventListener("mousedown", (ev) => {
      const r = cv.getBoundingClientRect(), x = ev.clientX - r.left;
      if (x < M.l || x > r.width - M.r) return;
      ev.preventDefault();
      drag = { active: true, x0: x, x1: x };
      draw();
    });

    document.addEventListener("pointermove", (ev) => {
      const r = cv.getBoundingClientRect(), x = ev.clientX - r.left;
      if (drag.active) {
        drag.x1 = Math.max(M.l, Math.min(r.width - M.r, x));
        draw();
        return;
      }
      const inside = ev.target === cv && x >= M.l && x <= r.width - M.r;
      const nh = inside ? x : null;
      if (nh !== hoverPx) { hoverPx = nh; draw(); }
      if (inside) updateReadout(x);
    });

    document.addEventListener("mouseup", () => {
      if (!drag.active) return;
      drag.active = false;
      if (Math.abs(drag.x1 - drag.x0) < 6) {
        draw();
        return;
      }
      const r = cv.getBoundingClientRect(), [xlo, xhi] = currentRange();
      const toE = (px) => xlo + (px - M.l) / (r.width - M.r - M.l) * (xhi - xlo);
      const lo = Math.max(20, toE(Math.min(drag.x0, drag.x1))), hi = Math.min(10000, toE(Math.max(drag.x0, drag.x1)));
      if (hi - lo < 1) {
        draw();
        return;
      }
      zoom = { lo, hi };
      markPreset(null);
      draw();
    });

    cv.addEventListener("dblclick", () => {
      zoom = null;
      setPreset(20, 10000);
    });

    const ro = new ResizeObserver(entries => {
      for (const entry of entries) {
        if (entry.contentRect.width !== entry.target.previousWidth) {
          entry.target.previousWidth = entry.contentRect.width;
          draw();
        }
      }
    });
    ro.observe(cv.parentNode);

    function humanOrigin(t) {   // «gamma/primary/World 58.1%» -> человекочитаемая строка (RU/EN по LANG)
      return String(t).split('; ').map((p) => {
        const m = /^(\S+)\/(\S+)\/(\S+) ([\d.]+)%$/.exec(p.trim());
        const pct = LANG === 'en' ? m && m[4] : m && m[4].replace('.', ',');
        return m ? russianOrigin(m[1] + ' / ' + m[2] + ' / ' + m[3], LANG) + ' — ' + pct + ' %' : p;
      }).join('\n');
    }
    function buildLegend() {
      const [xlo, xhi] = currentRange();
      const st = buildStack(layersNow(true), idx(xlo), idx(xhi));
      const rows = legendRows(series, st, { dE, hidden: hiddenLayer, colorOf, nameOf });
      const box = document.getElementById('series-box');
      box.innerHTML = legendHtml(rows, shown.total, getComputedStyle(document.documentElement).getPropertyValue('--ink').trim());
      box.querySelectorAll('input').forEach((inp) => inp.addEventListener('change', () => {
        if (inp.dataset.id === 'total') shown.total = inp.checked; else hiddenLayer[inp.dataset.id] = !inp.checked;
        draw();
      }));
    }
    // здесь была старая ручная сборка легенды (const groups=[] пуст, цикл никогда не
    // выполнялся, seriesContainer никуда не добавлялся в DOM) — мёртвый код, живой путь ниже через
    // buildLegend()/stackui.js. Удалено.

    const noteText = { comp: STR.modeComp, orig: STR.modeOrig };
    const modeBtns = document.querySelectorAll('#mode-seg .btn');
    function setMode(m) {
      mode = m; legendKey = '';
      modeBtns.forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.mode === m)));
      document.getElementById('mode-note').textContent = noteText[m];
      draw();
    }
    modeBtns.forEach((b) => b.addEventListener('click', () => setMode(b.dataset.mode)));
    document.getElementById('reset').addEventListener('click', () => {   // сброс: все слои, полный спектр, весь диапазон
      for (const k in hiddenLayer) delete hiddenLayer[k];
      shown.total = true; zoom = null; logY = true;
      const lg = document.getElementById('logy-toggle'); if (lg) lg.checked = true;
      markPreset({ lo: 20, hi: 10000 }); setMode('comp');
    });
    markPreset({ lo: 20, hi: 10000 });
    setMode('comp');

    buildTables();
    draw();
  