Напиши ОДИН файл JavaScript (ES5: var, function, без модулей, без сборки,
без стрелочных функций — как остальной код страницы). Верни ТОЛЬКО код, без
объяснений и без markdown-заборов. Комментарии — на РУССКОМ языке. Никаких
сетевых запросов, eval, innerHTML с неэкранированными данными (все строки из
данных пропускать через esc).

ИМЯ ФАЙЛА: src/scripts/g1s-amticseu-rerun.js

КОНТЕКСТ. Страница уже загрузила window.AMTICSEU (данные, схема ниже) и
scripts/g1s-amticseu.js, который выставил window.G1SA — набор помощников:
  num(x, d) — число с d знаками после запятой, десятичная запятая;
  cnt(x) — целое с разделителем тысяч;
  esc(s) — экранирование HTML;
  fit(canvas) — подготовка canvas под devicePixelRatio, возвращает {g, w, h};
  pal() — палитра {ink, faint, rule, grid, sum};
  mapX(v, lo, hi, x0, x1) — линейное отображение.
Файл — самовызывающаяся функция с "use strict". Если нет window.AMTICSEU или
window.G1SA — console.error("…") и return.
Экспортирует window.G1SA_RERUN = { showNet: …, showCmp: …, showCal: …, redraw: … }.

СХЕМА ДАННЫХ D = window.AMTICSEU (только то, что нужно этому файлу):
  D.nuclides: [{key, label_ru, color}]
  D.passport[key]: {A_Bq, dA_Bq}
  D.method1: {groups[key]: {A_Bq, dA_Bq (может быть null), A_over_passport},
    chi2, ndof, chi2_ndof, birge, bands: [{lo, hi, chi2, share_pct}],
    zone_50_70: {ratio, peak_meas_keV, peak_model_keV}}
  D.method1_bg_by_channel: {chi2_ndof, zone_ratio, groups[key]: {A_over_passport}}
  D.netarea: {definition: {win_fwhm, side_fwhm: [a, b], grid_keV, blend_fwhm},
    lines: [{nuclide, E_keV, main (bool), device: {E_keV, fwhm_keV, area, d_area} или null,
             fwhm_curve_keV, fwhm_ratio, eff_meas_pct, eff_model_pct, ratio,
             d_ratio_stat, blends: [{E_keV, nuclide, kind}]}]}
  D.cross_geometry: [{record, nuclide, E_keV, A_Bq, device: {E_keV, fwhm_keV, area, d_area},
    fwhm_curve_keV, fwhm_ratio, eff_meas_pct, eff_model_pct, ratio, template}]
  D.fwhm_cal: {points: [{E_keV, fwhm_keV, d_fwhm_keV, source}], curve: [[E, F], …],
    sqrt_law: [[E, F], …], f662_keV, device_peaks: [{E_keV, fwhm_keV, fwhm_curve_keV, ratio}]}
  D.meta: {model: {decays: [n, …]}, tail_T, live_s, real_s, bg_live_s, bg_real_s,
    bg_scale_time, cal_sample: {coefs_file: [c0, …], refs: [{E_true, ch, E_file,
    shift_keV, shift_fwhm}]}, cal_bg: {coefs_file: [c0, c1], coefs_refit: [c0, c1],
    rms_keV, anchors: [{E_keV, ch, netsum, resid_keV}]}}
Любое из этих полей может отсутствовать — тогда соответствующий блок пропускается
с console.warn, страница не должна падать.

ФУНКЦИИ:
 1. labelRu(key), colorOf(key) — по D.nuclides. sw(color) — строка
    "<span class='sw' style='background:" + color + "'></span>".
 2. mainNet(key) — строка D.netarea.lines с nuclide == key и main == true (или null).
 3. interp(arr, x) — линейная интерполяция по массиву пар [[x, y], …], возрастающему по x.
 4. fillValues() — каждому элементу document.querySelectorAll("[data-v]") выставить
    textContent по значению атрибута:
      n_refs → D.meta.cal_sample.refs.length
      n_fwhm_points → D.fwhm_cal.points.length
      zone_ratio → num(D.method1.zone_50_70.ratio, 3)
      decays → cnt(D.meta.model.decays[0]), если все элементы равны; иначе через « / »
      tail_T → num(D.meta.tail_T, 2)
      net_win → num(D.netarea.definition.win_fwhm, 2)
      net_side → num(side[0], 1) + "…" + num(side[1], 1)
      blend_k → num(D.netarea.definition.blend_fwhm, 1)
      am_fw_marinelli / am_fw_point / am_fw_petri → num(fwhm_ratio, 3) из
        D.cross_geometry с nuclide == "Am241" и record == "маринелли" /
        "точечный 5 см" / "Петри-60"
      am_point_ratio / am_petri_ratio → num(ratio, 3) из тех же строк
        («точечный 5 см» / «Петри-60», nuclide == "Am241")
      am_fw_curve_marinelli → num(fwhm_curve_keV, 2) строки Am241 «маринелли»
      am_marinelli_ratio → num(ratio, 3) строки Am241 «маринелли»
      cs_point_ratio / cs_petri_ratio → num(ratio, 3) строк nuclide == "Cs137chain"
        с record == "точечный 5 см" / "Петри-60"
      geo_spread → num(100 × max |ratio − 1| по ВСЕМ строкам D.cross_geometry,
        кроме строки Am241 «маринелли», 0)   (целое число процентов)
      bg_rms → num(D.meta.cal_bg.rms_keV, 2)
      n_bg_nat → D.meta.cal_bg.anchors.length − 1
      am_shift_fwhm → num(Math.abs(shift_fwhm), 2) репера с |E_true − 59.541| < 0.01
      sqrt59 → num(interp(D.fwhm_cal.sqrt_law, 59.541), 2)
      meas59 → num(fwhm_keV точки с |E_keV − 59.541| < 0.01, 2)
    Неизвестный ключ или отсутствующие данные — console.warn и текст «?».
 5. renderTwoPaths() → в #twoPathsBox таблица class='big': нуклид (sw + метка) |
    паспорт, Бк (cnt(A) + " ± " + cnt(dA)) | путь 1: отношение (num 3) |
    путь 2: отношение (num 3) и в скобках «по E кэВ» главной линии (num 1) |
    расхождение путей, % = |p2 − p1| / p1 · 100 (num 1) | пометка: если у главной
    линии есть бленды — «бленд: » + список «E (метка нуклида или «сумма»)»;
    иначе если расхождение ≤ 5 — «сходятся»; иначе «—».
 6. renderBands() → в #tblBands таблица: заголовок «вклады полос в χ²»: полоса
    «lo–hi кэВ» | χ² (num 0) | доля, % (num 1). Затем строки:
    «зона 50–70 кэВ: модель / измерение» | num(ratio, 3) | «максимумы: изм. X / мод. Y кэВ» (num 1);
    «χ²/ν, √(χ²/ν)» | num(chi2_ndof, 1) | num(birge, 2);
    «фон канал в канал (как в отчёте 11.09): зона 50–70, χ²/ν» |
      num(D.method1_bg_by_channel.zone_ratio, 3) | num(D.method1_bg_by_channel.chi2_ndof, 1).
 7. renderNet() → #sumM2: несколько блоков вида
    "<div><span class='lab'>ЯРЛЫК</span><span class='val'>ЗНАЧЕНИЕ</span></div>":
    окно (± num(win,2) ПШПВ), боковые полосы, сетка модели (num(grid,1) кэВ),
    порог бленда (num(blend,1) ПШПВ).
    #tblM2: все строки D.netarea.lines, сгруппированные по нуклиду в порядке
    D.nuclides, главная линия первой, её ячейка E — жирным (<b>). Столбцы:
    нуклид (sw + метка; только в первой строке группы) | E, кэВ (num 3) |
    пик прибора, кэВ (num 1) | ПШПВ прибора / по кривой, кэВ («num 2 / num 2») |
    отношение ширин (num 3) | площадь ± (cnt ± cnt) | эфф. измер., % (num 4) |
    эфф. модели, % (num 4) | изм./модель ± (num 3 ± num 3) |
    путь 1 (num 3, только в строке главной линии) | бленды («E метка» через запятую
    или «—»). Если device == null — в ячейках пика «нет пика в таблице прибора».
    Строке с отношением ширин вне [0,9; 1,1] добавить class='row-dirty'.
 8. renderAmCheck() → #tblAmCheck: строки D.cross_geometry. Столбцы: запись |
    нуклид (sw + метка) | E, кэВ (num 1) | пик прибора, кэВ (num 2) | ПШПВ прибора
    (num 2) | ПШПВ по кривой (num 2) | отношение ширин (num 3) | A, Бк (cnt) |
    эфф. измер., % (num 4) | эфф. модели, % (num 4) | изм./модель (num 3).
    Строка с отношением ширин вне [0,9; 1,1] — class='row-dirty'.
 9. coefsHtml(arr) — «c0 = …<br>c1 = …»: |c| от 0,01 до 10000 — num(c, 6),
    иначе c.toExponential(4) с заменой точки на запятую; каждое в <span class='mono'>.
10. renderCalTables() → #tblCal: параметр | образец | фон. Строки: живое время, с
    (num 2); реальное время, с (num 2); мёртвое время, % ((real − live)/real·100,
    num 3); шкала файла (coefsHtml coefs_file); шкала после проверки (образец:
    «кусочно-линейная по N реперам»; фон: coefsHtml(coefs_refit) + «<br>СКО num(rms,2) кэВ»);
    масштаб фона по времени (colspan 2, num 4).
    #tblRefs: подзаголовок-строка «реперы шкалы пробы»; столбцы E истинная (num 3) |
    канал (num 2) | E по шкале файла (num 2) | сдвиг, кэВ (num 2, со знаком) |
    сдвиг, доли ПШПВ (num 2, со знаком). Затем подзаголовок «якоря шкалы фона»;
    столбцы E (num 2) | канал (num 2) | чистый счёт (num 1) | невязка, кэВ (num 2).
11. renderFwhm() → #tblFwhm: E, кэВ (num 1) | источник (esc) | ПШПВ измер., кэВ
    (num 2 ± num 2) | закон √E, кэВ (interp(sqrt_law, E), num 2) | отклонение √E, %
    ((√E − измер.)/измер.·100, num 1, со знаком «+»/«−»).
    #tblFwhmDev: подзаголовок «пики этой записи (таблица прибора)»; E (num 1) |
    ПШПВ прибора (num 2) | по кривой (num 2) | отношение (num 3; строка вне
    [0,9; 1,1] — class='row-dirty').
    Canvas #cvFwhm через G1SA.fit; поля m = {l: 62, r: 16, t: 14, b: 34};
    x от 0 до 3000 кэВ, y от 0 до 1,15·max(всех точек, кривой, √E и пиков прибора);
    сетка pal().grid по x шагом 500 и по y шагом 20, подписи 11px system-ui цветом
    pal().faint, подписи осей «энергия, кэВ» (снизу) и «ПШПВ, кэВ» (повёрнута слева),
    рамка pal().rule толщиной 2. Кривая D.fwhm_cal.curve — сплошная #0f5aa8
    толщиной 2; закон √E — пунктир [5, 4] цветом pal().faint толщиной 1,5; точки
    points — закрашенные круги r = 4 цвета #c8541c; пики прибора этой записи —
    полые круги r = 4 цвета pal().ink. Легенда 4 строки в левом верхнем углу
    шрифтом "600 11px system-ui, sans-serif", каждая своим цветом:
    «измеренные точки комплекта», «кривая (интерполяция в логарифмах)»,
    «прежний закон √E», «пики этой записи (прибор)».
12. showNet() → renderNet(); showCmp() → renderAmCheck(); showCal() →
    renderCalTables() и renderFwhm(); redraw() → если элемент #viewCal не hidden —
    renderFwhm().
13. При загрузке файла: fillValues(); renderTwoPaths(); renderBands();
    window.addEventListener("resize", redraw).
