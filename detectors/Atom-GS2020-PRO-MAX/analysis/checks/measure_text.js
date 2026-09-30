// #GS-16 / W-152 (оператор 27.09 «текст должен быть по ширине страницы»): замер ПЕРЕД публикацией. Запуск —
// javascript_tool в Browser pane на локальной копии (python -m http.server в dist/); страница во фрейме ширины W.
// Приёмка: bad=[] на 1400 и 375 px, sw<=W. Блок «узкий», если он уже 90 % ширины содержимого своего родителя.
async function measure(W, url) {
  const f = document.createElement('iframe');
  f.style.cssText = 'position:absolute;left:-99999px;top:0;height:1200px;width:' + W + 'px';
  f.src = (url || '/') + '?m=' + W + Math.random(); document.body.appendChild(f);
  await new Promise(r => f.onload = r); await new Promise(r => setTimeout(r, 1500));
  const d = f.contentDocument, w = f.contentWindow, bad = [];
  const els = [...d.querySelectorAll('p,li,figcaption,dd')].filter(e => e.offsetParent && e.innerText.trim());
  for (const e of els) {
    const cs = w.getComputedStyle(e), pr = e.parentElement, ps = w.getComputedStyle(pr);
    const inner = pr.clientWidth - parseFloat(ps.paddingLeft) - parseFloat(ps.paddingRight);
    const frac = e.getBoundingClientRect().width / inner;
    if (frac < 0.9 || cs.textAlign === 'justify') bad.push({ t: e.innerText.slice(0, 40), frac: +frac.toFixed(2), al: cs.textAlign });
  }
  const sw = d.documentElement.scrollWidth; f.remove();
  return { W, n: els.length, bad, sw };
}
JSON.stringify([await measure(1400), await measure(375)]);
