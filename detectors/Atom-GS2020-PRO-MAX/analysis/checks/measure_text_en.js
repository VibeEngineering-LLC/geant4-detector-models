// #GS-79: копия measure_text.js для адресов с путём и #hash (английская страница /en/#src=...): ?m= ставится ДО хеша.
// Запуск — javascript_tool в Browser pane на локальной копии dist/; приёмка как у оригинала: bad=[] на 1400 и 375, sw<=W.
async function measure(W, url) {
  const u = url || '/', i = u.indexOf('#'), path = i < 0 ? u : u.slice(0, i), hash = i < 0 ? '' : u.slice(i);
  const f = document.createElement('iframe');
  f.style.cssText = 'position:absolute;left:-99999px;top:0;height:1200px;width:' + W + 'px';
  f.src = path + '?m=' + W + Math.random() + hash; document.body.appendChild(f);
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
  return { W, url: u, n: els.length, bad, sw };
}
const U = window.MEASURE_URLS || ['/en/'];
const R = []; for (const u of U) { R.push(await measure(1400, u)); R.push(await measure(375, u)); }
JSON.stringify(R);
