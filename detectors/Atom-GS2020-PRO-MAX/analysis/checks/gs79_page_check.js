// #GS-79: приёмка страницы GS2020 в Browser pane (javascript_tool): вкладки источников, подвкладки, канвы,
// всплывающие окна, протяжка-зум и двойной клик (#CHART-1, toDataURL до/после). Сгенерирован по specs/SPEC-gs79_page_check.md.
async function gs79check() {
  const out = { src: [], err: [], langSwitch: null };
  window.addEventListener("error", e => out.err.push(String(e.message)));

  const wait = ms => new Promise(r => setTimeout(r, ms));
  const vis = el => !!(el && el.offsetParent);
  const blank = c => {
    const cv = document.createElement("canvas");
    cv.width = c.width;
    cv.height = c.height;
    return c.toDataURL() === cv.toDataURL();
  };
  const mouse = (type, target, x, y) => {
    target.dispatchEvent(new MouseEvent(type, {
      bubbles: true,
      cancelable: true,
      clientX: x,
      clientY: y,
      view: window
    }));
  };

  const buttons = Array.from(document.querySelectorAll("#srcbar .btn")).filter(b => !b.disabled);

  for (const b of buttons) {
    b.click();
    await wait(600);
    const rec = {
      src: b.textContent.trim(),
      hash: location.hash,
      tabs: [],
      zoom: [],
      pops: []
    };

    const tabs = Array.from(document.querySelectorAll(".tab[data-tab]")).filter(vis);
    for (const t of tabs) {
      t.click();
      await wait(500);
      const canvases = Array.from(document.querySelectorAll("canvas")).filter(vis);
      const blankCount = canvases.filter(blank).length;
      rec.tabs.push({
        tab: t.textContent.trim().slice(0, 40),
        canvases: canvases.length,
        blank: blankCount
      });

      const c = canvases.find(cv => {
        const r = cv.getBoundingClientRect();
        return r.width > 200 && r.height > 100;
      });

      if (c) {
        c.scrollIntoView({ block: "center" });
        await wait(200);
        const r = c.getBoundingClientRect();
        const y = r.top + r.height / 2;
        const before = c.toDataURL();

        mouse("mousedown", c, r.left + 0.35 * r.width, y);
        mouse("mousemove", document, r.left + 0.5 * r.width, y);
        mouse("mousemove", document, r.left + 0.6 * r.width, y);
        mouse("mouseup", document, r.left + 0.6 * r.width, y);

        await wait(300);
        const zoomed = c.toDataURL();

        mouse("dblclick", c, r.left + r.width / 2, y);
        await wait(300);

        rec.zoom.push({
          tab: t.dataset.tab,
          id: c.id || c.className,
          changed: zoomed !== before,
          restored: c.toDataURL() === before
        });
      }
    }

    const pops = Array.from(document.querySelectorAll("[data-pop]")).filter(vis);
    for (const p of pops) {
      p.click();
      await wait(300);
      const pop = Array.from(document.querySelectorAll(".pop")).find(el => !el.hidden);
      const h2 = pop ? pop.querySelector("h2") : null;
      const record = {
        pop: p.dataset.pop,
        open: !!pop,
        chars: pop ? pop.innerText.trim().length : 0,
        title: h2 ? h2.textContent.trim() : ""
      };
      document.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape", bubbles: true }));
      await wait(200);
      record.closed = !pop || pop.hidden;
      rec.pops.push(record);
    }

    out.src.push(rec);
  }

  const ls = document.querySelector(".lang-switch");
  if (ls) {
    const href = ls.getAttribute("href");
    out.langSwitch = {
      text: ls.textContent,
      href: href,
      resolved: new URL(href, location.href).href
    };
  } else {
    out.langSwitch = null;
  }

  return JSON.stringify(out);
}

await gs79check();
