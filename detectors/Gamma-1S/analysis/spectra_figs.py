# -*- coding: utf-8 -*-
"""Спектры комплекта для веб-страницы: калибровка, поиск пиков, деконволюция.

Три рисунка на каждую запись комплекта, и все три — из ОДНОГО прогона тех же
функций, которыми считаются публикуемые числа. Это не иллюстрации к результату,
а сам результат в другом виде: если картинка разойдётся с таблицей, значит одна
из них считана иначе, и это будет видно.

  1. СПЕКТР. Проба и её вложенный фон, приведённый к живому времени пробы, и
     разность. Логарифм по счёту: иначе кроме 2614,5 ничего не разглядеть.
     Вертикали — аналитические линии записи.

  2. КАЛИБРОВКА. Найденная центроида против табличной энергии по каждой линии
     (bm.peak_find — первый момент по вычтенной подложке), измеренная ПШПВ
     против калиброванного закона ПШПВ² = a + b·E (deconv.fwhm). Показывается
     ТОЛЬКО у чистых линий: у мультиплета «центроида» есть центр тяжести группы, а
     не калибровочный сдвиг, и путать эти две вещи нельзя — на этом уже один
     раз испортили результат (см. шапку kit_recalc).

  3. ДЕКОНВОЛЮЦИЯ. Область подгонки: точки измерения за вычетом фона, полная
     модель, континуум со ступенькой и вклад каждой линии группы по
     отдельности. Кривые берутся из тех же колонок матрицы плана, которыми
     решалась задача (deconv возвращает их в ключе fit).

Модуль ничего не печатает и не пишет — только отдаёт куски разметки; собирает
страницу build_web.py.
"""
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, HERE)
import becqmoni as bm  # noqa: E402
import deconv as dc  # noqa: E402
import kit_recalc as kr  # noqa: E402

W, H = 760, 250
PAD_L, PAD_R, PAD_T, PAD_B = 62, 12, 12, 34
DW, DH = 370, 210
DPAD_L, DPAD_R, DPAD_T, DPAD_B = 52, 8, 10, 30

E_LO, E_HI = 0.0, 3000.0


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def ru(x, nd=2):
    return ("%.*f" % (nd, x)).replace(".", ",")


# --- рисунок 1: спектр -------------------------------------------------------

def _sx(E, w=W, lo=E_LO, hi=E_HI, pl=PAD_L, pr=PAD_R):
    return pl + (E - lo) / (hi - lo) * (w - pl - pr)


def _sy(v, lo, hi, h=H, pt=PAD_T, pb=PAD_B):
    v = max(v, lo)
    a, b = math.log10(lo), math.log10(hi)
    return h - pb - (math.log10(v) - a) / (b - a) * (h - pt - pb)


def _bin(sp, n=760):
    """Спектр, сведённый к n точкам: (энергии, отсчёты/с)."""
    ch = np.arange(len(sp.n), dtype=float)
    en = sp.energy(ch)
    m = (en >= E_LO) & (en <= E_HI)
    en, y = en[m], sp.n[m] / sp.live
    if len(en) <= n:
        return en, y
    k = int(math.ceil(len(en) / n))
    cut = len(en) - len(en) % k
    return (en[:cut].reshape(-1, k).mean(1),
            y[:cut].reshape(-1, k).mean(1))


def spectrum_data(sp, bg, marks, found=None, n=900):
    """Данные спектра для интерактивного графика — JSON-совместимый словарь.

    Отдаём именно СВЕДЁННЫЙ спектр (n точек вместо тысяч каналов): страница
    самодостаточна и не тянет ничего извне, а лишний вес разметки тут ни к
    чему. Сведение — среднее по группе каналов, без сглаживания: сглаживать
    данные нельзя, это правило отображения.
    """
    ex, ey = _bin(sp, n)
    out = dict(E=[round(float(v), 1) for v in ex],
               smp=[round(float(v), 6) for v in ey],
               marks=[[round(float(E), 1), lab, 1 if clean else 0]
                      for E, lab, clean in marks])
    if bg is not None:
        bx, by = _bin(bg, n)
        b = np.interp(ex, bx, by)
        out["bg"] = [round(float(v), 6) for v in b]
        out["net"] = [round(float(v), 6) for v in (ey - b)]
    if found:
        out["peaks"] = [[round(float(f["E"]), 1), round(float(f["sig"]), 1)]
                        for f in found]
    return out


# Один скрипт на все графики страницы: подставляется в разметку однажды.
SPECTRA_JS = r"""
<script>
(function(){
 var W=760,H=260,PL=62,PR=12,PT=12,PB=34;
 function fmt(v){var t=v>=100?v.toFixed(0):v>=1?v.toFixed(2):v.toExponential(1);return t.replace('.',',');}
 function draw(box){
  var d=box._d,st=box._st,sv=box.querySelector('svg');
  var lo=st.lo,hi=st.hi,E=d.E;
  var i0=0,i1=E.length-1;
  while(i0<i1&&E[i0]<lo)i0++; while(i1>i0&&E[i1]>hi)i1--;
  i0=Math.max(0,i0-1); i1=Math.min(E.length-1,i1+1);
  // пределы по вертикали — по видимому куску, иначе зум бесполезен
  var ymax=1e-12,ymin=Infinity,k,ser=[];
  if(st.smp)ser.push('smp'); if(st.bg&&d.bg)ser.push('bg'); if(st.net&&d.net)ser.push('net');
  for(var s=0;s<ser.length;s++){var a=d[ser[s]];
   for(k=i0;k<=i1;k++){var v=a[k]; if(v>ymax)ymax=v; if(v>0&&v<ymin)ymin=v;}}
  if(!isFinite(ymin)||ymin<=0)ymin=ymax*1e-5;
  if(st.log){ymin=Math.max(ymin,ymax*1e-6);}else{ymin=0;}
  // запас сверху, чтобы пик и метки найденных пиков не упирались в верхнюю
  // подпись оси (замечание оператора): на лог-шкале — почти декада, на
  // линейной — четверть высоты
  ymax*=st.log?3.2:1.35;
  function X(e){return PL+(e-lo)/(hi-lo)*(W-PL-PR);}
  function Y(v){
   if(!st.log)return H-PB-(v-ymin)/(ymax-ymin)*(H-PT-PB);
   v=Math.max(v,ymin);
   return H-PB-(Math.log(v)-Math.log(ymin))/(Math.log(ymax)-Math.log(ymin))*(H-PT-PB);}
  // сетка по энергии — вдвое плотнее прежнего (замечание оператора: крупная):
  // цель 8…16 крупных делений с подписями + мелкие штрихи между ними
  var o=[],step=Math.pow(10,Math.floor(Math.log10(hi-lo)))/2;
  if((hi-lo)/step>16)step*=2; if((hi-lo)/step<8)step/=2;
  var minor=step/5;  // мелкие деления между подписанными
  for(var em=Math.ceil(lo/minor)*minor;em<=hi;em+=minor){
   var xm=X(em);
   o.push('<line class="tick" x1="'+xm.toFixed(1)+'" y1="'+(H-PB)+'" x2="'+xm.toFixed(1)+'" y2="'+(H-PB+4)+'"/>');}
  for(var e=Math.ceil(lo/step)*step;e<=hi;e+=step){
   o.push('<line class="grid" x1="'+X(e).toFixed(1)+'" y1="'+PT+'" x2="'+X(e).toFixed(1)+'" y2="'+(H-PB)+'"/>');
   o.push('<line class="tick" x1="'+X(e).toFixed(1)+'" y1="'+(H-PB)+'" x2="'+X(e).toFixed(1)+'" y2="'+(H-PB+6)+'"/>');
   o.push('<text class="ax" x="'+X(e).toFixed(1)+'" y="'+(H-PB+16)+'" text-anchor="middle">'+Math.round(e)+'</text>');}
  if(st.log){for(var dec=Math.floor(Math.log10(ymin));dec<=Math.ceil(Math.log10(ymax));dec++){
    var v=Math.pow(10,dec); if(v>=ymin&&v<=ymax){var yy=Y(v);
     o.push('<line class="grid" x1="'+PL+'" y1="'+yy.toFixed(1)+'" x2="'+(W-PR)+'" y2="'+yy.toFixed(1)+'"/>');
     o.push('<line class="tick" x1="'+(PL-6)+'" y1="'+yy.toFixed(1)+'" x2="'+PL+'" y2="'+yy.toFixed(1)+'"/>');
     o.push('<text class="ax" x="'+(PL-9)+'" y="'+(yy+4).toFixed(1)+'" text-anchor="end">10<tspan dy="-4" font-size="7">'+dec+'</tspan></text>');}
    for(var mm=2;mm<=9;mm++){var vv2=mm*Math.pow(10,dec); if(vv2<ymin||vv2>ymax)continue; var y3=Y(vv2);
     o.push('<line class="tick" x1="'+(PL-3)+'" y1="'+y3.toFixed(1)+'" x2="'+PL+'" y2="'+y3.toFixed(1)+'"/>');}}}
  else{for(var q=0;q<=8;q++){var v2=ymin+(ymax-ymin)*q/8,y2=Y(v2);
    if(q%2===0){o.push('<line class="grid" x1="'+PL+'" y1="'+y2.toFixed(1)+'" x2="'+(W-PR)+'" y2="'+y2.toFixed(1)+'"/>');
     o.push('<text class="ax" x="'+(PL-9)+'" y="'+(y2+4).toFixed(1)+'" text-anchor="end">'+fmt(v2)+'</text>');}
    o.push('<line class="tick" x1="'+(PL-(q%2?3:6))+'" y1="'+y2.toFixed(1)+'" x2="'+PL+'" y2="'+y2.toFixed(1)+'"/>');}}
  // сами оси — сплошными линиями, темнее сетки
  o.push('<line class="axis" x1="'+PL+'" y1="'+(H-PB)+'" x2="'+(W-PR)+'" y2="'+(H-PB)+'"/>');
  o.push('<line class="axis" x1="'+PL+'" y1="'+PT+'" x2="'+PL+'" y2="'+(H-PB)+'"/>');
  o.push('<text class="axt" x="'+((PL+W-PR)/2)+'" y="'+(H-4)+'" text-anchor="middle">энергия, кэВ</text>');
  o.push('<text class="axt" transform="translate(12,'+((PT+H-PB)/2)+') rotate(-90)" text-anchor="middle">имп/с на канал</text>');
  // аналитические линии и найденные пики
  for(k=0;k<d.marks.length;k++){var m=d.marks[k]; if(m[0]<lo||m[0]>hi)continue;
   o.push('<line class="mark'+(m[2]?'':' dirty')+'" x1="'+X(m[0]).toFixed(1)+'" y1="'+PT+'" x2="'+X(m[0]).toFixed(1)+'" y2="'+(H-PB)+'"/>');
   o.push('<text class="mk" x="'+X(m[0]).toFixed(1)+'" y="'+(PT+9)+'" text-anchor="middle">'+m[1]+'</text>');}
  if(st.pk&&d.peaks)for(k=0;k<d.peaks.length;k++){var p=d.peaks[k]; if(p[0]<lo||p[0]>hi)continue;
   o.push('<polygon class="fpk" points="'+X(p[0]).toFixed(1)+','+(PT+13)+' '+(X(p[0])-4).toFixed(1)+','+(PT+5)+' '+(X(p[0])+4).toFixed(1)+','+(PT+5)+'"><title>найден пик '+String(p[0]).replace('.',',')+' кэВ, значимость '+String(p[1]).replace('.',',')+'</title></polygon>');}
  var cls={smp:'ssmp',bg:'sbg',net:'snet'};
  for(s=0;s<ser.length;s++){var arr=d[ser[s]],dd='',pen=false;
   for(k=i0;k<=i1;k++){var vv=arr[k];
    if(st.log&&vv<=0){pen=false;continue;}
    dd+=(pen?'L':'M')+X(E[k]).toFixed(1)+','+Y(vv).toFixed(1)+' ';pen=true;}
   if(dd)o.push('<path class="'+cls[ser[s]]+'" d="'+dd+'"/>');}
  sv.innerHTML=o.join('');
  box.querySelector('.si-range').textContent=Math.round(lo)+'…'+Math.round(hi)+' кэВ';
 }
 function init(box){
  var d=JSON.parse(box.querySelector('script[type="application/json"]').textContent);
  box._d=d; box._st={lo:d.E[0],hi:d.E[d.E.length-1],log:true,smp:true,
                     bg:!!d.bg,net:false,pk:!!d.peaks};
  var full=[d.E[0],d.E[d.E.length-1]];
  box.addEventListener('click',function(ev){
   var b=ev.target.closest('button'); if(!b)return;
   var a=b.dataset.act,st=box._st;
   if(a==='reset'){st.lo=full[0];st.hi=full[1];}
   else if(a==='log'){st.log=!st.log;}
   else {st[a]=!st[a];}
   box.querySelectorAll('button[data-act]').forEach(function(x){
    var k=x.dataset.act; if(k==='reset')return;
    x.classList.toggle('on',k==='log'?st.log:!!st[k]);});
   draw(box);});
  var sv=box.querySelector('svg');
  sv.addEventListener('wheel',function(ev){
   ev.preventDefault(); var st=box._st,r=sv.getBoundingClientRect();
   var f=(ev.clientX-r.left)/r.width*W; if(f<PL)f=PL; if(f>W-PR)f=W-PR;
   var e=st.lo+(f-PL)/(W-PL-PR)*(st.hi-st.lo);
   var z=ev.deltaY<0?0.8:1.25, w=(st.hi-st.lo)*z;
   if(w>full[1]-full[0])w=full[1]-full[0]; if(w<20)w=20;
   var t=(e-st.lo)/(st.hi-st.lo);
   st.lo=e-t*w; st.hi=st.lo+w;
   if(st.lo<full[0]){st.lo=full[0];st.hi=st.lo+w;}
   if(st.hi>full[1]){st.hi=full[1];st.lo=st.hi-w;}
   draw(box);},{passive:false});
  var drag=null;
  sv.addEventListener('mousedown',function(ev){drag={x:ev.clientX,lo:box._st.lo,hi:box._st.hi};});
  window.addEventListener('mouseup',function(){drag=null;});
  sv.addEventListener('mousemove',function(ev){
   var st=box._st,r=sv.getBoundingClientRect();
   if(drag){var dx=(ev.clientX-drag.x)/r.width*W;
    var de=-dx/(W-PL-PR)*(drag.hi-drag.lo);
    var lo=drag.lo+de,hi=drag.hi+de;
    if(lo<full[0]){hi+=full[0]-lo;lo=full[0];}
    if(hi>full[1]){lo-=hi-full[1];hi=full[1];}
    st.lo=lo;st.hi=hi;draw(box);return;}
   var f=(ev.clientX-r.left)/r.width*W;
   var e=st.lo+(f-PL)/(W-PL-PR)*(st.hi-st.lo);
   var d=box._d,i=0,best=1e18;
   for(var k=0;k<d.E.length;k++){var q=Math.abs(d.E[k]-e); if(q<best){best=q;i=k;}}
   var t=d.E[i]+' кэВ: проба '+fmt(d.smp[i]);
   if(d.bg)t+=', фон '+fmt(d.bg[i])+', разность '+fmt(d.net[i]);
   box.querySelector('.si-read').textContent=t;});
  draw(box);
 }
 document.querySelectorAll('.si').forEach(init);
})();
</script>"""


def spectrum_svg(sp, bg, marks):
    """Проба, фон и разность в одних осях. marks — [(E, подпись, чистая?)]."""
    ex, ey = _bin(sp)
    lo = max(1e-5, float(np.percentile(ey[ey > 0], 2)) if (ey > 0).any() else 1e-5)
    hi = float(ey.max()) * 1.6
    s = ['<svg viewBox="0 0 %d %d" class="plot" role="img" aria-label='
         '"Спектр записи, фон и разность">' % (W, H)]
    for E in range(0, 3001, 500):
        x = _sx(E)
        s.append('<line class="grid" x1="%.1f" y1="%d" x2="%.1f" y2="%d"/>'
                 % (x, PAD_T, x, H - PAD_B))
        s.append('<text class="ax" x="%.1f" y="%d" text-anchor="middle">%d'
                 '</text>' % (x, H - PAD_B + 15, E))
    d0, d1 = math.floor(math.log10(lo)), math.ceil(math.log10(hi))
    dec = d0
    while dec <= d1:
        v = 10.0 ** dec
        if lo <= v <= hi:
            y = _sy(v, lo, hi)
            s.append('<line class="grid" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                     % (PAD_L, y, W - PAD_R, y))
            s.append('<text class="ax" x="%d" y="%.1f" text-anchor="end">'
                     '10<tspan dy="-4" font-size="7">%d</tspan></text>'
                     % (PAD_L - 6, y + 4, dec))
        dec += 1
    s.append('<text class="axt" x="%.1f" y="%d" text-anchor="middle">энергия, '
             'кэВ</text>' % ((PAD_L + W - PAD_R) / 2, H - 4))
    s.append('<text class="axt" transform="translate(12,%.1f) rotate(-90)" '
             'text-anchor="middle">имп/с на канал</text>'
             % ((PAD_T + H - PAD_B) / 2))

    # аналитические линии записи
    for E, lab, clean in marks:
        if not (E_LO <= E <= E_HI):
            continue
        x = _sx(E)
        s.append('<line class="mark%s" x1="%.1f" y1="%d" x2="%.1f" y2="%d"/>'
                 % ("" if clean else " dirty", x, PAD_T, x, H - PAD_B))
        s.append('<text class="mk" x="%.1f" y="%d" text-anchor="middle">%s'
                 '</text>' % (x, PAD_T + 9, esc(lab)))

    def path(xs, ys, cls):
        d, pen = [], False
        for xv, yv in zip(xs, ys):
            if yv <= 0:
                pen = False
                continue
            d.append("%s%.1f,%.1f" % ("L" if pen else "M", _sx(xv),
                                      _sy(yv, lo, hi)))
            pen = True
        return ('<path class="%s" d="%s"/>' % (cls, " ".join(d))) if d else ""

    if bg is not None:
        bx, by = _bin(bg)
        s.append(path(bx, by, "sbg"))
        net = ey - np.interp(ex, bx, by)
        s.append(path(ex, net, "snet"))
    s.append(path(ex, ey, "ssmp"))
    s.append("</svg>")
    return "\n".join(s)


# --- рисунок 3: врезка деконволюции -----------------------------------------

def deconv_svg(res, E0, title):
    """Область подгонки: данные, полная модель, континуум, вклады линий."""
    f = res["fit"]
    x, y = f["x"], f["y"]
    lo, hi = float(x[0]), float(x[-1])
    ymax = max(float(y.max()), float(f["model"].max())) * 1.12
    ymin = min(0.0, float(y.min()))

    def px(E):
        return DPAD_L + (E - lo) / (hi - lo) * (DW - DPAD_L - DPAD_R)

    def py(v):
        return (DH - DPAD_B - (v - ymin) / (ymax - ymin)
                * (DH - DPAD_T - DPAD_B))

    s = ['<svg viewBox="0 0 %d %d" class="dplot" role="img" aria-label='
         '"Деконволюция группы %s">' % (DW, DH, esc(title))]
    step = 50 if (hi - lo) < 260 else 100
    e = step * math.ceil(lo / step)
    while e <= hi:
        s.append('<line class="grid" x1="%.1f" y1="%d" x2="%.1f" y2="%d"/>'
                 % (px(e), DPAD_T, px(e), DH - DPAD_B))
        s.append('<text class="ax" x="%.1f" y="%d" text-anchor="middle">%d'
                 '</text>' % (px(e), DH - DPAD_B + 14, e))
        e += step
    for k in range(5):
        v = ymin + (ymax - ymin) * k / 4
        s.append('<line class="grid" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (DPAD_L, py(v), DW - DPAD_R, py(v)))
        s.append('<text class="ax" x="%d" y="%.1f" text-anchor="end">%s</text>'
                 % (DPAD_L - 5, py(v) + 4,
                    (("%.0f" % v) if abs(v) >= 100
                     else ("%.3g" % v)).replace(".", ",")))
    s.append('<text class="axt" x="%.1f" y="%d" text-anchor="middle">энергия, '
             'кэВ</text>' % ((DPAD_L + DW - DPAD_R) / 2, DH - 3))

    def poly(vals, cls):
        d = " ".join("%s%.1f,%.1f" % ("L" if i else "M", px(xv), py(vv))
                     for i, (xv, vv) in enumerate(zip(x, vals)))
        return '<path class="%s" d="%s"/>' % (cls, d)

    # вклады отдельных линий — поверх континуума, чтобы читалась их доля
    cont = f["cont"]
    for E, I, g in res["fit"]["parts"]:
        s.append(poly(cont + g, "dcomp"))
    s.append(poly(cont, "dcont"))
    s.append(poly(f["model"], "dmodel"))
    for xv, vv in zip(x[::2], y[::2]):
        s.append('<circle class="ddata" cx="%.1f" cy="%.1f" r="1.5"/>'
                 % (px(xv), py(vv)))
    # подписи линий группы
    for E, I, _g in res["fit"]["parts"]:
        if not (lo <= E + f["shift"] <= hi):
            continue
        s.append('<text class="mk" x="%.1f" y="%d" text-anchor="middle">%s'
                 '</text>' % (px(E + f["shift"]), DPAD_T + 8, ru(E, 0)))
    s.append("</svg>")
    return "\n".join(s)


# --- сборка блока по записи --------------------------------------------------

def _record_files(geom, mask):
    kd = paths.kit_dir(geom)
    return sorted(str(p) for p in kd.rglob(mask)) if kd else []


def record_block(geom, mask, nuc, aspec, dpct, d0, mass, vol, geom_title):
    """Разметка одной записи: спектр, таблица калибровки, врезки деконволюции."""
    files = _record_files(geom, mask)
    if not files:
        return "нет спектра записи (маска %s)" % mask
    sp, bg, calib = bm.read_checked(files[0])
    txt = open(files[0], encoding="utf-8", errors="replace").read()
    md = re.search(r"<StartTime>(\d{4}-\d{2}-\d{2})", txt)
    md = md.group(1) if md else None
    A0 = aspec * mass / 1000.0 * kr.decay_factor(nuc, d0, md)
    rho = mass / vol
    R = float(sp.n.sum()) / sp.live
    pile = math.exp(2 * kr.TAU_SHAPE * R)
    lines, ckey = kr.VLINES[nuc]
    base = kr.RUNBASE.get((geom, ckey))
    if not base:
        # Спектр записи есть, а прогона распада в этой геометрии нет — считать
        # эффективность на распад нечем. Разница с «нет спектра» существенная,
        # и путать их нельзя: одно чинится поиском файла, другое — расчётом.
        return ("прогона распада %s в геометрии %s нет — эффективность на "
                "распад брать неоткуда" % (ckey, geom))

    # Поиск пиков — до всякого опознания и независимо от списка линий.
    # Именно он отвечает на вопрос «чего мы не видим», который по прежней
    # схеме нельзя было даже задать (см. шапку peaksearch.py).
    try:
        import peaksearch as ps
        psr = ps.analyze(sp, bg, base)
        found = psr["found"]
    except Exception as exc:                       # noqa: BLE001
        psr, found = None, None
        print("!! поиск пиков не отработал (%s): %s" % (geom, exc))

    marks, rows, panels = [], [], []
    for E in lines:
        # ВАЖНО, ДВА ЗАКОНА ПШПВ, И ПУТАТЬ ИХ НЕЛЬЗЯ.
        # Оконный съём (kit_recalc, публикуемые числа) работает по закону
        # корня от одной опорной точки 662 кэВ; деконволюция — по калибровке
        # ПШПВ² = a + b·E, потому что подгонке точность ширины важна. Здесь
        # столбец «окном» обязан ВОСПРОИЗВОДИТЬ публикуемое число, поэтому
        # окно и чистота считаются законом корня, а не калибровкой. Сначала
        # я взял здесь калибровку, и таблица страницы разошлась с
        # results/deconv_lines.csv — на Ra-226 351,9 до 12 % (причина ниже,
        # в тексте страницы: левая полка фона стоит на пике 295,2 кэВ).
        fw_win = kr.FWHM662 * math.sqrt(E / 661.657)
        fw_cal = dc.fwhm(E)
        frac, dirt = kr.purity(base, E, fw_win)
        clean = frac is not None and frac >= kr.CLEAN_FRAC
        marks.append((E, ru(E, 0), clean))
        pf = bm.peak_find(sp, E)
        cen = shift = fwm = None
        if pf:
            cen, fwm = pf[0], pf[1]
            shift = cen - E
        res = dc.deconvolve(sp, bg, base, E, geom=geom, rho_src=rho)
        A = res["A"] * pile if res else None
        nr = bm.net_rate(sp, bg, E, fw_win, roi=1.0, side=1.0)
        # Матрица источника — по правилу kit_recalc: у лёгких засыпок
        # (ρ ≤ 1,3) состав в файлах не записан, берётся вода.
        _k = min(kr.MU_O, key=lambda k: abs(k - E))
        eps = kr.eps_per_decay(geom, ckey, E, fw_win, rho,
                               kr.MU_O[_k] if rho > 1.3 else kr.MU_W[_k])
        win = (nr[0] * pile / eps) if (nr and eps and nr[0] > 0) else None
        rows.append(dict(E=E, frac=frac, clean=clean, cen=cen, shift=shift,
                         fwm=fwm, fwl=fw_win, fwc=fw_cal, A=A, win=win, A0=A0,
                         chi2=res["chi2"] if res else None,
                         nl=res["n_lines"] if res else 0,
                         dsh=res["shift"] if res else None,
                         dirt=dirt))
        if res:
            panels.append((E, res))

    h = ['<h4>%s — %s</h4>' % (esc(geom_title), esc(nuc))]
    h.append('<p class="cap">Живое время %s с, полный счёт %s имп/с, '
             'поправка на наложения %s. Паспорт на дату измерения %s Бк.</p>'
             % (ru(sp.live, 0), ru(R, 0), ru(pile, 3), ru(A0, 0)))
    # Калибровка ФОНА проверяется отдельно от калибровки пробы — правило ЛСРМ.
    if calib and calib.get("reason"):
        sh = calib.get("shifts") or []
        h.append('<p class="cap">Калибровка фона (проверяется независимо от '
                 'пробы): %s.%s</p>'
                 % (esc(calib["reason"]).replace(".", ","),
                    "" if not sh else
                    " Невязки якорных линий (ярких одиночных пиков фона, по "
                    "которым проверяется калибровка): " + ", ".join(
                        ("%s кэВ %+.2f ПШПВ" % (ru(E, 1), s)).replace(".", ",")
                        for E, s in sh)
                    + "."))
    import json
    data = spectrum_data(sp, bg, marks, found)
    uid = "sp%d" % (abs(hash((geom, nuc))) % 100000)
    h.append(
        '<div class="si" id="%s">'
        '<div class="si-bar">'
        '<button data-act="log" class="on">лог</button>'
        '<button data-act="smp" class="on">проба</button>'
        '<button data-act="bg" class="on">фон</button>'
        '<button data-act="net">разность</button>'
        '<button data-act="pk" class="on">найденные пики</button>'
        '<button data-act="reset">весь диапазон</button>'
        '<span class="si-range"></span></div>'
        '<svg viewBox="0 0 760 260" class="plot" role="img" aria-label='
        '"Спектр записи: проба, фон и разность"></svg>'
        '<div class="si-read">наведите курсор на график</div>'
        '<script type="application/json">%s</script></div>'
        % (uid, json.dumps(data, ensure_ascii=False)))
    h.append('<p class="cap"><b>Рисунок 2 — спектр поверочной записи: проба, '
             'фон и их разность (для выбранной записи).</b> Синяя линия — проба, '
             'оранжевая — фон, приведённый к живому времени пробы, розовая — '
             'разность; штриховые вертикали — линии в мультиплете, треугольники — '
             'найденные пики. Оси: энергия (кэВ) и скорость счёта на канал '
             '(имп/с). Управление: колесо — масштаб, перетаскивание — сдвиг.</p>')

    h.append('<p class="cap"><b>Таблица 1 — разметка линий записи: чистота, '
             'центроида, ПШПВ и активность (для выбранной записи).</b> «чистота» — '
             'доля выхода окна, приходящаяся на саму линию (порог годности для '
             'оконного съёма 0,95); «центроида» — найденный центр линии; «сдвиг» — '
             'центроида минус табличная энергия; «ПШПВ изм./&radic;E/калибр.» — '
             'полуширина измеренная, теоретическая (&prop;&radic;E) и по калибровке '
             'прибора; «A/пасп деконв./окном» — активность к паспорту двумя '
             'способами съёма площади.</p>')
    h.append('<div class="tw"><table><thead><tr>'
             '<th class="n">линия, кэВ</th><th class="n">чистота</th>'
             '<th class="n">центроида</th><th class="n">сдвиг, кэВ</th>'
             '<th class="n">ПШПВ изм.</th><th class="n">ПШПВ √E</th>'
             '<th class="n">ПШПВ калибр.</th>'
             '<th class="n">линий<br>в группе</th><th class="n">χ²/dof</th>'
             '<th class="n">A/пасп<br>деконв.</th>'
             '<th class="n">A/пасп<br>окном</th></tr></thead><tbody>')
    for r in rows:
        cl = "" if r["clean"] else ' class="dim"'
        h.append('<tr%s><td class="n">%s</td><td class="n">%s</td>'
                 '<td class="n">%s</td><td class="n">%s</td>'
                 '<td class="n">%s</td><td class="n">%s</td>'
                 '<td class="n">%s</td>'
                 '<td class="n">%d</td><td class="n">%s</td>'
                 '<td class="n">%s</td><td class="n">%s</td></tr>'
                 % (cl, ru(r["E"], 1),
                    ru(r["frac"], 2) if r["frac"] is not None else "—",
                    ru(r["cen"], 1) if r["cen"] else "—",
                    ("%+.1f" % r["shift"]).replace(".", ",")
                    if r["shift"] is not None else "—",
                    ru(r["fwm"], 1) if r["fwm"] else "—",
                    ru(r["fwl"], 1), ru(r["fwc"], 1), r["nl"],
                    ru(r["chi2"], 2) if r["chi2"] is not None else "—",
                    ru(r["A"] / r["A0"], 3) if r["A"] else "—",
                    ru(r["win"] / r["A0"], 3) if r["win"] else "—"))
    h.append("</tbody></table></div>")
    bad = [r for r in rows if not r["clean"]]
    if bad:
        h.append('<p class="cap">Строки серым — линии в мультиплете (чистота ниже '
                 '%s): их «центроида» есть центр тяжести группы, а не '
                 'калибровочный сдвиг. Для восстановления активности '
                 'оконный съём такую линию не использует, связанная '
                 'деконволюция — использует; на этом и основано сравнение '
                 'двух правых столбцов.</p>'
                 % ru(kr.CLEAN_FRAC, 2))
    if psr:
        h.append('<p class="cap">Поиск пиков (фильтр Марискотти, порог '
                 '%.0f&sigma;): найдено <b>%d</b>, опознано по спектру '
                 'испускания прогона <b>%d</b>; <b>%d</b> линий слито с '
                 'соседями из-за конечного разрешения детектора — такие '
                 'группы разделяются деконволюцией; <b>%d</b> не видно вовсе '
                 'в рабочем диапазоне '
                 '%.0f…%.0f кэВ. Худшая невязка калибровки %s ПШПВ — %s.</p>'
                 % (ps.SIGMA_THR, len(psr["found"]), len(psr["pairs"]),
                    len(psr["merged"]), len(psr["unseen"]),
                    psr["span"][0], psr["span"][1], ru(psr["cal_worst"], 2),
                    "пересчёт по правилу ЛСРМ НЕ требуется"
                    if not psr["cal_need"] else
                    "ПРЕВЫШЕН порог 0,3 ПШПВ, калибровку надо пересчитать"))
        if psr["unseen"]:
            h.append('<p class="cap">Сильнейшие из невидимых: %s.</p>'
                     % ", ".join("%s кэВ (выход %s)" % (ru(E, 1), ru(y, 3))
                                 for E, y in psr["unseen"][:5]))
        if psr["extra"]:
            h.append('<p class="cap">Пики без линии в библиотеке прогона '
                     '(%d): %s. Это не обязательно чужой нуклид — вторичные '
                     'структуры (вылеты, обратное рассеяние, сумм-пики) дают '
                     'то же самое, и разбирать их надо отдельно.</p>'
                     % (len(psr["extra"]),
                        ", ".join("%s кэВ" % ru(f["E"], 0) for f in
                                  sorted(psr["extra"],
                                         key=lambda f: -f["sig"])[:5])))
    if panels:
        h.append('<div class="dgrid">')
        for E, res in panels:
            nlg = res["n_lines"]
            wl = ("линия" if nlg % 10 == 1 and nlg % 100 != 11
                  else "линии" if 2 <= nlg % 10 <= 4 and not 12 <= nlg % 100 <= 14
                  else "линий")
            kind = ("мультиплета" if nlg > 1 else "одиночной линии")
            h.append('<figure>%s<figcaption class="cap"><b>Рисунок 3 — '
                     'деконволюция %s %s кэВ (для выбранной записи).</b> '
                     '%d %s в группе, сдвиг %s кэВ, χ²/dof %s.</figcaption></figure>'
                     % (deconv_svg(res, E, "%s %s" % (nuc, ru(E, 0))),
                        kind, ru(E, 1), nlg, wl,
                        ("%+.1f" % res["shift"]).replace(".", ","),
                        ru(res["chi2"], 2)))
        h.append("</div>")
    return "\n".join(h), rows


DECONV_LEGEND = (
    '<p class="leg"><span class="k"><i style="border-color:var(--ink)">'
    '</i>измерение за вычетом фона</span>'
    '<span class="k"><i style="border-color:var(--mc)"></i>полная модель'
    '</span><span class="k"><i style="border-color:var(--exp)"></i>континуум '
    'со ступенькой</span><span class="k"><i style="border-color:var(--corr);'
    'border-top-style:dotted"></i>вклад отдельной линии</span></p>')
