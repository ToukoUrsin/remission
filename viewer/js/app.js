/* Remission replay viewer. Plain DOM + canvas; data arrives through <script> tags (see inflate.js). */
(function () {
  'use strict';

  var D = window.REMISSION_DATA, IDX = D && D.idx;
  var $ = function (id) { return document.getElementById(id); };
  var DAYS = 364, BASE_SPEED = 14;               // days per second at 1x
  var SPEEDS = [0.5, 1, 2, 4];
  var NNBSP = ' ', NBSP = ' ';

  if (!IDX || !IDX.patients || !IDX.patients.length) {
    $('replay').hidden = true;
    var f = $('fatal');
    f.hidden = false;
    f.innerHTML = 'No run data found. Build it first:<br><code>python3 scripts/build_viewer.py "Reference scientist=JOB_DIR" "Agent name=JOB_DIR"</code><br>then reload this page.';
    return;
  }
  var LAB = IDX.lab, DRUGS = IDX.drugs;

  // ---------------------------------------------------------------- helpers
  function esc(s) {
    return String(s == null ? '' : s).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }
  function el(tag, cls, text) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text != null) e.textContent = text;
    return e;
  }
  function clamp(v, a, b) { return v < a ? a : v > b ? b : v; }
  function int(n) { return String(Math.round(n)).replace(/\B(?=(\d{3})+(?!\d))/g, NNBSP); }
  function cells(n) { return n >= 10 ? int(n) : n >= 0.05 ? n.toFixed(1) : '0'; }
  function pct(x, d) { var v = 100 * x; return (d != null ? v.toFixed(d) : v > 0.05 && v < 10 ? v.toFixed(1) : v.toFixed(0)) + NBSP + '%'; }
  function g(n) { return String(+(+n).toPrecision(4)); }
  function cssPx(name) { return parseFloat(getComputedStyle(document.documentElement).getPropertyValue(name)) || 0; }
  function shortClone(label) { return label.replace(/ clone$/, ''); }
  var ICON = {
    bad: '<svg class="mark mark-bad" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="none" stroke="currentColor" stroke-width="1.6"/><path d="M5.2 5.2l5.6 5.6M10.8 5.2l-5.6 5.6" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"/></svg>',
    good: '<svg class="mark mark-good" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="none" stroke="currentColor" stroke-width="1.6"/><path d="M4.6 8.3l2.3 2.3 4.5-4.9" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/></svg>',
    note: '<svg class="mark mark-note" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="none" stroke="currentColor" stroke-width="1.6"/><circle cx="8" cy="8" r="2" fill="currentColor"/></svg>'
  };

  // ---------------------------------------------------------------- state
  var S = { view: 'replay', patient: null, left: null, right: null, day: 0, playing: false, speed: 1 };
  var SIDE = {
    L: { key: 'L', run: null, fi: 0 },
    R: { key: 'R', run: null, fi: 0 }
  };
  var loadToken = 0, scale = 8;

  function patient(id) { return IDX.patients.filter(function (p) { return p.id === id; })[0]; }
  function options(p) {
    var o = p.runs.filter(function (r) { return !r.missing; });
    return p.best ? o.concat([p.best]) : o;
  }
  function entry(p, id) { return options(p).filter(function (r) { return r.id === id; })[0]; }
  function reference(p) {
    return p.runs.filter(function (r) { return r.kind === 'reference' && !r.missing; })[0] || p.best || null;
  }
  function kindText(r) {
    return r.kind === 'agent' ? 'AI agent' + (r.model ? ' · ' + r.model.split('/').pop() : '')
      : r.kind === 'reference' ? 'reference run: a scripted scientist using the same tools'
      : r.kind === 'best' ? 'searched with the hidden parameters' : 'baseline policy';
  }

  /* The default view: the most instructive agent failure next to the reference on the same patient. */
  function defaults() {
    var best = null;
    IDX.patients.forEach(function (p) {
      var ref = reference(p);
      p.runs.forEach(function (r) {
        if (r.missing || r.kind === 'reference') return;
        var gap = (ref ? ref.reward : 0.7) - r.reward;
        var score = (r.kind === 'agent' ? 100 : 0) + gap * 50 + (r.submitted ? 10 : 0) + (r.n_experiments > 0 ? 10 : 0)
          + (/regrew/.test(r.tag || '') ? 8 : 0) - (/toxic/.test(r.tag || '') ? 4 : 0);
        if (!best || score > best.score) best = { score: score, p: p, r: r, ref: ref };
      });
    });
    if (!best) {
      var p0 = IDX.patients[0], o = options(p0);
      return { patient: p0.id, left: o[0].id, right: (reference(p0) || o[o.length - 1]).id };
    }
    return { patient: best.p.id, left: best.r.id, right: (best.ref || best.r).id };
  }

  function readHash() {
    var out = {};
    location.hash.replace(/^#/, '').split('&').forEach(function (kv) {
      var i = kv.indexOf('=');
      if (i > 0) out[decodeURIComponent(kv.slice(0, i))] = decodeURIComponent(kv.slice(i + 1));
    });
    return out;
  }
  var lastHash = '';
  function writeHash() {
    var h = S.view === 'overview' ? '#view=overview'
      : '#patient=' + encodeURIComponent(S.patient) + '&left=' + encodeURIComponent(S.left) + '&right=' + encodeURIComponent(S.right)
        + '&day=' + Math.round(S.day);
    if (h === lastHash) return;
    lastHash = h;
    try { history.replaceState(null, '', h); } catch (e) { location.hash = h; }
  }

  // ---------------------------------------------------------------- lattice
  var RGB = {
    0: [10, 10, 9], 1: [42, 148, 145], 2: [240, 67, 122], 3: [245, 183, 53], 4: [194, 162, 242]
  };
  var JIT = { 1: [0.80, 0.90, 1.0, 1.09], 2: [0.90, 0.96, 1.0, 1.06], 3: [0.92, 0.97, 1.0, 1.05], 4: [0.92, 0.97, 1.0, 1.05] };
  var SHADE = [];                                  // SHADE[state][variant] -> uint32 (little-endian RGBA)
  (function () {
    for (var s = 0; s < 5; s++) {
      SHADE[s] = [];
      for (var v = 0; v < 4; v++) {
        var m = s ? JIT[s][v] : 1, c = RGB[s];
        var r = Math.min(255, Math.round(c[0] * m)), gg = Math.min(255, Math.round(c[1] * m)), b = Math.min(255, Math.round(c[2] * m));
        SHADE[s][v] = (255 << 24 | b << 16 | gg << 8 | r) >>> 0;
      }
    }
  })();
  function hash2(i, k) {
    var h = Math.imul(i + 1, 73856093) ^ Math.imul(k + 1, 19349663);
    h = Math.imul(h ^ (h >>> 13), 1274126177);
    return ((h ^ (h >>> 16)) >>> 0) / 4294967296;
  }

  function initLattice(side) {
    var cv = side.canvas, dpr = window.devicePixelRatio || 1;
    var P = Math.max(1, Math.round(scale * dpr)), n = 64 * P;
    side.P = P;
    side.gap = P >= 5 ? Math.max(1, Math.round(P / 8)) : 0;
    cv.width = cv.height = n;
    side.ctx = cv.getContext('2d');
    side.img = side.ctx.createImageData(n, n);
    side.px = new Uint32Array(side.img.data.buffer);
    side.px.fill(SHADE[0][0]);
    side.shown = new Uint8Array(4096).fill(255);
    side.count = [0, 0, 0, 0, 0];
  }

  function paint(side, day) {
    var run = side.run;
    if (!run || !side.px) return;
    var days = run.replay.days, F = run.frames, k = side.fi, last = days.length - 1;
    if (k > last - 1) k = Math.max(0, last - 1);
    while (k > 0 && days[k] > day) k--;
    while (k < last - 1 && days[k + 1] <= day) k++;
    side.fi = k;
    var a = F[k], b = F[Math.min(k + 1, last)];
    var frac = last > 0 ? clamp((day - days[k]) / (days[k + 1] - days[k]), 0, 1) : 0;
    var px = side.px, shown = side.shown, P = side.P, w = P - side.gap, N = 64 * P;
    var count = [0, 0, 0, 0, 0], dirty = false;
    for (var i = 0; i < 4096; i++) {
      var v = a[i];
      if (v !== b[i] && hash2(i, k) < frac) v = b[i];
      count[v]++;
      if (v === shown[i]) continue;
      shown[i] = v; dirty = true;
      var col = SHADE[v][(Math.imul(i, 2654435761) >>> 28) & 3];
      var x0 = (i & 63) * P, row = (i >> 6) * P * N;
      for (var y = 0; y < w; y++) {
        var o = row + y * N + x0;
        for (var x = 0; x < w; x++) px[o + x] = col;
      }
    }
    side.count = count;
    if (dirty) side.ctx.putImageData(side.img, 0, 0);
  }

  // ---------------------------------------------------------------- clock
  var lastTs = 0, rafId = 0;
  function setDay(d, quiet) {
    S.day = clamp(d, 0, DAYS);
    paint(SIDE.L, S.day);
    paint(SIDE.R, S.day);
    var di = Math.floor(S.day + 1e-6);
    if (di !== setDay.shown) {
      setDay.shown = di;
      $('day-num').textContent = di;
      $('day-phase').textContent = di < LAB.treat_days ? 'treatment window' : 'drug-free follow-up';
      $('year-hit').setAttribute('aria-valuenow', di);
    }
    liveCounts(SIDE.L); liveCounts(SIDE.R);
    var hit = $('year-hit');
    $('playhead').style.transform = 'translateX(' + (S.day / DAYS * hit.clientWidth) + 'px)';
    if (!quiet) writeHash();
  }
  function frame(ts) {
    if (!S.playing) return;
    var dt = Math.min(ts - lastTs, 100) / 1000;
    lastTs = ts;
    var d = S.day + dt * BASE_SPEED * S.speed;
    if (d >= DAYS) { setDay(DAYS, true); pause(); return; }
    setDay(d, true);
    rafId = requestAnimationFrame(frame);
  }
  function play() {
    if (S.playing) return;
    if (S.day >= DAYS - 0.5) setDay(0, true);
    S.playing = true;
    $('play').classList.add('on');
    $('play').setAttribute('aria-label', 'Pause');
    lastTs = performance.now();
    rafId = requestAnimationFrame(frame);
  }
  function pause() {
    S.playing = false;
    cancelAnimationFrame(rafId);
    $('play').classList.remove('on');
    $('play').setAttribute('aria-label', 'Play');
    writeHash();
  }
  function liveCounts(side) {
    if (!side.countEls || !side.run) return;
    var c = side.count;
    for (var i = 0; i < side.countEls.length; i++) {
      var v = c[side.countEls[i].state];
      if (side.countEls[i].last !== v) { side.countEls[i].last = v; side.countEls[i].el.textContent = int(v); }
    }
  }

  // ---------------------------------------------------------------- top bar, patient line, pickers
  function renderTop() {
    var box = $('patients');
    box.textContent = '';
    if (IDX.patients.length > 8) {
      var sel = el('select');
      sel.setAttribute('aria-label', 'Patient');
      sel.name = 'patient';
      IDX.patients.forEach(function (p) {
        var o = el('option', null, p.id + ' · ' + p.chart.tissue);
        o.value = p.id; o.selected = p.id === S.patient;
        sel.appendChild(o);
      });
      sel.onchange = function () { choosePatient(sel.value); };
      box.appendChild(sel);
    } else {
      IDX.patients.forEach(function (p) {
        var b = el('button', null, 'Patient ' + p.id);
        b.type = 'button';
        b.setAttribute('aria-pressed', String(p.id === S.patient && S.view === 'replay'));
        b.title = p.chart.tissue + ' tumour';
        b.onclick = function () { choosePatient(p.id); };
        box.appendChild(b);
      });
    }
    [].forEach.call(document.querySelectorAll('.views button'), function (b) {
      b.setAttribute('aria-pressed', String(b.dataset.view === S.view));
    });
    var p = patient(S.patient), leg = document.querySelector('.legend');
    leg.textContent = '';
    function item(cls, text) { var li = el('li'); li.appendChild(el('i', 'sw ' + cls)); li.appendChild(document.createTextNode(text)); leg.appendChild(li); }
    item('sw-h', 'healthy cell');
    item('sw-c0', 'tumour: dominant clone');
    if (S.view === 'overview') item('sw-c1', 'tumour: hidden resistant clone');
    else p.clones.slice(1).forEach(function (c, i) { item('sw-c' + (i + 1), 'tumour: ' + c.label); });
    item('sw-e', 'empty site');
  }

  function renderPatientLine() {
    var p = patient(S.patient), line = $('patient-line');
    var chart = p.chart.bulk_sequencing.alterations.join(', ');
    var html = '<b>Patient ' + esc(p.id) + '</b> · ' + esc(p.chart.tissue) + ' tumour'
      + '<span class="sep">|</span>On the chart the agent gets: <b>' + esc(chart) + '</b>';
    var minors = p.clones.slice(1);
    if (minors.length) {
      html += '<span class="sep">|</span>Not on the chart: ' + minors.map(function (c, i) {
        return '<i class="sw sw-c' + (i + 1) + '"></i><span class="hid">a ' + esc(c.label) + ', ' + esc(pct(c.share)) + ' of the tumour, already resistant</span>';
      }).join(' and ');
    }
    line.innerHTML = html;
  }

  function renderPick(side) {
    var p = patient(S.patient), box = $('pick-' + side.key), cur = side.key === 'L' ? S.left : S.right;
    box.textContent = '';
    var sel = el('select');
    sel.setAttribute('aria-label', side.key === 'L' ? 'Run shown on the left' : 'Run shown on the right');
    sel.name = 'run-' + side.key;
    var groups = [['agent', 'AI agents'], ['reference', 'Reference'], ['baseline', 'Baselines'], ['best', 'Yardstick']];
    groups.forEach(function (gr) {
      var rs = options(p).filter(function (r) { return r.kind === gr[0]; });
      if (!rs.length) return;
      var og = el('optgroup');
      og.label = gr[1];
      rs.forEach(function (r) {
        var o = el('option', null, r.label + '  ·  ' + r.reward.toFixed(2));
        o.value = r.id; o.selected = r.id === cur;
        og.appendChild(o);
      });
      sel.appendChild(og);
    });
    sel.onchange = function () {
      if (side.key === 'L') S.left = sel.value; else S.right = sel.value;
      loadSides(false);
    };
    box.appendChild(sel);
    var e = entry(p, cur);
    box.appendChild(el('small', null, e ? kindText(e) : ''));
    // the visible option text carries the reward; the closed control should read as a title
    if (e) sel.options[sel.selectedIndex].textContent = e.label;
    sel.onfocus = sel.onmousedown = function () {
      [].forEach.call(sel.options, function (o) { var r = entry(p, o.value); o.textContent = r.label + '  ·  ' + r.reward.toFixed(2); });
    };
    sel.onblur = function () { var r = entry(p, sel.value); if (r) sel.options[sel.selectedIndex].textContent = r.label; };
  }

  // ---------------------------------------------------------------- caption under each lattice
  function renderCaption(side) {
    var run = side.run, p = patient(S.patient), dish = $('dish-' + side.key), cap = dish.querySelector('figcaption');
    var sc = run.score, ok = run.diag.success;
    dish.classList.toggle('fail', !ok);
    cap.textContent = '';
    var counts = el('div', 'counts');
    side.countEls = [];
    function cnt(state, cls, text) {
      var s = el('span'), b = el('b', null, '0');
      s.appendChild(el('i', 'sw ' + cls)); s.appendChild(document.createTextNode(text)); s.appendChild(b);
      counts.appendChild(s);
      side.countEls.push({ state: state, el: b, last: -1 });
    }
    cnt(1, 'sw-h', 'healthy');
    cnt(2, 'sw-c0', 'dominant');
    p.clones.slice(1).forEach(function (c, i) { cnt(3 + i, 'sw-c' + (i + 1), scale < 7 ? shortClone(c.label).replace('-amplified', '') : shortClone(c.label)); });
    cap.appendChild(counts);

    var out = el('div', 'outcome');
    var rw = el('div', 'reward', sc.reward.toFixed(2));
    rw.appendChild(el('small', null, 'reward'));
    out.appendChild(rw);
    var right = el('div');
    var v = el('div', 'verdict' + (run.diag.verdict.length > 62 ? ' long' : ''));
    v.innerHTML = ICON[ok ? 'good' : 'bad'] + esc(run.diag.verdict);
    right.appendChild(v);
    var eq = el('div', 'eq');
    var reps = (run.replay.reward_reps || []).length;
    eq.innerHTML = 'healed <b>' + sc.healed.toFixed(2) + '</b> − harm <b>' + sc.harm.toFixed(2) + '</b>'
      + (reps && sc.cured != null ? ' · tumour-free in <b>' + Math.round(sc.cured * reps) + '</b> of ' + reps + ' lesions' : '');
    right.appendChild(eq);
    out.appendChild(right);
    cap.appendChild(out);
  }

  // ---------------------------------------------------------------- flank: lab timeline, plate, diagnosis
  function srcClass(run, cellsName) {
    var ln = run.lines && run.lines[cellsName];
    return !ln ? 's-tumour' : ln.kind === 'healthy' ? 's-healthy' : ln.kind === 'derived' ? 's-line' : 's-tumour';
  }
  function srcName(run, cellsName) {
    var ln = run.lines && run.lines[cellsName];
    return !ln ? cellsName : ln.kind === 'healthy' ? 'healthy cells' : ln.kind === 'derived' ? 'line ' + cellsName : 'tumour biopsy';
  }
  var ROWS = [['sequencing', 'Sequencing'], ['phospho', 'Signalling'], ['viability', 'Viability'], ['derive', 'New line']];

  function labEvents(run) {
    var ev = (run.experiments || []).map(function (e, i) {
      return { type: e.kind, id: e.id, start: e.start, end: e.end, cells: e.cells, cls: srcClass(run, e.cells),
               n: e.wells ? e.wells.length : (e.n_wells || 16), exp: e, order: i };
    });
    (run.derives || []).forEach(function (d, i) {
      ev.push({ type: 'derive', id: d.line || 'no line', start: d.day, end: d.day + (d.ok ? 3 : 0), cells: d.line, cls: 's-line' + (d.ok ? '' : ' fail'), n: 8,
                derive: d, order: 1000 + i });
    });
    ev.sort(function (a, b) { return a.start - b.start || a.order - b.order; });
    return ev;
  }

  function renderFlank(side) {
    var run = side.run, box = $('flank-' + side.key);
    box.textContent = '';
    box.classList.remove('dense');
    closePop();
    if (run.kind === 'best') {
      var h = el('h2', null, 'No lab work');
      box.appendChild(h);
      box.appendChild(el('p', 'nolab', 'This regimen was searched with the patient’s hidden parameters, which no run ever sees. It shows what a cure looks like on this tissue.'));
    } else {
      renderLab(side, box);
    }
    var diag = el('div', 'diag');
    var dh = el('h2', null, run.kind === 'best' ? 'What it is' : run.diag.success ? 'What made it work' : 'Where it went wrong');
    diag.appendChild(dh);
    if (run.diag.source === 'analysis' && run.diag.summary && (run.diag.callouts || []).length) diag.appendChild(el('p', 'nolab', run.diag.summary));
    (run.diag.callouts || []).forEach(function (c) {
      var jump = c.day != null;
      var b = el(jump ? 'button' : 'div', 'callout');
      if (jump) { b.type = 'button'; b.title = 'Go to day ' + c.day; }
      b.innerHTML = ICON[c.tone] || ICON.bad;
      var t = el('div', 'c-title', c.title);
      b.appendChild(t);
      var d = el('div', 'c-detail', c.detail);
      if (jump) { d.appendChild(document.createTextNode(' ')); d.appendChild(el('span', 'c-jump', '→ day ' + c.day)); }
      b.appendChild(d);
      if (jump) b.onclick = function () { pause(); setDay(c.day); };
      if (c.exp && c.exp.length) {
        b.onmouseenter = b.onfocus = function () { highlight(side, c.exp); };
        b.onmouseleave = b.onblur = function () { highlight(side, null); };
      }
      diag.appendChild(b);
    });
    if (!(run.diag.callouts || []).length) diag.appendChild(el('p', 'nolab', run.diag.summary || ''));
    box.appendChild(diag);
  }

  function highlight(side, ids) {
    var box = $('flank-' + side.key), gantt = box.querySelector('.gantt'), plate = box.querySelector('.plate');
    if (!gantt) return;
    gantt.classList.toggle('dim', !!ids);
    if (plate) plate.classList.toggle('dim', !!ids);
    [].forEach.call(box.querySelectorAll('.blk, .plate i'), function (b) {
      b.classList.toggle('hl', !!ids && ids.indexOf(b.dataset.id) >= 0);
    });
  }

  function renderLab(side, box) {
    var run = side.run, sc = run.score, ev = labEvents(run);
    var h = el('h2', null, 'In the lab');
    var stat = el('span');
    stat.innerHTML = '<b>' + g(sc.lab_days_used) + '</b> of ' + LAB.days + ' days · <b>' + ev.filter(function (e) { return e.type !== 'derive'; }).length
      + '</b> experiments · <b>' + (run.n_calls || 0) + '</b> tool calls';
    h.appendChild(stat);
    box.appendChild(h);

    var gantt = el('div', 'gantt'), nLanes = 0;
    ROWS.forEach(function (rw) {                    // greedy lane packing, per row
      var ends = [];
      ev.filter(function (e) { return e.type === rw[0]; }).forEach(function (e) {
        var vis = Math.max(e.end, e.start + 0.9), k = 0;
        while (k < ends.length && ends[k] > e.start + 1e-6) k++;
        ends[k] = vis;
        e.lane = k;
      });
      rw.lanes = Math.max(1, ends.length);
      nLanes += rw.lanes;
    });
    var dense = nLanes > 9;                         // busy labs get thinner lanes so the panel keeps its height
    box.classList.toggle('dense', dense);
    ROWS.forEach(function (rw) {
      var LH = dense ? (rw[0] === 'phospho' && rw.lanes > 1 ? 9 : 14) : nLanes > 7 ? 15 : 19;
      var items = ev.filter(function (e) { return e.type === rw[0]; });
      var row = el('div', 'g-row' + (items.length ? '' : ' none'));
      row.appendChild(el('span', null, rw[1]));
      var lane = el('div', 'g-lane');
      lane.style.height = rw.lanes * LH + 3 + 'px';
      items.forEach(function (e) {
        var b = el('button', 'blk ' + e.cls + (e.type === 'derive' ? ' derive' : ''));
        b.type = 'button';
        b.dataset.id = e.id;
        b.style.left = (e.start / LAB.days * 100) + '%';
        b.style.width = 'max(7px, ' + ((e.end - e.start) / LAB.days * 100) + '%)';
        b.style.top = (3 + e.lane * LH) + 'px';
        if (LH < 19) { b.style.height = b.style.lineHeight = (LH - 3) + 'px'; b.classList.add('sm'); }
        if ((e.end - e.start) >= 2.5 && LH >= 12) b.textContent = e.id + (e.end - e.start >= 5 && e.type !== 'derive' ? ' · ' + e.n : '');
        b.setAttribute('aria-label', e.id + ', ' + rw[1] + ' on ' + srcName(run, e.cells) + ', lab days ' + g(e.start) + ' to ' + g(e.end));
        b.onmouseenter = function () { if (!popPinned) showPop(side, e, b, false); };
        b.onmouseleave = function () { if (!popPinned) closePop(); };
        b.onclick = function (evn) { evn.stopPropagation(); showPop(side, e, b, true); };
        lane.appendChild(b);
      });
      row.appendChild(lane);
      gantt.appendChild(row);
    });
    var axis = el('div', 'g-axis');
    axis.appendChild(el('span', null, 'lab day'));
    var ticks = el('div');
    [0, 7, 14, 21, 28].forEach(function (d) { var t = el('i', null, String(d)); t.style.left = (d / LAB.days * 100) + '%'; ticks.appendChild(t); });
    axis.appendChild(ticks);
    gantt.appendChild(axis);
    var over = el('div', 'g-over');
    var sub = run.submitted_day, stop = sub != null ? sub : sc.lab_days_used;
    if (LAB.days - stop >= 0.25) {
      var un = el('div', 'g-unused');
      un.style.left = (stop / LAB.days * 100) + '%';
      if (LAB.days - stop >= 5) un.appendChild(el('b', null, g(LAB.days - stop) + ' days unused'));
      over.appendChild(un);
    } else if (sub != null) {
      var sm = el('div', 'g-submit');
      sm.style.left = (sub / LAB.days * 100) + '%';
      over.appendChild(sm);
    }
    gantt.appendChild(over);
    box.appendChild(gantt);

    // one tick per tool call, on the lab clock
    var calls = el('div', 'g-calls'), seen = {};
    (run.calls || []).forEach(function (c) {
      var k = String(c[1]);
      seen[k] = (seen[k] || 0) + 1;
      var t = el('i', c[3] ? '' : 'fail');
      t.style.left = 'min(100% - 1px, calc(' + (c[1] / LAB.days * 100) + '% + ' + (seen[k] - 1) * 2 + 'px))';
      calls.appendChild(t);
    });
    calls.onmousemove = function (evn) {
      var r = calls.getBoundingClientRect(), day = (evn.clientX - r.left) / r.width * LAB.days, bestD = null;
      (run.calls || []).forEach(function (c) { if (c[1] <= day + 0.4 && (bestD == null || c[1] > bestD)) bestD = c[1]; });
      if (bestD == null) return hideTip();
      var list = run.calls.filter(function (c) { return c[1] === bestD; });
      var html = '<b>Tool calls on lab day ' + esc(g(bestD)) + '</b><table>' + list.slice(0, 12).map(function (c) {
        return '<tr><td>' + esc(c[0]) + (c[3] ? '' : ' (refused)') + '</td><td>' + esc(c[4] || '') + '</td></tr>';
      }).join('') + '</table>' + (list.length > 12 ? '… ' + (list.length - 12) + ' more' : '');
      showTip(html, evn.clientX, r.bottom + 8);
    };
    calls.onmouseleave = hideTip;
    box.appendChild(calls);

    // the well budget as one 384-well plate
    var budget = el('div', 'budget'), plate = el('div', 'plate'), used = 0;
    plate.setAttribute('role', 'img');
    plate.setAttribute('aria-label', sc.wells_used + ' of ' + LAB.wells + ' wells used');
    ev.forEach(function (e) {
      for (var i = 0; i < e.n && used < LAB.wells; i++, used++) { var d = el('i', e.cls); d.dataset.id = e.id; plate.appendChild(d); }
    });
    for (; used < LAB.wells; used++) plate.appendChild(el('i'));
    budget.appendChild(plate);
    var bt = el('div', 'budget-text');
    var big = el('span', 'big', sc.wells_used + ' ');
    big.appendChild(el('small', null, 'of ' + LAB.wells + ' wells used'));
    bt.appendChild(big);
    var wall = run.wall_seconds ? ', after ' + (run.wall_seconds >= 90 ? Math.round(run.wall_seconds / 60) + ' min' : Math.round(run.wall_seconds) + ' s') + ' of wall-clock time' : '';
    bt.appendChild(el('div', null, sub != null ? 'Regimen submitted on lab day ' + g(sub) + wall + '.' : 'No regimen submitted.'));
    var srcs = el('div', 'srcs');
    [['s-tumour', 'sw-c0', 'tumour biopsy'], ['s-healthy', 'sw-h', 'healthy cells'], ['s-line', 'sw-c1', 'derived line']].forEach(function (s) {
      var sp = el('span'); sp.appendChild(el('i', 'sw ' + s[1])); sp.appendChild(document.createTextNode(s[2])); srcs.appendChild(sp);
    });
    bt.appendChild(srcs);
    budget.appendChild(bt);
    box.appendChild(budget);
  }

  // ---------------------------------------------------------------- experiment popover
  var popPinned = false, popAnchor = null;
  function closePop() {
    $('pop').hidden = true;
    popPinned = false;
    if (popAnchor) popAnchor.classList.remove('pinned');
    popAnchor = null;
  }
  function comboText(w) {
    var ks = Object.keys(w).sort(function (a, b) { return DRUGS.indexOf(a) - DRUGS.indexOf(b); });
    return ks.length ? ks.map(function (k) { return k + ' ' + g(w[k]); }).join(' + ') : 'no drug';
  }
  function isFinal(run, w) {
    var fc = (run.diag.facts || {}).final_combo || [], ks = Object.keys(w);
    return fc.length > 0 && ks.length === fc.length && fc.every(function (d) { return d in w; });
  }

  function showPop(side, e, anchor, pin) {
    var run = side.run, p = patient(S.patient), pop = $('pop'), html = '';
    var kindName = { sequencing: 'Deep sequencing', phospho: 'Signalling: pERK and pAKT', viability: 'Viability', derive: 'New cell line' }[e.type];
    if (e.type === 'derive') {
      var d = e.derive, ln = d.line && run.lines[d.line];
      html += '<h3>' + esc(d.line || 'No line') + '<small>' + kindName + '</small></h3>';
      html += '<p class="sub">Expanded from the survivors of ' + esc(d.source) + ' well ' + esc(d.well) + ' · lab day' + (d.ok ? 's ' + g(d.day) + '–' + g(d.day + 3) : ' ' + g(d.day)) + ' · 8 wells</p>';
      if (!d.ok) html += '<p>Fewer than 30 live cells in that well; nothing grew out.</p>';
      else if (ln) {
        html += '<p>' + esc(ln.note) + '</p>';
        if (ln.frac) {
          html += '<div class="mix">' + ln.frac.map(function (fr, i) { return fr > 0 ? '<i style="flex:' + fr + ';background:var(--c' + i + ')"></i>' : ''; }).join('') + '</div>';
          html += '<p class="note">True mix, which the agent cannot see: ' + ln.frac.map(function (fr, i) { return pct(fr) + ' ' + esc(p.clones[i].label); }).join(', ') + '.</p>';
        }
      }
    } else {
      var x = e.exp;
      html += '<h3>' + esc(x.id) + '<small>' + kindName + '</small></h3>';
      html += '<p class="sub">' + esc(srcName(run, x.cells)) + ' · lab days ' + g(x.start) + '–' + g(x.end) + ' · ' + e.n + ' wells</p>';
      if (x.kind === 'sequencing') {
        var known = p.chart.bulk_sequencing.alterations;
        html += '<div class="wells"><table><tr><th>alteration</th><th>share of cells</th></tr>' + (x.alterations.length ? x.alterations.map(function (a) {
          var isNew = known.indexOf(a.alteration) < 0;
          return '<tr><td' + (isNew ? ' class="found"' : '') + '>' + esc(a.alteration) + (isNew ? ' · not on the chart' : '') + '</td><td class="n">' + pct(a.cell_fraction) + '</td></tr>';
        }).join('') : '<tr><td colspan="2">nothing reported</td></tr>') + '</table></div>';
        html += '<p class="note">Reports alterations carried by at least ' + pct(x.lod) + ' of the cells.</p>';
      } else if (x.kind === 'phospho') {
        html += '<div class="key"><span><i style="border-color:var(--ink)"></i>pERK</span><span><i style="border-color:var(--ink-3)"></i>pAKT</span><span>mid line = untreated (1.0)</span></div>';
        html += '<div class="wells"><table><tr><th>well, µM</th><th>pERK</th><th>pAKT</th></tr>' + x.wells.map(function (w, i) {
          function bar(v, c) { return '<span class="barbox"><span class="bar ' + c + '" style="width:' + clamp(v / 2, 0, 1) * 84 + 'px"></span></span>'; }
          return '<tr><td>' + esc(comboText(w)) + '</td><td class="n">' + bar(x.pERK[i], 'erk') + ' ' + x.pERK[i].toFixed(2) + '</td><td class="n">' + bar(x.pAKT[i], 'akt') + ' ' + x.pAKT[i].toFixed(2) + '</td></tr>';
        }).join('') + '</table></div>';
      } else {
        html += viabilityHtml(run, x);
      }
    }
    pop.innerHTML = html;
    pop.hidden = false;
    var r = anchor.getBoundingClientRect(), pw = pop.offsetWidth, ph = pop.offsetHeight;
    var left = clamp(r.left, 8, window.innerWidth - pw - 8);
    var top = r.bottom + 8 + ph > window.innerHeight - 8 ? Math.max(8, r.top - 8 - ph) : r.bottom + 8;
    pop.style.left = left + 'px';
    pop.style.top = top + 'px';
    if (popAnchor) popAnchor.classList.remove('pinned');
    popAnchor = anchor;
    popPinned = pin;
    if (pin) anchor.classList.add('pinned');
    // row hover lights its curve
    [].forEach.call(pop.querySelectorAll('tr[data-w]'), function (tr) {
      tr.onmouseenter = function () { var c = pop.querySelector('path[data-w="' + tr.dataset.w + '"]'); if (c) { c.classList.add('hot'); c.parentNode.appendChild(c); } };
      tr.onmouseleave = function () { var c = pop.querySelector('path[data-w="' + tr.dataset.w + '"]'); if (c) c.classList.remove('hot'); };
    });
  }

  function viabilityHtml(run, x) {
    var W = 390, H = 128, ml = 34, mr = 8, mt = 6, mb = 18, dur = x.end - x.start;
    var maxC = 10;
    x.counts.forEach(function (row) { row.forEach(function (v) { if (v > maxC) maxC = v; }); });
    var lo = 1, hi = Math.log10(maxC * 1.15);
    var d0 = x.read_days.length > 1 ? x.read_days[0] : 0, d1 = x.read_days.length > 1 ? x.read_days[x.read_days.length - 1] : dur;
    function X(d) { return ml + 6 + (d - d0) / ((d1 - d0) || 1) * (W - ml - mr - 12); }
    function Y(v) { return mt + (1 - (Math.log10(Math.max(v, 10)) - lo) / (hi - lo)) * (H - mt - mb); }
    var svg = '<svg width="' + W + '" height="' + H + '" role="img" aria-label="Live cells per well over time">';
    [10, 100, 1000, 10000, 100000].forEach(function (t) {
      if (Math.log10(t) > hi) return;
      svg += '<line class="ax" x1="' + ml + '" x2="' + (W - mr) + '" y1="' + Y(t) + '" y2="' + Y(t) + '"/><text x="' + (ml - 5) + '" y="' + (Y(t) + 3) + '" text-anchor="end">' + (t >= 1000 ? t / 1000 + 'k' : t) + '</text>';
    });
    x.read_days.forEach(function (d) { svg += '<text x="' + X(d) + '" y="' + (H - 4) + '" text-anchor="middle">' + g(d) + '</text>'; });
    svg += '<text x="' + (ml - 5) + '" y="' + (H - 4) + '" text-anchor="end">day</text>';
    var order = x.wells.map(function (w, i) { return i; });
    function cls(i) { return Object.keys(x.wells[i]).length === 0 ? ' ctl' : isFinal(run, x.wells[i]) ? ' final' : ''; }
    order.sort(function (a, b) { return (cls(a) ? 1 : 0) - (cls(b) ? 1 : 0); });
    order.forEach(function (i) {
      var d = x.counts[i].map(function (v, k) { return (k ? 'L' : 'M') + X(x.read_days[k]).toFixed(1) + ' ' + Y(v).toFixed(1); }).join('');
      if (x.read_days.length === 1) d += 'h0.1';
      svg += '<path class="w' + cls(i) + '" data-w="' + i + '" d="' + d + '"' + (x.read_days.length === 1 ? ' stroke-linecap="round" stroke-width="4"' : '') + '/>';
    });
    svg += '</svg>';
    var anyFinal = x.wells.some(function (w) { return isFinal(run, w); });
    var key = '<div class="key"><span><i></i>one well, live cells (log)</span><span><i class="ctl"></i>no drug</span>'
      + (anyFinal ? '<span><i class="final"></i>the submitted combination</span>' : '') + '</div>';
    var lastIdx = x.read_days.length - 1;
    var ctl = x.wells.map(function (w, i) { return Object.keys(w).length ? null : x.counts[i][lastIdx]; }).filter(function (v) { return v != null; });
    var ref = ctl.length ? ctl.reduce(function (a, b) { return a + b; }, 0) / ctl.length : maxC;
    var rows = x.wells.map(function (w, i) { return i; }).sort(function (a, b) { return x.counts[a][lastIdx] - x.counts[b][lastIdx]; });
    var table = '<div class="wells"><table><tr><th>well, µM (fewest survivors first)</th><th>cells, day ' + g(x.read_days[lastIdx]) + '</th></tr>' + rows.map(function (i) {
      var v = x.counts[i][lastIdx];
      return '<tr data-w="' + i + '"' + (isFinal(run, x.wells[i]) ? ' class="final"' : '') + '><td>' + esc(comboText(x.wells[i])) + '</td><td class="n"><span class="bar" style="width:' + clamp(v / ref, 0, 1.2) * 60 + 'px"></span> ' + int(v) + '</td></tr>';
    }).join('') + '</table></div>';
    return svg + key + table;
  }

  // ---------------------------------------------------------------- tooltip
  function showTip(html, x, y) {
    var t = $('tip');
    t.innerHTML = html;
    t.hidden = false;
    var w = t.offsetWidth, h = t.offsetHeight;
    t.style.left = clamp(x + 12, 8, window.innerWidth - w - 8) + 'px';
    t.style.top = clamp(y, 8, window.innerHeight - h - 8) + 'px';
  }
  function hideTip() { $('tip').hidden = true; }

  // ---------------------------------------------------------------- the year: drugs, cells, normal tissue
  var TR = { drug: 30, cells: 58, norm: 27, gap: 6 };
  function yearGeom() {
    var W = $('year').clientWidth, gl = cssPx('--gl'), gr = cssPx('--gr');
    return { W: W, gl: gl, gr: gr, pw: W - gl - gr, x: function (d) { return gl + d / DAYS * (W - gl - gr); } };
  }

  function renderYear(side) {
    var run = side.run, p = patient(S.patient), G = yearGeom(), box = $('year-' + side.key);
    if (G.pw <= 0) return;                          // hidden (overview): drawn again by show()
    var H = TR.drug + TR.cells + TR.norm + 2 * TR.gap, x = G.x, rp = run.replay;
    var s = '<svg width="' + G.W + '" height="' + H + '" role="img" aria-label="' + esc(run.label) + ': drug levels, cells in the lesion and normal tissue over 364 days">';
    var te = run.treat_end || 0;
    if (te > 0) {
      s += '<rect class="treat" x="' + x(0) + '" y="0" width="' + (x(te) - x(0)) + '" height="' + H + '"/>';
      s += '<line class="stop" x1="' + x(te) + '" x2="' + x(te) + '" y1="0" y2="' + H + '"/>';
    }
    // ---- drugs
    var nm = run.label, room = G.gl - 10, need = nm.length * 7.3, fit = '';
    if (need > room * 1.12) { nm = nm.slice(0, Math.floor(room * 1.12 / 7.3) - 1) + '…'; need = nm.length * 7.3; }
    if (need > room) fit = ' textLength="' + room + '" lengthAdjust="spacingAndGlyphs"';
    s += '<text class="t-run" x="0" y="' + (TR.drug / 2 + 5) + '"' + fit + '><title>' + esc(run.label) + '</title>' + esc(nm) + '</text>';
    var drugs = DRUGS.filter(function (d) { return run.regimen && run.regimen.phases.some(function (ph) { return ph.doses[d] > 0; }); });
    if (!run.regimen) {
      s += '<text class="t-drug" x="' + (x(0) + 8) + '" y="' + (TR.drug / 2 + 4) + '">no drug: nothing was submitted, the patient is untreated</text>';
    } else {
      var rh = Math.min(9, (TR.drug - 2) / Math.max(1, drugs.length));
      var top = (TR.drug - rh * drugs.length) / 2;
      drugs.forEach(function (d, i) {
        var d0 = 0, y1 = top + i * rh;
        run.regimen.phases.forEach(function (ph) {
          var v = ph.doses[d] || 0, d1 = Math.min(d0 + ph.days, LAB.treat_days);
          if (v > 0 && d1 > d0) s += '<rect class="a-d" x="' + x(d0).toFixed(1) + '" y="' + y1.toFixed(1) + '" width="' + (x(d1) - x(d0)).toFixed(1) + '" height="' + Math.max(1.5, rh - 1.5).toFixed(1) + '" opacity="' + (0.3 + 0.7 * Math.min(v, 2) / 2).toFixed(2) + '"/>';
          d0 += ph.days;
        });
      });
      var list = drugs.map(function (d) {
        var seq = [];
        run.regimen.phases.forEach(function (ph) { var v = ph.doses[d] || 0; if (!seq.length || seq[seq.length - 1] !== v) seq.push(v); });
        return esc(d) + ' <tspan>' + seq.map(g).join('→') + '×</tspan>';
      }).join('   ');
      var nph = run.regimen.phases.length;
      s += '<text class="t-drug" x="' + (x(te) + 12) + '" y="' + (TR.drug / 2 + 4) + '">' + list
        + '<tspan dx="14" style="fill:var(--ink-3);font-family:var(--sans)">× = multiple of the labelled daily dose · ' + (nph > 1 ? nph + ' phases · ' : '') + 'dosed to day ' + te + ', then no drug</tspan></text>';
    }
    // ---- cells in the lesion (log)
    var cy0 = TR.drug + TR.gap, ch = TR.cells, lo = Math.log10(0.5), hi = Math.log10(4096);
    function Y(v) { return cy0 + (1 - (Math.log10(Math.max(v, 0.5)) - lo) / (hi - lo)) * ch; }
    s += '<text class="t-trk" x="0" y="' + (cy0 + ch / 2 + 4) + '">cells in lesion</text>';
    [10, 1000].forEach(function (t) {
      s += '<line class="grid" x1="' + G.gl + '" x2="' + (G.gl + G.pw) + '" y1="' + Y(t) + '" y2="' + Y(t) + '"/>';
      s += '<text class="t-num" text-anchor="end" x="' + (G.gl - 6) + '" y="' + (Y(t) + 3.5) + '">' + int(t) + '</text>';
    });
    s += '<line class="grid" x1="' + G.gl + '" x2="' + (G.gl + G.pw) + '" y1="' + (cy0 + ch) + '" y2="' + (cy0 + ch) + '"/>';
    var facts = run.diag.facts || null;
    var series = [{ cls: 'l-h', name: 'healthy', v: rp.healthy, known: facts ? facts.healthy_gone_day : undefined }];
    rp.tumour.forEach(function (col, i) {
      series.push({ cls: 'l-c' + i, name: i ? shortClone(p.clones[i].label) : 'dominant', v: col, tumour: true,
                    known: facts && facts.clones[i] ? facts.clones[i].gone_day : undefined });
    });
    var ends = [], notes = [];
    series.forEach(function (sr) {
      var v = sr.v, n = v.length, gone = n;
      if (sr.known !== undefined) gone = sr.known == null ? n : sr.known;   // exact: no cell left in any replicate lesion
      else while (gone > 0 && v[gone - 1] < 0.5) gone--;
      var stop = Math.min(n - 1, gone), dpath = '';
      for (var k = 0; k <= stop; k++) dpath += (k ? 'L' : 'M') + x(k).toFixed(1) + ' ' + Y(v[k]).toFixed(1);
      s += '<path class="ln ' + sr.cls + '" d="' + dpath + '"/>';
      if (gone < n) { sr.gone = gone; s += '<path class="ln ' + sr.cls + '" d="M' + x(gone) + ' ' + (cy0 + ch - 4) + 'v8"/>'; }
      else ends.push({ y: Y(v[n - 1]), name: sr.name, val: v[n - 1] });
    });
    var tum = series.filter(function (sr) { return sr.tumour; });
    if (tum.every(function (sr) { return sr.gone != null; })) {
      var lastGone = Math.max.apply(null, tum.map(function (sr) { return sr.gone; }));
      notes.push({ d: lastGone, text: 'last tumour cell gone, day ' + lastGone });
    }
    if (series[0].gone != null) notes.push({ d: series[0].gone, text: 'no healthy cells left, day ' + series[0].gone });
    notes.forEach(function (nt) {
      s += '<text class="t-end" x="' + (x(nt.d) + 7) + '" y="' + (cy0 + ch - 6) + '">' + esc(nt.text) + '</text>';
    });
    ends.sort(function (a, b) { return a.y - b.y; });
    for (var i = 1; i < ends.length; i++) if (ends[i].y - ends[i - 1].y < 13) ends[i].y = ends[i - 1].y + 13;
    ends.forEach(function (e) {
      s += '<text class="t-end" x="' + (G.gl + G.pw + 8) + '" y="' + (e.y + 4) + '">' + esc(e.name) + ' <tspan>' + esc(cells(e.val)) + '</tspan></text>';
    });
    // ---- normal tissue (linear, 0 to 100 %)
    var ny0 = cy0 + ch + TR.gap, nh = TR.norm - 3;
    function NY(v) { return ny0 + 2 + (1 - clamp(v, 0, 1)) * nh; }
    s += '<text class="t-trk" x="0" y="' + (ny0 + TR.norm / 2 + 4) + '">normal tissue</text>';
    s += '<text class="t-num" text-anchor="end" x="' + (G.gl - 6) + '" y="' + (NY(1) + 3.5) + '">100' + NBSP + '%</text>';
    s += '<text class="t-num" text-anchor="end" x="' + (G.gl - 6) + '" y="' + (NY(0) + 1) + '">0</text>';
    var np = '', nv = rp.normal, mn = 0;
    for (var k2 = 0; k2 < nv.length; k2++) { np += (k2 ? 'L' : 'M') + x(k2).toFixed(1) + ' ' + NY(nv[k2]).toFixed(1); if (nv[k2] < nv[mn]) mn = k2; }
    s += '<path class="a-n" d="' + np + 'L' + x(nv.length - 1) + ' ' + NY(0) + 'L' + x(0) + ' ' + NY(0) + 'Z"/>';
    s += '<path class="ln l-n" d="' + np + '"/>';
    if (facts && facts.harm_day != null) mn = facts.harm_day;
    var harm = run.score && run.score.harm != null ? run.score.harm : 1 - nv[mn];
    if (harm > 0.005) {
      s += '<line class="harm-drop" x1="' + x(mn) + '" x2="' + x(mn) + '" y1="' + NY(1) + '" y2="' + NY(nv[mn]) + '"/>';
      s += '<circle class="harm-dot" cx="' + x(mn) + '" cy="' + NY(nv[mn]) + '" r="3.5"/>';
    }
    s += '<text class="t-end" x="' + (G.gl + G.pw + 8) + '" y="' + (ny0 + TR.norm / 2 + 4) + '">harm <tspan>' + harm.toFixed(2) + (harm > 0.005 ? ' · day ' + mn : '') + '</tspan></text>';
    s += '</svg>';
    box.innerHTML = s;
  }

  function renderAxis() {
    var G = yearGeom(), s = '<svg width="' + G.W + '" height="22">';
    for (var d = 0; d <= DAYS; d += 28) {
      s += '<line class="grid" x1="' + G.x(d) + '" x2="' + G.x(d) + '" y1="0" y2="5"/>';
      s += '<text class="t-num" text-anchor="' + (d === 0 ? 'start' : d === DAYS ? 'end' : 'middle') + '" x="' + G.x(d) + '" y="17">' + (d === 0 ? 'day 0' : d) + '</text>';
    }
    s += '<text class="t-trk" x="0" y="17">after submission</text>';
    $('year-axis').innerHTML = s + '</svg>';
  }

  function yearTip(day, evn) {
    var p = patient(S.patient), runs = [SIDE.L.run, SIDE.R.run];
    if (!runs[0] || !runs[1]) return;
    var html = '<b>Day ' + day + '</b><table><tr><th></th><th>' + esc(runs[0].label) + '</th><th>' + esc(runs[1].label) + '</th></tr>';
    function row(name, color, fn) {
      html += '<tr><td><i class="k" style="border-color:' + color + '"></i>' + esc(name) + '</td>' + runs.map(function (r) { return '<td class="n">' + esc(fn(r)) + '</td>'; }).join('') + '</tr>';
    }
    row('healthy cells', 'var(--healthy)', function (r) { return cells(r.replay.healthy[day]); });
    p.clones.forEach(function (c, i) { row(c.label, 'var(--c' + i + ')', function (r) { return cells(r.replay.tumour[i][day]); }); });
    row('normal tissue', 'var(--ink-2)', function (r) { return pct(r.replay.normal[day], 0); });
    showTip(html + '</table>', evn.clientX, evn.clientY + 16);
  }

  function bindYear() {
    var hit = $('year-hit'), drag = false;
    function dayAt(evn) { var r = hit.getBoundingClientRect(); return clamp((evn.clientX - r.left) / r.width * DAYS, 0, DAYS); }
    hit.addEventListener('pointerdown', function (evn) {
      drag = true; hit.setPointerCapture(evn.pointerId);
      pause(); setDay(Math.round(dayAt(evn)), true);
    });
    hit.addEventListener('pointermove', function (evn) {
      var d = Math.round(dayAt(evn));
      if (drag) setDay(d, true);
      $('hoverline').style.opacity = drag ? 0 : 1;
      $('hoverline').style.transform = 'translateX(' + (d / DAYS * hit.clientWidth) + 'px)';
      yearTip(d, evn);
    });
    function end() { if (drag) { drag = false; writeHash(); } }
    hit.addEventListener('pointerup', end);
    hit.addEventListener('pointercancel', end);
    hit.addEventListener('pointerleave', function () { $('hoverline').style.opacity = 0; hideTip(); });
  }

  // ---------------------------------------------------------------- overview
  function renderOverview() {
    var t = $('ov-table'), labels = [];
    IDX.patients.forEach(function (p) {
      p.runs.forEach(function (r) {
        var base = r.label.replace(/ · run \d+$/, '');
        if (!labels.some(function (l) { return l.label === r.label; })) labels.push({ label: r.label, kind: r.kind, order: r.order, base: base });
      });
    });
    var rank = { reference: 0, agent: 1, baseline: 2 };
    labels.sort(function (a, b) { return rank[a.kind] - rank[b.kind] || a.order - b.order || (a.label < b.label ? -1 : 1); });
    var html = '<thead><tr><th></th>' + IDX.patients.map(function (p) {
      var minors = p.clones.slice(1).map(function (c) { return shortClone(c.label) + ' ' + pct(c.share); }).join(', ');
      return '<th><b>Patient ' + esc(p.id) + '</b><span>' + esc(p.chart.tissue) + ' · hidden: ' + esc(minors || 'no minor clone') + '</span></th>';
    }).join('') + '<th>mean</th></tr></thead><tbody>';
    var prevKind = null;
    function rowHtml(label, kindLabel, get, sep) {
      var sum = 0, n = 0;
      var tds = IDX.patients.map(function (p) {
        var r = get(p);
        if (!r) return '<td><div class="cell">not run</div></td>';
        if (r.missing) return '<td><div class="cell" title="' + esc(r.reason) + '">' + (r.running ? 'still running' : 'run failed: no episode') + '</div></td>';
        sum += r.reward; n++;
        return '<td><button type="button" data-p="' + esc(p.id) + '" data-r="' + esc(r.id) + '" style="--r:' + clamp(r.reward, 0, 1).toFixed(3) + '"><b>' + r.reward.toFixed(2) + '</b><span>' + esc(r.tag || '') + '</span></button></td>';
      }).join('');
      return '<tr' + (sep ? ' class="sep"' : '') + '><th>' + esc(label) + '<small>' + esc(kindLabel) + '</small></th>' + tds + '<td class="mean">' + (n ? (sum / n).toFixed(2) : '') + '</td></tr>';
    }
    var kn = { reference: 'reference run', agent: 'AI agent', baseline: 'baseline' };
    labels.forEach(function (l) {
      html += rowHtml(l.label, kn[l.kind], function (p) { return p.runs.filter(function (r) { return r.label === l.label; })[0]; }, prevKind && prevKind !== l.kind);
      prevKind = l.kind;
    });
    if (IDX.patients.some(function (p) { return p.best; })) html += rowHtml('Best-known regimen', 'searched with the hidden parameters', function (p) { return p.best; }, true);
    t.innerHTML = html + '</tbody>';
    t.style.width = 'min(100%, ' + (250 + IDX.patients.length * 300 + 90) + 'px)';
    if (!$('ov-scale')) {
      var sc = el('div', 'ov-scale');
      sc.id = 'ov-scale';
      sc.innerHTML = 'reward 0<i></i>1 · reward = healed − harm, mean of 16 replicate lesions per patient';
      t.parentNode.parentNode.appendChild(sc);
    }
    [].forEach.call(t.querySelectorAll('button'), function (b) {
      b.onclick = function () {
        var p = patient(b.dataset.p), r = entry(p, b.dataset.r), ref = reference(p);
        S.patient = p.id; S.left = r.id;
        S.right = ref && ref.id !== r.id ? ref.id : p.best && r.id !== 'best' ? 'best' : r.id;
        S.view = 'replay';
        show(); loadSides(true);
      };
    });
  }

  // ---------------------------------------------------------------- layout and wiring
  function layout() {
    var W = document.documentElement.clientWidth, H = window.innerHeight, narrow = W <= 1240, tight = W <= 1500;
    var pad = tight ? 16 : 24, gap = tight ? 16 : 20, flankMin = W >= 1800 ? 380 : 320;
    var availW = narrow ? (W - 2 * pad - 16) / 2 : (W - 2 * pad - 16 - 2 * gap - 2 * flankMin) / 2;
    var fixedV = 46 + 34 + 4 + (tight ? 104 : 58) + 100 + (2 * (TR.drug + TR.cells + TR.norm + 2 * TR.gap) + 22 + 30) + 8;
    var availH = narrow ? 1e9 : H - fixedV;
    var s = Math.floor(Math.min(availW, availH) / 64), fits = !narrow && s >= 6;
    // too short to show everything at a legible size: let the year band fall below the fold
    if (s < 6 && !narrow) s = Math.floor(Math.min(availW, availH + 2 * (TR.drug + TR.cells + TR.norm + 2 * TR.gap) + 52) / 64);
    s = clamp(s, 4, 12);
    var changed = s !== scale || !SIDE.L.px;
    scale = s;
    document.documentElement.style.setProperty('--lat', s * 64 + 'px');
    // when the whole page fits, a very busy lab panel scrolls inside itself instead of pushing the year band out
    document.documentElement.style.setProperty('--flank-max', fits ? ((tight ? 104 : 58) + clamp(s, 4, 12) * 64 + 108) + 'px' : 'none');
    return changed;
  }

  function show() {
    $('replay').hidden = S.view !== 'replay';
    $('overview').hidden = S.view !== 'overview';
    if (S.view === 'overview') { pause(); closePop(); renderOverview(); }
    else if (SIDE.L.run && SIDE.R.run) {            // the year band cannot be measured while hidden: draw it now
      if (layout()) { initLattice(SIDE.L); initLattice(SIDE.R); }
      renderYear(SIDE.L); renderYear(SIDE.R); renderAxis();
      setDay(S.day, true);
    }
    renderTop();
    writeHash();
  }

  function choosePatient(id) {
    var p = patient(id), ref = reference(p), keepL = entry(p, S.left), keepR = entry(p, S.right);
    S.patient = id;
    var fails = p.runs.filter(function (r) { return !r.missing && r.kind !== 'reference'; })
      .sort(function (a, b) { return (b.kind === 'agent') - (a.kind === 'agent') || a.reward - b.reward; });
    S.left = keepL ? keepL.id : (fails[0] || options(p)[0]).id;
    S.right = keepR ? keepR.id : (ref || options(p)[0]).id;
    S.view = 'replay';
    show(); loadSides(true);
  }

  function loadSides(restart) {
    var p = patient(S.patient), token = ++loadToken;
    var eL = entry(p, S.left) || options(p)[0], eR = entry(p, S.right) || reference(p) || options(p)[0];
    S.left = eL.id; S.right = eR.id;
    renderTop(); renderPatientLine(); renderPick(SIDE.L); renderPick(SIDE.R);
    return Promise.all([D.load(p.id, eL.id, eL.file), D.load(p.id, eR.id, eR.file)]).then(function (runs) {
      if (token !== loadToken) return;
      SIDE.L.run = runs[0]; SIDE.R.run = runs[1];
      layout();
      ['L', 'R'].forEach(function (k) {
        var side = SIDE[k];
        side.fi = 0;
        initLattice(side);
        renderCaption(side); renderFlank(side); renderYear(side);
      });
      renderAxis();
      if (restart) { setDay(0, true); if (!matchMedia('(prefers-reduced-motion: reduce)').matches) play(); }
      else setDay(S.day, true);
      writeHash();
    }).catch(function (err) {
      var f = $('fatal'); f.hidden = false; f.textContent = String(err.message || err);
    });
  }

  function applyHash(first) {
    var h = readHash(), def = defaults();
    var p = patient(h.patient) || patient(def.patient);
    var sameP = p.id === def.patient;
    S.patient = p.id;
    S.left = entry(p, h.left) ? h.left : sameP ? def.left : null;
    S.right = entry(p, h.right) ? h.right : sameP ? def.right : null;
    if (!S.left) {
      var c = p.runs.filter(function (r) { return !r.missing && r.kind !== 'reference'; }).sort(function (a, b) { return (b.kind === 'agent') - (a.kind === 'agent') || a.reward - b.reward; });
      S.left = (c[0] || options(p)[0]).id;
    }
    if (!S.right) S.right = (reference(p) || options(p)[0]).id;
    S.view = h.view === 'overview' ? 'overview' : 'replay';
    var hasDay = h.day != null && !isNaN(+h.day);
    S.day = hasDay ? clamp(+h.day, 0, DAYS) : 0;
    show();
    return loadSides(!hasDay && first).then(function () { if (hasDay) { pause(); setDay(S.day, true); } });
  }

  SIDE.L.canvas = $('dish-L').querySelector('canvas');
  SIDE.R.canvas = $('dish-R').querySelector('canvas');

  var sp = $('speed');
  SPEEDS.forEach(function (v) {
    var b = el('button', null, (v === 0.5 ? '½' : v) + '×');
    b.type = 'button';
    b.setAttribute('aria-pressed', String(v === S.speed));
    b.onclick = function () {
      S.speed = v;
      [].forEach.call(sp.children, function (c) { c.setAttribute('aria-pressed', String(c === b)); });
    };
    sp.appendChild(b);
  });
  $('play').onclick = function () { S.playing ? pause() : play(); };
  [].forEach.call(document.querySelectorAll('.views button'), function (b) {
    b.onclick = function () { S.view = b.dataset.view; show(); if (S.view === 'replay') loadSides(false); };
  });
  bindYear();

  document.addEventListener('keydown', function (evn) {
    if (evn.target && /^(SELECT|INPUT|TEXTAREA)$/.test(evn.target.tagName)) return;
    if (evn.key === 'Escape') { closePop(); return; }
    if (S.view !== 'replay' || evn.metaKey || evn.ctrlKey || evn.altKey) return;
    var step = evn.shiftKey ? 7 : 1;
    if (evn.key === ' ' || evn.code === 'Space') {
      if (evn.target && evn.target.tagName === 'BUTTON' && evn.target.id !== 'play') return;
      evn.preventDefault(); S.playing ? pause() : play();
    } else if (evn.key === 'ArrowRight') { evn.preventDefault(); pause(); setDay(Math.floor(S.day) + step); }
    else if (evn.key === 'ArrowLeft') { evn.preventDefault(); pause(); setDay(Math.ceil(S.day) - step); }
    else if (evn.key === 'Home') { evn.preventDefault(); pause(); setDay(0); }
    else if (evn.key === 'End') { evn.preventDefault(); pause(); setDay(DAYS); }
  });
  document.addEventListener('click', function (evn) { if (popPinned && !$('pop').contains(evn.target)) closePop(); });
  window.addEventListener('hashchange', function () { if (location.hash !== lastHash) applyHash(false); });
  var rz = 0;
  window.addEventListener('resize', function () {
    clearTimeout(rz);
    rz = setTimeout(function () {
      if (S.view !== 'replay' || !SIDE.L.run) return;
      var changed = layout();
      ['L', 'R'].forEach(function (k) { if (changed) initLattice(SIDE[k]); renderYear(SIDE[k]); });
      renderAxis();
      setDay(S.day, true);
    }, 80);
  });

  applyHash(true);
})();
