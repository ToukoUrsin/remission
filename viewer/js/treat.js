/* Remission: "Treat a patient". A visitor doses a simulated patient and watches one copy of the lesion for a year.
   Talks to viewer/serve.py (GET api/patients, POST api/treat); needs nothing from viewer/data/. */
(function (global) {
  'use strict';

  var $ = function (id) { return document.getElementById(id); };
  var DAYS = 364, BASE_SPEED = 28, SPEEDS = [0.5, 1, 2, 4];
  var NNBSP = ' ', NBSP = ' ', STORE = 'remission.treat.v1';
  var TR = { drug: 30, cells: 58, norm: 27, gap: 6 };

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
  function pct(x) { var v = 100 * x; return (v > 0.05 && v < 10 ? v.toFixed(1).replace(/\.0$/, '') : v.toFixed(0)) + NBSP + '%'; }
  function g(n) { return String(+(+n).toPrecision(4)); }
  function shortClone(label) { return label.replace(/ clone$/, ''); }
  function own(c) { return (c.own || []).map(function (a) { return a.replace(/^BCL2L1 \(BCL-xL\)/, 'BCL-xL').replace(/^MAP2K1 \(MEK1\)/, 'MEK1'); }).join(' + '); }
  var ICON = {
    bad: '<svg class="mark mark-bad" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="none" stroke="currentColor" stroke-width="1.6"/><path d="M5.2 5.2l5.6 5.6M10.8 5.2l-5.6 5.6" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"/></svg>',
    good: '<svg class="mark mark-good" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="none" stroke="currentColor" stroke-width="1.6"/><path d="M4.6 8.3l2.3 2.3 4.5-4.9" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/></svg>'
  };
  function b64bytes(text) {
    var bin = atob(text), n = bin.length, u8 = new Uint8Array(n), i;
    for (i = 0; i < n; i++) u8[i] = bin.charCodeAt(i);
    return u8;
  }
  function inflate(b64, n) { return global.REMISSION_DATA.inflate(b64bytes(b64), n); }

  // ---------------------------------------------------------------- storage (per visitor, per patient)
  function loadStore() { try { return JSON.parse(localStorage.getItem(STORE) || '{}') || {}; } catch (e) { return {}; } }
  function saveStore(o) { try { localStorage.setItem(STORE, JSON.stringify(o)); } catch (e) { /* private window: keep it in memory */ } }
  var MEM = loadStore();
  function rec(pid) { return MEM[pid] || (MEM[pid] = { attempts: [] }); }

  // ---------------------------------------------------------------- state
  var API = null, DRUGS = [], LIM = null;
  var S = {
    mounted: false, visible: false, pid: null, phases: [], sel: 0,
    result: null, revealed: false, busy: false, day: 0, playing: false, speed: 1, fi: 0,
    frames: null, days: null, day0: null, scale: 7, offline: false
  };
  var onHash = function () {};

  function patient(id) { return API && API.patients.filter(function (p) { return p.id === id; })[0]; }
  function regimenFromState() {
    return { phases: S.phases.map(function (ph) {
      var d = {};
      DRUGS.forEach(function (n) { if (ph.doses[n] > 0) d[n] = ph.doses[n]; });
      return { days: ph.days, doses: d };
    }) };
  }
  function setRegimen(reg) {
    var phases = (reg && reg.phases || []).slice(0, LIM.max_phases).map(function (ph) {
      var d = {};
      DRUGS.forEach(function (n) { d[n] = +(ph.doses && ph.doses[n]) || 0; });
      return { days: clamp(Math.round(ph.days), 1, LIM.treat_days), doses: d };
    });
    if (!phases.length) phases = [{ days: LIM.treat_days, doses: {} }];
    S.phases = phases;
    S.sel = 0;
    renderPhases(); renderDrugs();
  }
  function shortReg(reg) {
    return reg.phases.map(function (ph) {
      var ks = DRUGS.filter(function (d) { return ph.doses[d] > 0; });
      return (ks.length ? ks.map(function (d) { return d + ' ' + g(ph.doses[d]); }).join(' + ') : 'no drug') + ' · ' + ph.days + NBSP + 'd';
    }).join('  →  ');
  }
  function sameReg(a, b) { return JSON.stringify(canon(a)) === JSON.stringify(canon(b)); }
  function canon(reg) {
    return reg.phases.map(function (ph) {
      return [ph.days, DRUGS.map(function (d) { return +(ph.doses[d] || 0); })];
    });
  }

  // ---------------------------------------------------------------- lattice
  var RGB = { 0: [10, 10, 9], 1: [42, 148, 145], 2: [240, 67, 122], 3: [245, 183, 53], 4: [194, 162, 242] };
  var JIT = { 1: [0.80, 0.90, 1.0, 1.09], 2: [0.90, 0.96, 1.0, 1.06], 3: [0.92, 0.97, 1.0, 1.05], 4: [0.92, 0.97, 1.0, 1.05] };
  var SHADE = [];
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
  var LT = {};
  function initLattice() {
    var cv = $('t-canvas'), dpr = window.devicePixelRatio || 1;
    var P = Math.max(1, Math.round(S.scale * dpr)), n = 64 * P;
    LT.P = P;
    LT.gap = P >= 5 ? Math.max(1, Math.round(P / 8)) : 0;
    cv.width = cv.height = n;
    LT.ctx = cv.getContext('2d');
    LT.img = LT.ctx.createImageData(n, n);
    LT.px = new Uint32Array(LT.img.data.buffer);
    LT.px.fill(SHADE[0][0]);
    LT.shown = new Uint8Array(4096).fill(255);
    LT.count = [0, 0, 0, 0, 0];
  }
  function paintGrid(a, b, frac, k) {
    var px = LT.px, shown = LT.shown, P = LT.P, w = P - LT.gap, N = 64 * P;
    var count = [0, 0, 0, 0, 0], dirty = false;
    for (var i = 0; i < 4096; i++) {
      var v = a[i];
      if (b && v !== b[i] && hash2(i, k) < frac) v = b[i];
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
    LT.count = count;
    if (dirty) LT.ctx.putImageData(LT.img, 0, 0);
  }
  function paint(day) {
    if (!LT.px) return;
    if (!S.frames) { if (S.day0) paintGrid(S.day0, null, 0, 0); return; }
    var days = S.days, F = S.frames, k = S.fi, last = days.length - 1;
    if (k > last - 1) k = Math.max(0, last - 1);
    while (k > 0 && days[k] > day) k--;
    while (k < last - 1 && days[k + 1] <= day) k++;
    S.fi = k;
    var frac = last > 0 ? clamp((day - days[k]) / (days[k + 1] - days[k]), 0, 1) : 0;
    paintGrid(F[k], F[Math.min(k + 1, last)], frac, k);
  }

  // ---------------------------------------------------------------- layout
  function layout() {
    var W = document.documentElement.clientWidth, H = window.innerHeight;
    var mode = W <= 760 ? 'phone' : W <= 1180 ? 'mid' : 'wide';
    var pad = W <= 1500 ? 16 : 24;
    var availW = mode === 'phone' ? W - 32 : mode === 'mid' ? W - 2 * pad - 24 - 330 : W - 2 * pad - 2 * 28 - 2 * 340;
    var availH = mode === 'phone' ? 1e9 : H - 420;
    var s = clamp(Math.floor(Math.min(availW, Math.max(availH, 320)) / 64), 4, 10);
    var changed = s !== S.scale || !LT.px;
    S.scale = s;
    document.documentElement.style.setProperty('--tlat', s * 64 + 'px');
    if (changed) initLattice();
    return changed;
  }

  // ---------------------------------------------------------------- clock
  var lastTs = 0, rafId = 0;
  function setDay(d) {
    S.day = clamp(d, 0, DAYS);
    paint(S.day);
    var di = Math.floor(S.day + 1e-6), te = S.result ? S.result.treat_end : 0;
    if (di !== setDay.shown || setDay.res !== S.result) {
      setDay.shown = di; setDay.res = S.result;
      $('t-day').textContent = di;
      $('t-phase').textContent = !S.result ? 'before treatment' : di < te ? 'on treatment' : 'no drug, follow-up';
      $('t-year-hit').setAttribute('aria-valuenow', di);
    }
    liveCounts();
    var hit = $('t-year-hit'), xpx = S.day / DAYS * hit.clientWidth;
    $('t-playhead').style.transform = 'translateX(' + xpx + 'px)';
    $('t-veil').style.left = xpx + 'px';
    if (S.result && !S.revealed && S.day >= DAYS - 1e-6) reveal();
  }
  function frame(ts) {
    if (!S.playing) return;
    var dt = Math.min(ts - lastTs, 100) / 1000;
    lastTs = ts;
    var d = S.day + dt * BASE_SPEED * S.speed;
    if (d >= DAYS) { setDay(DAYS); pause(); return; }
    setDay(d);
    rafId = requestAnimationFrame(frame);
  }
  function play() {
    if (S.playing || !S.result) return;
    if (S.day >= DAYS - 0.5) setDay(0);
    S.playing = true;
    $('t-play').classList.add('on');
    $('t-play').setAttribute('aria-label', 'Pause');
    lastTs = performance.now();
    rafId = requestAnimationFrame(frame);
  }
  function pause() {
    S.playing = false;
    cancelAnimationFrame(rafId);
    $('t-play').classList.remove('on');
    $('t-play').setAttribute('aria-label', 'Play');
  }

  // ---------------------------------------------------------------- counts and outcome under the lattice
  function liveCounts() {
    var c = LT.count, box = $('t-counts');
    if (!box.firstChild) return;
    [].forEach.call(box.querySelectorAll('b[data-s]'), function (b) {
      var s = b.dataset.s.split(',').map(Number), v = 0;
      s.forEach(function (k) { v += c[k] || 0; });
      if (b._v !== v) { b._v = v; b.textContent = int(v); }
    });
  }
  function renderCounts() {
    var box = $('t-counts'), R = S.result;
    box.textContent = '';
    function cnt(states, sw, text) {
      var s = el('span'), b = el('b', null, '0');
      b.dataset.s = states.join(',');
      s.appendChild(el('i', 'sw ' + sw)); s.appendChild(document.createTextNode(text)); s.appendChild(b);
      box.appendChild(s);
    }
    cnt([1], 'sw-h', 'healthy');
    if (!R) { cnt([2, 3, 4], 'sw-c0', 'tumour'); return; }
    cnt([2], 'sw-c0', 'dominant');
    R.truth.clones.slice(1).forEach(function (c, i) {
      cnt([3 + i], 'sw-c' + (i + 1), S.scale < 6 || R.truth.clones.length > 2 ? shortClone(c.label).replace(/-amplified|-loss/, '') : shortClone(c.label));
    });
  }

  function verdict(R) {
    var f = R.facts, cl = R.truth.clones, start = f.start || 1, rw = R.reward.toFixed(2);
    var best = R.best && R.best.reward != null ? R.best.reward : null;
    var score = 'Healed ' + R.healed.toFixed(2) + ', harm ' + R.harm.toFixed(2) + ': reward ' + rw
      + (best != null ? ' (best known ' + best.toFixed(2) + ')' : '') + '.';
    var noDrug = R.regimen.phases.every(function (ph) { return !Object.keys(ph.doses).length; });
    var tox = '';
    if (R.harm >= 0.9 && f.healthy_gone_day != null) tox = ' The drugs killed every healthy cell in the lesion by day ' + f.healthy_gone_day + '.';
    else if (R.harm >= 0.5) tox = ' Normal tissue fell to ' + pct(f.min_normal) + ' on day ' + f.harm_day + '.';
    if (noDrug) return { ok: false, text: 'No drug given: the tumour took the lesion over. Reward 0.00.' };
    if (R.cured >= 1) {
      return { ok: R.reward > 0, text: (f.tumour_gone_day != null ? 'Tumour gone by day ' + f.tumour_gone_day + ' in all ' + R.reps + ' copies of the lesion.' : 'Tumour gone.') + tox + ' ' + score };
    }
    var main = 0;
    cl.forEach(function (c, i) { if (c.end > cl[main].end) main = i; });
    var T = R.replay.tumour, total = function (d) { var v = 0; for (var k = 0; k < T.length; k++) v += T[k][d]; return v; };
    var back = null;
    for (var d = f.nadir_day; d <= DAYS; d++) if (total(d) >= 0.5 * start) { back = d; break; }
    var low;
    if (f.nadir >= 0.9 * start) low = 'The tumour never shrank.';
    else if (f.nadir < 0.05 * start) low = 'Tumour down to ' + cells(f.nadir) + ' cells by day ' + f.nadir_day + ', and the lesion looked healed.';
    else low = 'Tumour down to ' + pct(f.nadir / start) + ' of its size by day ' + f.nadir_day + '.';
    var who, c = cl[main];
    if (main === 0) {
      who = (back != null && f.nadir < 0.9 * start ? ' Back by day ' + back + ': ' : ' ') + 'the dominant clone, the one on the chart, was never driven out.';
    } else {
      var seen = c.share >= R.truth.sequencing_lod
        ? 'too small for the chart, though the lab’s deep sequencing (' + pct(R.truth.sequencing_lod) + ') would have found it'
        : 'below what sequencing can see';
      who = (back != null ? ' Back by day ' + back + ': ' : ' It came back: ') + 'a ' + c.label + ' that was ' + pct(c.share) + ' of the biopsy, ' + seen + '.';
    }
    var part = R.cured > 0 ? ' Cleared in ' + Math.round(R.cured * R.reps) + ' of ' + R.reps + ' copies of the lesion.' : '';
    return { ok: false, text: low + who + part + tox + (R.reward > 0 ? ' ' + score : ' Reward ' + rw + '.') };
  }

  function renderOutcome() {
    var box = $('t-outcome'), R = S.result;
    box.textContent = '';
    box.className = 't-outcome';
    if (!R) {
      box.appendChild(el('p', 't-wait', 'Press Treat. The simulator doses the patient and follows ' + (API ? API.reps : 4) + ' copies of the lesion for 364 days; this one is shown.'));
      return;
    }
    if (!S.revealed) {
      var w = el('div', 'outcome');
      var rw = el('div', 'reward t-pending', '–');
      rw.appendChild(el('small', null, 'reward'));
      w.appendChild(rw);
      w.appendChild(el('div', 'eq', 'The outcome is read on day 364. Press End to skip there.'));
      box.appendChild(w);
      return;
    }
    var v = verdict(R), out = el('div', 'outcome' + (v.ok ? '' : ' fail'));
    var r = el('div', 'reward', R.reward.toFixed(2));
    r.appendChild(el('small', null, 'reward'));
    out.appendChild(r);
    var right = el('div');
    var vt = el('div', 'verdict t-verdict');
    vt.innerHTML = ICON[v.ok ? 'good' : 'bad'] + esc(v.text);
    right.appendChild(vt);
    var eq = el('div', 'eq');
    eq.innerHTML = 'healed <b>' + R.healed.toFixed(2) + '</b> − harm <b>' + R.harm.toFixed(2) + '</b> · tumour-free on day 364 in <b>'
      + Math.round(R.cured * R.reps) + '</b> of ' + R.reps + ' copies';
    right.appendChild(eq);
    out.appendChild(right);
    box.appendChild(out);
  }

  function reveal() {
    S.revealed = true;
    rec(S.pid).seen = true;
    saveStore(MEM);
    $('t-year').classList.add('revealed');
    renderOutcome(); renderHiding(); renderPresets(); renderAttempts(); renderYear();
  }

  // ---------------------------------------------------------------- the year chart (same marks as the replay)
  function yearGeom() {
    var W = $('t-year').clientWidth, narrow = W < 640;
    var gl = narrow ? 56 : W <= 1500 ? 128 : 148, gr = narrow ? 8 : 140;
    return { W: W, gl: gl, gr: gr, pw: W - gl - gr, narrow: narrow, x: function (d) { return gl + d / DAYS * (W - gl - gr); } };
  }
  function renderYear() {
    var G = yearGeom(), R = S.result, x = G.x, box = $('t-year-run');
    if (G.pw <= 0) return;
    var hit = $('t-year-hit');
    hit.style.left = G.gl + 'px'; hit.style.right = G.gr + 'px';
    var H = TR.drug + TR.cells + TR.norm + 2 * TR.gap;
    var s = '<svg width="' + G.W + '" height="' + H + '" role="img" aria-label="Drug levels, cells in the lesion and normal tissue over 364 days">';
    var te = R ? R.treat_end : 0;
    if (te > 0) {
      s += '<rect class="treat" x="' + x(0) + '" y="0" width="' + (x(te) - x(0)) + '" height="' + H + '"/>';
      s += '<line class="stop" x1="' + x(te) + '" x2="' + x(te) + '" y1="0" y2="' + H + '"/>';
    }
    s += '<text class="t-trk" x="0" y="' + (TR.drug / 2 + 4) + '">' + (G.narrow ? 'drugs' : 'drugs given') + '</text>';
    var cy0 = TR.drug + TR.gap, ch = TR.cells, lo = Math.log10(0.5), hi = Math.log10(4096);
    function Y(v) { return cy0 + (1 - (Math.log10(Math.max(v, 0.5)) - lo) / (hi - lo)) * ch; }
    s += '<text class="t-trk" x="0" y="' + (cy0 + ch / 2 + 4) + '">' + (G.narrow ? 'cells' : 'cells in lesion') + '</text>';
    [10, 1000].forEach(function (t) {
      s += '<line class="grid" x1="' + G.gl + '" x2="' + (G.gl + G.pw) + '" y1="' + Y(t) + '" y2="' + Y(t) + '"/>';
      s += '<text class="t-num" text-anchor="end" x="' + (G.gl - 6) + '" y="' + (Y(t) + 3.5) + '">' + int(t) + '</text>';
    });
    s += '<line class="grid" x1="' + G.gl + '" x2="' + (G.gl + G.pw) + '" y1="' + (cy0 + ch) + '" y2="' + (cy0 + ch) + '"/>';
    var ny0 = cy0 + ch + TR.gap, nh = TR.norm - 3;
    function NY(v) { return ny0 + 2 + (1 - clamp(v, 0, 1)) * nh; }
    s += '<text class="t-trk" x="0" y="' + (ny0 + TR.norm / 2 + 4) + '">' + (G.narrow ? 'normal' : 'normal tissue') + '</text>';
    s += '<text class="t-num" text-anchor="end" x="' + (G.gl - 6) + '" y="' + (NY(1) + 3.5) + '">100' + NBSP + '%</text>';
    s += '<text class="t-num" text-anchor="end" x="' + (G.gl - 6) + '" y="' + (NY(0) + 1) + '">0</text>';
    if (!R) {
      s += '<text class="t-trk" x="' + (G.gl + G.pw / 2) + '" y="' + (cy0 + ch / 2 + 4) + '" text-anchor="middle">the year appears here after you press Treat</text>';
      box.innerHTML = s + '</svg>';
      renderAxis(G);
      return;
    }
    var reg = R.regimen, rp = R.replay;
    // ---- drugs
    var drugs = DRUGS.filter(function (d) { return reg.phases.some(function (ph) { return ph.doses[d] > 0; }); });
    if (!drugs.length) {
      s += '<text class="t-drug" x="' + (x(0) + 8) + '" y="' + (TR.drug / 2 + 4) + '">no drug</text>';
    } else {
      var rh = Math.min(9, (TR.drug - 2) / Math.max(1, drugs.length)), top = (TR.drug - rh * drugs.length) / 2;
      drugs.forEach(function (d, i) {
        var d0 = 0, y1 = top + i * rh;
        reg.phases.forEach(function (ph) {
          var v = ph.doses[d] || 0, d1 = Math.min(d0 + ph.days, LIM.treat_days);
          if (v > 0 && d1 > d0) s += '<rect class="a-d" x="' + x(d0).toFixed(1) + '" y="' + y1.toFixed(1) + '" width="' + (x(d1) - x(d0)).toFixed(1) + '" height="' + Math.max(1.5, rh - 1.5).toFixed(1) + '" opacity="' + (0.3 + 0.7 * Math.min(v, 2) / 2).toFixed(2) + '"/>';
          d0 += ph.days;
        });
      });
      var list = drugs.map(function (d) {
        var seq = [];
        reg.phases.forEach(function (ph) { var v = ph.doses[d] || 0; if (!seq.length || seq[seq.length - 1] !== v) seq.push(v); });
        return esc(d) + ' <tspan>' + seq.map(g).join('→') + '×</tspan>';
      }).join('   ');
      s += '<text class="t-drug" x="' + (x(te) + 10) + '" y="' + (TR.drug / 2 + 4) + '">' + list
        + (G.narrow ? '' : '<tspan dx="14" style="fill:var(--ink-3);font-family:var(--sans)">dosed to day ' + te + ', then no drug</tspan>') + '</text>';
    }
    // ---- cells in the lesion (mean of the copies, log scale)
    var series = [{ cls: 'l-h', name: 'healthy', v: rp.healthy, known: R.facts.healthy_gone_day }];
    rp.tumour.forEach(function (col, i) {
      series.push({ cls: 'l-c' + i, name: i ? shortClone(R.truth.clones[i].label) : 'dominant', v: col, tumour: true, known: R.truth.clones[i].gone_day });
    });
    var ends = [], notes = [];
    series.forEach(function (sr) {
      var v = sr.v, n = v.length, gone = sr.known == null ? n : sr.known;
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
    if (!S.revealed || G.narrow) notes = [];
    notes.forEach(function (nt) { s += '<text class="t-end" x="' + (x(nt.d) + 7) + '" y="' + (cy0 + ch - 6) + '">' + esc(nt.text) + '</text>'; });
    if (!G.narrow && S.revealed) {
      ends.sort(function (a, b) { return a.y - b.y; });
      for (var i = 1; i < ends.length; i++) if (ends[i].y - ends[i - 1].y < 13) ends[i].y = ends[i - 1].y + 13;
      ends.forEach(function (e) {
        s += '<text class="t-end" x="' + (G.gl + G.pw + 8) + '" y="' + (e.y + 4) + '">' + esc(e.name) + ' <tspan>' + esc(cells(e.val)) + '</tspan></text>';
      });
    }
    // ---- normal tissue
    var np = '', nv = rp.normal, mn = R.facts.harm_day;
    for (var k2 = 0; k2 < nv.length; k2++) np += (k2 ? 'L' : 'M') + x(k2).toFixed(1) + ' ' + NY(nv[k2]).toFixed(1);
    s += '<path class="a-n" d="' + np + 'L' + x(nv.length - 1) + ' ' + NY(0) + 'L' + x(0) + ' ' + NY(0) + 'Z"/>';
    s += '<path class="ln l-n" d="' + np + '"/>';
    if (R.harm > 0.005) {
      s += '<line class="harm-drop" x1="' + x(mn) + '" x2="' + x(mn) + '" y1="' + NY(1) + '" y2="' + NY(nv[mn]) + '"/>';
      s += '<circle class="harm-dot" cx="' + x(mn) + '" cy="' + NY(nv[mn]) + '" r="3.5"/>';
    }
    if (!S.revealed) { /* read on day 364 */ }
    else if (G.narrow) s += '<text class="t-end" text-anchor="end" x="' + (G.gl + G.pw) + '" y="' + (ny0 + 10) + '">harm <tspan>' + R.harm.toFixed(2) + '</tspan></text>';
    else s += '<text class="t-end" x="' + (G.gl + G.pw + 8) + '" y="' + (ny0 + TR.norm / 2 + 4) + '">harm <tspan>' + R.harm.toFixed(2) + '</tspan></text>';
    box.innerHTML = s + '</svg>';
    renderAxis(G);
  }
  function renderAxis(G) {
    var s = '<svg width="' + G.W + '" height="22">', step = G.narrow ? 91 : 28;
    for (var d = 0; d <= DAYS; d += step) {
      s += '<line class="grid" x1="' + G.x(d) + '" x2="' + G.x(d) + '" y1="0" y2="5"/>';
      s += '<text class="t-num" text-anchor="' + (d === 0 ? 'start' : d === DAYS ? 'end' : 'middle') + '" x="' + G.x(d) + '" y="17">' + (d === 0 ? 'day 0' : d) + '</text>';
    }
    if (!G.narrow) s += '<text class="t-trk" x="0" y="17">from the first dose</text>';
    $('t-year-axis').innerHTML = s + '</svg>';
  }
  function yearTip(day, evn) {
    var R = S.result;
    if (!R || !S.revealed && day > S.day) return hideTip();
    var html = '<b>Day ' + day + '</b><table>';
    function row(name, color, val) { html += '<tr><td><i class="k" style="border-color:' + color + '"></i>' + esc(name) + '</td><td class="n">' + esc(val) + '</td></tr>'; }
    row('healthy cells', 'var(--healthy)', cells(R.replay.healthy[day]));
    R.truth.clones.forEach(function (c, i) { row(c.label, 'var(--c' + i + ')', cells(R.replay.tumour[i][day])); });
    row('normal tissue', 'var(--ink-2)', pct(R.replay.normal[day]));
    showTip(html + '</table><div class="t-tipnote">mean of ' + R.reps + ' copies of the lesion</div>', evn.clientX, evn.clientY + 16);
  }
  function showTip(html, x, y) {
    var t = $('tip');
    t.innerHTML = html;
    t.hidden = false;
    var w = t.offsetWidth, h = t.offsetHeight;
    t.style.left = clamp(x + 12, 8, window.innerWidth - w - 8) + 'px';
    t.style.top = clamp(y, 8, window.innerHeight - h - 8) + 'px';
  }
  function hideTip() { $('tip').hidden = true; }
  function bindYear() {
    var hit = $('t-year-hit'), drag = false;
    function dayAt(evn) { var r = hit.getBoundingClientRect(); return clamp((evn.clientX - r.left) / r.width * DAYS, 0, DAYS); }
    hit.addEventListener('pointerdown', function (evn) {
      if (!S.result) return;
      drag = true; hit.setPointerCapture(evn.pointerId);
      pause(); setDay(Math.round(dayAt(evn)));
    });
    hit.addEventListener('pointermove', function (evn) {
      var d = Math.round(dayAt(evn));
      if (drag) setDay(d);
      $('t-hoverline').style.opacity = drag || !S.result ? 0 : 1;
      $('t-hoverline').style.transform = 'translateX(' + (d / DAYS * hit.clientWidth) + 'px)';
      if (evn.pointerType === 'mouse') yearTip(d, evn);
    });
    function end() { drag = false; }
    hit.addEventListener('pointerup', end);
    hit.addEventListener('pointercancel', end);
    hit.addEventListener('pointerleave', function () { $('t-hoverline').style.opacity = 0; hideTip(); });
    hit.addEventListener('keydown', function (evn) {
      if (!S.result) return;
      var step = evn.shiftKey ? 7 : 1;
      if (evn.key === 'ArrowRight') { evn.preventDefault(); pause(); setDay(Math.floor(S.day) + step); }
      else if (evn.key === 'ArrowLeft') { evn.preventDefault(); pause(); setDay(Math.ceil(S.day) - step); }
    });
  }

  // ---------------------------------------------------------------- the plan: patient, presets, phases, doses
  function renderPatient() {
    var p = patient(S.pid), box = $('t-chart'), seq = p.chart.bulk_sequencing;
    box.innerHTML = '<p class="t-type">' + esc(p.tumour_type) + '</p>'
      + '<dl><dt>Biopsy sequencing</dt><dd><b>' + esc(seq.alterations.join(', ')) + '</b></dd>'
      + '<dt>Wild type</dt><dd>' + esc(seq.wild_type.join(', ')) + '</dd></dl>'
      + '<p class="t-note">The panel reports clones above about 10' + NBSP + '% of the tumour.</p>';
  }
  function presets() {
    var p = patient(S.pid), out = (p.guideline || []).map(function (gd) { return { label: gd.name, reg: gd.regimen }; });
    out.push({ label: 'Every drug at full dose', reg: API.full_dose });
    var r = rec(S.pid), open = !!(r.best && r.attempts.length && r.seen);
    out.push({ label: 'Best regimen we know', reg: open ? r.best.regimen : null, best: true, locked: !open });
    return out;
  }
  function renderPresets() {
    var box = $('t-presets'), cur = regimenFromState();
    box.textContent = '';
    presets().forEach(function (ps) {
      var b = el('button', 't-chip' + (ps.best ? ' t-best' : ''), ps.label);
      b.type = 'button';
      if (ps.locked) { b.classList.add('t-locked-chip'); b.disabled = true; b.setAttribute('aria-hidden', 'true'); b.tabIndex = -1; box.appendChild(b); return; }
      b.setAttribute('aria-pressed', String(sameReg(cur, ps.reg)));
      b.onclick = function () { setRegimen(ps.reg); renderPresets(); };
      box.appendChild(b);
    });
  }
  function usedDays(except) {
    return S.phases.reduce(function (a, ph, i) { return a + (i === except ? 0 : ph.days); }, 0);
  }
  function renderPhases() {
    var box = $('t-phases');
    box.textContent = '';
    var tabs = el('div', 't-tabs');
    tabs.setAttribute('role', 'group');
    tabs.setAttribute('aria-label', 'Phases');
    var d0 = 0;
    S.phases.forEach(function (ph, i) {
      var b = el('button', 't-tab');
      b.type = 'button';
      b.setAttribute('aria-pressed', String(i === S.sel));
      b.innerHTML = (S.phases.length > 1 ? 'Phase ' + (i + 1) + ' · ' : '') + 'days ' + d0 + '–' + (d0 + ph.days);
      b.onclick = function () { S.sel = i; renderPhases(); renderDrugs(); };
      tabs.appendChild(b);
      d0 += ph.days;
    });
    if (S.phases.length < LIM.max_phases) {
      var add = el('button', 't-tab t-add', '+ phase');
      add.type = 'button';
      add.title = 'Add a phase after the last one';
      add.onclick = function () {
        var last = S.phases[S.phases.length - 1], free = LIM.treat_days - usedDays(-1), days;
        if (free >= 7) days = Math.min(free, 28);
        else { days = Math.max(1, Math.floor(last.days / 2)); last.days -= days; if (last.days < 1) { last.days = 1; } }
        S.phases.push({ days: days, doses: JSON.parse(JSON.stringify(last.doses)) });
        S.sel = S.phases.length - 1;
        renderPhases(); renderDrugs(); renderPresets();
      };
      tabs.appendChild(add);
    }
    box.appendChild(tabs);
    var ph = S.phases[S.sel], max = LIM.treat_days - usedDays(S.sel);
    var row = el('div', 't-dayrow');
    var lab = el('label', null, S.phases.length > 1 ? 'Days, phase ' + (S.sel + 1) : 'Days of dosing');
    lab.htmlFor = 't-days-in';
    row.appendChild(lab);
    var rng = el('input');
    rng.type = 'range'; rng.id = 't-days-in'; rng.min = 1; rng.max = max; rng.step = 1; rng.value = ph.days;
    rng.style.setProperty('--f', ((ph.days - 1) / Math.max(1, max - 1) * 100) + '%');
    var val = el('output', 't-val', ph.days + ' d');
    val.htmlFor = 't-days-in';
    rng.oninput = function () {
      ph.days = +rng.value;
      val.textContent = ph.days + ' d';
      rng.style.setProperty('--f', ((ph.days - 1) / Math.max(1, max - 1) * 100) + '%');
      var d = 0;
      [].forEach.call(tabs.querySelectorAll('.t-tab:not(.t-add)'), function (b, i) {
        b.innerHTML = (S.phases.length > 1 ? 'Phase ' + (i + 1) + ' · ' : '') + 'days ' + d + '–' + (d + S.phases[i].days);
        d += S.phases[i].days;
      });
      renderPresets();
    };
    row.appendChild(rng);
    row.appendChild(val);
    if (S.phases.length > 1) {
      var rm = el('button', 't-rm', 'Remove');
      rm.type = 'button';
      rm.setAttribute('aria-label', 'Remove phase ' + (S.sel + 1));
      rm.onclick = function () { S.phases.splice(S.sel, 1); S.sel = Math.max(0, S.sel - 1); renderPhases(); renderDrugs(); renderPresets(); };
      row.appendChild(rm);
    }
    box.appendChild(row);
    var rest = LIM.treat_days - usedDays(-1);
    box.appendChild(el('p', 't-k t-after', 'No drug from day ' + usedDays(-1) + (rest > 0 ? ' (up to ' + rest + ' more days allowed)' : '') + ' to day ' + DAYS + '.'));
  }
  function renderDrugs() {
    var box = $('t-drugs'), ph = S.phases[S.sel];
    box.textContent = '';
    API.drugs.forEach(function (dr) {
      var id = 't-dose-' + dr.name, v = ph.doses[dr.name] || 0;
      var row = el('div', 't-drug' + (v > 0 ? ' on' : ''));
      var lab = el('label');
      lab.htmlFor = id;
      lab.appendChild(el('b', null, dr.name));
      lab.appendChild(el('small', null, dr['class']));
      row.appendChild(lab);
      var rng = el('input');
      rng.type = 'range'; rng.id = id; rng.min = LIM.dose_min; rng.max = LIM.dose_max; rng.step = LIM.dose_step; rng.value = v;
      rng.style.setProperty('--f', (v / LIM.dose_max * 100) + '%');
      var out = el('output', 't-val', v > 0 ? g(v) + '×' : 'off');
      out.htmlFor = id;
      rng.oninput = function () {
        var nv = +rng.value;
        ph.doses[dr.name] = nv;
        out.textContent = nv > 0 ? g(nv) + '×' : 'off';
        row.classList.toggle('on', nv > 0);
        rng.style.setProperty('--f', (nv / LIM.dose_max * 100) + '%');
        renderPresets();
      };
      row.appendChild(rng);
      row.appendChild(out);
      box.appendChild(row);
    });
  }

  // ---------------------------------------------------------------- what was hiding; attempts
  function renderHiding() {
    var box = $('t-hiding'), r = rec(S.pid), R = S.result && S.revealed ? S.result : null;
    var truth = R ? R.truth : r.truth, best = R ? R.best : r.best;
    box.textContent = '';
    if (!truth || !r.attempts.length || (S.result && !S.revealed && !r.seen)) {
      box.appendChild(el('h2', null, 'What was hiding'));
      box.appendChild(el('p', 't-k t-locked', 'Shown after your first treatment: every clone in this tumour and whether sequencing could see it.'));
      return;
    }
    box.appendChild(el('h2', null, 'What was hiding'));
    var t = '<table class="t-clones"><thead><tr><th>clone in the biopsy</th><th class="n">share</th><th>found by</th>'
      + (R ? '<th class="n">day 364</th>' : '') + '</tr></thead><tbody>';
    truth.clones.forEach(function (c, i) {
      var seen = c.on_chart ? 'the chart' : c.share >= truth.sequencing_lod ? 'lab sequencing' : '<span class="t-no">nothing</span>';
      t += '<tr><td><i class="sw sw-c' + i + '"></i>' + esc(c.label) + '<small>' + esc(i ? own(c) : c.alterations.join(', ')) + '</small></td><td class="n">'
        + esc(pct(c.share)) + '</td><td>' + seen + '</td>' + (R ? '<td class="n">' + esc(cells(R.truth.clones[i].end)) + '</td>' : '') + '</tr>';
    });
    t += '</tbody></table>';
    box.insertAdjacentHTML('beforeend', t);
    box.appendChild(el('p', 't-k', 'The chart reports clones above about 10' + NBSP + '% of the tumour; the lab’s deep sequencing reports those above '
      + pct(truth.sequencing_lod) + '. Clones below that are invisible to both, and each is already resistant to something.'));
    if (best && best.regimen) {
      var bb = el('div', 't-bestline');
      bb.innerHTML = '<span>Best regimen we know for this patient: <b>' + esc(shortReg(best.regimen)) + '</b>, reward <b class="num">' + best.reward.toFixed(2) + '</b>.</span>';
      var go = el('button', 't-chip t-best', 'Try it');
      go.type = 'button';
      go.onclick = function () { setRegimen(best.regimen); renderPresets(); treat(); };
      bb.appendChild(go);
      box.appendChild(bb);
    }
  }
  function renderAttempts() {
    var box = $('t-att'), r = rec(S.pid);
    box.textContent = '';
    var h = el('div', 't-h');
    h.appendChild(el('h2', null, 'Your attempts'));
    var so = el('button', 't-reset', 'Start over');
    so.type = 'button';
    so.disabled = !r.attempts.length && !S.result;
    so.onclick = startOver;
    h.appendChild(so);
    box.appendChild(h);
    if (!r.attempts.length) {
      box.appendChild(el('p', 't-k t-locked', 'None yet for this patient. They are kept in this browser only.'));
      return;
    }
    var t = el('table', 't-atts');
    t.innerHTML = '<thead><tr><th>regimen</th><th class="n">healed</th><th class="n">harm</th><th class="n">reward</th></tr></thead>';
    var tb = el('tbody'), list = r.attempts.slice();
    if (S.result && !S.revealed) list.pop();
    if (!list.length) { box.appendChild(el('p', 't-k t-locked', 'Your first attempt appears here on day 364.')); return; }
    list.reverse().forEach(function (a) {
      var tr = el('tr');
      var td = el('td');
      var b = el('button', 't-again', shortReg(a.regimen));
      b.type = 'button';
      b.title = 'Load this regimen and treat again';
      b.onclick = function () { setRegimen(a.regimen); renderPresets(); treat(); };
      td.appendChild(b);
      tr.appendChild(td);
      tr.insertAdjacentHTML('beforeend', '<td class="n">' + a.healed.toFixed(2) + '</td><td class="n">' + a.harm.toFixed(2) + '</td><td class="n t-rw' + (a.reward > 0 ? ' ok' : '') + '">' + a.reward.toFixed(2) + '</td>');
      tb.appendChild(tr);
    });
    t.appendChild(tb);
    box.appendChild(t);
  }
  function startOver() {
    pause();
    MEM[S.pid] = { attempts: [] };
    saveStore(MEM);
    clearResult();
    setRegimen(patient(S.pid).guideline[0].regimen);
    renderAll();
    $('t-msg').textContent = 'Cleared. The patient is untreated again.';
    $('t-go').focus();
  }

  // ---------------------------------------------------------------- treat
  function clearResult() {
    S.result = null; S.revealed = false; S.frames = null; S.days = null; S.fi = 0;
    $('t-year').classList.remove('revealed');
    $('t-play').disabled = true;
    if (LT.shown) LT.shown.fill(255);
    setDay(0);
  }
  var reqToken = 0;
  function treat() {
    if (S.busy) return;
    var reg = regimenFromState(), pid = S.pid, token = ++reqToken;
    pause();
    S.busy = true;
    var go = $('t-go'), msg = $('t-msg');
    go.disabled = true; go.textContent = 'Treating…';
    msg.textContent = 'Simulating 364 days in ' + API.reps + ' copies of the lesion.';
    fetch('api/treat', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ patient: pid, regimen: reg }) })
      .then(function (res) {
        return res.json().catch(function () { return {}; }).then(function (body) {
          if (!res.ok) throw new Error(body.error || (res.status === 429 ? 'The simulator is busy; try again in a few seconds.' : 'The simulator did not answer (' + res.status + ').'));
          return body;
        });
      })
      .then(function (R) {
        if (token !== reqToken || pid !== S.pid) return;
        var n = R.replay.days.length, raw = inflate(R.replay.grids, n * 4096);
        S.frames = [];
        for (var k = 0; k < n; k++) S.frames.push(raw.subarray(k * 4096, (k + 1) * 4096));
        S.days = R.replay.days;
        S.result = R; S.revealed = false; S.fi = 0;
        var r = rec(pid);
        r.attempts.push({ regimen: R.regimen, healed: R.healed, harm: R.harm, reward: R.reward, cured: R.cured, t: Date.now() });
        if (r.attempts.length > 30) r.attempts = r.attempts.slice(-30);
        r.truth = { clones: R.truth.clones.map(function (c) { return { label: c.label, alterations: c.alterations, own: c.own, share: c.share, on_chart: c.on_chart }; }), sequencing_lod: R.truth.sequencing_lod };
        r.best = R.best;
        saveStore(MEM);
        msg.textContent = '';
        $('t-play').disabled = false;
        if (LT.shown) LT.shown.fill(255);
        renderCounts(); renderYear(); renderOutcome(); renderHiding();
        setDay(0);
        var dish = document.querySelector('.t-dish').getBoundingClientRect();
        if (dish.top < 0 || dish.top > window.innerHeight * 0.5) {
          try { document.querySelector('.t-dish').scrollIntoView({ behavior: 'smooth', block: 'start' }); } catch (e) { /* old browsers */ }
        }
        if (matchMedia('(prefers-reduced-motion: reduce)').matches) setDay(DAYS); else play();
      })
      .catch(function (err) {
        if (token !== reqToken) return;
        msg.textContent = String(err && err.message || err).replace(/^TypeError: .*/, 'The simulator did not answer. Is python3 viewer/serve.py running?');
      })
      .then(function () {
        if (token !== reqToken) return;
        S.busy = false;
        go.disabled = false; go.textContent = 'Treat';
      });
  }

  // ---------------------------------------------------------------- patient switch, wiring
  function choose(pid) {
    var p = patient(pid) || API.patients[0];
    if (S.pid === p.id && S.phases.length) return;
    reqToken++; S.busy = false;
    $('t-go').disabled = false; $('t-go').textContent = 'Treat'; $('t-msg').textContent = '';
    pause();
    S.pid = p.id;
    $('t-patient').value = p.id;
    var raw = inflate(p.day0.grid, 4096);
    S.day0 = raw.subarray(0, 4096);
    clearResult();
    setRegimen(p.guideline[0].regimen);
    renderAll();
    onHash();
  }
  function renderAll() {
    renderPatient(); renderPresets(); renderCounts(); renderOutcome(); renderYear(); renderHiding(); renderAttempts();
    if (LT.shown) LT.shown.fill(255);
    setDay(S.day);
  }

  function offline(why) {
    S.offline = true;
    $('t-grid').hidden = true;
    var o = $('t-off');
    o.hidden = false;
    o.innerHTML = esc(why) + ' Start it with <code>python3 viewer/serve.py</code> and open <code>http://localhost:8377/#view=treat</code>.';
  }

  function mount() {
    if (S.mounted) return;
    S.mounted = true;
    if (location.protocol === 'file:') { offline('Treating a patient needs the simulator, which this page cannot reach from a file.'); return; }
    var sp = $('t-speed');
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
    $('t-play').onclick = function () { S.playing ? pause() : play(); };
    $('t-go').onclick = treat;
    $('t-patient').onchange = function () { choose($('t-patient').value); };
    bindYear();
    document.addEventListener('keydown', function (evn) {
      if (!S.visible || !S.result || evn.metaKey || evn.ctrlKey || evn.altKey) return;
      var tg = evn.target && evn.target.tagName;
      if (/^(SELECT|INPUT|TEXTAREA)$/.test(tg)) return;
      if (evn.key === ' ' || evn.code === 'Space') {
        if (tg === 'BUTTON' && evn.target.id !== 't-play') return;
        evn.preventDefault(); S.playing ? pause() : play();
      } else if (evn.key === 'Home') { evn.preventDefault(); pause(); setDay(0); }
      else if (evn.key === 'End') { evn.preventDefault(); pause(); setDay(DAYS); }
      else if (evn.target && evn.target.id === 't-year-hit') return;
      else if (evn.key === 'ArrowRight' && tg !== 'BUTTON') { evn.preventDefault(); pause(); setDay(Math.floor(S.day) + (evn.shiftKey ? 7 : 1)); }
      else if (evn.key === 'ArrowLeft' && tg !== 'BUTTON') { evn.preventDefault(); pause(); setDay(Math.ceil(S.day) - (evn.shiftKey ? 7 : 1)); }
    });
    var rz = 0;
    window.addEventListener('resize', function () {
      clearTimeout(rz);
      rz = setTimeout(function () {
        if (!S.visible || !API) return;
        if (layout()) renderCounts();
        renderYear();
        setDay(S.day);
      }, 80);
    });
    S.ready = fetch('api/patients', { headers: { Accept: 'application/json' } })
      .then(function (res) { if (!res.ok) throw new Error(res.status); return res.json(); })
      .then(function (data) {
        if (!data || !data.patients || !data.patients.length) throw new Error('empty');
        API = data; DRUGS = data.drugs.map(function (d) { return d.name; }); LIM = data.limits;
        var sel = $('t-patient');
        data.patients.forEach(function (p) {
          var o = el('option', null, 'Patient ' + p.id + ' · ' + p.tumour_type);
          o.value = p.id;
          sel.appendChild(o);
        });
        $('t-grid').hidden = false;
      })
      .catch(function () { offline('The simulator is not answering, so patients cannot be treated here.'); throw new Error('offline'); });
  }

  global.RemissionTreat = {
    onHash: function (fn) { onHash = fn; },
    patientId: function () { return S.pid || S.want || 'dev-01'; },
    show: function (params) {
      mount();
      S.visible = true;
      if (params && params.patient) S.want = params.patient;
      if (S.offline || !S.ready) return;
      S.ready.then(function () {
        if (!S.visible) return;
        layout();
        var want = patient(S.want) ? S.want : S.pid || API.patients[0].id;
        if (want !== S.pid) choose(want);
        else { renderYear(); setDay(S.day); }
      }, function () {});
    },
    hide: function () { S.visible = false; pause(); hideTip(); }
  };
})(window);
