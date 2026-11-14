/* Data registry and a small zlib inflate, so the viewer needs no network and no fetch():
   every data file is a <script> that calls REMISSION_DATA.index(...) or REMISSION_DATA.add(...). */
(function (global) {
  'use strict';

  var LBASE = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
  var LEXT = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
  var DBASE = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
  var DEXT = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
  var ORDER = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

  function huffman(lengths) {
    var count = new Uint16Array(16), offs = new Uint16Array(16), symbol = new Uint16Array(lengths.length), i;
    for (i = 0; i < lengths.length; i++) count[lengths[i]]++;
    count[0] = 0;
    for (i = 1; i < 16; i++) offs[i] = offs[i - 1] + count[i - 1];
    for (i = 0; i < lengths.length; i++) if (lengths[i]) symbol[offs[lengths[i]]++] = i;
    return { count: count, symbol: symbol };
  }

  /* zlib stream (as written by Python's zlib.compress) -> Uint8Array */
  function inflate(src, sizeHint) {
    var ip = 2, bitbuf = 0, bitcnt = 0, op = 0;
    var out = new Uint8Array(sizeHint || 1 << 17);

    function bits(n) {
      var v = bitbuf;
      while (bitcnt < n) { v |= src[ip++] << bitcnt; bitcnt += 8; }
      bitbuf = v >>> n; bitcnt -= n;
      return v & ((1 << n) - 1);
    }
    function decode(h) {
      var code = 0, first = 0, index = 0, len, c;
      for (len = 1; len < 16; len++) {
        code |= bits(1);
        c = h.count[len];
        if (code - c < first) return h.symbol[index + (code - first)];
        index += c; first += c; first <<= 1; code <<= 1;
      }
      throw new Error('inflate: bad code');
    }
    function room(n) {
      if (op + n <= out.length) return;
      var bigger = new Uint8Array(Math.max(out.length * 2, op + n));
      bigger.set(out);
      out = bigger;
    }

    var last, type, i, len, sym, dist, lcode, dcode, fixedL, fixedD;
    do {
      last = bits(1); type = bits(2);
      if (type === 0) {
        bitbuf = 0; bitcnt = 0;
        len = src[ip] | (src[ip + 1] << 8); ip += 4;
        room(len);
        out.set(src.subarray(ip, ip + len), op); op += len; ip += len;
        continue;
      }
      if (type === 1) {
        if (!fixedL) {
          var l = new Uint8Array(288);
          for (i = 0; i < 144; i++) l[i] = 8;
          for (; i < 256; i++) l[i] = 9;
          for (; i < 280; i++) l[i] = 7;
          for (; i < 288; i++) l[i] = 8;
          fixedL = huffman(l);
          var d = new Uint8Array(30);
          for (i = 0; i < 30; i++) d[i] = 5;
          fixedD = huffman(d);
        }
        lcode = fixedL; dcode = fixedD;
      } else if (type === 2) {
        var nlen = bits(5) + 257, ndist = bits(5) + 1, ncode = bits(4) + 4;
        var lengths = new Uint8Array(320), cl = new Uint8Array(19);
        for (i = 0; i < ncode; i++) cl[ORDER[i]] = bits(3);
        var lencode = huffman(cl);
        i = 0;
        while (i < nlen + ndist) {
          sym = decode(lencode);
          if (sym < 16) { lengths[i++] = sym; continue; }
          var prev = 0, rep;
          if (sym === 16) { prev = lengths[i - 1]; rep = 3 + bits(2); }
          else if (sym === 17) rep = 3 + bits(3);
          else rep = 11 + bits(7);
          while (rep--) lengths[i++] = prev;
        }
        lcode = huffman(lengths.subarray(0, nlen));
        dcode = huffman(lengths.subarray(nlen, nlen + ndist));
      } else {
        throw new Error('inflate: bad block type');
      }
      for (;;) {
        sym = decode(lcode);
        if (sym < 256) { room(1); out[op++] = sym; continue; }
        if (sym === 256) break;
        sym -= 257;
        len = LBASE[sym] + bits(LEXT[sym]);
        sym = decode(dcode);
        dist = DBASE[sym] + bits(DEXT[sym]);
        room(len);
        for (i = 0; i < len; i++, op++) out[op] = out[op - dist];
      }
    } while (!last);
    return out.subarray(0, op);
  }

  function b64bytes(text) {
    var bin = atob(text), n = bin.length, u8 = new Uint8Array(n), i;
    for (i = 0; i < n; i++) u8[i] = bin.charCodeAt(i);
    return u8;
  }

  var waiting = {};
  global.REMISSION_DATA = {
    idx: null,
    runs: {},
    index: function (o) { this.idx = o; },
    add: function (o) {
      var key = o.patient + '.' + o.id;
      var L = o.replay.L, n = o.replay.days.length;
      var raw = inflate(b64bytes(o.replay.grids), n * L * L);
      o.frames = [];
      for (var k = 0; k < n; k++) o.frames.push(raw.subarray(k * L * L, (k + 1) * L * L));
      delete o.replay.grids;
      this.runs[key] = o;
      (waiting[key] || []).forEach(function (fn) { fn(o); });
      delete waiting[key];
    },
    /* load one run's data file by <script>; resolves with the run */
    load: function (patient, id, file) {
      var self = this, key = patient + '.' + id;
      if (self.runs[key]) return Promise.resolve(self.runs[key]);
      return new Promise(function (resolve, reject) {
        var first = !waiting[key];
        (waiting[key] = waiting[key] || []).push(resolve);
        if (!first) return;
        var s = document.createElement('script');
        s.src = 'data/' + file;
        s.onerror = function () { delete waiting[key]; reject(new Error('could not load data/' + file)); };
        document.head.appendChild(s);
      });
    },
    inflate: inflate
  };
})(window);
