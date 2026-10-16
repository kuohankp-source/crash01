/* =========================================================
   正向彈性碰撞—微觀碰撞過程與質心座標系分析器
   物理模型：簡諧運動（SHM）解析解，力學能 100% 精確守恆
   ========================================================= */
'use strict';

(function () {

  // ---------- 錯誤顯示（方便課堂現場除錯） ----------
  window.addEventListener('error', (e) => {
    const b = document.getElementById('errBanner');
    if (b) { b.style.display = 'block'; b.textContent = '程式錯誤：' + e.message + '\n' + (e.filename || '') + ':' + (e.lineno || ''); }
  });

  // =====================================================
  // 一、物理參數設定
  // =====================================================
  const SCENARIOS = {
    chase: { v1: 3.0, v2: 1.0, name: '同向追撞' },     // v1 > v2 > 0
    still: { v1: 3.0, v2: 0.0, name: '一動撞一靜' },   // v1 > 0, v2 = 0
    head:  { v1: 2.5, v2: -1.5, name: '迎面對撞' },    // v1 > 0, v2 < 0
  };
  const MASSES = {
    equal:   { m1: 1, m2: 1, l1: 'm',  l2: 'm' },
    unequal: { m1: 2, m2: 1, l1: '2m', l2: 'm' },
  };

  const SQUASH = 0.34;          // 最大形變時，每顆球接觸面內縮量（以 R 為單位）≈ 34%
  const X_MAX = 2 * SQUASH;     // 系統最大壓縮量 x_max（兩球合計）
  const BULGE = 0.45;           // 麻糬上下外鼓係數
  const TOP_H = 1 + BULGE * SQUASH; // 球體最大半高（R）
  const P0 = 15, PM = 50, P1 = 85;  // 剛接觸、最大形變、剛分離
  const VK = 0.36;              // 速度箭頭比例：每 1 m/s 對應 0.36 R
  const FK = 0.9;               // 最大彈力箭頭長度（R）
  const BASE_DURATION = 5;      // 1.0x 時，0% → 100% 所需秒數

  const COL = {
    v: '#22d3ee', vp: '#ff4fd8', F: '#84ff57', vc: '#ffc940',
    grid: 'rgba(120,160,230,0.10)', gridMajor: 'rgba(120,160,230,0.22)',
  };
  const BALL_COLORS = {
    equal:   [[59, 130, 246], [249, 115, 22]],
    unequal: [[30, 64, 175],  [253, 186, 116]], // 重球顏色較深、輕球較淺
  };

  // =====================================================
  // 二、狀態
  // =====================================================
  const state = {
    scenario: 'chase',
    mass: 'equal',
    mode: 'A',
    layers: { 1: true, 2: false, 3: false },
    p: 0,
    playing: false,
    speed: 0.25,
  };
  let model = null;

  // =====================================================
  // 三、物理模型
  // =====================================================
  function buildModel() {
    const sc = SCENARIOS[state.scenario], ms = MASSES[state.mass];
    const m1 = ms.m1, m2 = ms.m2, M = m1 + m2;
    const v1i = sc.v1, v2i = sc.v2;
    const vc = (m1 * v1i + m2 * v2i) / M;
    const vrel0 = v1i - v2i;                 // 接近速率
    const omega = vrel0 / X_MAX;             // 自動調整等效彈性係數：k = μω²，確保 x_max 固定
    const Tc = Math.PI / omega;              // 接觸時間（半週期）
    const tPerP = Tc / (P1 - P0);            // 每 1% 進度對應的時間
    const E = 0.5 * m1 * v1i * v1i + 0.5 * m2 * v2i * v2i;
    const Kc = 0.5 * M * vc * vc;
    const mu = m1 * m2 / M;
    const KintMax = 0.5 * mu * vrel0 * vrel0;
    const v1f = vc - (m2 / M) * vrel0;
    const v2f = vc + (m1 / M) * vrel0;

    const m = {
      m1, m2, M, v1i, v2i, vc, vrel0, omega, Tc, tPerP, E, Kc, KintMax, v1f, v2f,
      l1: ms.l1, l2: ms.l2,
    };
    m.Xc50 = vc * (PM - P0) * tPerP;

    // 計算畫面所需水平範圍（以 R 為單位，含箭頭）
    let ext = 0;
    for (let p = 0; p <= 100; p += 0.5) {
      const S = stateAt(m, p);
      const xs = [
        S.x1 - 1.05, S.x2 + 1.05,
        S.x1 + S.v1 * VK, S.x2 + S.v2 * VK,
        S.x1 + S.v1p * VK, S.x2 + S.v2p * VK,
        S.Xc + vc * VK,
      ];
      for (const cam of [m.Xc50, S.Xc]) {
        for (const X of xs) ext = Math.max(ext, Math.abs(X - cam));
      }
    }
    m.ext = ext;
    return m;
  }

  // 依進度 p（0~100）回傳瞬時狀態；長度單位為 R，速度單位為 m/s
  function stateAt(m, p) {
    const t = (p - P0) * m.tPerP;          // t = 0 為剛接觸瞬間
    let x = 0, d, vrel, theta = null;
    if (p < P0) {                          // 接觸前：等速接近
      d = 2 - m.vrel0 * t;
      vrel = m.vrel0;
    } else if (p <= P1) {                  // 接觸中：SHM 解析解
      theta = (p - P0) / (P1 - P0) * Math.PI;
      x = X_MAX * Math.sin(theta);
      d = 2 - x;
      vrel = m.vrel0 * Math.cos(theta);
    } else {                               // 分離後：等速遠離
      d = 2 + m.vrel0 * (t - m.Tc);
      vrel = -m.vrel0;
    }
    const Xc = m.vc * t;
    const x1 = Xc - (m.m2 / m.M) * d;
    const x2 = Xc + (m.m1 / m.M) * d;
    const v1 = m.vc + (m.m2 / m.M) * vrel;
    const v2 = m.vc - (m.m1 / m.M) * vrel;
    const v1p = v1 - m.vc, v2p = v2 - m.vc;

    const K1 = 0.5 * m.m1 * v1 * v1;
    const K2 = 0.5 * m.m2 * v2 * v2;
    const Kt = K1 + K2;
    const Kint = 0.5 * m.m1 * v1p * v1p + 0.5 * m.m2 * v2p * v2p;
    const s = x / X_MAX;                       // = sin θ
    const Us = m.KintMax * s * s;              // = K_int,max · sin²θ
    const pc = (v) => v / m.E * 100;

    return {
      p, t, x, s, d, theta, vrel, Xc, x1, x2, v1, v2, v1p, v2p,
      contactX: (x1 + x2) / 2,
      inContact: p >= P0 && p <= P1,
      K1: pc(K1), K2: pc(K2), Kt: pc(Kt), Kint: pc(Kint), Kc: pc(m.Kc), Us: pc(Us), E: pc(Kt + Us),
    };
  }

  // =====================================================
  // 四、工具函式
  // =====================================================
  const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
  const fmtV = (v) => (Math.abs(v) < 0.005 ? '0' : (v > 0 ? '+' : '−') + Math.abs(v).toFixed(2));
  const fmtP = (v) => (Math.abs(v) < 0.05 ? '0.0' : v.toFixed(1)) + '%';
  const rgb = (c, a = 1) => `rgba(${c[0] | 0},${c[1] | 0},${c[2] | 0},${a})`;
  const mix = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
  const lighten = (c, t) => mix(c, [255, 255, 255], t);
  const darken = (c, t) => mix(c, [0, 0, 0], t);
  const lum = (c) => (0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]) / 255;
  const hex = (c) => '#' + c.map((v) => ('0' + Math.round(v).toString(16)).slice(-2)).join('');

  function roundRectPath(ctx, x, y, w, h, r) {
    r = Math.min(r, w / 2, h / 2);
    ctx.beginPath();
    ctx.moveTo(x + r, y);
    ctx.lineTo(x + w - r, y); ctx.quadraticCurveTo(x + w, y, x + w, y + r);
    ctx.lineTo(x + w, y + h - r); ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
    ctx.lineTo(x + r, y + h); ctx.quadraticCurveTo(x, y + h, x, y + h - r);
    ctx.lineTo(x, y + r); ctx.quadraticCurveTo(x, y, x + r, y);
    ctx.closePath();
  }

  // =====================================================
  // 五、Canvas 繪圖元件
  // =====================================================

  /**
   * 「麻糬式」形變球體路徑（四段三次貝茲曲線 + 平截接觸面）
   * delta：接觸面內縮量（px）；dir = +1 接觸面在右側（球1），-1 在左側（球2）
   * delta = 0 時精確退化為標準貝茲正圓。
   */
  function mochiPath(ctx, cx, cy, R, delta, dir) {
    const k = 0.5523;
    const t = delta / R;
    const a = R - delta;                                         // 球心到平截面距離
    const h = Math.sqrt(Math.max(0, R * R - a * a)) * (1 + 0.25 * t); // 平截面半高
    const H = R * (1 + BULGE * t);                               // 上下外鼓
    const Rb = R * (1 - 0.05 * t);                               // 背面微縮
    const bx = 0.25 * delta;                                     // 外鼓峰值略偏向接觸面
    const P = (X, Y) => [cx + dir * X, cy + Y];

    ctx.beginPath();
    ctx.moveTo(...P(a, -h));
    ctx.lineTo(...P(a, h));
    ctx.bezierCurveTo(...P(a, h + k * (H - h)), ...P(bx + k * (a - bx), H), ...P(bx, H));
    ctx.bezierCurveTo(...P(bx - k * (bx + Rb), H), ...P(-Rb, k * H), ...P(-Rb, 0));
    ctx.bezierCurveTo(...P(-Rb, -k * H), ...P(bx - k * (bx + Rb), -H), ...P(bx, -H));
    ctx.bezierCurveTo(...P(bx + k * (a - bx), -H), ...P(a, -h - k * (H - h)), ...P(a, -h));
    ctx.closePath();
    return { a, h, H };
  }

  function drawBall(ctx, cx, cy, R, delta, dir, base) {
    ctx.save();
    mochiPath(ctx, cx, cy, R, delta, dir);
    const g = ctx.createRadialGradient(cx - R * 0.35, cy - R * 0.42, R * 0.08, cx, cy, R * 1.25);
    g.addColorStop(0, rgb(lighten(base, 0.6)));
    g.addColorStop(0.45, rgb(base));
    g.addColorStop(1, rgb(darken(base, 0.5)));
    ctx.fillStyle = g;
    ctx.shadowColor = rgb(base, 0.45);
    ctx.shadowBlur = R * 0.25;
    ctx.fill();
    ctx.shadowBlur = 0;
    ctx.lineWidth = Math.max(1.5, R * 0.022);
    ctx.strokeStyle = rgb(lighten(base, 0.45), 0.95);
    ctx.stroke();
    ctx.restore();
  }

  // 應力發光熱區：強度 ∝ U_s（= sin²θ），顏色由黃漸變至亮紅
  function stressColor(s2) { return mix([255, 232, 80], [255, 35, 30], clamp(s2, 0, 1)); }

  function drawStressInside(ctx, cx, cy, R, delta, dir, s2) {
    if (s2 < 0.002) return;
    ctx.save();
    mochiPath(ctx, cx, cy, R, delta, dir);
    ctx.clip();
    const faceX = cx + dir * (R - delta);
    const col = stressColor(s2);
    const g = ctx.createLinearGradient(faceX, 0, faceX - dir * R * 0.85, 0);
    g.addColorStop(0, rgb(col, 0.95 * Math.sqrt(s2)));
    g.addColorStop(0.3, rgb(mix(col, [255, 200, 60], 0.4), 0.55 * s2));
    g.addColorStop(1, rgb(col, 0));
    ctx.fillStyle = g;
    ctx.fillRect(cx - R * 1.3, cy - R * 1.4, R * 2.6, R * 2.8);
    ctx.restore();
  }

  function drawContactGlow(ctx, faceX, cy, h, R, s2) {
    if (s2 < 0.002) return;
    const col = stressColor(s2);
    ctx.save();
    ctx.globalCompositeOperation = 'lighter';
    // 外圍橢圓光暈
    ctx.translate(faceX, cy);
    ctx.scale(0.45, 1);
    const rr = Math.max(h * 1.35, 4);
    const g = ctx.createRadialGradient(0, 0, 0, 0, 0, rr);
    g.addColorStop(0, rgb(col, 0.55 * s2));
    g.addColorStop(1, rgb(col, 0));
    ctx.fillStyle = g;
    ctx.beginPath(); ctx.arc(0, 0, rr, 0, Math.PI * 2); ctx.fill();
    ctx.restore();
    // 接觸平截面亮線
    ctx.save();
    ctx.strokeStyle = rgb(lighten(col, 0.35), 0.5 + 0.5 * s2);
    ctx.shadowColor = rgb(col, 1);
    ctx.shadowBlur = 8 + 22 * s2;
    ctx.lineWidth = 2 + 3.5 * s2 * (R / 100);
    ctx.lineCap = 'round';
    ctx.beginPath(); ctx.moveTo(faceX, cy - h); ctx.lineTo(faceX, cy + h); ctx.stroke();
    ctx.restore();
  }

  function drawArrow(ctx, x0, y, len, color, lw, alpha = 1) {
    const dir = Math.sign(len), L = Math.abs(len);
    ctx.save();
    ctx.globalAlpha = alpha;
    ctx.strokeStyle = color; ctx.fillStyle = color;
    ctx.lineCap = 'round'; ctx.lineJoin = 'round';
    ctx.shadowColor = color; ctx.shadowBlur = 8;
    if (L < 3) {
      ctx.beginPath(); ctx.arc(x0, y, lw * 1.2, 0, Math.PI * 2); ctx.fill();
      ctx.restore(); return;
    }
    const hl = Math.min(L * 0.6, lw * 3.4 + 6);
    const hw = lw * 2.0 + 3;
    ctx.lineWidth = lw;
    ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x0 + dir * (L - hl * 0.7), y); ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x0 + dir * L, y);
    ctx.lineTo(x0 + dir * (L - hl), y - hw);
    ctx.lineTo(x0 + dir * (L - hl), y + hw);
    ctx.closePath(); ctx.fill();
    // 起點小圓
    ctx.beginPath(); ctx.arc(x0, y, lw * 0.75, 0, Math.PI * 2); ctx.fill();
    ctx.restore();
  }

  /**
   * 物理符號標籤：spec = { pre, b, prime, sub, vec, val, post, bg }
   * 例：{b:'v', prime:true, sub:'1', vec:true, val:'+1.00'} → v⃗′₁ = +1.00
   */
  function drawSym(ctx, spec, x, y, align, color, fs, alpha = 1) {
    const fSym = `italic 700 ${fs}px "Times New Roman","Cambria Math",serif`;
    const fSub = `700 ${Math.round(fs * 0.62)}px "Times New Roman",serif`;
    const fTxt = `600 ${Math.round(fs * 0.78)}px "Microsoft JhengHei","Noto Sans TC","PingFang TC",sans-serif`;
    const fVal = `700 ${Math.round(fs * 0.8)}px "Segoe UI","Microsoft JhengHei",sans-serif`;
    const items = [];
    if (spec.pre) items.push({ t: spec.pre, f: fTxt, dy: 0 });
    if (spec.b) items.push({ t: spec.b, f: fSym, dy: 0, base: true });
    if (spec.prime) items.push({ t: '′', f: fSym, dy: 0 });
    if (spec.sub) items.push({ t: spec.sub, f: fSub, dy: fs * 0.26 });
    if (spec.val != null) items.push({ t: ' = ' + spec.val, f: fVal, dy: 0 });
    if (spec.post) items.push({ t: spec.post, f: fTxt, dy: 0 });

    ctx.save();
    ctx.globalAlpha = alpha;
    let total = 0;
    for (const it of items) { ctx.font = it.f; it.w = ctx.measureText(it.t).width; total += it.w; }
    let cx = align === 'right' ? x - total : align === 'center' ? x - total / 2 : x;
    const base = y + fs * 0.33;
    if (spec.bg) {
      ctx.fillStyle = 'rgba(5,9,18,0.72)';
      roundRectPath(ctx, cx - 5, y - fs * 0.9, total + 10, fs * 1.6, 6);
      ctx.fill();
    }
    ctx.fillStyle = color; ctx.strokeStyle = color;
    ctx.textBaseline = 'alphabetic';
    ctx.shadowColor = 'rgba(0,0,0,0.9)'; ctx.shadowBlur = 4;
    for (const it of items) {
      ctx.font = it.f;
      ctx.fillText(it.t, cx, base + it.dy);
      if (it.base && spec.vec) {
        const ax0 = cx + fs * 0.1, ax1 = cx + it.w + fs * 0.14, ay = base - fs * 0.84, hs = fs * 0.16;
        ctx.lineWidth = Math.max(1.2, fs * 0.075);
        ctx.beginPath(); ctx.moveTo(ax0, ay); ctx.lineTo(ax1 - hs * 0.4, ay); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(ax1, ay); ctx.lineTo(ax1 - hs, ay - hs * 0.6); ctx.lineTo(ax1 - hs, ay + hs * 0.6); ctx.closePath(); ctx.fill();
      }
      cx += it.w;
    }
    ctx.restore();
    return total;
  }

  // 質心符號：黑白（金）四象限圓 + 十字線
  function drawCM(ctx, x, y, R, fs) {
    const r = Math.max(9, R * 0.11);
    ctx.save();
    ctx.shadowColor = 'rgba(255,201,64,0.9)'; ctx.shadowBlur = 10;
    ctx.strokeStyle = '#fff'; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x - r * 1.8, y); ctx.lineTo(x + r * 1.8, y);
    ctx.moveTo(x, y - r * 1.8); ctx.lineTo(x, y + r * 1.8); ctx.stroke();
    ctx.shadowBlur = 0;
    for (let i = 0; i < 4; i++) {
      ctx.beginPath(); ctx.moveTo(x, y);
      ctx.arc(x, y, r, -Math.PI / 2 + i * Math.PI / 2, i * Math.PI / 2);
      ctx.closePath();
      ctx.fillStyle = i % 2 === 0 ? COL.vc : '#111827';
      ctx.fill();
    }
    ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.stroke();
    ctx.restore();
    drawSym(ctx, { b: 'CM', bg: true }, x, y - r - fs * 0.85, 'center', COL.vc, fs * 0.95);
  }

  // =====================================================
  // 六、版面配置
  // =====================================================
  const canvas = document.getElementById('cv');
  const ctx = canvas.getContext('2d');
  const wrap = document.getElementById('canvasWrap');
  let cssW = 800, cssH = 400, dpr = 1;

  function resizeCanvas() {
    const r = wrap.getBoundingClientRect();
    dpr = window.devicePixelRatio || 1;
    cssW = Math.max(100, r.width); cssH = Math.max(100, r.height);
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
  }

  function computeLayout() {
    const W = cssW, H = cssH;
    const split = state.mode === 'AB';
    const panels = split
      ? [{ x: 0, y: 0, w: W, h: H / 2, frame: 'ground' }, { x: 0, y: H / 2, w: W, h: H / 2, frame: 'cm' }]
      : [{ x: 0, y: 0, w: W, h: H, frame: state.mode === 'A' ? 'ground' : 'cm' }];
    const ph = panels[0].h;
    let fs = clamp(Math.min(W, H * 1.8) * 0.014, 14, 24);
    const labelPx = fs * 7.6;
    const Rw = (W - 2 * labelPx - 24) / (2 * model.ext);
    const Rh = (ph - 1.8 * fs - 26) / 3.55;
    const R = clamp(Math.min(Rw, Rh), 26, 175);
    fs = clamp(R * 0.24, 13, 24);
    const g1 = Math.max(0.32 * R, fs * 1.15);
    const g2 = Math.max(0.46 * R, fs * 1.65);
    const g3 = Math.max(0.42 * R, fs * 1.5);
    const above = TOP_H * R + g1 + g2 + fs * 0.9;
    const below = TOP_H * R + g3 + fs * 0.9;
    for (const pan of panels) {
      pan.cy = pan.y + Math.max(0, (pan.h - above - below) / 2) + above;
      pan.row0 = pan.cy - TOP_H * R - g1;
      pan.row1 = pan.row0 - g2;
      pan.rowB = pan.cy + TOP_H * R + g3;
    }
    return { panels, R, fs, W, H };
  }

  // =====================================================
  // 七、單一視窗繪製
  // =====================================================
  function drawGrid(ctx, pan, toX, camX, cy, R) {
    const step = R * 0.5;
    ctx.save();
    // 垂直線（隨座標系移動）
    const Xl = camX - (pan.w / 2) / R - 1, Xr = camX + (pan.w / 2) / R + 1;
    const i0 = Math.floor(Xl * 2), i1 = Math.ceil(Xr * 2);
    for (let i = i0; i <= i1; i++) {
      const sx = toX(i / 2);
      ctx.strokeStyle = i % 2 === 0 ? COL.gridMajor : COL.grid;
      ctx.lineWidth = i % 2 === 0 ? 1.2 : 1;
      ctx.beginPath(); ctx.moveTo(sx, pan.y); ctx.lineTo(sx, pan.y + pan.h); ctx.stroke();
    }
    // 水平線
    const j0 = Math.floor((pan.y - cy) / step) - 1, j1 = Math.ceil((pan.y + pan.h - cy) / step) + 1;
    for (let j = j0; j <= j1; j++) {
      const sy = cy + j * step;
      ctx.strokeStyle = j % 2 === 0 ? COL.gridMajor : COL.grid;
      ctx.lineWidth = j % 2 === 0 ? 1.2 : 1;
      ctx.beginPath(); ctx.moveTo(pan.x, sy); ctx.lineTo(pan.x + pan.w, sy); ctx.stroke();
    }
    ctx.restore();
  }

  function renderPanel(ctx, pan, S, Lay) {
    const { R, fs } = Lay;
    const L1 = state.layers[1], L2 = state.layers[2];
    const isCM = pan.frame === 'cm';
    const camX = isCM ? S.Xc : model.Xc50;
    const cx0 = pan.x + pan.w / 2;
    const toX = (X) => cx0 + (X - camX) * R;
    const cy = pan.cy;
    const lw = clamp(R * 0.045, 3, 7);

    ctx.save();
    ctx.beginPath(); ctx.rect(pan.x, pan.y, pan.w, pan.h); ctx.clip();

    // 背景
    const bg = ctx.createLinearGradient(0, pan.y, 0, pan.y + pan.h);
    if (isCM) { bg.addColorStop(0, '#0d0f26'); bg.addColorStop(1, '#090b1c'); }
    else { bg.addColorStop(0, '#08132a'); bg.addColorStop(1, '#060e1f'); }
    ctx.fillStyle = bg; ctx.fillRect(pan.x, pan.y, pan.w, pan.h);

    drawGrid(ctx, pan, toX, camX, cy, R);

    // 碰撞軸
    ctx.save();
    ctx.strokeStyle = 'rgba(200,215,240,0.28)'; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(pan.x, cy); ctx.lineTo(pan.x + pan.w, cy); ctx.stroke();
    ctx.fillStyle = 'rgba(200,215,240,0.6)';
    ctx.font = `italic 700 ${Math.round(fs * 0.9)}px "Times New Roman",serif`;
    ctx.textAlign = 'right';
    ctx.fillText('+x →', pan.x + pan.w - 10, cy - 8);
    ctx.restore();

    // 質心座標系：中央十字線
    if (isCM) {
      ctx.save();
      ctx.strokeStyle = 'rgba(255,201,64,0.45)'; ctx.lineWidth = 1.5; ctx.setLineDash([8, 6]);
      ctx.beginPath(); ctx.moveTo(cx0, pan.y); ctx.lineTo(cx0, pan.y + pan.h); ctx.stroke();
      ctx.restore();
    } else if (L2) {
      // 對地座標系：質心追蹤線
      ctx.save();
      ctx.strokeStyle = 'rgba(255,201,64,0.35)'; ctx.lineWidth = 1.5; ctx.setLineDash([6, 6]);
      const sx = toX(S.Xc);
      ctx.beginPath(); ctx.moveTo(sx, pan.y); ctx.lineTo(sx, pan.y + pan.h); ctx.stroke();
      ctx.restore();
    }

    // ---- 球體 ----
    const deltaPx = (S.x / 2) * R;         // 每球接觸面內縮量
    const s2 = S.s * S.s;                  // ∝ U_s
    const c1 = toX(S.x1), c2 = toX(S.x2);
    const cols = BALL_COLORS[state.mass];
    drawBall(ctx, c1, cy, R, deltaPx, +1, cols[0]);
    drawBall(ctx, c2, cy, R, deltaPx, -1, cols[1]);
    if (S.inContact) {
      drawStressInside(ctx, c1, cy, R, deltaPx, +1, s2);
      drawStressInside(ctx, c2, cy, R, deltaPx, -1, s2);
      const a = R - deltaPx;
      const h = Math.sqrt(Math.max(0, R * R - a * a)) * (1 + 0.25 * deltaPx / R);
      drawContactGlow(ctx, toX(S.contactX), cy, h, R, s2);
    }

    // 質量標示（永遠顯示）
    const txt1 = lum(cols[0]) > 0.6 ? '#0b1020' : '#ffffff';
    const txt2 = lum(cols[1]) > 0.6 ? '#0b1020' : '#ffffff';
    drawSym(ctx, { b: 'm', sub: '1', val: model.l1 }, c1 - deltaPx * 0.4, cy - R * 0.45, 'center', txt1, fs * 1.0);
    drawSym(ctx, { b: 'm', sub: '2', val: model.l2 }, c2 + deltaPx * 0.4, cy - R * 0.45, 'center', txt2, fs * 1.0);

    // ---- 階段 1：球心、連線、彈力 ----
    if (L1) {
      ctx.save();
      ctx.strokeStyle = 'rgba(255,255,255,0.75)'; ctx.lineWidth = 1.6; ctx.setLineDash([5, 5]);
      ctx.beginPath(); ctx.moveTo(c1, cy); ctx.lineTo(c2, cy); ctx.stroke();
      ctx.setLineDash([]);
      for (const c of [c1, c2]) {
        ctx.beginPath(); ctx.arc(c, cy, Math.max(3.5, R * 0.045), 0, Math.PI * 2);
        ctx.fillStyle = '#fff'; ctx.shadowColor = '#000'; ctx.shadowBlur = 4; ctx.fill();
      }
      ctx.restore();
      drawSym(ctx, { b: 'C', sub: '1' }, c1 - R * 0.08, cy + R * 0.22, 'right', txt1, fs * 0.95);
      drawSym(ctx, { b: 'C', sub: '2' }, c2 + R * 0.08, cy + R * 0.22, 'left', txt2, fs * 0.95);

      if (S.inContact && S.s > 0.004) {
        const faceX = toX(S.contactX);
        const FL = FK * R * S.s;              // F = kx ∝ sinθ
        const yF = cy + R * 0.56;
        const flw = clamp(R * 0.05, 3, 8);
        drawArrow(ctx, faceX - 2, yF, -FL, COL.F, flw);
        drawArrow(ctx, faceX + 2, yF, +FL, COL.F, flw);
        if (S.s > 0.15) {
          drawSym(ctx, { b: 'F', sub: '12', vec: true }, faceX - FL / 2, yF + fs * 1.05, 'center', COL.F, fs * 0.95);
          drawSym(ctx, { b: 'F', sub: '21', vec: true }, faceX + FL / 2, yF + fs * 1.05, 'center', COL.F, fs * 0.95);
        }
      }
    }

    // ---- 質心標記 ----
    if (L2 || isCM) drawCM(ctx, toX(S.Xc), cy, R, fs);

    // ---- 速度向量列 ----
    const rows = [];
    if (!isCM) {
      if (L1) rows.push({ type: 'v' });
      if (L2) rows.push({ type: 'vp' });
    } else {
      if (L1 || L2) rows.push({ type: 'vp' });
      if (L2) rows.push({ type: 'v', faded: true });
    }
    rows.forEach((row, i) => {
      const y = i === 0 ? pan.row0 : pan.row1;
      const isV = row.type === 'v';
      const col = isV ? COL.v : COL.vp;
      const alpha = row.faded ? 0.6 : 1;
      const va = isV ? S.v1 : S.v1p, vb = isV ? S.v2 : S.v2p;
      const la = va * VK * R, lb = vb * VK * R;
      drawArrow(ctx, c1, y, la, col, lw, alpha);
      drawArrow(ctx, c2, y, lb, col, lw, alpha);
      const pad = lw * 2 + 6;
      const pre = row.faded ? '對地 ' : '';
      drawSym(ctx, { pre, b: 'v', prime: !isV, sub: '1', vec: true, val: fmtV(va) }, Math.min(c1, c1 + la) - pad, y, 'right', col, fs, alpha);
      drawSym(ctx, { pre, b: 'v', prime: !isV, sub: '2', vec: true, val: fmtV(vb) }, Math.max(c2, c2 + lb) + pad, y, 'left', col, fs, alpha);
    });

    // ---- 質心速度 / 方格網移動提示 ----
    if (!isCM && L2) {
      const sx = toX(S.Xc), lv = model.vc * VK * R;
      drawArrow(ctx, sx, pan.rowB, lv, COL.vc, lw);
      const right = lv >= 0;
      drawSym(ctx, { b: 'v', sub: 'c', vec: true, val: fmtV(model.vc), post: '（全程恆定）' },
        right ? Math.max(sx, sx + lv) + lw * 2 + 6 : Math.min(sx, sx + lv) - lw * 2 - 6, pan.rowB, right ? 'left' : 'right', COL.vc, fs);
    }
    if (isCM) {
      const lv = -model.vc * VK * R;
      const gcol = 'rgba(190,205,235,0.8)';
      drawArrow(ctx, cx0, pan.rowB, lv, gcol, Math.max(2, lw * 0.6));
      const txt = Math.abs(model.vc) < 0.005
        ? { pre: '此情境 ', b: 'v', sub: 'c', vec: true, val: '0', post: '：質心系與地面系重合' }
        : { pre: '地面方格以 −', b: 'v', sub: 'c', vec: true, val: fmtV(-model.vc), post: ' 移動' };
      drawSym(ctx, txt, lv <= 0 ? cx0 + 14 : cx0 - 14, pan.rowB, lv <= 0 ? 'left' : 'right', gcol, fs * 0.9);
    }

    // ---- 視窗標籤 ----
    const badge = isCM ? 'B｜質心座標系：CM 固定於中央，方格網以質心速度反向滑動' : 'A｜對地座標系：方格網靜止（地面觀察者）';
    ctx.save();
    ctx.font = `700 ${Math.round(fs * 0.78)}px "Microsoft JhengHei","Noto Sans TC",sans-serif`;
    const bw = ctx.measureText(badge).width + 22, bh = fs * 1.55;
    const bxp = pan.x + 10, byp = pan.y + pan.h - bh - 8;
    roundRectPath(ctx, bxp, byp, bw, bh, bh / 2);
    ctx.fillStyle = isCM ? 'rgba(80,40,120,0.75)' : 'rgba(20,70,110,0.75)';
    ctx.fill();
    ctx.strokeStyle = isCM ? 'rgba(255,79,216,0.7)' : 'rgba(34,211,238,0.7)'; ctx.lineWidth = 1.5; ctx.stroke();
    ctx.fillStyle = '#fff'; ctx.textBaseline = 'middle';
    ctx.fillText(badge, bxp + 11, byp + bh / 2 + 1);
    ctx.restore();

    // 圖例（右下）
    drawLegend(ctx, pan, fs, isCM);

    ctx.restore();
  }

  function drawLegend(ctx, pan, fs, isCM) {
    const L1 = state.layers[1], L2 = state.layers[2];
    const items = [];
    if (!isCM && L1) items.push({ c: COL.v, s: { b: 'v', vec: true, post: ' 對地速度' } });
    if ((isCM && (L1 || L2)) || (!isCM && L2)) items.push({ c: COL.vp, s: { b: 'v', prime: true, vec: true, post: ' 相對質心速度' } });
    if (L1) items.push({ c: COL.F, s: { b: 'F', sub: '12', vec: true, post: '：球2 對球1 的彈力' } });
    if (!items.length) return;
    const f = fs * 0.8;
    let y = pan.y + pan.h - 10 - f * 0.9;
    for (let i = items.length - 1; i >= 0; i--) {
      const it = items[i];
      ctx.save();
      ctx.font = `600 ${Math.round(f * 0.78)}px "Microsoft JhengHei",sans-serif`;
      ctx.restore();
      const w = drawSymMeasure(ctx, it.s, f);
      const x = pan.x + pan.w - 14 - w;
      ctx.save();
      ctx.fillStyle = it.c; ctx.fillRect(x - 26, y - 2, 18, 4);
      ctx.restore();
      drawSym(ctx, it.s, x, y, 'left', it.c, f, 0.95);
      y -= f * 1.55;
    }
  }

  function drawSymMeasure(ctx, spec, fs) {
    // 以全透明繪製取得寬度
    return drawSym(ctx, spec, -9999, -9999, 'left', 'rgba(0,0,0,0)', fs, 0);
  }

  function drawSplitOverlay(ctx, Lay) {
    const { W, H, fs } = Lay;
    ctx.save();
    // 分隔線
    ctx.fillStyle = '#22314f'; ctx.fillRect(0, H / 2 - 1.5, W, 3);
    // 質心對齊線（貫穿上下）
    ctx.strokeStyle = 'rgba(255,201,64,0.85)'; ctx.lineWidth = 2; ctx.setLineDash([10, 7]);
    ctx.shadowColor = 'rgba(255,201,64,0.8)'; ctx.shadowBlur = 6;
    ctx.beginPath(); ctx.moveTo(W / 2, 0); ctx.lineTo(W / 2, H); ctx.stroke();
    ctx.restore();
    drawSym(ctx, { pre: '▲▼ 質心對齊線（下：CM 固定於此；上：CM 於 50% 最大形變瞬間通過此線）', bg: true },
      W / 2 + 12, H / 2, 'left', COL.vc, fs * 0.82);
  }

  // =====================================================
  // 八、主繪製
  // =====================================================
  function draw(S) {
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    const Lay = computeLayout();
    for (const pan of Lay.panels) renderPanel(ctx, pan, S, Lay);
    if (state.mode === 'AB') drawSplitOverlay(ctx, Lay);
  }

  // =====================================================
  // 九、DOM 介面更新
  // =====================================================
  const $ = (s) => document.querySelector(s);
  const $$ = (s) => Array.from(document.querySelectorAll(s));
  const el = {
    play: $('#btnPlay'), reset: $('#btnReset'),
    fill: $('#fill'), thumb: $('#thumb'), track: $('#track'), scrubArea: $('#scrubArea'),
    pLabel: $('#pLabel'), thetaLabel: $('#thetaLabel'),
    phaseTag: $('#phaseTag'), phaseDesc: $('#phaseDesc'), readout: $('#readout'),
    grpBasic: $('#grpBasic'), grpKoenig: $('#grpKoenig'), eEmpty: $('#eEmpty'),
  };

  function colRefs(id) {
    const c = document.getElementById(id);
    return { c, val: c.querySelector('.bar-val'), fills: c.querySelectorAll('.seg-fill') };
  }
  const bars = {
    K1: colRefs('colK1'), K2: colRefs('colK2'), Kt: colRefs('colKt'), Us: colRefs('colUs'), E: colRefs('colE'),
    Kc: colRefs('colKc'), Kint: colRefs('colKint'), St: colRefs('colStack'),
  };
  const ktBracket = $('#ktBracket');

  const H = (v) => clamp(v, 0, 100).toFixed(3) + '%';
  let lastHTML = { tag: '', desc: '', ro: '' };

  function keyAt(p) {
    for (const k of [P0, PM, P1]) if (Math.abs(p - k) < 0.06) return k;
    return null;
  }

  const SYM = {
    v1: 'v<sub>1</sub>', v2: 'v<sub>2</sub>', vc: 'v<sub>c</sub>',
    v1p: 'v′<sub>1</sub>', v2p: 'v′<sub>2</sub>',
    Kt: 'K<sub>t</sub>', Kc: 'K<sub>c</sub>', Ki: 'K<sub>int</sub>', Us: 'U<sub>s</sub>',
  };

  function phaseContent(S) {
    const k = keyAt(S.p), L2 = state.layers[2] || state.mode !== 'A';
    if (k === P0) return {
      key: true, tag: '⏮ 碰撞起點（剛接觸 15%）',
      desc: `兩球仍為<b>正圓</b>、表面剛接觸：形變 <b>x = 0</b>、<span class="F">F = 0</span>、<span class="hl">${SYM.Us} = 0%</span><br>` +
        `系統內動能 <span class="hl">${SYM.Ki} 最大 = ${fmtP(S.Kint)}</span>，系統總動能 ${SYM.Kt} = ${fmtP(S.Kt)}`,
    };
    if (k === PM) return {
      key: true, tag: '🎯 極限狀態（最大形變／最接近 50%）',
      desc: `麻糬壓最扁、C<sub>1</sub>C<sub>2</sub> 最接近（d = ${S.d.toFixed(2)}R）；<span class="v">${SYM.v1} = ${SYM.v2} = ${SYM.vc} = ${fmtV(model.vc)}</span>，<span class="vp">${SYM.v1p} = ${SYM.v2p} = 0</span>（相對質心靜止）<br>` +
        `<span class="hl">${SYM.Ki} = 0%</span>、<span class="hl">${SYM.Us} 達極大 = ${fmtP(S.Us)}</span>、${SYM.Kt} 降至最低 = ${SYM.Kc} = ${fmtP(S.Kc)}`,
    };
    if (k === P1) return {
      key: true, tag: '⏭ 碰撞終點（剛分離 85%）',
      desc: `恢復<b>正圓</b>：x = 0、<span class="F">F = 0</span>；末速 <span class="v">${SYM.v1} = ${fmtV(model.v1f)}、${SYM.v2} = ${fmtV(model.v2f)}</span> m/s<br>` +
        `${SYM.Us} 全數釋放 → <span class="hl">${SYM.Ki} = ${fmtP(S.Kint)}</span>、${SYM.Kt} = ${fmtP(S.Kt)}（回到碰撞前）`,
    };
    if (S.p < P0) return {
      tag: '① 接觸前：等速接近',
      desc: `兩球維持正圓、尚未接觸：x = 0、無交互作用力，各自以初速 <span class="v">${SYM.v1} = ${fmtV(model.v1i)}、${SYM.v2} = ${fmtV(model.v2i)}</span> m/s 等速接近。`,
    };
    if (S.p < PM) return {
      tag: '② 壓縮期：擠壓形變',
      desc: `球心距 d 縮小、形變 x 增大；<span class="F">F</span> 與 <span class="v">v</span> 反向者減速、同向者加速，${SYM.Kt} 轉為 ${SYM.Us}。` +
        (L2 ? `<br>質心系：<span class="F">F</span> 與 <span class="vp">v′</span> <b>反向</b> → ${SYM.v1p}、${SYM.v2p} 同步縮短，${SYM.Ki} → ${SYM.Us}。` : ''),
    };
    if (S.p <= P1) return {
      tag: '③ 恢復期：恢復形變',
      desc: `形變恢復、彈力持續推開兩球，${SYM.Us} 釋放回動能 ${SYM.Kt}。` +
        (L2 ? `<br>質心系：<span class="F">F</span> 與 <span class="vp">v′</span> <b>同向</b> → ${SYM.v1p}、${SYM.v2p} 反向變長，${SYM.Us} → ${SYM.Ki}。` : ''),
    };
    return {
      tag: '④ 分離後：等速遠離',
      desc: `兩球恢復正圓，以末速 <span class="v">${SYM.v1} = ${fmtV(model.v1f)}、${SYM.v2} = ${fmtV(model.v2f)}</span> m/s 等速遠離（彈性碰撞：動能完全恢復）。`,
    };
  }

  function readoutHTML(S) {
    const L1 = state.layers[1], L2 = state.layers[2] || state.mode !== 'A';
    const ro = (c, sym, val, unit = '') => `<div class="ro"><i style="color:${c}">${sym}</i> <b>${val}</b><span class="unit">${unit}</span></div>`;
    let h = '';
    h += ro(COL.v, SYM.v1, fmtV(S.v1)) + ro(COL.v, SYM.v2, fmtV(S.v2));
    if (L2) {
      h += ro(COL.vc, SYM.vc, fmtV(model.vc), ' 🔒') + ro('#cbd5e1', '', '', '<span style="opacity:0">.</span>');
      h += ro(COL.vp, SYM.v1p, fmtV(S.v1p)) + ro(COL.vp, SYM.v2p, fmtV(S.v2p));
    }
    h += ro('#e2e8f0', 'd', S.d.toFixed(2), ' R') + ro(COL.F, 'x', (S.s * 100).toFixed(0) + '%', ' x<sub>max</sub>');
    h += `<div class="ro unit">速度單位</div><div class="ro unit">m/s</div>`;
    void L1;
    return h;
  }

  function updateUI(S) {
    // 時間軸
    el.fill.style.width = S.p + '%';
    el.thumb.style.left = S.p + '%';
    el.pLabel.textContent = `進度 ${S.p.toFixed(1)}%`;
    el.thetaLabel.innerHTML = S.theta != null
      ? `碰撞相位 θ = ${(S.theta * 180 / Math.PI).toFixed(0)}°　x = x<sub>max</sub> sin θ = ${(S.s * 100).toFixed(0)}% x<sub>max</sub>`
      : (S.p < P0 ? '尚未接觸' : '已分離');

    const k = keyAt(S.p);
    $$('.btn.key').forEach((b) => b.classList.toggle('active', Number(b.dataset.key) === k));

    // 階段說明
    const pc = phaseContent(S);
    if (pc.tag !== lastHTML.tag) { el.phaseTag.textContent = pc.tag; lastHTML.tag = pc.tag; }
    el.phaseTag.classList.toggle('key', !!pc.key);
    if (pc.desc !== lastHTML.desc) { el.phaseDesc.innerHTML = pc.desc; lastHTML.desc = pc.desc; }
    const ro = readoutHTML(S);
    if (ro !== lastHTML.ro) { el.readout.innerHTML = ro; lastHTML.ro = ro; }

    // 能量儀表板
    const L1 = state.layers[1], L3 = state.layers[3];
    el.grpBasic.classList.toggle('hidden', !L1);
    el.grpKoenig.classList.toggle('hidden', !L3);
    el.eEmpty.classList.toggle('show', !L1 && !L3);

    const setOne = (b, v) => { b.fills[0].style.height = H(v); b.val.textContent = fmtP(v); };
    setOne(bars.K1, S.K1);
    setOne(bars.K2, S.K2);
    setOne(bars.Us, S.Us);
    setOne(bars.E, S.E);
    setOne(bars.Kc, S.Kc);
    setOne(bars.Kint, S.Kint);
    // K_t：開啟階段 3 時以 K_c（下）+ K_int（上）雙色顯示
    bars.Kt.fills[0].style.height = H(L3 ? S.Kc : S.Kt);
    bars.Kt.fills[1].style.height = H(L3 ? S.Kint : 0);
    bars.Kt.fills[0].classList.toggle('split', L3);
    bars.Kt.fills[1].classList.toggle('split', L3);
    bars.Kt.c.classList.toggle('split', L3);
    bars.Kt.val.textContent = fmtP(S.Kt);
    // 堆疊條
    bars.St.fills[0].style.height = H(S.Kc);
    bars.St.fills[1].style.height = H(S.Kint);
    bars.St.fills[2].style.height = H(S.Us);
    bars.St.fills[0].querySelector('span').style.opacity = S.Kc > 12 ? 1 : 0;
    bars.St.fills[1].querySelector('span').style.opacity = S.Kint > 12 ? 1 : 0;
    bars.St.fills[2].querySelector('span').style.opacity = S.Us > 12 ? 1 : 0;
    bars.St.val.textContent = fmtP(S.Kc + S.Kint + S.Us);
    ktBracket.style.height = H(S.Kt);

    // 關鍵瞬間高亮對應能量柱
    const hl = { K1: false, K2: false, Kt: false, Us: false, E: false, Kc: false, Kint: false, St: false };
    if (k === P0) { hl.Us = true; hl.Kint = true; }
    if (k === PM) { hl.Us = true; hl.Kint = true; hl.Kt = true; hl.Kc = true; }
    if (k === P1) { hl.Us = true; hl.Kint = true; hl.Kt = true; }
    for (const key in bars) bars[key].c.classList.toggle('hl', hl[key]);
  }

  function syncControls() {
    $$('[data-scenario]').forEach((b) => b.classList.toggle('active', b.dataset.scenario === state.scenario));
    $$('[data-mass]').forEach((b) => b.classList.toggle('active', b.dataset.mass === state.mass));
    $$('[data-mode]').forEach((b) => b.classList.toggle('active', b.dataset.mode === state.mode));
    $$('[data-speed]').forEach((b) => b.classList.toggle('active', Number(b.dataset.speed) === state.speed));
    $$('[data-layer]').forEach((b) => b.classList.toggle('on', !!state.layers[b.dataset.layer]));
    el.play.textContent = state.playing ? '⏸ 暫停' : '▶ 播放';
    el.play.classList.toggle('playing', state.playing);
    const cols = BALL_COLORS[state.mass];
    document.documentElement.style.setProperty('--b1', hex(cols[0]));
    document.documentElement.style.setProperty('--b2', hex(cols[1]));
  }

  function setPlaying(v) {
    state.playing = v;
    syncControls();
  }

  // =====================================================
  // 十、事件綁定
  // =====================================================
  $$('[data-scenario]').forEach((b) => b.addEventListener('click', () => {
    state.scenario = b.dataset.scenario; model = buildModel(); syncControls();
  }));
  $$('[data-mass]').forEach((b) => b.addEventListener('click', () => {
    state.mass = b.dataset.mass; model = buildModel(); syncControls();
  }));
  $$('[data-mode]').forEach((b) => b.addEventListener('click', () => {
    state.mode = b.dataset.mode; syncControls();         // 保留目前時間軸進度
  }));
  $$('[data-layer]').forEach((b) => b.addEventListener('click', () => {
    const k = b.dataset.layer; state.layers[k] = !state.layers[k]; syncControls();
  }));
  $$('[data-speed]').forEach((b) => b.addEventListener('click', () => {
    state.speed = Number(b.dataset.speed); syncControls();
  }));
  $$('.btn.key').forEach((b) => b.addEventListener('click', () => {
    state.p = Number(b.dataset.key); setPlaying(false);
  }));
  el.play.addEventListener('click', () => {
    if (!state.playing && state.p >= 100) state.p = 0;
    setPlaying(!state.playing);
  });
  el.reset.addEventListener('click', () => { state.p = 0; setPlaying(false); });

  // 時間軸拖曳（支援觸控／滑鼠／觸控筆）
  let dragging = false;
  function seekFromEvent(e) {
    const r = el.track.getBoundingClientRect();
    let p = clamp((e.clientX - r.left) / r.width * 100, 0, 100);
    for (const k of [P0, PM, P1]) if (Math.abs(p - k) < 0.9) p = k;   // 關鍵節點吸附
    state.p = p;
  }
  el.scrubArea.addEventListener('pointerdown', (e) => {
    dragging = true;
    el.scrubArea.setPointerCapture(e.pointerId);
    if (state.playing) setPlaying(false);
    seekFromEvent(e);
    e.preventDefault();
  });
  el.scrubArea.addEventListener('pointermove', (e) => { if (dragging) seekFromEvent(e); });
  const endDrag = () => { dragging = false; };
  el.scrubArea.addEventListener('pointerup', endDrag);
  el.scrubArea.addEventListener('pointercancel', endDrag);

  // 鍵盤快捷鍵（空白鍵播放、←→ 逐格、1/2/3 關鍵瞬間）
  window.addEventListener('keydown', (e) => {
    if (e.code === 'Space') { e.preventDefault(); el.play.click(); }
    else if (e.code === 'ArrowRight') { setPlaying(false); state.p = clamp(state.p + (e.shiftKey ? 5 : 0.5), 0, 100); }
    else if (e.code === 'ArrowLeft') { setPlaying(false); state.p = clamp(state.p - (e.shiftKey ? 5 : 0.5), 0, 100); }
    else if (e.key === '1') { state.p = P0; setPlaying(false); }
    else if (e.key === '2') { state.p = PM; setPlaying(false); }
    else if (e.key === '3') { state.p = P1; setPlaying(false); }
    else if (e.key === 'r' || e.key === 'R') { state.p = 0; setPlaying(false); }
  });

  // =====================================================
  // 十一、主迴圈
  // =====================================================
  let last = performance.now();
  function tick(now) {
    const dt = Math.min(0.05, (now - last) / 1000);
    last = now;
    if (state.playing && !dragging) {
      state.p += dt * 100 / BASE_DURATION * state.speed;
      if (state.p >= 100) { state.p = 100; setPlaying(false); }
    }
    const S = stateAt(model, state.p);
    draw(S);
    updateUI(S);
    requestAnimationFrame(tick);
  }

  // 初始化
  model = buildModel();
  resizeCanvas();
  new ResizeObserver(resizeCanvas).observe(wrap);
  window.addEventListener('resize', resizeCanvas);
  syncControls();

  // 供自動測試使用：?p=50&mode=AB&layers=123&scenario=chase&mass=unequal
  try {
    const q = new URLSearchParams(location.search);
    if (q.has('scenario') && SCENARIOS[q.get('scenario')]) state.scenario = q.get('scenario');
    if (q.has('mass') && MASSES[q.get('mass')]) state.mass = q.get('mass');
    if (q.has('mode')) state.mode = q.get('mode');
    if (q.has('layers')) { const s = q.get('layers'); state.layers = { 1: s.includes('1'), 2: s.includes('2'), 3: s.includes('3') }; }
    if (q.has('p')) state.p = clamp(Number(q.get('p')) || 0, 0, 100);
    model = buildModel(); syncControls();
  } catch (_) { /* 忽略 */ }

  // 能量守恆自我檢查（主控台輸出）
  (function selfCheck() {
    let maxErr = 0;
    for (const sc in SCENARIOS) for (const ms in MASSES) {
      const save = [state.scenario, state.mass];
      state.scenario = sc; state.mass = ms;
      const m = buildModel();
      for (let p = 0; p <= 100; p += 0.25) {
        const S = stateAt(m, p);
        maxErr = Math.max(maxErr, Math.abs(S.E - 100), Math.abs(S.Kc + S.Kint - S.Kt));
      }
      [state.scenario, state.mass] = save;
    }
    console.log('[能量守恆自檢] 最大誤差 =', maxErr.toExponential(2), '%');
    model = buildModel();
  })();

  requestAnimationFrame(tick);
})();
