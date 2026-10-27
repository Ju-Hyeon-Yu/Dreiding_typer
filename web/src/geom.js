// SMILES만 있을 때 RDKit.js 2D 좌표에서 3D 구조를 만드는 간단한 최소화 (DREIDING 결합·각 + 반발)
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory(require("./typer.js"));
  else root.DreidingGeom = factory(root.DreidingTyper);
})(typeof self !== "undefined" ? self : this, function (T) {
  "use strict";

  const RVDW = { H: 1.2, C: 1.7, N: 1.55, O: 1.52, F: 1.47, P: 1.8, S: 1.8, Cl: 1.75, Br: 1.85, I: 1.98,
    Si: 2.1, B: 1.92, Se: 1.9, Na: 2.27, K: 2.75, Li: 1.82, Ca: 2.31, Mg: 1.73 };
  const RCOV = { H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57, P: 1.07, S: 1.05, Cl: 1.02, Br: 1.2, I: 1.39,
    Si: 1.11, B: 0.84, Se: 1.2 };
  const DEG = Math.PI / 180;

  function mulberry32(seed) {
    let a = seed >>> 0;
    return () => {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /** RDKit.js 분자(명시적 H 포함)의 2D 좌표 [[x,y], ...]. */
  function coords2D(mol) {
    mol.set_new_coords();
    const lines = mol.get_molblock().split(/\r?\n/);
    const xy = [];
    if (/V3000/.test(lines[3] || "")) {
      let inAtoms = false;
      for (const l of lines) {
        if (/M {2}V30 BEGIN ATOM/.test(l)) { inAtoms = true; continue; }
        if (/M {2}V30 END ATOM/.test(l)) break;
        if (inAtoms) { const t = l.trim().split(/\s+/); xy.push([parseFloat(t[4]), parseFloat(t[5])]); }
      }
    } else {
      const n = parseInt(lines[3].slice(0, 3), 10);
      for (let k = 0; k < n; k++) {
        const l = lines[4 + k];
        xy.push([parseFloat(l.slice(0, 10)), parseFloat(l.slice(10, 20))]);
      }
    }
    return xy;
  }

  /** 중심 원자의 목표 결합각(도). 4배위·3배위 이상 S 는 사면체로 둔다(좌표 생성 전용). */
  function targetAngle(label, deg, ff) {
    if (deg >= 4) return 109.47;
    if (label.startsWith("S_3") && deg >= 3) return 109.47;
    if (label.startsWith("P_3") && deg >= 3) return 109.47;
    const t = ff.theta0[label];
    return t === undefined ? 109.47 : t;
  }

  /** 한 조각을 최소화. idx: 전역 원자 번호 목록, x: Float64Array(3n) 시작 좌표. */
  function optimizeFragment(g, labels, ff, idx, x, maxIter = 6000) {
    const n = idx.length, loc = new Map(idx.map((a, k) => [a, k]));
    const sym = idx.map(a => g.atoms[a].sym);
    const bonds = [], angles = [], planes = [], pairs = [];
    const adj = idx.map(() => []);
    for (const b of g.bonds) {
      if (!loc.has(b.a)) continue;
      const i = loc.get(b.a), j = loc.get(b.b);
      adj[i].push(j); adj[j].push(i);
      let r0 = ff.bondR0(labels[b.a], labels[b.b]);
      if (r0 === null) r0 = (RCOV[sym[i]] || 1.0) + (RCOV[sym[j]] || 1.0);
      bonds.push([i, j, r0]);
    }
    for (let c = 0; c < n; c++) {
      const nb = adj[c], lab = labels[idx[c]];
      const cos0 = Math.cos(targetAngle(lab, nb.length, ff) * DEG);
      for (let p = 0; p < nb.length; p++) for (let q = p + 1; q < nb.length; q++) angles.push([nb[p], c, nb[q], cos0]);
      if (nb.length === 3 && T.sp2Class(lab)) planes.push([c, nb[0], nb[1], nb[2]]);
    }
    // 위상 거리 (1-2, 1-3 은 반발 제외, 1-4 는 약하게)
    const topo = idx.map((_, s) => {
      const d = new Array(n).fill(99); d[s] = 0; const q = [s];
      while (q.length) { const u = q.shift(); if (d[u] >= 3) continue; for (const v of adj[u]) if (d[v] === 99) { d[v] = d[u] + 1; q.push(v); } }
      return d;
    });
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (topo[i][j] <= 2) continue;
        const s = (RVDW[sym[i]] || 1.8) + (RVDW[sym[j]] || 1.8);
        pairs.push([i, j, topo[i][j] === 3 ? 0.62 * s : 0.8 * s]);
      }
    }

    const KB = 300, KA = 100, KP = 100, KR = 30;
    const f = new Float64Array(3 * n);
    const force = () => {
      f.fill(0);
      let fmax = 0;
      for (const [i, j, r0] of bonds) {
        const dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1], dz = x[3 * i + 2] - x[3 * j + 2];
        const r = Math.hypot(dx, dy, dz) || 1e-6, s = -2 * KB * (r - r0) / r;
        f[3 * i] += s * dx; f[3 * i + 1] += s * dy; f[3 * i + 2] += s * dz;
        f[3 * j] -= s * dx; f[3 * j + 1] -= s * dy; f[3 * j + 2] -= s * dz;
      }
      for (const [a, c, b, cos0] of angles) {
        const u = [x[3 * a] - x[3 * c], x[3 * a + 1] - x[3 * c + 1], x[3 * a + 2] - x[3 * c + 2]];
        const v = [x[3 * b] - x[3 * c], x[3 * b + 1] - x[3 * c + 1], x[3 * b + 2] - x[3 * c + 2]];
        const ru = Math.hypot(...u) || 1e-6, rv = Math.hypot(...v) || 1e-6;
        const cs = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (ru * rv);
        const pre = -2 * KA * (cs - cos0);
        for (let k = 0; k < 3; k++) {
          const ga = v[k] / (ru * rv) - cs * u[k] / (ru * ru);   // d cos / d a
          const gb = u[k] / (ru * rv) - cs * v[k] / (rv * rv);   // d cos / d b
          f[3 * a + k] += pre * ga; f[3 * b + k] += pre * gb; f[3 * c + k] -= pre * (ga + gb);
        }
      }
      for (const [c, p, q, r] of planes) {
        const P = [p, q, r].map(k => [x[3 * k], x[3 * k + 1], x[3 * k + 2]]);
        const e1 = [0, 1, 2].map(k => P[1][k] - P[0][k]), e2 = [0, 1, 2].map(k => P[2][k] - P[0][k]);
        const nrm = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
        const ln = Math.hypot(...nrm);
        if (ln < 1e-8) continue;
        const nh = nrm.map(v => v / ln);
        const h = [0, 1, 2].reduce((s, k) => s + nh[k] * (x[3 * c + k] - (P[0][k] + P[1][k] + P[2][k]) / 3), 0);
        for (let k = 0; k < 3; k++) {
          const g0 = -2 * KP * h * nh[k];
          f[3 * c + k] += g0;
          for (const m of [p, q, r]) f[3 * m + k] -= g0 / 3;
        }
      }
      for (const [i, j, d0] of pairs) {
        const dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1], dz = x[3 * i + 2] - x[3 * j + 2];
        const r = Math.hypot(dx, dy, dz);
        if (r >= d0) continue;
        const rr = r || 1e-6, s = 2 * KR * (d0 - rr) / rr;
        f[3 * i] += s * dx; f[3 * i + 1] += s * dy; f[3 * i + 2] += s * dz;
        f[3 * j] -= s * dx; f[3 * j + 1] -= s * dy; f[3 * j + 2] -= s * dz;
      }
      for (let k = 0; k < n; k++) fmax = Math.max(fmax, Math.hypot(f[3 * k], f[3 * k + 1], f[3 * k + 2]));
      return fmax;
    };

    // FIRE 최소화
    const v = new Float64Array(3 * n);
    let dt = 0.002, alpha = 0.1, npos = 0, it = 0;
    const DTMAX = 0.01, MAXSTEP = 0.1;
    for (; it < maxIter; it++) {
      const fmax = force();
      if (fmax < 0.05) break;
      let P = 0, vn = 0, fn = 0;
      for (let k = 0; k < 3 * n; k++) { P += f[k] * v[k]; vn += v[k] * v[k]; fn += f[k] * f[k]; }
      vn = Math.sqrt(vn); fn = Math.sqrt(fn) || 1e-12;
      if (P > 0) {
        for (let k = 0; k < 3 * n; k++) v[k] = (1 - alpha) * v[k] + alpha * f[k] / fn * vn;
        if (++npos > 5) { dt = Math.min(dt * 1.1, DTMAX); alpha *= 0.99; }
      } else {
        v.fill(0); dt *= 0.5; alpha = 0.1; npos = 0;
      }
      for (let k = 0; k < n; k++) {
        let sx = v[3 * k] * dt + f[3 * k] * dt * dt, sy = v[3 * k + 1] * dt + f[3 * k + 1] * dt * dt,
          sz = v[3 * k + 2] * dt + f[3 * k + 2] * dt * dt;
        const s = Math.hypot(sx, sy, sz);
        if (s > MAXSTEP) { sx *= MAXSTEP / s; sy *= MAXSTEP / s; sz *= MAXSTEP / s; }
        x[3 * k] += sx; x[3 * k + 1] += sy; x[3 * k + 2] += sz;
        v[3 * k] += f[3 * k] * dt; v[3 * k + 1] += f[3 * k + 1] * dt; v[3 * k + 2] += f[3 * k + 2] * dt;
      }
    }
    return it;
  }

  /** 조각마다 따로 최소화하고 x 축으로 gap Å 씩 띄워 놓는다 (Python embed 와 같은 배치). */
  function embed3D(g, labels, ff, xy, seed = 1, gap = 4.0) {
    const rnd = mulberry32(seed);
    const xyz = new Array(g.atoms.length);
    let xEnd = null;
    for (const fr of T.fragments(g)) {
      const x = new Float64Array(3 * fr.length);
      fr.forEach((a, k) => { x[3 * k] = xy[a][0]; x[3 * k + 1] = xy[a][1]; x[3 * k + 2] = (rnd() - 0.5) * 0.8; });
      if (fr.length > 1) optimizeFragment(g, labels, ff, fr, x);
      let cx = 0, cy = 0, cz = 0;
      fr.forEach((_, k) => { cx += x[3 * k]; cy += x[3 * k + 1]; cz += x[3 * k + 2]; });
      cx /= fr.length; cy /= fr.length; cz /= fr.length;
      const pts = fr.map((_, k) => [x[3 * k] - cx, x[3 * k + 1] - cy, x[3 * k + 2] - cz]);
      const minX = Math.min(...pts.map(p => p[0]));
      const shift = xEnd === null ? 0 : xEnd + gap - minX;
      fr.forEach((a, k) => { xyz[a] = [pts[k][0] + shift, pts[k][1], pts[k][2]]; });
      xEnd = Math.max(...pts.map(p => p[0] + shift));
    }
    return xyz;
  }

  return { coords2D, embed3D, targetAngle };
});
