// 캔버스 2D 로 그리는 볼-스틱 분자 뷰어 — 회전·확대·이동, 원자 타입 라벨, 강조, 마우스 정보
(function (root) {
  "use strict";

  const ELEMENT_COLOR = { H: "#E8EEF2", C: "#8C99A6", N: "#4C7BF0", O: "#E8513F", F: "#7FD15A", Cl: "#3DBA5E",
    Br: "#B5562F", I: "#9B4FD1", S: "#E8C23A", P: "#F08A2E", Si: "#C9A27E", B: "#F2A8A0", Na: "#9C6BE0",
    Li: "#C080FF", K: "#8F40D4", Ca: "#4FD07A", Se: "#F0A020", Fe: "#E06633", Zn: "#7D80B0" };
  const ATOM_R = { H: 0.22 };
  // 타입 색 (어두운 무대에서 서로 구분되는 순서)
  const TYPE_PALETTE = ["#5BC0EB", "#F4A259", "#9BC53D", "#E55934", "#C792EA", "#FDE74C", "#4ECDC4", "#FF8FAB",
    "#A3A3FF", "#B8E0D2", "#F7B267", "#7AE582", "#FF6B6B", "#6C9BCF", "#D4A5A5", "#E2C044"];

  function shade(hex, f) {
    const n = parseInt(hex.slice(1), 16);
    const c = [(n >> 16) & 255, (n >> 8) & 255, n & 255].map(v => Math.max(0, Math.min(255, Math.round(v * f))));
    return `rgb(${c[0]},${c[1]},${c[2]})`;
  }

  function principalFrame(pts) {
    const n = pts.length, c = [0, 0, 0];
    pts.forEach(p => { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; });
    c.forEach((_, k) => { c[k] /= n || 1; });
    const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    pts.forEach(p => {
      const d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) A[i][j] += d[i] * d[j];
    });
    // 야코비 고유값 분해
    const V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 30; sweep++) {
      for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
        if (Math.abs(A[p][q]) < 1e-12) continue;
        const th = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
        const cs = Math.cos(th), sn = Math.sin(th);
        for (let k = 0; k < 3; k++) {
          const akp = A[k][p], akq = A[k][q];
          A[k][p] = cs * akp - sn * akq; A[k][q] = sn * akp + cs * akq;
        }
        for (let k = 0; k < 3; k++) {
          const apk = A[p][k], aqk = A[q][k];
          A[p][k] = cs * apk - sn * aqk; A[q][k] = sn * apk + cs * aqk;
        }
        for (let k = 0; k < 3; k++) {
          const vkp = V[k][p], vkq = V[k][q];
          V[k][p] = cs * vkp - sn * vkq; V[k][q] = sn * vkp + cs * vkq;
        }
      }
    }
    const order = [0, 1, 2].sort((a, b) => A[b][b] - A[a][a]);
    const rows = order.map(k => [V[0][k], V[1][k], V[2][k]]);
    const z = [rows[0][1] * rows[1][2] - rows[0][2] * rows[1][1], rows[0][2] * rows[1][0] - rows[0][0] * rows[1][2],
      rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]];
    return { center: c, rot: [rows[0], rows[1], z] };
  }

  const mul = (A, B) => A.map(r => [0, 1, 2].map(j => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
  const rotX = a => [[1, 0, 0], [0, Math.cos(a), -Math.sin(a)], [0, Math.sin(a), Math.cos(a)]];
  const rotY = a => [[Math.cos(a), 0, Math.sin(a)], [0, 1, 0], [-Math.sin(a), 0, Math.cos(a)]];

  class MolViewer {
    constructor(canvas, tooltip) {
      this.canvas = canvas;
      this.ctx = canvas.getContext("2d");
      this.tip = tooltip;
      this.mol = null;
      this.labelMode = "main";      // main | all | none
      this.colorMode = "type";      // type | element
      this.selected = new Set();    // 강조할 타입
      this.hover = -1;
      this.typeColor = new Map();
      this._bind();
      new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
      this.resize();
    }

    /** mol: {atoms:[{sym, label, name, q, xyz, frag}], bonds:[[i,j]], flagged:Set} */
    setMolecule(mol, typeOrder) {
      this.mol = mol;
      this.nbrs = mol.atoms.map(() => []);
      mol.bonds.forEach(([i, j]) => { this.nbrs[i].push(j); this.nbrs[j].push(i); });
      this.typeColor = new Map(typeOrder.map((t, k) => [t, TYPE_PALETTE[k % TYPE_PALETTE.length]]));
      this.selected.clear();
      this.resetView();
    }

    colorOf(i) {
      const a = this.mol.atoms[i];
      if (this.colorMode === "type") return this.typeColor.get(a.label) || "#999999";
      return ELEMENT_COLOR[a.sym] || "#D07FD0";
    }

    resetView() {
      if (!this.mol) { this.draw(); return; }
      const fr = principalFrame(this.mol.atoms.map(a => a.xyz));
      this.center = fr.center;
      this.rot = fr.rot;
      // 주축 기준 x·y 범위로 화면에 꽉 차게 맞춤
      const lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
      this.mol.atoms.forEach(a => {
        const d = [0, 1, 2].map(k => a.xyz[k] - this.center[k]);
        for (let r = 0; r < 2; r++) {
          const v = fr.rot[r][0] * d[0] + fr.rot[r][1] * d[1] + fr.rot[r][2] * d[2];
          lo[r] = Math.min(lo[r], v); hi[r] = Math.max(hi[r], v);
        }
      });
      this.span = [hi[0] - lo[0] + 2.4, hi[1] - lo[1] + 2.4];
      this.mid = [(hi[0] + lo[0]) / 2, (hi[1] + lo[1]) / 2];
      this.zoom = 1;
      this.pan = [0, 0];
      this.draw();
    }

    resize() {
      const box = this.canvas.parentElement.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      this.w = Math.max(200, box.width);
      this.h = Math.max(200, box.height);
      this.canvas.width = Math.round(this.w * dpr);
      this.canvas.height = Math.round(this.h * dpr);
      this.canvas.style.width = this.w + "px";
      this.canvas.style.height = this.h + "px";
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.draw();
    }

    _project() {
      // 위쪽 도구 막대(약 50px)와 아래 안내(약 30px)를 피해서 맞춘다
      const fit = Math.min((this.w - 48) / this.span[0], (this.h - 110) / this.span[1], 110);
      const s = Math.max(4, fit) * this.zoom, R = this.rot, c = this.center;
      this.scale = s;
      const ox = -s * this.mid[0], oy = s * this.mid[1] + 10;
      this.proj = this.mol.atoms.map(a => {
        const d = [a.xyz[0] - c[0], a.xyz[1] - c[1], a.xyz[2] - c[2]];
        const x = R[0][0] * d[0] + R[0][1] * d[1] + R[0][2] * d[2];
        const y = R[1][0] * d[0] + R[1][1] * d[1] + R[1][2] * d[2];
        const z = R[2][0] * d[0] + R[2][1] * d[1] + R[2][2] * d[2];
        return [this.w / 2 + ox + this.pan[0] + s * x, this.h / 2 + oy + this.pan[1] - s * y, z];
      });
    }

    draw() {
      const ctx = this.ctx;
      ctx.clearRect(0, 0, this.w, this.h);
      if (!this.mol) return;
      this._project();
      const P = this.proj, s = this.scale, atoms = this.mol.atoms;
      const zs = P.map(p => p[2]), zmin = Math.min(...zs), zmax = Math.max(...zs), zr = (zmax - zmin) || 1;
      const order = atoms.map((_, i) => i).sort((a, b) => P[a][2] - P[b][2]);
      const dim = i => this.selected.size > 0 && !this.selected.has(atoms[i].label);
      const radius = i => (ATOM_R[atoms[i].sym] || 0.34) * s;

      for (const i of order) {
        const depth = 0.62 + 0.38 * (P[i][2] - zmin) / zr;
        const col = this.colorOf(i);
        ctx.globalAlpha = dim(i) ? 0.16 : 1;
        // 반쪽 결합 (원자 색)
        ctx.lineCap = "round";
        for (const j of this.nbrs[i]) {
          const mx = (P[i][0] + P[j][0]) / 2, my = (P[i][1] + P[j][1]) / 2;
          ctx.strokeStyle = shade(col.startsWith("#") ? col : "#999999", depth * 0.82);
          ctx.lineWidth = Math.max(2, 0.16 * s);
          ctx.beginPath(); ctx.moveTo(P[i][0], P[i][1]); ctx.lineTo(mx, my); ctx.stroke();
        }
        const r = radius(i);
        const g = ctx.createRadialGradient(P[i][0] - r * 0.35, P[i][1] - r * 0.35, r * 0.1, P[i][0], P[i][1], r);
        g.addColorStop(0, shade(col, Math.min(1.35, depth * 1.35)));
        g.addColorStop(1, shade(col, depth * 0.7));
        ctx.fillStyle = g;
        ctx.beginPath(); ctx.arc(P[i][0], P[i][1], r, 0, 2 * Math.PI); ctx.fill();
        if (this.mol.flagged.has(i)) {
          ctx.strokeStyle = "#FFB020"; ctx.lineWidth = 2.5;
          ctx.beginPath(); ctx.arc(P[i][0], P[i][1], r + 4, 0, 2 * Math.PI); ctx.stroke();
        }
        if (i === this.hover) {
          ctx.strokeStyle = "#FFFFFF"; ctx.lineWidth = 2;
          ctx.beginPath(); ctx.arc(P[i][0], P[i][1], r + 2, 0, 2 * Math.PI); ctx.stroke();
        }
      }
      ctx.globalAlpha = 1;

      // 라벨 (맨 위에). 앞쪽 원자부터 자리를 잡고, 겹치면 다른 방향으로 옮긴다
      if (this.labelMode !== "none") {
        const fs = Math.max(10, Math.min(13, 0.32 * s));
        ctx.font = `600 ${fs}px "Cascadia Mono", Consolas, "SFMono-Regular", Menlo, monospace`;
        ctx.textBaseline = "middle";
        const placed = [], boxes = [];
        const hit = b => boxes.some(o => b[0] < o[0] + o[2] && o[0] < b[0] + b[2] && b[1] < o[1] + o[3] && o[1] < b[1] + b[3]);
        for (const i of order.slice().reverse()) {
          const a = atoms[i];
          if (dim(i)) continue;
          if (this.labelMode === "main" && a.sym === "H" && a.label !== "H_HB") continue;
          const r = radius(i), tw = ctx.measureText(a.label).width + 6, th = fs * 1.24;
          const cx = P[i][0], cy = P[i][1], o = r * 0.7 + 3;
          const spots = [[cx + o, cy - o - th], [cx + o, cy + o], [cx - o - tw, cy - o - th], [cx - o - tw, cy + o],
            [cx + r + 3, cy - th / 2], [cx - r - 3 - tw, cy - th / 2], [cx - tw / 2, cy - r - 3 - th], [cx - tw / 2, cy + r + 3]];
          let b = spots.map(([x, y]) => [x, y, tw, th]).find(bb => !hit(bb)) || [spots[0][0], spots[0][1], tw, th];
          boxes.push(b);
          placed.push([i, b]);
        }
        for (const [i, b] of placed.reverse()) {
          const a = atoms[i];
          ctx.fillStyle = "rgba(12,18,25,0.8)";
          ctx.beginPath();
          if (ctx.roundRect) ctx.roundRect(b[0], b[1], b[2], b[3], 3); else ctx.rect(b[0], b[1], b[2], b[3]);
          ctx.fill();
          ctx.fillStyle = this.colorMode === "type" ? (this.typeColor.get(a.label) || "#FFFFFF") : "#F2F6FA";
          ctx.fillText(a.label, b[0] + 3, b[1] + b[3] / 2);
        }
      }
    }

    _pick(mx, my) {
      if (!this.mol || !this.proj) return -1;
      let best = -1, bz = -Infinity;
      this.proj.forEach((p, i) => {
        const r = (ATOM_R[this.mol.atoms[i].sym] || 0.34) * this.scale + 3;
        if ((p[0] - mx) ** 2 + (p[1] - my) ** 2 <= r * r && p[2] > bz) { best = i; bz = p[2]; }
      });
      return best;
    }

    _bind() {
      const cv = this.canvas, pts = new Map();
      let last = null, pinch = null;
      cv.addEventListener("contextmenu", e => e.preventDefault());
      cv.addEventListener("pointerdown", e => {
        cv.setPointerCapture(e.pointerId);
        pts.set(e.pointerId, [e.offsetX, e.offsetY]);
        last = { x: e.offsetX, y: e.offsetY, pan: e.shiftKey || e.button === 1 || e.button === 2 };
        if (pts.size === 2) { const [a, b] = [...pts.values()]; pinch = { d: Math.hypot(a[0] - b[0], a[1] - b[1]), z: this.zoom }; }
      });
      cv.addEventListener("pointermove", e => {
        if (!this.mol) return;
        if (pts.has(e.pointerId)) pts.set(e.pointerId, [e.offsetX, e.offsetY]);
        if (pinch && pts.size === 2) {
          const [a, b] = [...pts.values()];
          this.zoom = Math.max(0.2, Math.min(8, pinch.z * Math.hypot(a[0] - b[0], a[1] - b[1]) / (pinch.d || 1)));
          this.draw();
          return;
        }
        if (last && pts.size === 1) {
          const dx = e.offsetX - last.x, dy = e.offsetY - last.y;
          if (last.pan) { this.pan[0] += dx; this.pan[1] += dy; }
          else this.rot = mul(mul(rotY(dx * 0.01), rotX(dy * 0.01)), this.rot);
          last.x = e.offsetX; last.y = e.offsetY;
          this._setHover(-1);
          this.draw();
          return;
        }
        this._setHover(this._pick(e.offsetX, e.offsetY), e.offsetX, e.offsetY);
      });
      const end = e => { pts.delete(e.pointerId); if (pts.size < 2) pinch = null; if (!pts.size) last = null; };
      cv.addEventListener("pointerup", end);
      cv.addEventListener("pointercancel", end);
      cv.addEventListener("pointerleave", () => this._setHover(-1));
      cv.addEventListener("wheel", e => {
        if (!this.mol) return;
        e.preventDefault();
        this.zoom = Math.max(0.2, Math.min(8, this.zoom * Math.exp(-e.deltaY * 0.0015)));
        this.draw();
      }, { passive: false });
      cv.addEventListener("dblclick", () => this.resetView());
      cv.addEventListener("keydown", e => {
        if (!this.mol) return;
        const k = { ArrowLeft: [0, -0.15], ArrowRight: [0, 0.15], ArrowUp: [-0.15, 0], ArrowDown: [0.15, 0] }[e.key];
        if (k) { e.preventDefault(); this.rot = mul(mul(rotY(k[1]), rotX(k[0])), this.rot); this.draw(); }
        else if (e.key === "+" || e.key === "=") { this.zoom = Math.min(8, this.zoom * 1.15); this.draw(); }
        else if (e.key === "-") { this.zoom = Math.max(0.2, this.zoom / 1.15); this.draw(); }
        else if (e.key === "0") this.resetView();
      });
    }

    _setHover(i, mx, my) {
      if (i !== this.hover) { this.hover = i; this.draw(); }
      if (i < 0) { this.tip.hidden = true; return; }
      const a = this.mol.atoms[i];
      this.tip.textContent = "";
      const head = document.createElement("strong");
      head.textContent = a.label;
      const meta = document.createElement("span");
      meta.textContent = `원자 ${i + 1} (${a.name}), 조각 ${a.frag}, 전하 ${a.q.toFixed(4)}`;
      this.tip.append(head, meta);
      if (this.mol.flagged.has(i)) {
        const w = document.createElement("em");
        w.textContent = "경고 대상 원자";
        this.tip.append(w);
      }
      this.tip.hidden = false;
      const tw = this.tip.offsetWidth, th = this.tip.offsetHeight;
      this.tip.style.left = Math.min(this.w - tw - 6, mx + 14) + "px";
      this.tip.style.top = Math.max(6, my - th - 10) + "px";
    }

    toPNG() { return this.canvas.toDataURL("image/png"); }
  }

  root.MolViewer = MolViewer;
  root.MolViewer.TYPE_PALETTE = TYPE_PALETTE;
})(typeof self !== "undefined" ? self : this);
