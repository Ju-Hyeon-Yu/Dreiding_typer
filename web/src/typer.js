// dreiding_typer.py 를 그대로 옮긴 브라우저/Node 용 DREIDING 라벨 타이퍼
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.DreidingTyper = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const SYMBOLS = ["", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",
    "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
    "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"];

  class TypingError extends Error {}

  // ---------------------------------------------------------------- 규칙표 (Python DreidingFF 대응)
  class DreidingFF {
    constructor(tables) {
      this.labels = tables.labels;
      this.pairs = tables.pairs;
      this.theta0 = tables.theta0;
    }
    _pair(a, b) { return this.pairs[a < b ? a + "|" + b : b + "|" + a]; }
    bondOrder(a, b) { const p = this._pair(a, b); return p ? p[0] : null; }
    bondR0(a, b) { const p = this._pair(a, b); return p ? p[1] : null; }
    hasDihedral(a, b) { const p = this._pair(a, b); return !!(p && p[2]); }
  }

  // ---------------------------------------------------------------- 분자 그래프
  // graph = {atoms: [{sym, chg, aromatic}], bonds: [{a, b, order}]}  (order: 1, 2, 3, 1.5=방향족)
  function graphFromRDKit(mol) {
    const j = JSON.parse(mol.get_json());
    const m = j.molecules[0];
    const def = Object.assign({ atom: { z: 6, chg: 0 }, bond: { bo: 1 } }, j.defaults || {});
    const ext = (m.extensions || []).find(e => e.name === "rdkitRepresentation") || {};
    const aroA = new Set(ext.aromaticAtoms || []);
    const aroB = new Set(ext.aromaticBonds || []);
    const atoms = m.atoms.map((a, i) => {
      const z = a.z !== undefined ? a.z : def.atom.z;
      return { sym: SYMBOLS[z], chg: a.chg !== undefined ? a.chg : (def.atom.chg || 0), aromatic: aroA.has(i) };
    });
    const bonds = (m.bonds || []).map((b, k) => ({
      a: b.atoms[0], b: b.atoms[1],
      order: aroB.has(k) ? 1.5 : (b.bo !== undefined ? b.bo : def.bond.bo),
    }));
    return finishGraph({ atoms, bonds });
  }

  function finishGraph(g) {
    g.nbrs = g.atoms.map(() => []);      // [이웃 원자, 결합 차수] (결합 순서대로)
    g.bonds.forEach(b => { g.nbrs[b.a].push([b.b, b.order]); g.nbrs[b.b].push([b.a, b.order]); });
    g.deg = g.nbrs.map(x => x.length);
    return g;
  }

  /** perm[new] = old 로 원자 순서를 바꾼 그래프 (RDKit RenumberAtoms 와 같이 결합 순서는 유지). */
  function renumber(g, perm) {
    const inv = new Array(perm.length);
    perm.forEach((old, nw) => { inv[old] = nw; });
    return finishGraph({
      atoms: perm.map(old => Object.assign({}, g.atoms[old])),
      bonds: g.bonds.map(b => ({ a: inv[b.a], b: inv[b.b], order: b.order })),
    });
  }

  /** 연결 성분. RDKit GetMolFrags 처럼 가장 작은 원자 번호 순으로 정렬. */
  function fragments(g) {
    const seen = new Array(g.atoms.length).fill(false), frags = [];
    for (let s = 0; s < g.atoms.length; s++) {
      if (seen[s]) continue;
      const fr = [], st = [s];
      seen[s] = true;
      while (st.length) {
        const u = st.pop();
        fr.push(u);
        for (const [v] of g.nbrs[u]) if (!seen[v]) { seen[v] = true; st.push(v); }
      }
      frags.push(fr.sort((x, y) => x - y));
    }
    return frags;
  }

  function formalCharge(g) { return g.atoms.reduce((s, a) => s + a.chg, 0); }

  // ---------------------------------------------------------------- 타입 규칙
  const sp2Class = lab => lab.includes("_2") || lab.includes("_R") || lab.includes("_1");

  function resonanceGroups(g) {
    const groups = [];
    g.atoms.forEach((atom, i) => {
      if ((atom.sym !== "C" && atom.sym !== "N") || atom.aromatic) return;
      const term = g.nbrs[i].filter(([o]) => g.atoms[o].sym === "O" && g.deg[o] === 1);
      const hasDouble = term.some(([, od]) => od === 2);
      const hasAnion = term.some(([o, od]) => od === 1 && g.atoms[o].chg === -1);
      if (hasDouble && hasAnion) groups.push([i, term.map(([o]) => o)]);
    });
    return groups;
  }

  function isSp2Atom(g, i, resAtoms) {
    if (g.atoms[i].aromatic || resAtoms.has(i)) return true;
    return g.nbrs[i].some(([, o]) => o > 1.0);
  }

  function fail(g, i) {
    throw new TypingError(`원자 ${i + 1} (${g.atoms[i].sym}, 이웃 ${g.deg[i]}개): DREIDING 라벨 규칙이 없습니다`);
  }

  function baseLabel(g, i, resAtoms, warnings) {
    const sym = g.atoms[i].sym, deg = g.deg[i];
    const orders = g.nbrs[i].map(([, o]) => o);
    const has2 = orders.includes(2), has3 = orders.includes(3);
    const resonant = g.atoms[i].aromatic || resAtoms.has(i);
    const anySp2Nbr = () => g.nbrs[i].some(([n]) => isSp2Atom(g, n, resAtoms));

    if (sym === "H") return ["H", ""];
    if (sym === "C") {
      if (resonant) return ["C_R", ""];
      const m = { 4: "C_3", 3: "C_2", 2: "C_1" }[deg];
      if (!m) fail(g, i);
      return [m, ""];
    }
    if (sym === "N") {
      if (resonant) return ["N_R", `_d${deg - 1}`];
      if (deg === 4) return ["N_3", ""];
      if (has3) return ["N_1", ""];
      if (has2) return ["N_2", `_d${deg - 1}`];
      if (deg === 3) return anySp2Nbr() ? ["N_R", "_d2"] : ["N_3", ""];
      if (deg === 2) return ["N_R", "_d1"];
      fail(g, i);
    }
    if (sym === "O") {
      if (resonant) return ["O_R", ""];
      if (deg === 1) {
        const nb = g.atoms[g.nbrs[i][0][0]].sym;
        if (nb === "S" || nb === "P" || has2) return ["O_2", ""];
        if (has3) return ["O_1", ""];
        return ["O_3", ""];
      }
      if (deg === 3) warnings.push(`원자 ${i + 1} O: 3배위 O 는 DREIDING 에 없어 O_3 로 둡니다`);
      return ["O_3", ""];
    }
    if (sym === "F" || sym === "Cl" || sym === "Br" || sym === "I") {
      if (sym === "Cl") warnings.push("dreiding.lt 의 Cl 질량이 34.45 로 들어가 있습니다(실제 35.45)");
      return [sym, ""];
    }
    if (sym === "S") {
      if (deg > 2 || has2) {
        warnings.push(`원자 ${i + 1} S: ${deg}배위/다중결합 S 는 DREIDING 원 논문에 없어 S_3 로 둡니다 ` +
          "(평형각 92.1°, 이면각 장벽 나눗셈 가정이 실제와 다름)");
      }
      return ["S_3", ""];
    }
    if (sym === "P") return ["P_3", `_d${deg - 1}`];
    if (sym === "Si") return ["Si_3", ""];
    if (deg === 0) return [sym, ""];
    fail(g, i);
  }

  function hFlag(g, i, resAtoms) {
    const sym = g.atoms[i].sym;
    if (sym !== "N" && sym !== "O" && sym !== "F") return "";
    if (g.nbrs[i].some(([n]) => g.atoms[n].sym === "H")) return "_hd";
    if (sym === "O" || sym === "F") return "_ha";
    const deg = g.deg[i];
    if (g.atoms[i].chg > 0) return "";
    if (g.atoms[i].aromatic) return deg === 2 ? "_ha" : "";
    if (resAtoms.has(i)) return "";
    const orders = g.nbrs[i].map(([, o]) => o);
    if (orders.includes(3) || orders.includes(2) || deg === 2) return "_ha";
    if (deg === 3 && !g.nbrs[i].some(([n]) => isSp2Atom(g, n, resAtoms))) return "_ha";
    return "";
  }

  const fmtG = x => String(Number(x.toPrecision(6)));   // Python "{:g}"

  /** graph → {labels, warnings}. Python assign_labels 와 같은 순서·동률 처리. */
  function assignLabels(g, ff, maxSteps = 1000000) {
    const warnings = [];
    const resAtoms = new Set(), resBonds = new Set();
    for (const [c, os] of resonanceGroups(g)) {
      resAtoms.add(c);
      os.forEach(o => { resAtoms.add(o); resBonds.add(c < o ? `${c},${o}` : `${o},${c}`); });
    }
    const n = g.atoms.length;
    const cands = [], wantedH = [];
    for (let i = 0; i < n; i++) {
      const [eh, d] = baseLabel(g, i, resAtoms, warnings);
      const h = hFlag(g, i, resAtoms);
      wantedH.push(h);
      const bflags = sp2Class(eh) ? ["", "_b1", "_b2"] : [""];
      const hs = h ? [h, ""] : [""];
      const c = [];
      for (const hh of hs) for (const b of bflags) if ((eh + b + d + hh) in ff.labels) c.push(eh + b + d + hh);
      if (!c.length) throw new TypingError(`원자 ${i + 1} (${g.atoms[i].sym}): '${eh}${d}' 계열 라벨이 dreiding.lt 에 없습니다`);
      cands.push(c);
    }
    for (let i = 0; i < n; i++) {
      if (g.atoms[i].sym === "H" && g.deg[i] === 1 && ["N", "O", "F"].includes(g.atoms[g.nbrs[i][0][0]].sym)) {
        cands[i] = ["H_HB"];
      }
    }

    // 기대 결합 차수 (결합 순서대로, 키는 작은 번호 먼저)
    const exp = new Map();
    for (const b of g.bonds) {
      const i = Math.min(b.a, b.b), j = Math.max(b.a, b.b), key = `${i},${j}`;
      exp.set(key, { i, j, e: (b.order === 1.5 || resBonds.has(key)) ? 1.5 : b.order });
    }
    const adj = Array.from({ length: n }, () => []);
    for (const { i, j } of exp.values()) { adj[i].push(j); adj[j].push(i); }

    const multi = [];
    for (let i = 0; i < n; i++) if (cands[i].length > 1) multi.push(i);
    const order = [], seen = new Set();
    for (const s of multi) {
      if (seen.has(s)) continue;
      const q = [s];
      seen.add(s);
      while (q.length) {
        const u = q.shift();
        order.push(u);
        for (const v of adj[u]) if (!seen.has(v) && cands[v].length > 1) { seen.add(v); q.push(v); }
      }
    }

    let chosen = cands.map(c => (c.length === 1 ? c[0] : null));
    const cost = (i, lab) => {
      let c = lab.includes("_b1") ? 1.0 : lab.includes("_b2") ? 1.5 : 0.0;
      if (wantedH[i] && !lab.endsWith(wantedH[i])) c += 10.0;
      return c;
    };
    const deg = g.deg;
    const bondOk = (i, j) => {
      const a = chosen[i], b = chosen[j];
      if (deg[i] > 1 && deg[j] > 1 && !ff.hasDihedral(a, b)) return false;
      if (sp2Class(a) && sp2Class(b)) {
        const o = ff.bondOrder(a, b);
        return o !== null && Math.abs(o - exp.get(i < j ? `${i},${j}` : `${j},${i}`).e) < 1e-6;
      }
      return true;
    };
    const ok = i => adj[i].every(j => chosen[j] === null || bondOk(i, j));

    let best = null, bestCost = Infinity, steps = 0;
    const solve = (k, acc) => {
      if (acc >= bestCost) return;
      if (k === order.length) { best = chosen.slice(); bestCost = acc; return; }
      const i = order[k];
      const sorted = cands[i].map((lab, idx) => [cost(i, lab), idx, lab]).sort((x, y) => x[0] - y[0] || x[1] - y[1]);
      for (const [c, , lab] of sorted) {
        steps += 1;
        if (steps > maxSteps) break;
        chosen[i] = lab;
        if (ok(i)) solve(k + 1, acc + c);
      }
      chosen[i] = null;
    };
    solve(0, 0.0);
    if (best === null) {
      const bad = [];
      for (const { i, j } of exp.values()) {
        if (cands[i].length === 1 && cands[j].length === 1 && !bondOk(i, j)) bad.push(`${i + 1}-${j + 1}`);
      }
      throw new TypingError("결합 차수와 이면각 규칙을 모두 만족하는 _b1/_b2 플래그 조합이 없습니다 " +
        "(DREIDING 라벨 체계로 표현할 수 없는 결합 패턴)" + (bad.length ? ` — 고정 라벨끼리 맞지 않는 결합: ${bad.join(", ")}` : ""));
    }
    if (steps > maxSteps) warnings.push("플래그 탐색 단계 한도에 걸려 최소 비용 해가 아닐 수 있습니다(결합 차수는 만족)");
    chosen = best;

    const flagged = new Set();   // 경고와 관련된 원자 (뷰어 강조용)
    for (const { i, j, e } of exp.values()) {
      const imp = ff.bondOrder(chosen[i], chosen[j]);
      if (imp === null) throw new TypingError(`결합 ${i + 1}-${j + 1} (${chosen[i]}–${chosen[j]}): dreiding.lt 에 결합 규칙이 없습니다`);
      if (Math.abs(imp - e) > 1e-6) {
        warnings.push(`결합 ${i + 1}-${j + 1} (${chosen[i]}–${chosen[j]}): 실제 차수 ${fmtG(e)} 를 DREIDING 에서 ${fmtG(imp)} 로 둡니다`);
        flagged.add(i); flagged.add(j);
      }
    }
    chosen.forEach((lab, i) => {
      if (lab.startsWith("N_3") && deg[i] === 4) {
        warnings.push(`원자 ${i + 1} N: 4배위 N(${lab})은 dreiding.lt 가 이웃 3개를 가정해 이 N 을 지나는 이면각 장벽이 1.5배로 들어갑니다`);
        flagged.add(i);
      }
    });
    wantedH.forEach((h, i) => {
      if (h && !chosen[i].endsWith(h)) {
        warnings.push(`원자 ${i + 1}: '${chosen[i]}${h}' 라벨이 없어 수소결합 플래그 ${h} 를 뺐습니다`);
        flagged.add(i);
      }
    });
    g.atoms.forEach((a, i) => {
      if (a.sym === "S" && (deg[i] > 2 || g.nbrs[i].some(([, o]) => o === 2))) flagged.add(i);
      if (a.sym === "O" && deg[i] === 3) flagged.add(i);
    });
    return { labels: chosen, warnings, flagged: [...flagged].sort((x, y) => x - y) };
  }

  return { SYMBOLS, TypingError, DreidingFF, graphFromRDKit, finishGraph, renumber, fragments, formalCharge,
    assignLabels, sp2Class };
});
