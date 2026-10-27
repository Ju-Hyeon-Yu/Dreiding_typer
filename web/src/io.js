// mol2·pdb·xyz 읽기, SMILES 그래프와 원자 대응, moltemplate .lt 쓰기 (dreiding_lt.py 대응)
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory(require("./typer.js"));
  else root.DreidingIO = factory(root.DreidingTyper);
})(typeof self !== "undefined" ? self : this, function (T) {
  "use strict";

  const ELEMENTS = new Set(T.SYMBOLS.filter(Boolean));
  const GAFF_PREFIX = [["cl", "Cl"], ["br", "Br"], ["i", "I"], ["f", "F"], ["c", "C"], ["n", "N"],
    ["o", "O"], ["s", "S"], ["p", "P"], ["h", "H"]];
  const isLowerWord = s => /[a-z]/.test(s) && s === s.toLowerCase();
  const capitalize = s => s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();

  class InputError extends Error {}

  /** mol2 원자 타입(SYBYL 또는 GAFF) → 원소, 안 되면 원자 이름 앞 글자 (Python _element 와 같음). */
  function mol2Element(name, atype) {
    const base = atype.split(".")[0];
    if (ELEMENTS.has(base)) return base;
    if (isLowerWord(atype) || /^[a-z]/.test(atype)) {
      for (const [pre, el] of GAFF_PREFIX) if (atype.startsWith(pre)) return el;
    }
    const m = /^[A-Za-z]+/.exec(name);
    if (m) {
      for (const cand of [capitalize(m[0].slice(0, 2)), m[0].slice(0, 1).toUpperCase()]) {
        if (ELEMENTS.has(cand)) return cand;
      }
    }
    throw new InputError(`mol2 원자 '${name} ${atype}' 의 원소를 알 수 없습니다`);
  }

  // 공유결합 반지름 (Å, Cordero 2008) — pdb(CONECT 없음)·xyz 결합 판정용
  const RCOV = { H: 0.31, B: 0.84, C: 0.76, N: 0.71, O: 0.66, F: 0.57, Si: 1.11, P: 1.07, S: 1.05, Cl: 1.02,
    Se: 1.20, Br: 1.20, I: 1.39, Li: 1.28, Na: 1.66, K: 2.03, Mg: 1.41, Ca: 1.76, Al: 1.21, Zn: 1.22, Fe: 1.32,
    Ga: 1.22, Ge: 1.20, As: 1.19, In: 1.42, Sn: 1.39, Sb: 1.39, Te: 1.38 };

  /** 반환: {kind, atoms:[{name, el, xyz, q}], bonds:[[i,j]] | null, hasCharge} — 인덱스는 0부터. */
  function parseMol2(text) {
    const atoms = [], rawBonds = [], ids = new Map();
    let sec = null, allQ = true;
    for (const line of text.split(/\r?\n/)) {
      const s = line.trim();
      if (s.startsWith("@<TRIPOS>")) { sec = s.slice(9); continue; }
      if (!s || s.startsWith("#")) continue;
      const t = s.split(/\s+/);
      if (sec === "ATOM") {
        if (t.length < 6) throw new InputError(`mol2 ATOM 줄 형식 오류: ${s}`);
        ids.set(parseInt(t[0], 10), atoms.length);
        const hasQ = t.length > 8;
        if (!hasQ) allQ = false;
        atoms.push({ name: t[1], el: mol2Element(t[1], t[5]), xyz: t.slice(2, 5).map(Number),
          q: hasQ ? parseFloat(t[8]) : 0.0 });
      } else if (sec === "BOND") {
        rawBonds.push([parseInt(t[1], 10), parseInt(t[2], 10)]);
      }
    }
    if (!atoms.length) throw new InputError("mol2 의 ATOM 섹션이 비어 있습니다");
    const bonds = rawBonds.map(([a, b]) => {
      if (!ids.has(a) || !ids.has(b)) throw new InputError(`mol2 BOND 가 없는 원자 번호를 가리킵니다: ${a}-${b}`);
      return [ids.get(a), ids.get(b)];
    });
    return { kind: "mol2", atoms, bonds, hasCharge: true, partialCharge: !allQ };
  }

  function pdbElement(line) {
    const el = line.slice(76, 78).trim();
    if (el) return capitalize(el);
    const name = line.slice(12, 16).trim().replace(/[^A-Za-z]/g, "");
    for (const cand of [capitalize(name.slice(0, 2)), name.slice(0, 1).toUpperCase()]) {
      if (ELEMENTS.has(cand)) return cand;
    }
    throw new InputError(`pdb 원자 '${line.slice(12, 16).trim()}' 의 원소를 알 수 없습니다`);
  }

  function parsePdb(text) {
    const atoms = [], serial = new Map(), pairs = new Set(), bonds = [];
    let sawConect = false;
    for (const line of text.split(/\r?\n/)) {
      const rec = line.slice(0, 6).trim();
      if (rec === "ATOM" || rec === "HETATM") {
        serial.set(parseInt(line.slice(6, 11), 10), atoms.length);
        atoms.push({ name: line.slice(12, 16).trim(), el: pdbElement(line),
          xyz: [line.slice(30, 38), line.slice(38, 46), line.slice(46, 54)].map(Number), q: 0.0 });
      } else if (rec === "CONECT") {
        sawConect = true;
        const nums = [];
        for (const k of [6, 11, 16, 21, 26]) {
          const v = parseInt(line.slice(k, k + 5), 10);
          if (!Number.isNaN(v)) nums.push(v);
        }
        const a = nums[0];
        for (const b of nums.slice(1)) {
          if (!serial.has(a) || !serial.has(b)) continue;
          const i = serial.get(a), j = serial.get(b), key = i < j ? `${i},${j}` : `${j},${i}`;
          if (i !== j && !pairs.has(key)) { pairs.add(key); bonds.push([Math.min(i, j), Math.max(i, j)]); }
        }
      }
    }
    if (!atoms.length) throw new InputError("pdb 에 ATOM/HETATM 줄이 없습니다");
    if (atoms.some(a => a.xyz.some(Number.isNaN))) throw new InputError("pdb 좌표를 읽지 못했습니다");
    return { kind: "pdb", atoms, bonds: sawConect ? bonds : null, hasCharge: false };
  }

  function parseXyz(text) {
    const lines = text.split(/\r?\n/);
    const n = parseInt(lines[0], 10);
    if (!(n > 0)) throw new InputError("xyz 첫 줄에 원자 수가 없습니다");
    const atoms = [];
    for (let k = 2; k < 2 + n; k++) {
      const t = (lines[k] || "").trim().split(/\s+/);
      if (t.length < 4) throw new InputError(`xyz ${k + 1}번째 줄 형식 오류`);
      let el = /^\d+$/.test(t[0]) ? T.SYMBOLS[parseInt(t[0], 10)] : capitalize(t[0].replace(/[^A-Za-z]/g, ""));
      if (!ELEMENTS.has(el)) throw new InputError(`xyz 원소 '${t[0]}' 를 알 수 없습니다`);
      atoms.push({ name: `${el}${k - 1}`, el, xyz: t.slice(1, 4).map(Number), q: 0.0 });
    }
    if (atoms.some(a => a.xyz.some(Number.isNaN))) throw new InputError("xyz 좌표를 읽지 못했습니다");
    return { kind: "xyz", atoms, bonds: null, hasCharge: false };
  }

  function parseCoordFile(name, text) {
    const ext = (name.split(".").pop() || "").toLowerCase();
    if (ext === "mol2") return parseMol2(text);
    if (ext === "pdb" || ext === "ent") return parsePdb(text);
    if (ext === "xyz") return parseXyz(text);
    throw new InputError(`지원하지 않는 파일 형식입니다: .${ext} (mol2, pdb, xyz)`);
  }

  /** 결합 정보가 없는 파일: 템플릿에 실제로 있는 원소 쌍만, 공유결합 반지름 합 × 1.2 이내를 결합으로 본다. */
  function bondsByDistance(file, template) {
    const allowed = new Set(template.bonds.map(b => {
      const x = template.atoms[b.a].sym, y = template.atoms[b.b].sym;
      return x < y ? `${x}|${y}` : `${y}|${x}`;
    }));
    const A = file.atoms, bonds = [];
    for (let i = 0; i < A.length; i++) {
      for (let j = i + 1; j < A.length; j++) {
        const x = A[i].el, y = A[j].el;
        if (!allowed.has(x < y ? `${x}|${y}` : `${y}|${x}`)) continue;
        const d = Math.hypot(A[i].xyz[0] - A[j].xyz[0], A[i].xyz[1] - A[j].xyz[1], A[i].xyz[2] - A[j].xyz[2]);
        const cut = 1.2 * ((RCOV[x] || 1.5) + (RCOV[y] || 1.5));
        if (d > 0.4 && d <= cut) bonds.push([i, j]);
      }
    }
    return bonds;
  }

  /** 원소+연결 그래프 동형 대응: 반환 t2f[템플릿 원자] = 파일 원자. 없으면 null. */
  function matchGraphs(tEl, tAdj, fEl, fAdj) {
    const n = tEl.length;
    if (fEl.length !== n) return null;
    // 두 그래프를 합쳐 같은 코드북으로 원자 불변량을 정제 (원소·이웃 수 → 이웃 라벨 반복)
    const allEl = tEl.concat(fEl);
    const allAdj = tAdj.concat(fAdj.map(a => a.map(j => j + n)));
    let lab = allEl.map((e, i) => `${e}:${allAdj[i].length}`);
    for (let it = 0; it < 6; it++) {
      const next = lab.map((l, i) => l + "(" + allAdj[i].map(j => lab[j]).sort().join(",") + ")");
      const ids = new Map();
      lab = next.map(l => { if (!ids.has(l)) ids.set(l, String(ids.size)); return ids.get(l); });
    }
    const tl = lab.slice(0, n), fl = lab.slice(n);
    const count = arr => arr.reduce((m, x) => m.set(x, (m.get(x) || 0) + 1), new Map());
    const ct = count(tl), cf = count(fl);
    for (const [k, v] of ct) if (cf.get(k) !== v) return null;

    // 탐색 순서: 연결을 따라가며 (이미 정한 원자와 이웃인 원자 우선)
    const order = [], seen = new Array(n).fill(false);
    const byRarity = [...Array(n).keys()].sort((a, b) => ct.get(tl[a]) - ct.get(tl[b]) || a - b);
    for (const s of byRarity) {
      if (seen[s]) continue;
      const q = [s];
      seen[s] = true;
      while (q.length) {
        const u = q.shift();
        order.push(u);
        for (const v of tAdj[u]) if (!seen[v]) { seen[v] = true; q.push(v); }
      }
    }
    const fAdjSet = fAdj.map(a => new Set(a));
    const byLabel = new Map();
    fl.forEach((l, j) => { if (!byLabel.has(l)) byLabel.set(l, []); byLabel.get(l).push(j); });
    const t2f = new Array(n).fill(-1), used = new Array(n).fill(false);
    let steps = 0;
    const go = k => {
      if (k === n) return true;
      if (++steps > 2000000) throw new InputError("원자 대응 탐색이 너무 길어 중단했습니다");
      const u = order[k];
      const mappedNbrs = tAdj[u].filter(v => t2f[v] >= 0);
      let cands = byLabel.get(tl[u]);
      if (mappedNbrs.length) cands = cands.filter(j => fAdjSet[t2f[mappedNbrs[0]]].has(j));
      for (const j of cands) {
        if (used[j]) continue;
        if (!mappedNbrs.every(v => fAdjSet[t2f[v]].has(j))) continue;
        t2f[u] = j; used[j] = true;
        if (go(k + 1)) return true;
        t2f[u] = -1; used[j] = false;
      }
      return false;
    };
    return go(0) ? t2f : null;
  }

  /** SMILES 그래프와 좌표 파일 대응 → perm[파일 순번] = 템플릿 원자. */
  function matchFile(template, file) {
    const n = template.atoms.length;
    let bonds = file.bonds;
    if (!bonds) bonds = bondsByDistance(file, template);
    const fAdj = file.atoms.map(() => []);
    for (const [i, j] of bonds) { fAdj[i].push(j); fAdj[j].push(i); }
    if (file.atoms.length !== n || bonds.length !== template.bonds.length) {
      const how = file.bonds ? "" : " (결합은 원자 간 거리로 판정)";
      throw new InputError(`원자/결합 수가 다릅니다: SMILES ${n}/${template.bonds.length}, ` +
        `${file.kind} ${file.atoms.length}/${bonds.length}${how}. 수소가 모두 들어 있는지 확인하십시오`);
    }
    const tAdj = template.atoms.map((_, i) => template.nbrs[i].map(([j]) => j));
    const t2f = matchGraphs(template.atoms.map(a => a.sym), tAdj, file.atoms.map(a => a.el), fAdj);
    if (!t2f) throw new InputError(`SMILES 와 ${file.kind} 의 연결 구조가 일치하지 않습니다 (원소 판정이나 결합 확인 필요)`);
    const perm = new Array(n);
    t2f.forEach((f, t) => { perm[f] = t; });
    return perm;
  }

  // ---------------------------------------------------------------- .lt 쓰기 (Python write_lt 와 같은 형식)
  const pad = (s, w) => (s.length >= w ? s : s + " ".repeat(w - s.length));
  const fix = (x, d, w) => {
    let s = x.toFixed(d);
    if (Object.is(x, -0)) s = "-" + s;
    return s.length >= w ? s : " ".repeat(w - s.length) + s;
  };

  /** rows: [{name, frag, label, q, xyz}], bonds: [[i,j]] (i<j, 정렬됨). */
  function writeLt(name, smiles, rows, bonds, chargeNote, tool = "dreiding_typer.html") {
    const out = [];
    out.push(`# ${tool} 생성 — SMILES: ${smiles}\n# 전하: ${chargeNote}\n\n`);
    out.push('import "dreiding.lt"\n\n');
    out.push(`${name} inherits DREIDING {\n\n`);
    out.push('  write("Data Atoms") {\n');
    for (const r of rows) {
      out.push(`    $atom:${pad(r.name, 6)} $mol:f${r.frag} @atom:${pad(r.label, 14)} ${fix(r.q, 6, 10)} ` +
        `${fix(r.xyz[0], 5, 11)} ${fix(r.xyz[1], 5, 11)} ${fix(r.xyz[2], 5, 11)}\n`);
    }
    out.push("  }\n\n");
    out.push('  write("Data Bond List") {\n');
    bonds.forEach(([i, j], k) => {
      out.push(`    $bond:b${pad(String(k + 1), 4)} $atom:${pad(rows[i].name, 6)} $atom:${rows[j].name}\n`);
    });
    out.push("  }\n\n}\n");
    return out.join("");
  }

  return { InputError, mol2Element, parseMol2, parsePdb, parseXyz, parseCoordFile, bondsByDistance,
    matchGraphs, matchFile, writeLt };
});
