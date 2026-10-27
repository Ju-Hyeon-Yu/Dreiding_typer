// SMILES(+좌표 파일) → 라벨·좌표·전하·.lt 까지 한 번에 만드는 처리 흐름 (dreiding_lt.py main 대응)
(function (root, factory) {
  if (typeof module === "object" && module.exports) {
    module.exports = factory(require("./typer.js"), require("./io.js"), require("./geom.js"));
  } else root.DreidingPipeline = factory(root.DreidingTyper, root.DreidingIO, root.DreidingGeom);
})(typeof self !== "undefined" ? self : this, function (T, IO, G) {
  "use strict";

  const STEREO = /[@\\/]/;

  /**
   * opts: {smiles, name, file: {name, text} | null, seed}
   * 반환: {name, rows, bonds, labels, warnings, flagged, lt, graph, source}
   * 실패하면 T.TypingError 또는 IO.InputError 를 던진다.
   */
  function run(R, ff, opts) {
    const smiles = (opts.smiles || "").trim();
    if (!smiles) throw new IO.InputError("SMILES 를 입력하십시오");
    if (/\s/.test(smiles)) throw new IO.InputError("SMILES 에 공백이 있습니다 (이온쌍은 '.' 으로 연결)");
    const name = (opts.name || "MOL").trim();
    if (!/^[A-Za-z_][A-Za-z0-9_]*$/.test(name)) {
      throw new IO.InputError("분자 이름은 영문자·숫자·밑줄만 쓸 수 있고 숫자로 시작할 수 없습니다");
    }
    const mol = R.get_mol(smiles);
    if (!mol || !mol.is_valid()) { if (mol) mol.delete(); throw new IO.InputError("SMILES 를 읽지 못했습니다"); }
    const warnings = [];
    let graph, xyz, q, chargeNote, source, xy = null;
    try {
      mol.add_hs_in_place();
      graph = T.graphFromRDKit(mol);
      const n = graph.atoms.length;
      if (opts.file) {
        const file = IO.parseCoordFile(opts.file.name, opts.file.text);
        const perm = IO.matchFile(graph, file);
        graph = T.renumber(graph, perm);         // 이후 번호는 파일 순서 (Python 과 같음)
        xyz = file.atoms.map(a => a.xyz);
        q = file.atoms.map(a => a.q);
        source = file.kind;
        if (file.kind === "mol2") {
          if (file.partialCharge) warnings.push("mol2 일부 원자에 전하 열이 없어 0 으로 둡니다");
          chargeNote = `mol2 (${opts.file.name})`;
          const qsum = q.reduce((s, v) => s + v, 0), formal = T.formalCharge(graph);
          if (Math.abs(qsum - formal) > 1e-3) {
            warnings.push(`mol2 전하 합 ${qsum.toFixed(4)} 이 SMILES 형식 전하 ${formal} 와 다릅니다`);
          }
        } else {
          chargeNote = `없음 — 모두 0 (${file.kind} 에는 전하 정보가 없음)`;
          warnings.push(`${file.kind} 에는 전하 정보가 없어 전하를 모두 0 으로 썼습니다`);
        }
      } else {
        source = "smiles";
        q = new Array(n).fill(0.0);
        chargeNote = "없음 — 모두 0 (RESP 는 mol2 로 넣으세요)";
        warnings.push("좌표 파일이 없어 전하를 모두 0 으로 썼습니다");
      }
      if (source === "smiles") xy = G.coords2D(mol);
    } finally {
      mol.delete();
    }

    const typed = T.assignLabels(graph, ff);
    warnings.push(...typed.warnings);
    if (source === "smiles") {
      xyz = G.embed3D(graph, typed.labels, ff, xy, opts.seed || 1);
      warnings.push("3D 좌표는 간이 최소화로 만든 시작 구조입니다 — MD 전에 에너지 최소화를 하십시오");
      if (STEREO.test(smiles)) warnings.push("SMILES 의 입체 정보(@, /, \\)는 3D 좌표 생성에 반영되지 않습니다");
    }

    const frags = T.fragments(graph);
    const fragOf = new Array(graph.atoms.length);
    frags.forEach((fr, k) => fr.forEach(i => { fragOf[i] = k + 1; }));
    const rows = graph.atoms.map((a, i) => ({ name: `${a.sym}${i + 1}`, frag: fragOf[i], label: typed.labels[i],
      q: q[i], xyz: xyz[i], sym: a.sym }));
    const bonds = graph.bonds.map(b => [Math.min(b.a, b.b), Math.max(b.a, b.b)])
      .sort((u, v) => u[0] - v[0] || u[1] - v[1]);
    const lt = IO.writeLt(name, smiles, rows, bonds, chargeNote);
    return { name, rows, bonds, labels: typed.labels, warnings: [...new Set(warnings)], flagged: typed.flagged,
      lt, graph, source, nFrags: frags.length, ff };
  }

  return { run };
});
